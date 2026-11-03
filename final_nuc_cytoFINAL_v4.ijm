var NUCLEI_CHANNEL  = 1; // numbers defined by .czi
var CELLS_CHANNEL   = 1;
var PROTEIN_CHANNEL = 4;

var BORDER_WIDTH = 3;      // telling which channel to use to define cell border
var BORDER_VALUE = 65535;  // pixel width of the line dawn around mask

// Resolution embedded in every saved TIFF. ImageJ writes no DPI info by
// default, which makes Illustrator assume 72 DPI when PLACING an image —
// making it appear huge relative to a normal artboard. This is display/
// print-size metadata only — it doesn't touch pixel data, masks, ROIs,
// or any of the intensity measurements below.
var OUTPUT_DPI = 300;


list = getFileList(inputDir);
setBatchMode(false);

for (f = 0; f < list.length; f++) {
    if (endsWith(list[f], ".czi")) {
        
        // Open the raw multi-channel image
        open(inputDir + list[f]);
        rawID = getImageID();
        // rawTitle = getTitle();
        baseName = File.nameWithoutExtension;
        
        // Folder for per-image output
        imageOutDir = outputDir + baseName + File.separator;
        File.makeDirectory(imageOutDir);
        File.copy(inputDir + list[f], imageOutDir + list[f]);
        
        // Extract nuclei channel
        selectImage(rawID);
        run("Duplicate...", "title=nuclei_chan duplicate channels=" + NUCLEI_CHANNEL);
        nucleiChanID = getImageID();
        nucleiChanTitle = getTitle();
        selectImage(nucleiChanID);
        run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
		saveAs("Tiff", imageOutDir + baseName + "_DAPI.tif");
        nucleiChanID = getImageID();   // saveAs renames the window; ID stays valid, refresh title too
        nucleiChanTitle = getTitle();
        
        // Extract whole cell channel
        selectImage(rawID);
        run("Duplicate...", "title=cells_chan duplicate channels=1-4");
        cellsChanID = getImageID();
        cellsChanTitle = getTitle();
        
        // Extract protein channel (in line 46, variables cell_diameter and min_size need to be defined empirically)
        // first line is to mask nuclei, second is to maske the whole cell (therfore cyt = whole cell - nuclei)
        selectImage(rawID);
        run("Duplicate...", "title=protein_chan duplicate channels=" + PROTEIN_CHANNEL);
        proteinChanID = getImageID();
        proteinChanTitle = getTitle();
        run("Subtract Background...", "rolling=250");
        
        titlesBefore = getList("image.titles");
        selectImage(nucleiChanID);
        run("Cellpose...", "cp_model=nuclei custom_model= cell_diameter=35 cyto_channel=None nuclei_channel=1 min_size=25 normalize=true resample=true return_rois=false cellprob_threshold=0.0 flow_threshold=0.60 tile_overlap=0.1 niter=0 compute_flows=false mode_3d=None stitch_threshold=0.0 flow3d_smooth=0 torchversion=cpu usegpu=false");
        nucleiMaskTitle = nucleiChanTitle;
        titlesAfter = getList("image.titles");
        for (i = 0; i < titlesAfter.length; i++) {
            isOld = false;
            for (j = 0; j < titlesBefore.length; j++) {
                if (titlesAfter[i] == titlesBefore[j]) isOld = true;
            }
            if (!isOld) nucleiMaskTitle = titlesAfter[i];
        }
        selectWindow(nucleiMaskTitle);
        nucleiMaskID = getImageID();
        // SAVE RAW CELLPOSE NUCLEI MASK
        run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
		saveAs("Tiff", imageOutDir + baseName + "_Cellpose_NucleiRaw.tif");
        
        titlesBefore = getList("image.titles");
        selectImage(cellsChanID);
        run("Cellpose...", "cp_model=cyto3 custom_model= cell_diameter=90 cyto_channel=4 nuclei_channel=1 min_size=80 normalize=true resample=true return_rois=false cellprob_threshold=0.0 flow_threshold=0.4 tile_overlap=0.1 niter=0 compute_flows=false mode_3d=None stitch_threshold=0.0 flow3d_smooth=0 torchversion=cpu usegpu=false");
        cellsMaskTitle = cellsChanTitle;
        titlesAfter = getList("image.titles");
        for (i = 0; i < titlesAfter.length; i++) {
            isOld = false;
            for (j = 0; j < titlesBefore.length; j++) {
                if (titlesAfter[i] == titlesBefore[j]) isOld = true;
            }
            if (!isOld) cellsMaskTitle = titlesAfter[i];
        }
        selectWindow(cellsMaskTitle);
        cellsMaskID = getImageID();
        // SAVE RAW CELLPOSE Cell MASK
        run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
		saveAs("Tiff", imageOutDir + baseName + "_Cellpose_CellsRaw.tif");
        
        // nucleus matching (ensures that both masks match)
        roiManager("reset");
        selectImage(nucleiMaskID);
        selectImage(nucleiMaskID);
		print("nucleiMaskID = " + nucleiMaskID);
		print("Current image = " + getTitle());
        run("Label image to ROIs");
        nNuclei = roiManager("count");

        selectImage(cellsMaskID);
        cellsW = getWidth();
        cellsH = getHeight();
        newImage("nuclei_relabeled", "16-bit black", cellsW, cellsH, 1);
        nucleiRelabeledID = getImageID();

        orphanNuclei = 0;
        dupCollisions = 0;
        seenCellID = newArray(0); 

        for (i = 0; i < nNuclei; i++) {
            roiManager("select", i);
            Roi.getContainedPoints(xpts, ypts);
            px = xpts[0]; py = ypts[0];

            selectImage(cellsMaskID);
            parentCellID = getPixel(px, py);

            if (parentCellID == 0) {
                orphanNuclei++; // nucleus with no underlying cell segmentation; excluded
            } else {
                isDup = false;
                for (k = 0; k < seenCellID.length; k++) {
                    if (seenCellID[k] == parentCellID) isDup = true;
                }
                if (isDup) {
                    dupCollisions++; // two nuclei mapped to the same cell; only first kept
                } else {
                    seenCellID = Array.concat(seenCellID, parentCellID);
                    roiManager("select", i);
                    selectImage(nucleiRelabeledID);
                    setColor(parentCellID);
                    fill();
                }
            }
        }
        roiManager("reset");

        if (orphanNuclei > 0 || dupCollisions > 0) {
            print(baseName + ": " + orphanNuclei + " nuclei with no matching cell, " +
                  dupCollisions + " duplicate nucleus/cell collisions (excluded from analysis).");
        }

        // cytoplasmic mask subtraction
        imageCalculator("Subtract create", cellsMaskID, "nuclei_relabeled");
        cytoMaskTitle = getTitle();
        cytoMaskID = getImageID();

        // nuclear mask saved to folder
        selectImage(nucleiRelabeledID);
        run("Select None");
        run("Duplicate...", "title=nuclei_relabeled_bordered");
        borderedID = getImageID();

        roiManager("reset");
        run("Label image to ROIs");
        nBorderROIs = roiManager("count");
        setColor(BORDER_VALUE);
        setLineWidth(BORDER_WIDTH);
        for (b = 0; b < nBorderROIs; b++) {
            roiManager("select", b);
            run("Draw", "slice");
        }
        roiManager("reset");

        selectImage(borderedID);
        run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
        saveAs("Tiff", imageOutDir + baseName + "_NucleiMask.tif");
        close();
        selectImage(cellsMaskID);
        run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
        saveAs("Tiff", imageOutDir + baseName + "_CellsMask.tif");
        cellsMaskID = getImageID();
        selectImage(cytoMaskID);
        run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
        saveAs("Tiff", imageOutDir + baseName + "_CytoplasmMask.tif");
        cytoMaskID = getImageID();

        // set analysis measurements (if want to do more, edit line 157 and within the for-loop 171)
        run("Set Measurements...", "mean integrated redirect=None decimal=3");

        // nuclear measurement
        selectImage(nucleiRelabeledID);
        run("Label image to ROIs");

        nucCellIDs = newArray(0);
        nucMeansArr = newArray(0);
        nNucROIs = roiManager("count");
        if (nNucROIs > 0) {
            selectImage(proteinChanID);
            roiManager("deselect");
            roiManager("Measure");

            for (i = 0; i < nNucROIs; i++) {
                roiManager("select", i);
                Roi.getContainedPoints(xpts, ypts);
                selectImage(nucleiRelabeledID);
                cellID = getPixel(xpts[0], ypts[0]);

                nucCellIDs = Array.concat(nucCellIDs, cellID);
                nucMeansArr = Array.concat(nucMeansArr, getResult("Mean", i));

                setResult("Label", i, baseName + "_Cell_" + cellID);
                setResult("Compartment", i, "Nucleus");
            }
            updateResults();
            saveAs("Results", imageOutDir + baseName + "_Nuclei_Results.csv");
            run("Clear Results");
        }
        roiManager("reset");

        // cytoplasmic measurement (see nuclei above)
        selectImage(cytoMaskID);
        run("Label image to ROIs");

        cytoCellIDs = newArray(0);
        cytoMeansArr = newArray(0);
        nCytoROIs = roiManager("count");
        if (nCytoROIs > 0) {
            selectImage(proteinChanID);
            roiManager("deselect");
            roiManager("Measure");

            for (i = 0; i < nCytoROIs; i++) {
                roiManager("select", i);
                Roi.getContainedPoints(xpts, ypts);
                selectImage(cytoMaskID);
                cellID = getPixel(xpts[0], ypts[0]);

                cytoCellIDs = Array.concat(cytoCellIDs, cellID);
                cytoMeansArr = Array.concat(cytoMeansArr, getResult("Mean", i));

                setResult("Label", i, baseName + "_Cell_" + cellID);
                setResult("Compartment", i, "Cytoplasm");
            }
            updateResults();
            saveAs("Results", imageOutDir + baseName + "_Cytoplasm_Results.csv");
            run("Clear Results");
        }
        roiManager("reset");

        // calculate and save n/c ratio
        pairRow = 0;
        unmatchedCount = 0;
        for (i = 0; i < nucCellIDs.length; i++) {
            matchIdx = -1;
            for (j = 0; j < cytoCellIDs.length; j++) {
                if (cytoCellIDs[j] == nucCellIDs[i]) { matchIdx = j; }
            }
            if (matchIdx >= 0) {
                setResult("Label", pairRow, baseName + "_Cell_" + nucCellIDs[i]);
                setResult("Nucleus_Mean", pairRow, nucMeansArr[i]);
                setResult("Cytoplasm_Mean", pairRow, cytoMeansArr[matchIdx]);
                setResult("NC_Ratio", pairRow, nucMeansArr[i] / cytoMeansArr[matchIdx]);
                pairRow++;
            } else {
                unmatchedCount++;
            }
        }
        if (unmatchedCount > 0) {
            print(baseName + ": " + unmatchedCount + " cells had a nucleus but no matched cytoplasm ROI (excluded from N:C table).");
        }
        if (pairRow > 0) {
            updateResults();
            saveAs("Results", imageOutDir + baseName + "_NC_Ratio_Results.csv");
            run("Clear Results");
        }

        run("Close All");
    }
}

setBatchMode(false);
waitForUser("Batch Finished!", "All images have been processed! Check your output folder for CSV files.");