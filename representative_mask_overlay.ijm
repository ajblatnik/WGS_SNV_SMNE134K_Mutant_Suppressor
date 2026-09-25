// =====================================================================
// representative_mask_overlay.ijm
//
// Generates a representative figure showing where intensity values are
// sampled from: the protein channel (with the same levels/gamma you use
// in the batch export script) with the nuclear and/or cytoplasmic mask
// regions overlaid as semi-transparent color fills.
//
// Produces three files for ONE chosen image:
//   <baseName>_nuclearOverlay.tif    - protein channel + nuclear fill only
//   <baseName>_cytoplasmOverlay.tif  - protein channel + cytoplasmic fill only
//   <baseName>_combinedOverlay.tif   - protein channel + both fills together
//
// Run this on a single representative image (it prompts you to pick the
// raw .czi and its matching _NucleiMask.tif / _CytoplasmMask.tif from
// the segmentation macro's output folder).
// =====================================================================

var PROTEIN_CHANNEL = 4;

// Match these to whatever you used in the batch export script so this
// figure looks consistent with your other figures.
var PROTEIN_MIN   = 2400;
var PROTEIN_MAX   = 4500;
var PROTEIN_GAMMA = 1.0;

// Fill colors (standard ImageJ color names: red, green, blue, cyan,
// magenta, yellow, orange, pink, white, black, gray) and opacity (0-100).
var NUCLEI_FILL_COLOR = "cyan";
var CYTO_FILL_COLOR   = "magenta";
var FILL_OPACITY      = 45;

var OUTPUT_DPI = 300;
// ======================================================

rawPath      = File.openDialog("Choose the representative .czi image");
nucMaskPath  = File.openDialog("Choose the matching *_NucleiMask.tif");
cytoMaskPath = File.openDialog("Choose the matching *_CytoplasmMask.tif");
outputDir    = getDirectory("Choose an output folder");
baseName     = File.getNameWithoutExtension(rawPath);

setBatchMode(true);

// ---- Build the background: protein channel with levels+gamma, grayscale ----
open(rawPath);
rawID = getImageID();
selectImage(rawID);
run("Duplicate...", "title=protein_bg duplicate channels=" + PROTEIN_CHANNEL);
bgID = getImageID();
setMinAndMax(PROTEIN_MIN, PROTEIN_MAX);
run("8-bit");
run("Gamma...", "value=" + PROTEIN_GAMMA);
run("Grays");
selectImage(rawID);
close();

// ---- Nuclear selection, stored in ROI Manager slot 0 ----
open(nucMaskPath);
nucMaskID = getImageID();
selectImage(nucMaskID);
setThreshold(1, 65535);
setOption("BlackBackground", true);
run("Convert to Mask");
nucMaskID = getImageID();
selectImage(nucMaskID);
run("Create Selection");
roiManager("reset");
roiManager("add");
selectImage(nucMaskID);
close();

// ---- Cytoplasmic selection, stored in ROI Manager slot 1 ----
open(cytoMaskPath);
cytoMaskID = getImageID();
selectImage(cytoMaskID);
setThreshold(1, 65535);
setOption("BlackBackground", true);
run("Convert to Mask");
cytoMaskID = getImageID();
selectImage(cytoMaskID);
run("Create Selection");
roiManager("add");
selectImage(cytoMaskID);
close();

// ---- Nuclear-only overlay ----
selectImage(bgID);
run("Duplicate...", "title=bg_nuclear");
nuclearOutID = getImageID();
run("RGB Color");
nuclearOutID = getImageID();
roiManager("select", 0);
run("Overlay Options...", "stroke=none fill=" + NUCLEI_FILL_COLOR + " opacity=" + FILL_OPACITY);
run("Add Selection...");
run("Flatten");
flatNuclearID = getImageID();
run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
saveAs("Tiff", outputDir + baseName + "_nuclearOverlay.tif");

// ---- Cytoplasmic-only overlay ----
selectImage(bgID);
run("Duplicate...", "title=bg_cyto");
cytoOutID = getImageID();
run("RGB Color");
cytoOutID = getImageID();
roiManager("select", 1);
run("Overlay Options...", "stroke=none fill=" + CYTO_FILL_COLOR + " opacity=" + FILL_OPACITY);
run("Add Selection...");
run("Flatten");
flatCytoID = getImageID();
run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
saveAs("Tiff", outputDir + baseName + "_cytoplasmOverlay.tif");

// ---- Combined overlay (both regions together) ----
selectImage(bgID);
run("Duplicate...", "title=bg_combined");
combinedID = getImageID();
run("RGB Color");
combinedID = getImageID();
roiManager("select", 0);
run("Overlay Options...", "stroke=none fill=" + NUCLEI_FILL_COLOR + " opacity=" + FILL_OPACITY);
run("Add Selection...");
roiManager("select", 1);
run("Overlay Options...", "stroke=none fill=" + CYTO_FILL_COLOR + " opacity=" + FILL_OPACITY);
run("Add Selection...");
run("Flatten");
flatCombinedID = getImageID();
run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
saveAs("Tiff", outputDir + baseName + "_combinedOverlay.tif");

roiManager("reset");
close("*");
setBatchMode(false);

print("Representative mask overlay images saved to: " + outputDir);
