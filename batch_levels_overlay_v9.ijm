// =====================================================================
// batch_levels_overlay.ijm
//
// Batch-processes all .czi images in a folder:
//   1. Applies FIXED levels (min/max) and gamma (contrast) to channel 1
//      (nuclei) and channel 4 (protein), converts each to 8-bit, saves
//      as a single-channel grayscale TIFF (ready to drop into
//      Illustrator).
//   2. Colorizes each channel with its named LUT, then additively
//      combines them into one flattened RGB pseudocolor overlay TIFF.
//
// HOW TO SET LEVELS (MIN/MAX):
//   Open one representative image in Fiji, select the nuclei channel,
//   go to Image > Adjust > Brightness/Contrast, and use "Set" to try
//   values until it looks right. Read off the Minimum/Maximum you
//   land on and enter them below. Do the same for the protein channel.
//   These are applied identically to every image in the batch, so all
//   images end up on the same scale and are directly comparable.
//
// HOW TO SET GAMMA (CONTRAST):
//   Gamma reshapes the tone curve on top of the min/max stretch above.
//   1.0 = no change. Values BELOW 1.0 brighten midtones (lower
//   contrast look). Values ABOVE 1.0 darken midtones / punch up
//   contrast. Try 0.7-1.5 as a starting range and adjust by eye on a
//   test image (Image > Adjust > Gamma) before committing a value here.
//
// HOW TO SET THE LUTS (PSEUDOCOLOR):
//   Any name from Image > Lookup Tables works here (not just primary
//   colors) — e.g. "Red", "Green", "Blue", "Magenta", "Cyan", "Grays",
//   "Fire", "Orange Hot", "Green Fire Blue", "Spectrum", etc. Check
//   that menu in Fiji for the exact spelling of the LUT you want.
// =====================================================================

var NUCLEI_CHANNEL  = 1;
var PROTEIN_CHANNEL = 4;

// ---- EDIT THESE for your dataset before running ----
// ---- Y12 MIN 400, MAX 5500, BD-SMN MIN , MAX 
var NUCLEI_MIN   = 400;
var NUCLEI_MAX   = 5500;
var NUCLEI_GAMMA = 1.0;

// ---- Y12 MIN 2600, MAX 7000, BD-SMN MIN 2400, MAX 4500
var PROTEIN_MIN   = 2400;
var PROTEIN_MAX   = 4500;
var PROTEIN_GAMMA = 1.0;

// Overlay pseudocolors — any Image > Lookup Tables name
var NUCLEI_LUT  = "Blue";
var PROTEIN_LUT = "Red";

// Resolution embedded in every saved TIFF. ImageJ writes no DPI info by
// default, which makes Illustrator assume 72 DPI when PLACING an image —
// making it appear huge relative to a normal artboard. This doesn't
// resample/change any pixel data, it only sets how big it displays.
var OUTPUT_DPI = 300;
// ======================================================

inputDir  = getDirectory("Choose the folder of .czi images");
outputDir = getDirectory("Choose an output folder");

list = getFileList(inputDir);
setBatchMode(true);

print("=== Batch levels/overlay export ===");
print("Nuclei channel " + NUCLEI_CHANNEL + ": min=" + NUCLEI_MIN + " max=" + NUCLEI_MAX + " gamma=" + NUCLEI_GAMMA + " lut=" + NUCLEI_LUT);
print("Protein channel " + PROTEIN_CHANNEL + ": min=" + PROTEIN_MIN + " max=" + PROTEIN_MAX + " gamma=" + PROTEIN_GAMMA + " lut=" + PROTEIN_LUT);

for (f = 0; f < list.length; f++) {
    if (endsWith(list[f], ".czi")) {

        print("Processing: " + list[f]);

        open(inputDir + list[f]);
        rawID = getImageID();
        baseName = File.nameWithoutExtension;

        imageOutDir = outputDir + baseName + File.separator;
        File.makeDirectory(imageOutDir);

        // ---- Nuclei channel: fixed levels + gamma -> 8-bit -> save grayscale ----
        selectImage(rawID);
        run("Duplicate...", "title=nuclei_ch duplicate channels=" + NUCLEI_CHANNEL);
        nucleiChID = getImageID();
        setMinAndMax(NUCLEI_MIN, NUCLEI_MAX);
        run("8-bit");
        run("Gamma...", "value=" + NUCLEI_GAMMA);
        nucleiChID = getImageID();
        run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
        saveAs("Tiff", imageOutDir + baseName + "_ch" + NUCLEI_CHANNEL + "_nuclei.tif");
        nucleiChID = getImageID();      // saveAs renames the window; refresh
        nucleiChTitle = getTitle();

        // ---- Protein channel: fixed levels + gamma -> 8-bit -> save grayscale ----
        selectImage(rawID);
        run("Duplicate...", "title=protein_ch duplicate channels=" + PROTEIN_CHANNEL);
        proteinChID = getImageID();
        setMinAndMax(PROTEIN_MIN, PROTEIN_MAX);
        run("8-bit");
        run("Gamma...", "value=" + PROTEIN_GAMMA);
        proteinChID = getImageID();
        run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
        saveAs("Tiff", imageOutDir + baseName + "_ch" + PROTEIN_CHANNEL + "_protein.tif");
        proteinChID = getImageID();     // saveAs renames the window; refresh
        proteinChTitle = getTitle();

        // ---- Pseudocolor overlay: colorize each grayscale channel with its
        //      LUT, bake into RGB, then additively combine ----
        selectImage(nucleiChID);
        run("Duplicate...", "title=nuclei_rgb");
        run(NUCLEI_LUT);
        run("RGB Color");
        nucleiRGBID = getImageID();

        selectImage(proteinChID);
        run("Duplicate...", "title=protein_rgb");
        run(PROTEIN_LUT);
        run("RGB Color");
        proteinRGBID = getImageID();

        imageCalculator("Add create", nucleiRGBID, proteinRGBID);
        overlayID = getImageID();
        run("Set Scale...", "distance=" + OUTPUT_DPI + " known=1 unit=inch");
        saveAs("Tiff", imageOutDir + baseName + "_overlay_ch" + NUCLEI_CHANNEL + "ch" + PROTEIN_CHANNEL + ".tif");

        // ---- Clean up before next image ----
        close("*");
    }
}

print("=== Batch complete ===");
setBatchMode(false);
