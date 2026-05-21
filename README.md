# Mutant Suppressor WGS Variant Analysis

R scripts for identifying and prioritizing candidate suppressor variants from whole-genome sequencing (WGS) of *Smn*^E134K^ mutant suppressor cell lines. The analysis integrates variant filtering, gene essentiality data, GO-term annotation, and differential gene expression to rank candidate genes that may suppress SMN loss-of-function.

---

## Scripts

### `mutsupp_wgs_variant.R`

Main variant processing and candidate prioritization pipeline.

**Parts:**

| Part | Description |
|------|-------------|
| 1 | Import BGI SNP annotation files; filter out parental-line variants; annotate cross-line overlap; generate Venn diagrams and Jaccard similarity heatmaps |
| 2 | Build GRCm38/mm10 ↔ GRCm39/mm11 cross-reference annotation table via biomaRt; collapse GO-terms (BP, CC, MF) to one row per gene |
| 3 | Map SNP coordinates to gene intervals (data.table `foverlaps`); assign ENU/O6-BG mutation priority scores; integrate RNA-seq variant and expression data |
| 4 | Multi-step candidate filtering (zygosity, dbSNP, shared-line count, variant consequence, essentiality, expression); GO enrichment with `clusterProfiler`; network analysis; SMN-associated GO-term filtering |
| 5 | Rank top-10 candidates per mutant suppressor line by SMN GO-term association frequency and DepMap essentiality score |
| 6 | Differential gene expression integration; volcano, density, and CDF plots stratified by variant consequence class; Wilcoxon rank-sum tests |

---

### `essentiality.R`

Builds the gene essentiality reference table used in `mutsupp_wgs_variant.R` (Part 3 onward).

**Parts:**

| Part | Description |
|------|-------------|
| 1 | Downloads IMPC DR15 viability data; maps lethal/subviable/viable phenotypes to gene symbols; joins to collapsed GO annotation |
| 2 | Downloads DepMap CRISPR gene-effect scores; summarises mean effect and likelihood of essentiality per gene; joins to IMPC data |
| 3 | Defines `core_go_id` — GO terms enriched in lethal or highly essential genes |
| 4 | Annotates GRCm39 genes with a `core_fxn` flag for GO-based essentiality |

> **Note:** Run `essentiality.R` first to produce `grcm3839_final_master.tsv`, which is required by `mutsupp_wgs_variant.R`.

---

## Requirements

### R packages

```r
# Data wrangling
library(readr); library(dplyr); library(tidyr); library(tidyverse)
library(stringr); library(data.table); library(fuzzyjoin)

# Genomics / annotation
library(biomaRt); library(GO.db); library(GenomicRanges)
library(clusterProfiler); library(DESeq2)

# Visualisation
library(ggplot2); library(ggpubr); library(gridExtra)
library(gtable); library(patchwork); library(ggpattern)
library(ggvenn); library(pheatmap); library(ggpointdensity)
library(igraph); library(ggraph)

# Other
library(MASS); library(jsonlite)
```

Install Bioconductor packages with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("biomaRt", "GO.db", "GenomicRanges", "clusterProfiler",
                        "DESeq2", "AnnotationDbi"))
```

### External data (required before running)

| File | Source | Used in |
|------|--------|---------|
| `bgi_variant_tables/*.snp.annot.csv` | BGI WGS SNP annotation output | `mutsupp_wgs_variant.R` Part 1 |
| `candidate_list/SNPs/SNPs-equalInParents-variableInTreated.csv` | RNA-seq variant analysis | `mutsupp_wgs_variant.R` Part 3 |
| `candidate_list/Expression/diffexp-T2_v_P2-results.csv` | Differential expression results | `mutsupp_wgs_variant.R` Parts 3 & 6 |
| IMPC viability data | Downloaded automatically from EBI FTP | `essentiality.R` Part 1 |
| DepMap CRISPR gene-effect scores | Downloaded automatically from DepMap portal API | `essentiality.R` Part 2 |

---

## Setup and Usage

### 1. Configure paths

At the top of each script, set `data_dir` and `output_dir` to match your local directory structure:

**`mutsupp_wgs_variant.R`**
```r
data_dir   <- "/path/to/your/input/data"
output_dir <- "/path/to/your/output/folder"
```

**`essentiality.R`**
```r
output_dir <- "/path/to/your/output/folder"
```

Output subdirectories (`biomaRt/`, `wgs_r_tables/`, `wgs_r_figures/`, etc.) are created automatically.

### 2. Run order

```
1. essentiality.R          → produces biomaRt/grcm3839_final_master.tsv
2. mutsupp_wgs_variant.R   → produces all variant tables, figures, and candidate lists
```

---

## Expected Output Structure

```
output_dir/
├── biomaRt/
│   ├── grcm38_gene_annotation.tsv
│   ├── grcm39_gene_annotation.tsv
│   ├── grcm3839_gene_annotation.tsv
│   ├── grcm3839_collapsed_gene_annotation.tsv
│   ├── grcm3839_final_master.tsv          ← from essentiality.R
│   └── impc_viability.tsv                  ← from essentiality.R
├── suppressor_mutant_unique_variant_master.tsv
├── wgs_r_tables/
│   ├── annotated_variant_master.tsv
│   ├── *_snv_freq.tsv                      ← per-line SNV frequency tables
│   ├── descriptives/
│   │   ├── snv_func_breakdown.tsv
│   │   ├── snv_biotype_breakdown.tsv
│   │   └── variant_gex_function.tsv
│   └── candidates/
│       ├── X622_805_candidates.tsv
│       ├── X752_799_814_candidates.tsv
│       ├── smn_associated_candidates.tsv
│       └── top10_*_candidates.tsv
└── wgs_r_figures/
    ├── mutsupp_variant_venn_raw.png
    ├── mutsupp_variant_jacc-heat_raw.png
    ├── mutsupp_gene_venn.png
    ├── mutsupp_gene_jacc-heat.png
    ├── mutsupp_gene_filt_venn.png
    ├── mutsupp_gene_smn_venn.png
    ├── go_enrichment/
    │   ├── GO_Top10_*.png
    │   └── mutsupp_go_jacc-heat.png
    └── gex/
        ├── GEX_volcano.png
        ├── GEX_*_volcano_density_cdf.png
        └── Top10_GO_DEGs_Dotplot.png
```

---

## Key Analysis Decisions

- **Parental filtering:** Variants present in either parental line (E1 or E3) are excluded from all suppressor-line comparisons.
- **ENU priority scoring:** SNVs are weighted by expected ENU/O6-BG mutagenic spectrum (A:T → G:C transitions highest; G:C transversions lowest); scores are normalised by gene length.
- **Essentiality threshold:** Candidates are retained if IMPC viability = lethal **or** DepMap `likelihood_essential > 0.1`.
- **SMN GO association:** Candidates are cross-referenced against GO terms annotated to *Smn1* to prioritise functionally relevant suppressors.
- **Candidate split:** Lines 622 and 805 share a large proportion of variants (likely clonal or co-derived) and are analysed as a group separate from 752, 799, and 814.
