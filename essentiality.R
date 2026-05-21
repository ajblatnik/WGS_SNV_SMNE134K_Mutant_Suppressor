
#----------------------------
# This script associates GO-terms with phenotypic data correlated with cell essentiality or viability
#----------------------------

# SET OUTPUT DIRECTORY
# Change this to the root directory where output subdirectories (biomaRt/, etc.) should be written.
output_dir <- "."   # default: current working directory

# PART 1: Download International Mouse Phenotyping Consortium (IMPC) Viability Data and Processing

  library(readr)
  library(dplyr)

  # specify the URL for the IMPC viability release (example: DR15)
    impc_url  <- "https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-15.0/results/viability.csv.gz"  # DR15 example :contentReference[oaicite:4]{index=4}
    impc_file <- tempfile(fileext = ".csv.gz")

    download.file(impc_url, impc_file)
    impc_viab <- read_csv(impc_file)

  # Quick check of key columns
    glimpse(impc_viab)

  # Filter & clean: pull the gene + viability status
    impc_viab2 <- impc_viab %>% 
      dplyr::select(`Gene Accession Id`, `Gene Symbol`, `Viability Phenotype HOMs/HEMIs`) %>%
      mutate(
        viability_status = case_when(
          `Viability Phenotype HOMs/HEMIs` == "lethal" ~ "Lethal",
          `Viability Phenotype HOMs/HEMIs` == "subviable" ~ "Subviable",
          `Viability Phenotype HOMs/HEMIs` == "viable" ~ "Viable",
          TRUE ~ NA_character_
        )
      ) %>%
      distinct()
  # rename impc_viab2 to allow annotation with GO terms
    names(impc_viab2)[ names(impc_viab2) == "Gene Accession ID" ] <- "accession_id"
    names(impc_viab2)[ names(impc_viab2) == "Gene Symbol" ] <- "external_gene_name"
    names(impc_viab2)[ names(impc_viab2) == "Viability Phenotype HOMs/HEMIs" ] <- "impc_viability"
    impc_viab2 <- impc_viab2[,c(2:3)]
    
  # now cross-reference the Gene Symbols with GO terms
    # access BioMart and find attributes of the database
      ensembl <- useMart("ensembl")
      ensembl_m <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
      attributes <- listAttributes(ensembl_m)
    
    # this pulls information for the current build (GRCm39, mm11) but the coordinates may not translate so also need to pull the mm10 GRCm38 build
      grcm39 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id"),
                      mart = ensembl_m,
                      useCache = F)
  # use 'collapsed' from mutsupp_wgs_variant.R, which is essentially a collapsed version of grcm39
    impc_essential <- left_join(collapsed, impc_viab2, by = "external_gene_name")
    
  # write to directory
    write_tsv(as.data.frame(impc_essential), file = file.path(output_dir, "biomaRt", "impc_viability.tsv"), quote = "none")
    
    


# PART 2: DepMap/CRISPR Knockout Essentiality Data 

library(readr)
library(dplyr)

# Browse DepMap data page for latest release :contentReference[oaicite:5]{index=5}
    
  files_df <- jsonlite::fromJSON("https://depmap.org/portal/api/download/files") 
  
  gene_eff_file <- files_df %>%
    dplyr::filter(grepl("gene_effect", file_name, ignore.case = TRUE),
                  release == max(release)) %>%
    dplyr::pull(download_url) %>%
    .[1]  # first match
  
  tmpfile <- tempfile(fileext = ".csv")
  download.file(gene_eff_file, tmpfile)
  gene_effects <- readr::read_csv(tmpfile)
  
    
# Example: gene effect file
  depmap_url  <- "https://storage.googleapis.com/depmap-external-downloads/downloads-by-canonical-id/25q3-public-6202.1/CRISPRGeneEffect.csv?response-content-disposition=attachment%3B+filename%3D%22CRISPRGeneEffect.csv%22&GoogleAccessId=depmap-external-downloads%40broad-achilles.iam.gserviceaccount.com&Expires=1764194432&Signature=FI87Hjiq5J3Y6ITK6tOOBfvyY2rfYJ2mhPIAJAk6j9zD2RumIus4%252FD6nWYLcXJCK6q02v3XLza9BuR9F8g7g1%252FoScEbYHUhuaTX5gXLi6Dg2ABUiXR%252FAoth%252BL%252B8s5xECwc5ERQbCEd2rY4CvmALGy%252BKhPFMUQ7AckdW9LFnSDVPmugwaYXNI%252BUXbgs4rJmvIuvSXuwaeEr0kHoNSLa62ZlPpvP8d3gTXqHcurfl4lfnLq%252BOjPPncz0eCZH4SF4As55KkwEkKrnRB6SDDpiGlij1YL3I6nhTg%252FhwgXyUuHyNU0ufvwUNJcf9DPiQcNvT5M1vIEninabrS8Z2dj4%252FCJg%3D%3D&userProject=broad-achilles"  # you'll need to find exact filename from DepMap portal
  depmap_file <- tempfile(fileext = ".csv")

  download.file(depmap_url, depmap_file)
  depmap_ce   <- read_csv(depmap_file)
  
  # Clean column names: remove Entrez IDs in parentheses
    clean_names <- names(depmap_ce)
    clean_names <- str_replace(clean_names, " \\(.*\\)", "")  # remove " (####)"
    names(depmap_ce) <- clean_names
    
  # rename first column
    names(depmap_ce)[1] <- "DepMap_ID"
    
    depmap_long <- depmap_ce %>%
      pivot_longer(
        cols = -DepMap_ID,
        names_to = "gene_symbol",
        values_to = "gene_effect"
      )
    
    gene_summary <- depmap_long %>%
      group_by(gene_symbol) %>%
      summarise(
        mean_effect_viability = mean(gene_effect, na.rm = TRUE),
        n_lines = n(),
        likelihood_essential = mean(gene_effect < -0.5, na.rm = TRUE)
      ) %>%
      arrange(mean_effect_viability)
    # Mean effect near –1 → strong dependency → likely essential.
    # Mean effect near 0 → non-essential.
    # prop_essential close to 1 means the gene is required for survival in nearly all tested lines.
    
    # Create lowercase columns *before* joining
      gene_summary2 <- gene_summary %>% mutate(gene_symbol_lower = tolower(gene_symbol))
    
      impc_essential2 <- impc_essential %>% mutate(external_gene_name_lower = tolower(external_gene_name))
    
    # Join on lowercase gene names
      essential <- full_join(gene_summary2, impc_essential2, by = c("gene_symbol_lower" = "external_gene_name_lower"))
    # remove those rows without a ensembl_gene_id
      essential <- essential %>% filter(!is.na(ensembl_gene_id))
      essential <- essential[-c(1,3,5)]
      
    # reorder then write to directory
      names(essential)[ names(essential) == "cc_go_desc" ] <- "cc_go_term"
      names(essential)[ names(essential) == "mf_go_desc" ] <- "mf_go_term"
      names(essential)[ names(essential) == "unspecified_goID" ] <- "unsp_goID"
      names(essential)[ names(essential) == "unspecified_go_desc" ] <- "unsp_go_term"
    
      essential <- essential %>% dplyr::select("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "gene_length", "external_gene_name", "gene_biotype", "description", "go_id", "bp_goID", "bp_go_term", "cc_goID", "cc_go_term", "mf_goID", "mf_go_term", "unsp_goID", "unsp_go_term", "impc_viability", "mean_effect_viability", "likelihood_essential")
      write_tsv(as.data.frame(essential), file = file.path(output_dir, "biomaRt", "grcm3839_final_master.tsv"), quote = "none")
      

# PART 3: create a list of essential go_ids
      
  # determine a threshold for essential cuttof using impc_viability = lethal
    lethal <- essential %>% filter(impc_viability == "lethal")
    median(lethal$likelihood_essential, na.rm = TRUE)
    core_go_id <- essential %>% filter(impc_viability == "lethal" | likelihood_essential > 0.34)
    list_go_id <- core_go_id %>% dplyr::select("go_id")
    list_go_id <- list_go_id %>% separate_rows(go_id, sep = ";")
    list_go_id <- distinct(list_go_id)
    list_go_id <- list_go_id %>% filter(!is.na(go_id), go_id != "")
    

    
      

