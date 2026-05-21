

#-----------------------------
# Perform filtering for variant analysis for suppression of SMNE134K loss of function
#-----------------------------

#-----------------------------
# PATH CONFIGURATION
# Set these two directories before running the script.
#   data_dir   : root folder containing input data subdirectories
#                (bgi_variant_tables/, candidate_list/SNPs/, candidate_list/Expression/)
#   output_dir : root folder where output subdirectories will be written
#                (biomaRt/, wgs_r_tables/, wgs_r_figures/, etc.)
# Output subdirectories are created automatically if they do not exist.
#-----------------------------
data_dir   <- "."   # e.g. "/path/to/data"
output_dir <- "."   # e.g. "/path/to/output"

# Create required output subdirectories
dirs_to_create <- c(
  file.path(output_dir, "biomaRt"),
  file.path(output_dir, "wgs_r_tables"),
  file.path(output_dir, "wgs_r_tables", "descriptives"),
  file.path(output_dir, "wgs_r_tables", "candidates"),
  file.path(output_dir, "wgs_r_figures"),
  file.path(output_dir, "wgs_r_figures", "go_enrichment"),
  file.path(output_dir, "wgs_r_figures", "gex")
)
invisible(lapply(dirs_to_create, dir.create, recursive = TRUE, showWarnings = FALSE))

#-----------------------------
# load libraries
#-----------------------------

library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(tidyverse)
library(fuzzyjoin) # only good for smaller datasets
library(data.table)

library(DESeq2)
library(ggpointdensity)
library(MASS)

library(biomaRt)
library(GO.db)
library(GenomicRanges)
library(clusterProfiler)

library(ggplot2)
library(ggpubr)
library(gridExtra)
library(gtable)
library(patchwork)
library(ggpattern)
library(ggvenn)
library(pheatmap)

#-----------------------------
# PART 1: Filter for Variants Unique to Mutant Suppressor Lines
#-----------------------------

# STEP 1: import snp_annot.csv for all cell lines sequenced by BGI
  
  # import
    X622 <- read_csv(file.path(data_dir, "bgi_variant_tables", "622.5.041.snp.annot.csv"))
    X752 <- read_csv(file.path(data_dir, "bgi_variant_tables", "752.33.snp.annot.csv"))
    X799 <- read_csv(file.path(data_dir, "bgi_variant_tables", "799.36.snp.annot.csv"))
    X805 <- read_csv(file.path(data_dir, "bgi_variant_tables", "805.044.snp.annot.csv"))
    X814 <- read_csv(file.path(data_dir, "bgi_variant_tables", "814.292.snp.annot.csv"))
    E1 <- read_csv(file.path(data_dir, "bgi_variant_tables", "E1.snp.annot.csv"))
    E3 <- read_csv(file.path(data_dir, "bgi_variant_tables", "E3.snp.annot.csv"))

  # remove placeholder columns 7-9, 11-14
    X622 <- X622[,c(1:6,10,15:22)]
    X752 <- X752[,c(1:6,10,15:22)]
    X799 <- X799[,c(1:6,10,15:22)]
    X805 <- X805[,c(1:6,10,15:22)]
    X814 <- X814[,c(1:6,10,15:22)]
    E1 <- E1[,c(1:6,10,15:22)]
    E3 <- E3[,c(1:6,10,15:22)]
    
    
# STEP 2: combine variants identified in both parental lines to make one file to filter against
  
  # mutate one parental file column to report "yes" if the same SNP is identified in both datasets
    parental_1 <- E1 %>% mutate(shared = if_else(paste(Chr, Start, End, Ref, Obs) %in% (E3 %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, "yes", "no"))  

    parental_3 <- E3 %>% mutate(shared = if_else(paste(Chr, Start, End, Ref, Obs) %in% (E1 %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, "yes", "no"))

  # remove shared "yes" SNPs from one file
    parental_3_unique <- parental_3 %>% filter(shared == "no")
  # append unique SNPs to the parental file
    parental <- bind_rows(parental_1, parental_3_unique)
    
    
# STEP 3: Filter out variants identically called in parent lines
    
  # annotate suppressor mutant files to indicate whether the identified SNPs are shared in parent lines
    X622 <- X622 %>% mutate(parental = if_else(paste(Chr, Start, End, Ref, Obs) %in% (parental %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, "yes", "no"))  

    X752 <- X752 %>% mutate(parental = if_else(paste(Chr, Start, End, Ref, Obs) %in% (parental %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, "yes", "no"))  

    X799 <- X799 %>% mutate(parental = if_else(paste(Chr, Start, End, Ref, Obs) %in% (parental %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, "yes", "no"))  

    X805 <- X805 %>% mutate(parental = if_else(paste(Chr, Start, End, Ref, Obs) %in% (parental %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, "yes", "no"))  

    X814 <- X814 %>% mutate(parental = if_else(paste(Chr, Start, End, Ref, Obs) %in% (parental %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, "yes", "no"))  

  # remove those SNPs shared with parent lines
    X622u <- X622 %>% filter(parental == "no")
    X752u <- X752 %>% filter(parental == "no")
    X799u <- X799 %>% filter(parental == "no")
    X805u <- X805 %>% filter(parental == "no")
    X814u <- X814 %>% filter(parental == "no")
    
    
# STEP 4: Determine if any SNPS are shared across the different lines
    
  # mutate column reporting 1 if the SNP is identified in another mutant suppressor line or 0 if not
    X622u <- X622u %>% mutate(in_622 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X622u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0:      0, 1: 198810
    X622u <- X622u %>% mutate(in_752 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X752u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 187864, 1:  10946 
    X622u <- X622u %>% mutate(in_799 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X799u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 183305, 1:  15505 
    X622u <- X622u %>% mutate(in_805 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X805u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0:  48877, 1: 149933
    X622u <- X622u %>% mutate(in_814 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X814u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 185679, 1:  13131
    
    X752u <- X752u %>% mutate(in_752 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X752u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0:      0, 1: 241251 
    X752u <- X752u %>% mutate(in_622 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X622u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 230305, 1:  10946 
    X752u <- X752u %>% mutate(in_799 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X799u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 229846, 1:  11405   
    X752u <- X752u %>% mutate(in_805 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X805u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 230273, 1:  10978
    X752u <- X752u %>% mutate(in_814 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X814u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 229909, 1:  11342 
    
    X799u <- X799u %>% mutate(in_799 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X799u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0:      0, 1: 159636 
    X799u <- X799u %>% mutate(in_622 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X622u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 144131, 1:  15505 
    X799u <- X799u %>% mutate(in_752 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X752u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 148231, 1:  11405  
    X799u <- X799u %>% mutate(in_805 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X805u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 143830, 1:  15806  
    X799u <- X799u %>% mutate(in_814 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X814u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 146456, 1:  13180
    
    X805u <- X805u %>% mutate(in_805 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X805u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0:      0, 1: 203126  
    X805u <- X805u %>% mutate(in_622 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X622u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0:  53193, 1: 149933  
    X805u <- X805u %>% mutate(in_752 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X752u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 192148, 1:  10978
    X805u <- X805u %>% mutate(in_799 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X799u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 187320, 1:  15806  
    X805u <- X805u %>% mutate(in_814 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X814u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 189885, 1:  13241
    
    X814u <- X814u %>% mutate(in_814 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X814u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0:      0, 1: 164722    
    X814u <- X814u %>% mutate(in_622 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X622u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 151591, 1:  13131    
    X814u <- X814u %>% mutate(in_752 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X752u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 153380, 1:  11342
    X814u <- X814u %>% mutate(in_799 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X799u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 151542, 1:  13180   
    X814u <- X814u %>% mutate(in_805 = if_else(paste(Chr, Start, End, Ref, Obs) %in% (X805u %>% transmute(k = paste(Chr, Start, End, Ref, Obs)))$k, 1, 0))  # 0: 151481, 1:  13241
    
  # mutate a column that sums the other mutant supressor columns
    X622u <- X622u %>% mutate(total_shared = in_622 + in_752 + in_799 + in_805 + in_814)

    X752u <- X752u %>% mutate(total_shared = in_622 + in_752 + in_799 + in_805 + in_814)

    X799u <- X799u %>% mutate(total_shared = in_622 + in_752 + in_799 + in_805 + in_814)

    X805u <- X805u %>% mutate(total_shared = in_622 + in_752 + in_799 + in_805 + in_814)

    X814u <- X814u %>% mutate(total_shared = in_622 + in_752 + in_799 + in_805 + in_814)

  # make a unique identifier column to compare variants across samples and generate a venn diagram
    X622u <- X622u %>% mutate(variant_id = paste0(Chr, ":", Start, "_", Ref, "->", Obs))
    X752u <- X752u %>% mutate(variant_id = paste0(Chr, ":", Start, "_", Ref, "->", Obs))
    X799u <- X799u %>% mutate(variant_id = paste0(Chr, ":", Start, "_", Ref, "->", Obs))
    X805u <- X805u %>% mutate(variant_id = paste0(Chr, ":", Start, "_", Ref, "->", Obs))
    X814u <- X814u %>% mutate(variant_id = paste0(Chr, ":", Start, "_", Ref, "->", Obs))
    
    
# STEP 5: plots for overlap and similarity
    
  # venn_diagram for overlap and independence
    mutsupvar <- list("622" = X622u$variant_id,
                      "752" = X752u$variant_id,
                      "799" = X799u$variant_id,
                      "805" = X805u$variant_id,
                      "814" = X814u$variant_id)
    
    mutsupp_var <- ggvenn(mutsupvar,
                          show_percentage = FALSE,
                          fill_color = c("#fa030375", "#e1005375", "#be007675", "#94009c75", "#040cd475"),
                          stroke_size = 0,
                          set_name_color =c("#fa0303", "#e10053", "#be0076", "#94009c", "#040cd4"),
                          set_name_size = 4,
                          text_color = "white",
                          text_size = 1.5)
    
    ggsave(filename = "mutsupp_variant_venn_raw.png", plot = mutsupp_var, device = "png", dpi = 300, height = 4, width = 4, path = file.path(output_dir, "wgs_r_figures"))
    
  # perform Fischer's exact test
    options(scipen = 999, digits = 22)
    
    pairs <- combn(names(mutsupvar),2, simplify=FALSE)
    
    pvals <- sapply(pairs, function(p){
      A <- mutsupvar[[p[1]]]
      B <- mutsupvar[[p[2]]]
      
      inter <- length(intersect(A,B))
      onlyA <- length(setdiff(A,B))
      onlyB <- length(setdiff(B,A))
      
      m <- matrix(c(inter, onlyA, onlyB, 0), nrow=2, byrow=TRUE)
      fisher.test(m)$p.value
    })
    
    padj <- p.adjust(pvals, method="BH")
    
    # Jaccard similarity
    jaccard <- sapply(pairs, function(p){
      A <- mutsupvar[[p[1]]]
      B <- mutsupvar[[p[2]]]
      length(intersect(A,B)) / length(union(A,B))
    })
    
    results <- data.frame(
      pair = sapply(pairs, paste, collapse=" vs "),
      fisher_p = pvals,
      fisher_p_adj = padj,
      jaccard = jaccard
    )
    
  # compute pairwise Jaccard matrix
    samples <- names(mutsupvar)
    jacmat <- matrix(0, nrow=length(samples), ncol=length(samples),
                     dimnames=list(samples, samples))
    
    for(i in samples){
      for(j in samples){
        A <- mutsupvar[[i]]
        B <- mutsupvar[[j]]
        jacmat[i,j] <- length(intersect(A,B)) / length(unique(c(A,B)))
      }
    }
    
    # heatmap
    my_breaks <- seq(0, 1, length.out = 101)  
    outpath <- path.expand(file.path(output_dir, "wgs_r_figures", "mutsupp_variant_jacc-heat_raw.png"))
      png(filename = outpath, width=2000, height=2000, res=300)
      pheatmap(jacmat,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               display_numbers = TRUE,
               number_format = "%.2f",
               fontsize = 12,
               cellwidth = 30,
               cellheight = 30,
               border_color = NA,
               treeheight_row = 20,
               treeheight_col =  20,
               fontsize_number = 10,
               gaps_row = NULL,
               gaps_col = NULL,
               angle_col = 0,
               number_color = "black",
               filename = outpath
               )
      
    
# STEP 6: generate a master table for all variants
    
    # use bind_rows to generate a table containing all variants and use distinct() to remove duplicate rows
      variants <- bind_rows(X622u, X752u, X799u, X805u, X814u) %>% distinct(variant_id, .keep_all = TRUE)

      # remove parental column
        variants <- variants %>% dplyr::select (-"parental")
      # reorder columns
        variants <- variants %>% dplyr::select("variant_id", everything())
    
    # write file to folder so don't have to download every time
      write_tsv(as.data.frame(variants), file = file.path(output_dir, "suppressor_mutant_unique_variant_master.tsv"), quote = "none")
      
      
      
      
      
#-----------------------------
# PART 2: Create Master Tables to Annotate with gene_ids and GO-terms
#-----------------------------  
      
# STEP 1: pull information from the biomaRt to annotate the variants with gene_ids and GO-terms

  # pull contempory mm11 build GO terms, but mm10 coordinates
    # access BioMart and find attributes of the database
      ensembl <- useMart("ensembl")
      ensembl_m <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
      attributes <- listAttributes(ensembl_m)
      
      # this pulls information for the current build (GRCm39, mm11) but the coordinates may not translate so also need to pull the mm10 GRCm38 build
        grcm39 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "description", "go_id", "name_1006", "namespace_1003"),
                        mart = ensembl_m,
                        useCache = F)
          x <- as.data.frame(table(grcm39$ensembl_gene_id))
          
        # provide a column for gene length
          grcm39 <- grcm39 %>% mutate(gene_length = abs(start_position - end_position))
      
    # connect to Nov 2020 Ensembl archive (GRCm38 / mm10)
      mart_grcm38 <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "https://nov2020.archive.ensembl.org")
      
      # pull coordinate information for the GRCm38/mm10 build
        grcm38 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "description", "go_id", "name_1006", "namespace_1003"),
                      mart = mart_grcm38,
                      useCache = F)
        # provide a column for gene length
          grcm38 <- grcm38 %>% mutate(gene_length = abs(start_position - end_position))
    
    # write file to folder so don't have to download every time
      write_tsv(as.data.frame(grcm39), file = file.path(output_dir, "biomaRt", "grcm39_gene_annotation.tsv"), quote = "none")
      write_tsv(as.data.frame(grcm38), file = file.path(output_dir, "biomaRt", "grcm38_gene_annotation.tsv"), quote = "none")
      

# STEP 2: make a master table for annotation
      
  # remove GO and descriptions from grcm38, remove coordinates from grcm39
    grcm38_nogo <- grcm38 %>% dplyr::select(-c(external_gene_name, gene_biotype, description, go_id, name_1006, namespace_1003))
    grcm39_nocoord <- grcm39 %>% dplyr::select(-c(chromosome_name, start_position, end_position, strand, gene_length))
    
    # isolate missing gene_ids and see if there are shared coordinates in grcm39
      lost_grcm38 <- grcm38_nogo %>% filter(!ensembl_gene_id %in% grcm39$ensembl_gene_id)
      # determine if any of these lost_grcm38 coordinates align with grcm39 annotated genes
        lost_grcm38 <- lost_grcm38 %>% mutate(in_grcm39 = if_else(paste(chromosome_name, start_position, end_position, strand) %in% (grcm39 %>% transmute(k = paste(chromosome_name, start_position, end_position, strand)))$k, "yes", "no"))  
        same_in_both <- lost_grcm38 %>% filter(in_grcm39 == "yes")  # use mm10 coordinates with mm11 gene_id and GO-terms
        # do same for grcm39
          lost_grcm39 <- grcm39 %>% mutate(in_grcm38 = if_else(paste(chromosome_name, start_position, end_position, strand) %in% (lost_grcm38 %>% transmute(k = paste(chromosome_name, start_position, end_position, strand)))$k, "yes", "no")) 
          same_in_both_39 <- lost_grcm39 %>% filter(in_grcm38 == "yes")
        x <- as.data.frame(table(same_in_both_39$ensembl_gene_id))
        
      # determine if any of these have altered coordinates, but same gene_names
        lost_grcm38 <- lost_grcm38 %>% mutate(name_in_grcm39 = if_else(external_gene_name %in% grcm39$external_gene_name, "yes", "no")) 
        same_name_in_both <- lost_grcm38 %>% filter(name_in_grcm39 == "yes")  # use mm10 coordinates with mm11 gene_id and GO-terms
        same_name_in_both_39 <- grcm39_nocoord %>% filter(external_gene_name %in% same_name_in_both$external_gene_name)
        
    # this table contains direct gene_id matches between mm10 (grcm38) and mm11 (grcm39) 
      master <- inner_join(grcm38_nogo, grcm39_nocoord, by = "ensembl_gene_id")
      master <- master %>% distinct() # this is in the format to build the final table
      
      # for coordinates with different gene_ids, want to update to mm11
        coord_38 <- same_in_both %>% mutate(coord_id = paste0(chromosome_name, "_", start_position, "-", end_position, "_", strand))
        coord_38 <- coord_38 %>% dplyr::select("chromosome_name", "start_position", "end_position", "strand", "coord_id", "gene_length")
        
        coord_39 <- same_in_both_39 %>% mutate(coord_id = paste0(chromosome_name, "_", start_position, "-", end_position, "_", strand))
        coord_39 <- coord_39 %>% dplyr::select("coord_id", "ensembl_gene_id", "external_gene_name", "gene_biotype", "description", "go_id", "name_1006", "namespace_1003")
      
        # inner_join() by coord_id
          coord <- inner_join(coord_39, coord_38, by = "coord_id")
          coord <- coord %>% dplyr::select(-c(1)) 
          coord <- coord %>% distinct()  # Formatted for use as input to bind_rows() final table construction.
      
      # for same gene_names but altered coordinates, use similar strategy as above
        name_38 <- same_name_in_both %>% dplyr::select("external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "gene_length")
        name_39 <- same_name_in_both %>% dplyr::select("ensembl_gene_id", "external_gene_name", "gene_biotype", "description", "go_id", "name_1006", "namespace_1003")
        
        # inner_join() by coord_id
          exname <- inner_join(name_39, name_38, by = "external_gene_name")
          exname <- exname %>% distinct()  # Formatted for use as input to bind_rows() final table construction.
          
          
# STEP 3: final build
          
  # use bind_rows() to make a single table      
    master <- bind_rows(master, coord, exname)  # Note: a small number of gene IDs present in GRCm38 have no direct equivalent in GRCm39.
    
    # for remaining, which have no clear correlates in mm11 build
      remaining <- grcm38 %>% filter(!ensembl_gene_id %in% master$ensembl_gene_id)
      remaining <- remaining %>% mutate(in_grcm39 = if_else(paste(chromosome_name, start_position, end_position, strand) %in% (grcm39 %>% transmute(k = paste(chromosome_name, start_position, end_position, strand)))$k, "yes", "no"))  
      remaining <- remaining %>% filter(in_grcm39 == "no")
      remaining <- remaining %>% mutate(name_in_grcm39 = if_else(external_gene_name %in% grcm39$external_gene_name, "yes", "no")) 
      remaining <- remaining %>% filter(name_in_grcm39 == "no")
      remaining <- remaining %>% dplyr::select(-c(13:14))
  
    master <- bind_rows(master, remaining)
  # write file to folder so don't have to download every time
    write_tsv(as.data.frame(master), file = file.path(output_dir, "biomaRt", "grcm3839_gene_annotation.tsv"), quote = "none")
    
    
# STEP 4: complete GO descriptions for go_ids
  
  # helper functions for retriving from the GO.db
    get_GO_term <- function(x) {Term(GOID(x))}
    
    get_GO_definition <- function(x) {Definition(GOID(x))}
    
    get_GO_ontology <- function(x) {Ontology(GOID(x))}
    
  # make a table of go_ids with missing descriptions
    missing_go_desc <- master %>% filter(is.na(name_1006) | name_1006 == "")
  # use above functions to retrieve missing information from the GO.db
    master <- master %>% mutate(name_1006 = case_when(!is.na(name_1006) & name_1006 != "" ~ name_1006, TRUE ~ map_chr(go_id, get_GO_term)))
    master <- master %>% mutate(namespace_1003 = case_when(!is.na(namespace_1003) & namespace_1003 != "" ~ namespace_1003, TRUE ~ map_chr(go_id, get_GO_ontology)))
    master <- master %>% mutate(namespace_1003 = case_when(namespace_1003 == "biological_process" ~ "BP", namespace_1003 == "molecular_function" ~ "MF", namespace_1003 == "cellular_component" ~ "CC", TRUE ~ namespace_1003))
    
    
# STEP 5: compress GO-terms to have a single row per gene_id 
    
  # make separate tables for each of the namespace_1003 column
    table(master$namespace_1003)
    # biological_process (BP)    cellular_component (CC)     molecular_function (MF) 

    master_bp <- master %>% filter(namespace_1003 == "BP")
    master_cc <- master %>% filter(namespace_1003 == "CC")
    master_mf <- master %>% filter(namespace_1003 == "MF")
    master_na <- master %>% filter(is.na(namespace_1003))
    # collapse tables so that multiple GO ids and descriptions are contained defined columns so one row per gene_id
      collapsed_bp <- master_bp %>% group_by(ensembl_gene_id, external_gene_name) %>% summarise(bp_goID = paste(unique(go_id), collapse = ";"), bp_go_term = paste(unique(name_1006), collapse = ";"), .groups = "drop")
      collapsed_cc <- master_cc %>% group_by(ensembl_gene_id, external_gene_name) %>% summarise(cc_goID = paste(unique(go_id), collapse = ";"), cc_go_desc = paste(unique(name_1006), collapse = ";"), .groups = "drop")
      collapsed_mf <- master_mf %>% group_by(ensembl_gene_id, external_gene_name) %>% summarise(mf_goID = paste(unique(go_id), collapse = ";"), mf_go_desc = paste(unique(name_1006), collapse = ";"), .groups = "drop")
      collapsed_na <- master_na %>% group_by(ensembl_gene_id, external_gene_name) %>% summarise(unspecified_goID = paste(unique(go_id), collapse = ";"), unspecified_go_desc = paste(unique(name_1006), collapse = ";"), .groups = "drop")
    # remove external_gene_name() to not duplicate column
      collapsed_bp <- collapsed_bp %>% dplyr::select(-c(2))
      collapsed_cc <- collapsed_cc %>% dplyr::select(-c(2))
      collapsed_mf <- collapsed_mf %>% dplyr::select(-c(2))
      collapsed_na <- collapsed_na %>% dplyr::select(-c(2))
      
    # collapse all associated go_ids from the master sheet
      collapsed_all <- master %>% dplyr::select(-c(11:12))
      collapsed_all <- collapsed_all %>% group_by(ensembl_gene_id, chromosome_name, start_position, end_position, strand, gene_length, external_gene_name, gene_biotype, description, ) %>% summarise(go_id = paste(unique(go_id), collapse = ";"), .groups = "drop")
  # combine using full_join()
      collapsed <- full_join(collapsed_bp, collapsed_cc, by = "ensembl_gene_id")
      collapsed <- full_join(collapsed, collapsed_mf, by = "ensembl_gene_id")
      collapsed <- full_join(collapsed, collapsed_na, by = "ensembl_gene_id")
      collapsed <- full_join(collapsed_all, collapsed, by = "ensembl_gene_id")
  # write file to folder so don't have to download every time
    write_tsv(as.data.frame(collapsed), file = file.path(output_dir, "biomaRt", "grcm3839_collapsed_gene_annotation.tsv"), quote = "none")
      
 
    
         
#-----------------------------
# PART 3: Combining variant table with reference tables for filtering
#-----------------------------     
    
# STEP 1: combine variant table with collapsed_master file
  
  # need to remove 'chr' preceding chromosome number
    variants <- variants %>% mutate(Chr = str_replace(Chr, "^chr", ""))
  # import grcm3839_final_master.tsv
    essential <- read_delim(file.path(output_dir, "biomaRt", "grcm3839_final_master.tsv"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
  # use data.table to map SNP coordinates to the appropriate gene_id
    setDT(variants)
    setDT(essential)
    
    setkey(essential, chromosome_name, start_position, end_position)
    
    annot <- foverlaps(
      variants[, .(Chr, Start, End = Start)],   # point overlap
      essential[, .(chromosome_name, start_position, end_position, ensembl_gene_id)],
      by.x = c("Chr", "Start", "End"),
      by.y = c("chromosome_name", "start_position", "end_position"),
      type = "within"
    )
  # compute width BEFORE removing columns
    annot[, gene_width := end_position - start_position]
    
  # choose best gene per variant (smallest gene span)
    annot <- annot[order(gene_width)][, .SD[1], by = .(Chr, Start)]
    
  # now mutate / attach gene_id to the original variants table
    variants[, ensembl_gene_id := annot$ensembl_gene_id[
      match(paste(Chr, Start), paste(annot$Chr, annot$Start))
    ]]
    
  # left_join() with collapsed
    annotated_variants <- left_join(variants, essential, by = "ensembl_gene_id")
  # clean up columns
    annotated_variants <- annotated_variants %>% dplyr::select(-c(3:5,9,24:27))
  
  # reorder for simplicity
    annotated_variants <- annotated_variants %>% dplyr::select("variant_id", "Chr", "Start", "End", "Ref", "Obs", "Func", "Codon_Change", "Impact", "dbSNP", "Otherinfo", "isoform", "in_622", "in_752", "in_799", "in_805", "in_814", "total_shared", "ensembl_gene_id", "external_gene_name", "gene_length", "gene_biotype", "description", "go_id", "bp_goID", "bp_go_term", "cc_goID", "cc_go_term", "mf_goID", "mf_go_term", "unsp_goID", "unsp_go_term", "impc_viability", "mean_effect_viability", "likelihood_essential")
    
    
# STEP 2: Edit table to include a priority scoring for ENU/O6-BG mutagenesis
    
  # ENU primarily mutates A/T bases, causing A->G/T->C transitions and A->T/T->A transversions
  # ENU strongly targets A:T base pairs. So any change where the REF is A or T should get higher prior weight.
  # Within A/T changes, A→G / T→C (A:T → G:C transitions) are especially common, and T→A / A→T (transversions) are also frequent.
  # G:C changes are less likely from ENU, so give them lower weight.
  # create:
  # enu_prior (0–1 weight)
  # correct weights for introduction of O6-BG, which increases the probability of G->A and C->T transitions.
    
  # first, calculate the frequency of specific nucleotide changes detected in mutagenized lines
    variants_dt <- copy(variants)  # just in case
      # should also do for individual lines
        X622u_dt <- as.data.table(X622u)
        X752u_dt <- as.data.table(X752u)
        X799u_dt <- as.data.table(X799u)
        X805u_dt <- as.data.table(X805u)
        X814u_dt <- as.data.table(X814u)
        
  # for simplicity
    x <- variants_dt
    x <- X622u_dt
    x <- X752u_dt
    x <- X799u_dt
    x <- X805u_dt
    x <- X814u_dt
    
    snv_freq <- x[
      nchar(Ref) == 1 & nchar(Obs) == 1 & Ref %in% c("A","T","C","G") & Obs %in% c("A","T","C","G"),
      .(count = .N),
      by = .(sub = paste0(Ref, "->", Obs))
    ][
      , total := sum(count)
    ][
      , freq := count / total
    ][
      order(-freq)
    ]
    
    variants_snv_freq <- snv_freq
    X622u_snv_freq <- snv_freq
    X752u_snv_freq <- snv_freq
    X799u_snv_freq <- snv_freq
    X805u_snv_freq <- snv_freq
    X814u_snv_freq <- snv_freq
    
    # write this to file to present
      write_tsv(as.data.frame(variants_snv_freq), file = file.path(output_dir, "wgs_r_tables", "variants_snv_freq.tsv"), quote = "none")
      write_tsv(as.data.frame(X622u_snv_freq), file = file.path(output_dir, "wgs_r_tables", "X622u_snv_freq.tsv"), quote = "none")
      write_tsv(as.data.frame(X752u_snv_freq), file = file.path(output_dir, "wgs_r_tables", "X752u_snv_freq.tsv"), quote = "none")
      write_tsv(as.data.frame(X799u_snv_freq), file = file.path(output_dir, "wgs_r_tables", "X799u_snv_freq.tsv"), quote = "none")
      write_tsv(as.data.frame(X805u_snv_freq), file = file.path(output_dir, "wgs_r_tables", "X805u_snv_freq.tsv"), quote = "none")
      write_tsv(as.data.frame(X814u_snv_freq), file = file.path(output_dir, "wgs_r_tables", "X814u_snv_freq.tsv"), quote = "none")
      
  # function for mutation priority scoring
    assign_enu_prior <- function(ref, obs) {
      # Assumes single-base REF/ALT uppercase strings
      if (is.na(ref) || is.na(obs)) return(0.01)
      ref <- toupper(ref); obs <- toupper(obs)
      if (nchar(ref) != 1 || nchar(obs) != 1) return(0.01)  # not a simple SNV
      
      # Define complementary mapping if you want strand-agnostic view (not needed here)
      # classify
      if (ref %in% c("A","T")) {
        # A/T reference = ENU-targeted base
        if ((ref == "A" && obs == "G") || (ref == "T" && obs == "C")) {
          return(0.097)   # A->G / T->C transition 
        } else if ((ref == "A" && obs == "T") || (ref == "T" && obs == "A")) {
          return(0.092)   # common transversion class seen with ENU
        } else {
          return(0.046)   # other A/T changes (e.g. A->C / T->G)
        }
      } else if ((ref == "G" && obs == "A") || (ref == "C" && obs == "T")) {
        # G/C changes are less likely from ENU, but weighted accordingly for O6-BG inhibition of MGMT
        return(0.19)
      } else if ((ref == "C" && obs == "A") || (ref == "G" && obs == "T")) {
        return(0.054)
      } else {
        return(0.019) # C->G and G->C changes
      }
    }
  
  # annotate column giving ENU priority
    annotated_variants <- annotated_variants %>% mutate(enu_prior = mapply(assign_enu_prior, Ref, Obs))
    
  # now normalize ENU priority by gene_length
    annotated_variants <- annotated_variants %>% mutate(enu_prior_adjusted = enu_prior / gene_length)
    
  # remove rows that are completely identical
    annotated_variants <- distinct(annotated_variants)
  # have duplicate variant_ids, the only difference is multiple gene_lengths
    annotated_variants_collapsed <- annotated_variants %>% group_by(variant_id) %>% slice_max(order_by = gene_length, with_ties = FALSE) %>% ungroup()
    
    
# STEP 3: Annotate with RNA sequencing data  
      
  # import RNA-seq variant analysis (from Pfaff lab)
    rna_var <- read_csv(file.path(data_dir, "candidate_list", "SNPs", "SNPs-equalInParents-variableInTreated.csv"))
    # make variant_id
      rna_var <- rna_var %>% mutate(variant_id = paste0(X.CHROM, ":", POS, "_", REF, "->", ALT))
      
      # label which variant_ids are in WGS
        rna_var <- rna_var %>% mutate(in_WGS = variant_id %in% annotated_variants_collapsed$variant_id)
        # filter for RNA-detected variants not in WGS
        missing_wgs <- rna_var %>% filter(in_WGS == FALSE)
  # mutate column in annotated_variants_collapsed for whether variant is detected in RNA-seq
    annotated_variants_collapsed <- annotated_variants_collapsed %>% mutate(variant_in_RNA = variant_id %in% rna_var$variant_id) 
    #  FALSE   TRUE 
    # Note: missing 414 detected variants in RNAseq, but these do not exist in the WGS
    
    # find missing variants
      missing_rna <- annotated_variants_collapsed %>% filter(var_in_RNA == "FALSE") 
      
  # import RNA-seq expression analysis (from Pfaff lab)
    gex <- read_csv(file.path(data_dir, "candidate_list", "Expression", "diffexp-T2_v_P2-results.csv"))
    # rename first column from 'id' to 'external_gene_name'
      gex <- gex %>% rename(external_gene_name = id)
    
    # combine with grcm3839 to provide gene_ids
      ann_gex <- left_join(gex, grcm3839, by = "external_gene_name")
  # annotate variant table to indicate gene is expressed
    annotated_variants_collapsed <- annotated_variants_collapsed %>% mutate(gene_expressed = ensembl_gene_id %in% ann_gex$ensembl_gene_id)
    
  # write file to folder so don't have to download every time
    write_tsv(as.data.frame(annotated_variants_collapsed), file = file.path(output_dir, "wgs_r_tables", "annotated_variant_master.tsv"), quote = "none")
   
    
    
    
#-----------------------------
# PART 4: Candidate Filtering
#-----------------------------         
    
# STEP 1: Filter out least likely candidates - Venn & Jaccard  
    
  # remove OtherInfo (hom), keep dbSNP is.na(), remove total_shared >= 4, remove Func (synonymous_variant)
    # remove all homozygous mutations
      candidates <- annotated_variants_collapsed %>% filter(Otherinfo == "het")
    # remove SNV in dbSNP
      candidates <- candidates %>% filter(is.na(dbSNP))
    # remove SNVs in <= 4 mutant suppressor lines
      candidates <- candidates %>% filter(total_shared < 3)
    # remove if silent mutations
      candidates <- candidates %>% filter(Func != "synonymous_variant")
      
      
# STEP 2: plots for overlap and similarity
    
  # venn_diagram for overlap and independence
    y622_805 <- candidates %>% filter(in_622 == 1 & in_805 == 1)
    y752 <- candidates %>% filter(in_752 == 1)
    y799 <- candidates %>% filter(in_799 == 1)
    y814 <- candidates %>% filter(in_814 == 1)
    var_gene <- list("622/805" = y622_805$ensembl_gene_id,
                      "752" = y752$ensembl_gene_id,
                      "799" = y799$ensembl_gene_id,
                      "814" = y814$ensembl_gene_id)
    
    gene_venn <- ggvenn(var_gene,
                        show_percentage = FALSE,
                        fill_color = c("#fa030375", "#e1005375", "#be007675", "#040cd475"),
                        stroke_size = 0,
                        set_name_color =c("#fa0303", "#e10053", "#be0076", "#040cd4"),
                        set_name_size = 4,
                        text_color = "white",
                        text_size = 3)
    
    ggsave(filename = "mutsupp_gene_venn.png", plot = gene_venn, device = "png", dpi = 300, height = 4, width = 4, path = file.path(output_dir, "wgs_r_figures"))
    
    # perform Fischer's exact test
      options(scipen = 999, digits = 22)
    
      pairs <- combn(names(var_gene),2, simplify=FALSE)
    
      pvals <- sapply(pairs, function(p){
        A <- var_gene[[p[1]]]
        B <- var_gene[[p[2]]]
      
        inter <- length(intersect(A,B))
        onlyA <- length(setdiff(A,B))
        onlyB <- length(setdiff(B,A))
      
      m <- matrix(c(inter, onlyA, onlyB, 0), nrow=2, byrow=TRUE)
      fisher.test(m)$p.value
    })
    
    padj <- p.adjust(pvals, method="BH")
    
    # Jaccard similarity
    jaccard <- sapply(pairs, function(p){
      A <- var_gene[[p[1]]]
      B <- var_gene[[p[2]]]
      length(intersect(A,B)) / length(union(A,B))
    })
    
    results_gene <- data.frame(
      pair = sapply(pairs, paste, collapse=" vs "),
      fisher_p = pvals,
      fisher_p_adj = padj,
      jaccard = jaccard
    )
    
  # compute pairwise Jaccard matrix
    samples_gene <- names(var_gene)
    jacmat_gene <- matrix(0, nrow=length(samples_gene), ncol=length(samples_gene),
                          dimnames=list(samples_gene, samples_gene))
    
    for(i in samples_gene){
      for(j in samples_gene){
        A <- var_gene[[i]]
        B <- var_gene[[j]]
        jacmat_gene[i,j] <- length(intersect(A,B)) / length(unique(c(A,B)))
      }
    }
    
    # heatmap
      my_breaks <- seq(0, 1, length.out = 101)
      outpath <- path.expand(file.path(output_dir, "wgs_r_figures", "mutsupp_gene_jacc-heat.png"))
      png(filename = outpath, width=2000, height=2000, res=300)
      pheatmap(jacmat_gene,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               display_numbers = TRUE,
               number_format = "%.2f",
               fontsize = 12,
               cellwidth = 40,
               cellheight = 40,
               border_color = NA,
               breaks = my_breaks,
               treeheight_row = 20,
               treeheight_col =  20,
               fontsize_number = 10,
               gaps_row = NULL,
               gaps_col = NULL,
               angle_col = 0,
               number_color = "black",
               filename = outpath
      )
    

# STEP 3: Split apart 622/805 from the rest, then provide descriptives
      
  # split into two groups 622/805 and 752/799/814
    X622_805 <- candidates %>% filter(in_622 == 1 & in_805 == 1)
    X752_799_814 <- candidates %>% filter(in_622 != 1 & in_805 != 1)
  # for full variant table list
    func <- as.data.frame(table(annotated_variants_collapsed$Func))
    biotype <- as.data.frame(table(annotated_variants_collapsed$gene_biotype))
    
  # for filtered/split groups
    func_622_805 <- as.data.frame(table(X622_805$Func))
    func_752_799_814 <- as.data.frame(table(X752_799_814$Func))
    biotype_622_805 <- as.data.frame(table(X622_805$gene_biotype))
    biotype_752_799_814 <- as.data.frame(table(X752_799_814$gene_biotype))
    
  # write file to folder so don't have to download every time
    write_tsv(as.data.frame(func), file = file.path(output_dir, "wgs_r_tables", "descriptives", "snv_func_breakdown.tsv"), quote = "none")
    write_tsv(as.data.frame(biotype), file = file.path(output_dir, "wgs_r_tables", "descriptives", "snv_biotype_breakdown.tsv"), quote = "none")
    write_tsv(as.data.frame(func_622_805), file = file.path(output_dir, "wgs_r_tables", "descriptives", "snv_func_622_805_breakdown.tsv"), quote = "none")
    write_tsv(as.data.frame(biotype_622_805), file = file.path(output_dir, "wgs_r_tables", "descriptives", "snv_biotype_622_805_breakdown.tsv"), quote = "none")
    write_tsv(as.data.frame(func_752_799_814), file = file.path(output_dir, "wgs_r_tables", "descriptives", "snv_func_752_799_814_breakdown.tsv"), quote = "none")
    write_tsv(as.data.frame(biotype_752_799_814), file = file.path(output_dir, "wgs_r_tables", "descriptives", "snv_biotype_752_799_814_breakdown.tsv"), quote = "none")
    
  
# STEP 4: Keep missense, filter by impc_viability, likelihood_essential, and expression
    
  # focus only on missense mutations
    X622_805_candidate <- X622_805 %>% filter(Func == "missense_variant")
    X752_799_814_candidate <- X752_799_814 %>% filter(Func == "missense_variant")
  # remove impc_viability and positive mean effect
    X622_805_candidate <- X622_805_candidate %>% filter(!(impc_viability %in% c("viable", "subviable")) & mean_effect_viability < 0)
    X752_799_814_candidate <- X752_799_814_candidate %>% filter(!(impc_viability %in% c("viable", "subviable")) & mean_effect_viability < 0)
  # determine a threshold for essential cutoff using impc_viability = lethal
    lethal <- candidates %>% filter(impc_viability == "lethal")
    median(lethal$likelihood_essential, na.rm = TRUE) # mean = 0.2, median = 0.013
      
    # keep likelihood_essential > 0.1 (dependency in 10% of cells tested) or impc = lethal
      X622_805_candidate <- X622_805_candidate %>% filter(impc_viability == "lethal" | likelihood_essential > 0.1)
      X752_799_814_candidate <- X752_799_814_candidate %>% filter(impc_viability == "lethal" | likelihood_essential > 0.1)
  # ensure retained genes are expressed
    X622_805_candidate <- X622_805_candidate %>% filter(gene_expressed == "TRUE")
    X752_799_814_candidate <- X752_799_814_candidate %>% filter(gene_expressed == "TRUE")
  # write file to directory
    write_tsv(as.data.frame(X622_805_candidate), file = file.path(output_dir, "wgs_r_tables", "candidates", "X622_805_candidates.tsv"), quote = "none")
    write_tsv(as.data.frame(X752_799_814_candidate), file = file.path(output_dir, "wgs_r_tables", "candidates", "X752_799_814_candidates.tsv"), quote = "none")
    
# STEP 5: test for non-random overlap      
      
  # venn_diagram for overlap and independence
    z752 <- X752_799_814_candidate %>% filter(in_752 == 1)
    z799 <- X752_799_814_candidate %>% filter(in_799 == 1)
    z814 <- X752_799_814_candidate %>% filter(in_814 == 1)
      var_gene_filt <- list("622/805" = X622_805_candidate$ensembl_gene_id,
                            "752" = z752$ensembl_gene_id,
                            "799" = z799$ensembl_gene_id,
                            "814" = z814$ensembl_gene_id)
      
      gene_filt_venn <- ggvenn(var_gene_filt,
                               show_percentage = FALSE,
                               fill_color = c("#fa030375", "#e1005375", "#be007675", "#040cd475"),
                               stroke_size = 0,
                               set_name_color =c("#fa0303", "#e10053", "#be0076", "#040cd4"),
                               set_name_size = 4,
                               text_color = "white",
                               text_size = 3)
      
      ggsave(filename = "mutsupp_gene_filt_venn.png", plot = gene_filt_venn, device = "png", dpi = 300, height = 4, width = 4, path = file.path(output_dir, "wgs_r_figures"))
      
      # perform Fischer's exact test
        options(scipen = 999, digits = 22)
      
        pairs <- combn(names(var_gene_filt),2, simplify=FALSE)
      
        pvals_filt <- sapply(pairs, function(p){
          A <- var_gene_filt[[p[1]]]
          B <- var_gene_filt[[p[2]]]
        
          inter <- length(intersect(A,B))
          onlyA <- length(setdiff(A,B))
          onlyB <- length(setdiff(B,A))
        
          m <- matrix(c(inter, onlyA, onlyB, 0), nrow=2, byrow=TRUE)
          fisher.test(m)$p.value
        })
      
        padj_filt <- p.adjust(pvals, method="BH")
        
      
# STEP 6:  GO-enrichment  
        
  # build reference gene/GO set
    grcm3839 <- read_delim(file.path(output_dir, "biomaRt", "grcm3839_collapsed_gene_annotation.tsv"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    grcm3839 <- grcm3839[,c(1,7,10)]
    gene2go_universe <- grcm3839 %>% separate_rows(go_id, sep = "[;,]") %>% filter(go_id != "" & !is.na(go_id)) %>% distinct()
    gene2go_universe <- gene2go_universe %>% dplyr::rename(term = go_id, gene = ensembl_gene_id, gene_symbol = external_gene_name)
    TERM2GENE <- gene2go_universe %>% dplyr::select(term, gene)
    
    unique_terms <- unique(gene2go_universe$term)
    
    go_names_universe <- AnnotationDbi::select(GO.db, keys = unique_terms, columns = c("GOID", "TERM", "ONTOLOGY"), keytype = "GOID") %>% rename(term = GOID, name = TERM, ontology = ONTOLOGY)
    
  # build long-form GO tables from SNV sets
    X622_805_snvc <- X622_805_candidate %>% dplyr::select("ensembl_gene_id", "go_id")
    X622_805_snvc <- X622_805_snvc %>% separate_rows(go_id, sep="[;,]") %>% filter(go_id != "" & !is.na(go_id)) %>% distinct(ensembl_gene_id, go_id)
    X752_snvc <- X752_799_814_candidate %>% filter(in_752 == 1)
    X752_snvc <- X752_snvc %>% dplyr::select("ensembl_gene_id", "go_id")
    X752_snvc <- X752_snvc %>% separate_rows(go_id, sep="[;,]") %>% filter(go_id != "" & !is.na(go_id)) %>% distinct(ensembl_gene_id, go_id)
    X799_snvc <- X752_799_814_candidate %>% filter(in_799 == 1)
    X799_snvc <- X799_snvc %>% dplyr::select("ensembl_gene_id", "go_id")
    X799_snvc <- X799_snvc %>% separate_rows(go_id, sep="[;,]") %>% filter(go_id != "" & !is.na(go_id)) %>% distinct(ensembl_gene_id, go_id)
    X814_snvc <- X752_799_814_candidate %>% filter(in_814 == 1)
    X814_snvc <- X814_snvc %>% dplyr::select("ensembl_gene_id", "go_id")
    X814_snvc <- X814_snvc %>% separate_rows(go_id, sep="[;,]") %>% filter(go_id != "" & !is.na(go_id)) %>% distinct(ensembl_gene_id, go_id)
  # extract unique gene lists for enrichment
    genes_snvc_622_805 <- unique(X622_805_snvc$ensembl_gene_id)
    genes_snvc_752 <- unique(X752_snvc$ensembl_gene_id)
    genes_snvc_799 <- unique(X799_snvc$ensembl_gene_id)
    genes_snvc_814 <- unique(X814_snvc$ensembl_gene_id)
    
  # perform GO-enrichment against reference
      
    ego_622_805 <- enricher(
      gene      = genes_snvc_622_805,
      TERM2GENE = TERM2GENE,
      TERM2NAME = go_names_universe,
      universe      = unique(gene2go_universe_for_enrich$gene)
    )
    
    ego_752 <- enricher(
      gene      = genes_snvc_752,
      TERM2GENE = TERM2GENE,
      TERM2NAME = go_names_universe,
      universe      = unique(gene2go_universe_for_enrich$gene)
    )
    
    ego_799 <- enricher(
      gene      = genes_snvc_799,
      TERM2GENE = TERM2GENE,
      TERM2NAME = go_names_universe,
      universe      = unique(gene2go_universe_for_enrich$gene)
    )
    
    ego_814 <- enricher(
      gene      = genes_snvc_814,
      TERM2GENE = TERM2GENE,
      TERM2NAME = go_names_universe,
      universe      = unique(gene2go_universe_for_enrich$gene)
    )
    
  # filter and dotplot
    
    compareCluster_result <- compareCluster(
      geneCluster = list(
        `622/805` = genes_snvc_622_805,
        `752`     = genes_snvc_752,
        `799`     = genes_snvc_799,
        `814`     = genes_snvc_814
      ),
      fun = "enricher",
      TERM2GENE = TERM2GENE,
      TERM2NAME = go_names_universe
    )
    
    df <- as.data.frame(compareCluster_result) %>% left_join(go_names_universe %>% select(term, ontology) %>% rename(ID = term), by = "ID")
    
  # Keep top 10 GO terms per cluster (by p.adjust)
    df_top10 <- df %>% group_by(Cluster) %>% slice_min(order_by = p.adjust, n = 10) %>% ungroup() %>% mutate(log10p = -log10(p.adjust))
    
  # Keep top 10 GO BP terms per cluster (by p.adjust)
    top_BP <- df %>% filter(ontology == "BP") %>% group_by(Cluster) %>% arrange(p.adjust, .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup() %>% mutate(log10p = -log10(p.adjust))
    
  # Keep top 10 GO CC terms per cluster (by p.adjust)
    top_CC <- df %>% filter(ontology == "CC") %>% group_by(Cluster) %>% arrange(p.adjust, .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup() %>% mutate(log10p = -log10(p.adjust))
    
  # Keep top 10 GO MF terms per cluster (by p.adjust)
    top_MF <- df %>% filter(ontology == "MF") %>% group_by(Cluster) %>% arrange(p.adjust, .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup() %>% mutate(log10p = -log10(p.adjust))
    
  # to simplify    
    o <- df_top10
    o <- top_BP
    o <- top_CC
    o <- top_MF
    
  # Dotplot
    p <- ggplot(o,
                aes(x = Cluster,
                    y = Description,
                    size = Count,
                    color = log10p)) +
      geom_point() +
      scale_color_gradient(low = "blue", high = "red") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 7),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
#            plot.title = element_text(hjust = 0.675)) +     # for top10
            plot.title = element_text(hjust = 0.675)) +
      labs(x = "Mutant Suppressor Line",
           color = "-log10(p.adjust)",
           size  = "Gene Count",
#           title = "Top 10 Enriched GO Terms per Mutant Suppressor Line"
#           title = "GO: Biological Process" 
#           title = "GO: Cellular Component"
           title = "GO: Molecular Function"
      ) 
  
    output_dir_go <- file.path(output_dir, "wgs_r_figures", "go_enrichment")
  
    top10 <- p
    BP <- p
    CC <- p
    MF <- p
      
    ggsave(filename = file.path(output_dir_go, "GO_Top10_all_Dotplot.png"), plot = top10, width = 6.5, height = 8, dpi = 300)
    ggsave(filename = file.path(output_dir_go, "GO_Top10_BP_Dotplot.png"), plot = BP, width = 6.5, height = 8, dpi = 300)
    ggsave(filename = file.path(output_dir_go, "GO_Top10_CC_Dotplot.png"), plot = CC, width = 6, height = 8, dpi = 300)
    ggsave(filename = file.path(output_dir_go, "GO_Top10_MF_Dotplot.png"), plot = MF, width = 5, height = 8, dpi = 300)
    
    
# STEP 7:  network analysis (exploratory; overlap between terms may be limited)
    
    # Helper function to get significant GO IDs per line
      get_sig_go_df <- function(ego_result, padj_cutoff = 0.0001, sample_name) {
        df <- as.data.frame(ego_result@result)
        df <- df %>% filter(p.adjust < padj_cutoff) %>%
          mutate(Sample = sample_name) %>%
          select(GO_ID = ID, p.adjust, Sample)
        return(df)
      }
    
      sig_622_805 <- get_sig_go_df(ego_622_805, padj_cutoff = 0.001, "622/805")
      sig_752     <- get_sig_go_df(ego_752, padj_cutoff = 0.001, "752")
      sig_799     <- get_sig_go_df(ego_799, padj_cutoff = 0.001, "799")
      sig_814     <- get_sig_go_df(ego_814, padj_cutoff = 0.001, "814")
    
    # Combine all significant GO IDs
      sig_go_all <- bind_rows(sig_622_805, sig_752, sig_799, sig_814)
    
    # Merge with gene mapping
      network_df_filtered <- network_df %>% inner_join(sig_go_all, by = c("GO_ID", "Sample"))  # Row count varies by p.adjust threshold: if padj 0.05, 1428 if 0.01, 942 if 0.001
    
    # build network
      library(igraph)
      library(ggraph)
    
    # Edge list
      edges <- network_df_filtered %>% select(Gene, GO_ID)
    
    # Graph
      g <- graph_from_data_frame(edges, directed = FALSE)
      V(g)$type <- ifelse(V(g)$name %in% network_df_filtered$GO_ID, "GO", "Gene")
    
    # Optionally, store padj as node attribute for GO nodes
      GO_padj <- sig_go_all %>% distinct(GO_ID, p.adjust)
      V(g)$padj <- ifelse(V(g)$type == "GO", GO_padj$p.adjust[match(V(g)$name, GO_padj$GO_ID)], NA)
    
    # plot filtered network
      ggraph(g, layout = "fr") +
        geom_edge_link(alpha = 0.4) +
        geom_node_point(aes(color = type, size = ifelse(type=="GO", -log10(padj), 3))) +
        geom_node_text(aes(label = name), repel = TRUE, size = 3) +
        theme_void() +
        scale_color_manual(values = c(Gene = "lightblue", GO = "salmon")) +
        ggtitle("Filtered Gene–GO Network (padj < 0.001)")
    

# STEP 8: Similarity in GO-terms and genes

  # Helper function to extract significant GO IDs
    get_sig_go <- function(ego_result, padj_cutoff = 0.05) {
      df <- as.data.frame(ego_result@result)
      sig_df <- df %>% filter(p.adjust < padj_cutoff)
      return(sig_df$ID)  # pull only GO IDs
    }
      
    go_622_805 <- get_sig_go(ego_622_805)
    go_752     <- get_sig_go(ego_752)
    go_799     <- get_sig_go(ego_799)
    go_814     <- get_sig_go(ego_814)
      
  # define go_universe
    go_universe <- gene2go_universe$term
      
  # pairwise jaccard indices
    go_sets <- list(`622/805` = go_622_805,
                    `752`     = go_752,
                    `799`     = go_799,
                    `814`     = go_814)
      
  # Jaccard function
    jaccard_index <- function(A, B) {
      length(intersect(A, B)) / length(union(A, B))
    }
      
  # Create pairwise matrix
    jaccard_mat <- outer(names(go_sets), names(go_sets), Vectorize(function(x,y) jaccard_index(go_sets[[x]], go_sets[[y]])))
    colnames(jaccard_mat) <- rownames(jaccard_mat) <- names(go_sets)
    jaccard_mat
      
    go_sets <- list(`622/805` = go_622_805,
                    `752`     = go_752,
                    `799`     = go_799,
                    `814`     = go_814)
      
  # Jaccard similarity heatmap
    outpath <- path.expand(file.path(output_dir, "wgs_r_figures", "go_enrichment", "mutsupp_go_jacc-heat.png"))
    png(filename = outpath, width=2000, height=2000, res=300)
    pheatmap(jaccard_mat, 
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             display_numbers = TRUE,
             number_format = "%.2f",
             fontsize = 12,
             cellwidth = 40,
             cellheight = 40,
             border_color = NA,
             treeheight_row = 20,
             treeheight_col =  20,
             fontsize_number = 10,
             gaps_row = NULL,
             gaps_col = NULL,
             angle_col = 0,
             number_color = "black",
             filename = outpath
             )
      
      
# STEP 9: Perform filtering for SMN-associated GO IDs 
      
  # grab all the go ids/terms associated with Smn
    smn_go <- gene2go_universe %>% filter(gene_symbol == "Smn1")
  # combine candidates for all lines
    candidates_combined <- bind_rows(X622_805_candidate, X752_799_814_candidate)
  # mutate column to flag 'yes' if go_id in smn_go$term is in the collapsed candidates_combined$go_id
    # provide number of SMN-associated GO terms to sort by
      smn_go_ids <- smn_go$term
      smn_associated <- candidates_combined %>% mutate(smn_association_freq = map_int(str_split(go_id, ";"), ~ sum(.x %in% smn_go_ids)), smn_associated   = if_else(smn_association_freq > 0, "yes", "no"))
      smn_associated <- smn_associated %>% filter(smn_associated == "yes")
  # venn diagram of mutated genes associated with SMN GO-terms
      v622_805 <- smn_associated %>% filter(in_622 == 1 & in_805 == 1)
      v752 <- smn_associated %>% filter(in_752 == 1)
      v799 <- smn_associated %>% filter(in_799 == 1)
      v814 <- smn_associated %>% filter(in_814 == 1)
      var_gene_smn <- list("622/805" = v622_805$ensembl_gene_id,
                            "752" = v752$ensembl_gene_id,
                            "799" = v799$ensembl_gene_id,
                            "814" = v814$ensembl_gene_id)
      
      gene_smn_venn <- ggvenn(var_gene_smn,
                               show_percentage = FALSE,
                               fill_color = c("#fa030375", "#e1005375", "#be007675", "#040cd475"),
                               stroke_size = 0,
                               set_name_color =c("#fa0303", "#e10053", "#be0076", "#040cd4"),
                               set_name_size = 4,
                               text_color = "white",
                               text_size = 3)
      
      ggsave(filename = "mutsupp_gene_smn_venn.png", plot = gene_smn_venn, device = "png", dpi = 300, height = 4, width = 4, path = file.path(output_dir, "wgs_r_figures"))

  # annotate smn_associated to include whether the gene is mutated in multiple cell lines
    smn_associated <- smn_associated %>% mutate(gene_in_622 = if_else(ensembl_gene_id %in% v622_805$ensembl_gene_id, 1, 0))
    smn_associated <- smn_associated %>% mutate(gene_in_752 = if_else(ensembl_gene_id %in% v752$ensembl_gene_id, 1, 0))
    smn_associated <- smn_associated %>% mutate(gene_in_799 = if_else(ensembl_gene_id %in% v799$ensembl_gene_id, 1, 0))
    smn_associated <- smn_associated %>% mutate(gene_in_814 = if_else(ensembl_gene_id %in% v814$ensembl_gene_id, 1, 0))
    smn_associated <- smn_associated %>% mutate(shared_mut_gene = gene_in_622 + gene_in_752 + gene_in_799 + gene_in_814)
  # this is the final combined table to publish
    write_tsv(as.data.frame(smn_associated), file = file.path(output_dir, "wgs_r_tables", "candidates", "smn_associated_candidates.tsv"), quote = "none")
    

#-----------------------------
# PART 5: Top 10 Variant Candidates per Mutant Suppressor Line
#-----------------------------   

# STEP 1: Make separate tables for each mutant suppressor line
    
  sa622_805 <- smn_associated %>% filter(in_622 == 1 & in_805 == 1)
  sa752 <- smn_associated %>% filter(in_752 == 1)
  sa799 <- smn_associated %>% filter(in_799 == 1)
  sa814 <- smn_associated %>% filter(in_814 == 1)
  
  
# STEP 2: Sort by smn_association_freq and likelihood_essential
  
  # First: sort by Smn-association
    sa622_805_top20 <- sa622_805 %>% arrange(desc(smn_association_freq)) %>% slice_head(n = 20)
    sa752_top20 <- sa752 %>% arrange(desc(smn_association_freq)) %>% slice_head(n = 20)
    sa799_top20 <- sa799 %>% arrange(desc(smn_association_freq)) %>% slice_head(n = 20)
    sa814_top20 <- sa814 %>% arrange(desc(smn_association_freq)) %>% slice_head(n = 20)
    
  # Second: sort by likelihood_essential
    sa622_805_top10 <- sa622_805_top20 %>% arrange(desc(likelihood_essential)) %>% slice_head(n = 10)
    sa752_top10 <- sa752_top20 %>% arrange(desc(likelihood_essential)) %>% slice_head(n = 10)
    sa799_top10 <- sa799_top20 %>% arrange(desc(likelihood_essential)) %>% slice_head(n = 10)
    sa814_top10 <- sa814_top20 %>% arrange(desc(likelihood_essential)) %>% slice_head(n = 10)
    
  # write these top10 tables to directory
    write_tsv(as.data.frame(sa622_805_top10), file = file.path(output_dir, "wgs_r_tables", "candidates", "top10_622_805_candidates.tsv"), quote = "none")
    write_tsv(as.data.frame(sa752_top10), file = file.path(output_dir, "wgs_r_tables", "candidates", "top10_752_candidates.tsv"), quote = "none")
    write_tsv(as.data.frame(sa799_top10), file = file.path(output_dir, "wgs_r_tables", "candidates", "top10_799_candidates.tsv"), quote = "none")
    write_tsv(as.data.frame(sa814_top10), file = file.path(output_dir, "wgs_r_tables", "candidates", "top10_814_candidates.tsv"), quote = "none")
    
    
    
#-----------------------------
# PART 6: GEX and variants
#-----------------------------      
    
# STEP 1: Import and process differential expression analysis
    
  # import RNA-seq expression analysis (from Pfaff lab)
    gex <- read_csv(file.path(data_dir, "candidate_list", "Expression", "diffexp-T2_v_P2-results.csv"))
  # rename first column from 'id' to 'external_gene_name'
    gex <- gex %>% rename(external_gene_name = id)
    
  # combine with grcm3839 to provide gene_ids
    ann_gex <- left_join(gex, grcm3839, by = "external_gene_name")
    
    
# STEP 2: Look at variant type (Func) in all, upregulated, and downregulated genes
  
  # create lists of gene_ids for all expressed, upregulated, and downregulated
    all_gex_ids <- ann_gex$ensembl_gene_id
  
    up_gex <- ann_gex %>% filter(logFC >= 0.5 & padj < 0.05)
    up_gex_ids <- up_gex$ensembl_gene_id
    
    down_gex <- ann_gex %>% filter(logFC <= -0.5 & padj < 0.05)
    down_gex_ids <- down_gex$ensembl_gene_id
    
  # use above to separate variants associated with up/down regulated genes
    all_gex_var <- annotated_variants_collapsed %>% filter(ensembl_gene_id %in% all_gex_ids)
    up_gex_var   <- annotated_variants_collapsed %>% filter(ensembl_gene_id %in% up_gex_ids)
    down_gex_var <- annotated_variants_collapsed %>% filter(ensembl_gene_id %in% down_gex_ids)

  # now variant function
    all_func <- as.data.frame(table(all_gex_var$Func))
      names(all_func) <- c("variant_function", "all")
      
    up_func <- as.data.frame(table(up_gex_var$Func))
      names(up_func) <- c("variant_function", "up")
      
    down_func <- as.data.frame(table(down_gex_var$Func))
      names(down_func) <- c("variant_function", "down")
      
  # combine into one table
    var_func_gex <- full_join(all_func, up_func, by = "variant_function")
    var_func_gex <- full_join(var_func_gex, down_func, by = "variant_function")
    
  # write to directory
    write_tsv(as.data.frame(var_func_gex), file = file.path(output_dir, "wgs_r_tables", "descriptives", "variant_gex_function.tsv"), quote = "none")
 
    
# STEP 3: Volcano plot for differential expression    
  
  # annotate ann_gex for significantly differentialy expressed genes
    ann_gex <- ann_gex %>% mutate(expr_change = case_when(external_gene_name %in% up_gex$external_gene_name ~ "Up", 
                                                          external_gene_name %in% down_gex$external_gene_name ~ "Down",
                                                          TRUE ~ "No Change"))
    
  # Volcano Plot
    vol <- ggplot(data = ann_gex, aes(x = logFC, y = -log10(padj), col = expr_change)) + 
        geom_point(alpha=0.50) + 
        labs(x = "log2(f∆) Suppressors vs Parent", y = "-log10(p-adj)") +
        theme_minimal() +
        scale_color_manual("Sig. DEGs", values = c("#040cd4", "grey", "#fa0303"), breaks = c("Down", "No Change", "Up" ), labels = c("Down", "No Change", "Up")) + 
        coord_cartesian(xlim = c(-10,10)) + 
        theme(panel.background = element_rect(fill = NA, color = "black", linewidth = rel(3)), 
              panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
              panel.grid.minor.x = element_blank(),
              axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
              axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
              #axis.title.x = element_blank(),
              #axis.text.x = element_blank(),
              legend.position = c(0.76, 0.80), 
              legend.box.background = element_rect(fill = "#FFFFFF80", color = "white", linewidth = 1), # element_blank(), 
              legend.key = element_blank(),
        legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
        legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"))+
        geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype = "dashed") + geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed")
    
    # write plot to directory
      ggsave(filename = "GEX_volcano.png", p = vol, device = "png", dpi = 300, height = 4, width = 4, path = file.path(output_dir, "wgs_r_figures", "gex"))
    
      
# STEP 4: Gene Ontology for Differentially Expressed Genes
      
  # filter and dotplot
      
      gex_go <- compareCluster(
        geneCluster = list(
          #`Up`    = up_gex_ids,
          `Down`  = down_gex_ids
        ),
        fun = "enricher",
        TERM2GENE = TERM2GENE,
        TERM2NAME = go_names_universe
      )
      
      df <- as.data.frame(gex_go) %>% left_join(go_names_universe %>% select(term, ontology) %>% rename(ID = term), by = "ID")
      
      # Keep top 10 GO terms per cluster (by p.adjust)
        top10 <- df %>% group_by(Cluster) %>% slice_min(order_by = p.adjust, n = 10) %>% ungroup() %>% mutate(log10p = -log10(p.adjust))
        
      # Keep top 10 GO BP terms per cluster (by p.adjust)
        top_BP <- df %>% filter(ontology == "BP") %>% group_by(Cluster) %>% arrange(p.adjust, .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup() %>% mutate(log10p = -log10(p.adjust))
        
      # Keep top 10 GO CC terms per cluster (by p.adjust)
        top_CC <- df %>% filter(ontology == "CC") %>% group_by(Cluster) %>% arrange(p.adjust, .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup() %>% mutate(log10p = -log10(p.adjust))
        
      # Keep top 10 GO MF terms per cluster (by p.adjust)
        top_MF <- df %>% filter(ontology == "MF") %>% group_by(Cluster) %>% arrange(p.adjust, .by_group = TRUE) %>% slice_head(n = 10) %>% ungroup() %>% mutate(log10p = -log10(p.adjust))
      
      df_top_combined <- bind_rows(top10, top_BP, top_CC, top_MF)
        
  # Dotplot
      p <- ggplot(df_top_combined, aes(x = Cluster, y = Description, size = Count, color = log10p)) +
        geom_point() +
        scale_color_gradient(low = "blue", high = "red") +
        facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 6),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text.y = element_text(size = 9, face = "bold"),
              plot.title = element_text(hjust = 0.5)) +
        labs(x = "Down-regulated Genes",
             color = "-log10(p.adjust)",
             size  = "Gene Count",
             title = "Top 10 GO Terms per Ontology (Down DEGs)")
      
        
        output_dir_gex <- file.path(output_dir, "wgs_r_figures", "gex")
        
        ggsave(filename = file.path(output_dir_gex, "Top10_GO_DEGs_Dotplot.png"), plot = p, width = 5, height = 8, dpi = 300)
        
     
# STEP 5: Label genes/transcripts that are associated with variants
    
  # make tables for specific Func categories for identified variants
    start_stop <- annotated_variants_collapsed %>% filter(Func == "start_lost" | Func == "stop_gained" | Func == "stop_lost" | Func == "5_prime_UTR_premature_start_codon_gain_variant")
    utr5 <- annotated_variants_collapsed %>% filter(Func == "5_prime_UTR_variant")
    utr3 <- annotated_variants_collapsed %>% filter(Func == "3_prime_UTR_variant")
    upstream <- annotated_variants_collapsed %>% filter(Func == "upstream_gene_variant")
    downstream <- annotated_variants_collapsed %>% filter(Func == "downstream_gene_variant")
    intergenic <- annotated_variants_collapsed %>% filter(Func == "intergenic_region")
    splicing <- annotated_variants_collapsed %>% filter(Func == "splice_acceptor_variant" | Func == "splice_donor_variant" | Func == "splice_region_variant")
  # label if gene/transcript is associated with a variant
    ann_gex <- ann_gex %>% mutate(variant = case_when(ensembl_gene_id %in% annotated_variants_collapsed$ensembl_gene_id ~ "variant", TRUE ~ "absent"))
    # absent variant 

    ann_gex <- ann_gex %>% mutate(start_stop_variant = case_when(ensembl_gene_id %in% start_stop$ensembl_gene_id ~ "start_stop", TRUE ~ "absent"))
    # absent start_stop 

    ann_gex <- ann_gex %>% mutate(utr_variant = case_when(ensembl_gene_id %in% utr5$ensembl_gene_id ~ "5UTR", ensembl_gene_id %in% utr3$ensembl_gene_id ~ "3UTR", TRUE ~ "absent"))
    # 3UTR   5UTR absent 

    ann_gex <- ann_gex %>% mutate(genic_variant = case_when(ensembl_gene_id %in% upstream$ensembl_gene_id ~ "upstream", ensembl_gene_id %in% downstream$ensembl_gene_id ~ "downstream", ensembl_gene_id %in% intergenic$ensembl_gene_id ~ "intergenic", TRUE ~ "absent"))
    # absent downstream intergenic   upstream 

    ann_gex <- ann_gex %>% mutate(splicing = case_when(ensembl_gene_id %in% splicing$ensembl_gene_id ~ "splicing", TRUE ~ "absent"))
    # absent splicing 

  # mutate a column to flag different variant classifications
    start_stop_ids  <- start_stop$ensembl_gene_id
    utr5_ids        <- utr5$ensembl_gene_id
    utr3_ids        <- utr3$ensembl_gene_id
    splicing_ids    <- splicing$ensembl_gene_id
    upstream_ids    <- upstream$ensembl_gene_id
    downstream_ids  <- downstream$ensembl_gene_id
    intergenic_ids  <- intergenic$ensembl_gene_id
    all_ids         <- annotated_variants_collapsed$ensembl_gene_id
    
    ann_gex <- ann_gex %>% mutate(Event_type = case_when(
      ensembl_gene_id %in% start_stop_ids  ~ "start_stop",
      ensembl_gene_id %in% utr5_ids        ~ "5UTR",
      ensembl_gene_id %in% utr3_ids        ~ "3UTR",
      ensembl_gene_id %in% splicing_ids    ~ "splicing",
      ensembl_gene_id %in% upstream_ids    ~ "upstream",
      ensembl_gene_id %in% downstream_ids  ~ "downstream",
      ensembl_gene_id %in% intergenic_ids  ~ "intergenic",
      ensembl_gene_id %in% all_ids         ~ "other",
      TRUE                                  ~ "none"))
    

# STEP 6: Visualizing variant effects on differential expression
  
  # Volcano Plot
#    all_vol <- ggplot(data = ann_gex, aes(x = logFC, y = -log10(padj), col = Event_type)) + 
#    stsp_vol <- ggplot(data = ann_gex, aes(x = logFC, y = -log10(padj), col = start_stop_variant)) + 
#    utr_vol <- ggplot(data = ann_gex, aes(x = logFC, y = -log10(padj), col = utr_variant)) + 
#    genic_vol <- ggplot(data = ann_gex, aes(x = logFC, y = -log10(padj), col = genic_variant)) + 
    splicing_vol <- ggplot(data = ann_gex, aes(x = logFC, y = -log10(padj), col = splicing)) + 
      geom_point(alpha=0.50) + 
      #labs(x = "log2(f∆) Mutant Suppressors vs Parent", y = "-log10(p-adj)") +
      theme_minimal() +
#      scale_color_manual("Variant:", values = c("#fa0303", "#f10038", "#e10053", "#d00065", "#be0076", "#ab0087", "#94009c","#7000bc", "#040cd4" ), breaks = c("start_stop", "5UTR", "3UTR", "splicing", "upstream", "downstream", "intergenic", "other", "none"), labels = c("start_stop", "5UTR", "3UTR", "splicing", "upstream", "downstream", "intergenic", "other", "none"))+
#      scale_color_manual("Start/Stop Variant", values = c("#040cd4", "#fa0303"), breaks = c("absent", "start_stop" ), labels = c("Absent", "Start/Stop")) + 
#      scale_color_manual("UTR Variant", values = c("#040cd4", "#be0076", "#fa0303"), breaks = c("absent", "3UTR", "5UTR" ), labels = c("Absent", "3'UTR", "5'UTR")) + 
#      scale_color_manual("Genic/Intergenic Variant", values = c("#040cd4", "#d60060", "#a4008d", "#fa0303"), breaks = c("absent", "intergenic", "downstream", "upstream" ), labels = c("Absent", "Intergenic", "Downstream", "Upstream")) + 
      scale_color_manual("Splicing Variant", values = c("#040cd4", "#fa0303"), breaks = c("absent", "splicing" ), labels = c("Absent", "Splicing")) + 
      coord_cartesian(xlim = c(-10,10)) + 
      theme(panel.background = element_rect(fill = NA, color = "black", linewidth = rel(3)), 
            panel.grid.major = element_line(colour = "grey90", linetype = "dashed"),
            panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
            panel.grid.minor.x = element_blank(),
            axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            legend.position = "none", #c(0.85, 0.70), 
            #legend.box.background = element_rect(fill = "#FFFFFF", color = "white", linewidth = 1), # element_blank(), 
            legend.key = element_blank()) + 
            #legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
            #legend.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"))+
      geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype = "dashed") + geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed")
    
    
  # Density Plot
#    all_den <- ggplot(data = ann_gex, aes(x = logFC, color = Event_type)) + 
#    stsp_den <- ggplot(data = ann_gex, aes(x = logFC, color = start_stop_variant)) + 
#    utr_den <- ggplot(data = ann_gex, aes(x = logFC, color = utr_variant)) + 
#    genic_den <- ggplot(data = ann_gex, aes(x = logFC, color = genic_variant)) + 
     splicing_den <- ggplot(data = ann_gex, aes(x = logFC, color = splicing)) + 
      geom_density(linewidth = 1) +
      labs(y = "Density") +
      coord_cartesian(xlim = c(-10,10), ylim = c(0,1)) +
#      scale_color_manual("Variant:", values = c("#fa0303", "#f10038", "#e10053", "#d00065", "#be0076", "#ab0087", "#94009c","#7000bc", "#040cd4" ), breaks = c("start_stop", "5UTR", "3UTR", "splicing", "upstream", "downstream", "intergenic", "other", "none"), labels = c("start_stop", "5UTR", "3UTR", "splicing", "upstream", "downstream", "intergenic", "other", "none"))+
#      scale_color_manual("Start/Stop Variant", values = c("#040cd4", "#fa0303"), breaks = c("absent", "start_stop" ), labels = c("Absent", "Start/Stop")) + 
#      scale_color_manual("UTR Variant", values = c("#040cd4", "#be0076", "#fa0303"), breaks = c("absent", "3UTR", "5UTR" ), labels = c("Absent", "3'UTR", "5'UTR")) + 
#      scale_color_manual("Genic/Intergenic Variant", values = c("#040cd4", "#d60060", "#a4008d", "#fa0303"), breaks = c("absent", "intergenic", "downstream", "upstream" ), labels = c("Absent", "Intergenic", "Downstream", "Upstream")) + 
      scale_color_manual("Splicing Variant", values = c("#040cd4", "#fa0303"), breaks = c("absent", "splicing" ), labels = c("Absent", "Splicing")) + 
      theme(legend.position = "top",
            legend.box.background = element_blank(), 
            legend.key = element_blank(), 
            legend.box = "horizontal",
            legend.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            legend.title = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
            axis.text = element_text(size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            axis.title.x = element_blank(),
            panel.grid.major = element_line(colour = "gray90", linetype = "dashed", linewidth = 0.25),
            panel.grid.minor = element_blank()) +
      geom_vline(xintercept = c(-0.6, 0.6), col="black", linetype = "dashed") +
      guides(color = guide_legend(nrow = 2))

    
  # CDF Plot 
#    all_cdf <- ggplot(ann_gex, aes(logFC, color = Event_type)) +
#    stsp_cdf <- ggplot(ann_gex, aes(logFC, color = start_stop_variant)) +
#    utr_cdf <- ggplot(ann_gex, aes(logFC, color = utr_variant)) +
#    genic_cdf <- ggplot(ann_gex, aes(logFC, color = genic_variant)) +
    splicing_cdf <- ggplot(ann_gex, aes(logFC, color = splicing)) +
      stat_ecdf(geom = "step", linewidth = 1.25, alpha = 1) +
      labs(x = "log2(f∆) Mutant Suppressors vs Parent", y = "Fraction of Data") +
      scale_y_continuous(position = "right") +
      coord_cartesian(xlim=c(-5,5)) +
#      scale_color_manual("Variant:", values = c("#fa0303", "#f10038", "#e10053", "#d00065", "#be0076", "#ab0087", "#94009c","#7000bc", "#040cd4" ), breaks = c("start_stop", "5UTR", "3UTR", "splicing", "upstream", "downstream", "intergenic", "other", "none"), labels = c("start_stop", "5UTR", "3UTR", "splicing", "upstream", "downstream", "intergenic", "other", "none"))+
#      scale_color_manual("Start/Stop Variant", values = c("#040cd4", "#fa0303"), breaks = c("absent", "start_stop" ), labels = c("Absent", "Start/Stop")) + 
#      scale_color_manual("UTR Variant", values = c("#040cd4", "#be0076", "#fa0303"), breaks = c("absent", "3UTR", "5UTR" ), labels = c("Absent", "3'UTR", "5'UTR")) + 
#      scale_color_manual("Genic/Intergenic Variant", values = c("#040cd4", "#d60060", "#a4008d", "#fa0303"), breaks = c("absent", "intergenic", "downstream", "upstream" ), labels = c("Absent", "Intergenic", "Downstream", "Upstream")) + 
      scale_color_manual("Splicing Variant", values = c("#040cd4", "#fa0303"), breaks = c("absent", "splicing" ), labels = c("Absent", "Splicing")) + 
      theme(panel.background = element_rect(fill = "NA", color = "black", linewidth = rel(3)), 
            panel.grid = element_blank(),
            axis.text = element_text( size = rel(1.2), family = "Arial", face = "bold", color = "black"),
            axis.title = element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            axis.title.x = element_blank(),
            legend.position = "none",
            #legend.box.background = element_rect("white"), 
            #legend.key = element_blank(), 
            #legend.direction = "horizontal",
            #legend.text = element_text(size = rel(1.25), family = "Arial", face = "bold", color = "black"),
            #legend.title = element_blank(), #element_text(size = rel(1.5), family = "Arial", face = "bold", color = "black"),
            panel.grid.major = element_line(colour = "gray90", linetype = "dashed"))
    
  # Perform Wilcoxon Rank Sum Test with Continuity Correction and Benjamini-Hochberg False Discovery Rate corrections for multiple testing
    
    table(ann_gex$start_stop_variant)
    a <- subset(ann_gex, start_stop_variant == "absent")$logFC
    s <- subset(ann_gex, start_stop_variant == "start_stop")$logFC
      sa <- wilcox.test(s, a, exact = FALSE, alternative = 'less')$p.value 
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(sa)
        p.adjust(raw_pvals, method = "BH") 
    
    table(ann_gex$utr_variant)
    a <- subset(ann_gex, utr_variant == "absent")$logFC
    f <- subset(ann_gex, utr_variant == "5UTR")$logFC
    t <- subset(ann_gex, utr_variant == "3UTR")$logFC
      fa <- wilcox.test(f, a, exact = FALSE, alternative = 'less')$p.value 
      ta <- wilcox.test(t, a, exact = FALSE, alternative = 'less')$p.value 
      tf <- wilcox.test(f, t, exact = FALSE, alternative = 'less')$p.value 
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(fa, ta, tf)
        p.adjust(raw_pvals, method = "BH") 
    
    table(ann_gex$genic_variant)
    a <- subset(ann_gex, genic_variant == "absent")$logFC
    i <- subset(ann_gex, genic_variant == "intergenic")$logFC
    d <- subset(ann_gex, genic_variant == "downstream")$logFC
    u <- subset(ann_gex, genic_variant == "upstream")$logFC
      ia <- wilcox.test(i, a, exact = FALSE, alternative = 'less')$p.value 
      da <- wilcox.test(d, a, exact = FALSE, alternative = 'less')$p.value 
      ua <- wilcox.test(u, a, exact = FALSE, alternative = 'less')$p.value 
      di <- wilcox.test(d, i, exact = FALSE, alternative = 'less')$p.value 
      ui <- wilcox.test(u, i, exact = FALSE, alternative = 'less')$p.value 
      ud <- wilcox.test(u, d, exact = FALSE, alternative = 'less')$p.value 
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(ia, da, ua, di, ui, ud)
        p.adjust(raw_pvals, method = "BH") 
    
    table(ann_gex$splicing)
    a <- subset(ann_gex, splicing == "absent")$logFC
    s <- subset(ann_gex, splicing == "splicing")$logFC
      sa <- wilcox.test(s, a, exact = FALSE, alternative = 'less')$p.value 
      # correction for multiple testing using Benjamini-Hochberg False Discovery Rate
        raw_pvals <- c(sa)
        p.adjust(raw_pvals, method = "BH")
    
  # Design for combining plots and writing to directory
    # layout of plot
        design <- "
      11
      24
      34
      "
    # all
      all_p <- guide_area() + all_vol + all_den + all_cdf + plot_layout(heights = c(1.5,6,1.5), design = design, guides = "collect")
      ggsave(filename = "GEX_volcano_density_cdf.png", p = all_p, device = "png", dpi = 300, height = 6, width = 8, path = file.path(output_dir, "wgs_r_figures", "gex"))
        
    # Start/Stop
      stsp_p <- guide_area() + stsp_vol + stsp_den + stsp_cdf + plot_layout(heights = c(1.5,6,1.5), design = design, guides = "collect")
      ggsave(filename = "GEX_StSp_volcano_density_cdf.png", p = stsp_p, device = "png", dpi = 300, height = 6, width = 8, path = file.path(output_dir, "wgs_r_figures", "gex"))
      
    # UTR
      utr_p <- guide_area() + utr_vol + utr_den + utr_cdf + plot_layout(heights = c(1.5,6,1.5), design = design, guides = "collect")
      ggsave(filename = "GEX_UTR_volcano_density_cdf.png", p = utr_p, device = "png", dpi = 300, height = 6, width = 8, path = file.path(output_dir, "wgs_r_figures", "gex"))
    
    # Genic/Intergenic
      genic_p <- guide_area() + genic_vol + genic_den + genic_cdf + plot_layout(heights = c(1.5,6,1.5), design = design, guides = "collect")
      ggsave(filename = "GEX_genic_volcano_density_cdf.png", p = genic_p, device = "png", dpi = 300, height = 6, width = 8, path = file.path(output_dir, "wgs_r_figures", "gex"))
      
    # splicing
      splicing_p <- guide_area() + splicing_vol + splicing_den + splicing_cdf + plot_layout(heights = c(1.5,6,1.5), design = design, guides = "collect")
      ggsave(filename = "GEX_splicing_volcano_density_cdf.png", p = splicing_p, device = "png", dpi = 300, height = 6, width = 8, path = file.path(output_dir, "wgs_r_figures", "gex"))
      
      
      
      
         
    