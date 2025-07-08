library(org.Hs.eg.db)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(coloc)

# Define the diseases and their corresponding files
# Define all diseases
diseases <- c(
  "Parkinson", "Osteoarthritis", "BC", "CAD", 
  "CKD", "Colorectal_cancer", "Gastric_cancer", 
  "Hyperlipidemia", "Osteoporosis"
)
gwas_files <- paste0(diseases, "_GWAS_summary.GCTA-COJO.tsv")

# Define base paths
gwas_base_path <- "methods/Colocalization/GWAS/"
gtex_base_path <- "methods/Colocalization/GTEx/"


# Complete function with gene-level analysis
run_coloc_analysis <- function(path, gwas, disease_name, remove_duplicate_snps = TRUE) {
  # Create disease-specific output directory if it doesn't exist
  output_dir <- file.path(gtex_base_path, disease_name)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # Step 1: Load the data
  df <- read.csv(path, sep = "\t", header = TRUE)
  
  # Filter data based on p-value threshold
  filtered_df <- df[df$pvalue < 1e-4, ]
    
  # Remove duplicates at GTEx level, keeping the most significant association
  filtered_df <- filtered_df %>%
    group_by(rsid, symbol) %>%
    slice_min(order_by = pvalue, n = 1) %>%
    ungroup()
  
  # Merge filtered data with provided GWAS data
  merged_data <- merge(gwas, 
                      filtered_df, by.x = "SNP", by.y = "rsid", 
                      suffixes = c("_GWAS", "_GTEX"), all = FALSE)
  
  # Clean merged data
  merged_data_cleaned <- merged_data[, 
                                   !names(merged_data) %in% c("molecular_trait_id", 
                                                            "molecular_trait_object_id")]
  merged_data_cleaned_no_na <- merged_data_cleaned[!is.na(merged_data_cleaned$symbol) & 
                                                   merged_data_cleaned$symbol != "", ]
  merged_data_cleaned_no_na$n <- merged_data_cleaned_no_na$an / 2
  
  # Initialize empty data frame for results
  coloc_results <- data.frame()
  
  # Loop through each unique gene symbol
  unique_genes <- unique(merged_data_cleaned_no_na$symbol)
  
  for (gene in unique_genes) {
    gene_data <- merged_data_cleaned_no_na %>% filter(symbol == gene)
    
    # Define dataset1 for GWAS
    dataset1 <- list(
      pvalues = gene_data$P, 
      type = "cc", 
      s = gene_data$s, 
      N = gene_data$N, 
      snp = gene_data$SNP, 
      beta = gene_data$b, 
      varbeta = (gene_data$se_GWAS)^2,
      MAF = gene_data$maf_GWAS
    )
    
    # Define dataset2 for eQTL
    dataset2 <- list(
      pvalues = gene_data$pvalue, 
      type = "quant", 
      N = gene_data$n, 
      snp = gene_data$SNP, 
      beta = gene_data$beta, 
      varbeta = (gene_data$se_GTEX)^2,
      MAF = gene_data$maf_GTEX
    )
    
    # Run coloc analysis for the gene
    coloc_result <- coloc.abf(dataset1, dataset2)
    
    # Extract SNP-level H4 values
    snp_h4 <- coloc_result$results %>% select(snp, SNP.PP.H4)
    snp_h4 <- snp_h4 %>% mutate(GeneSymbol = gene)
    gene_h4 <- coloc_result$summary["PP.H4.abf"]
    snp_h4 <- snp_h4 %>% mutate(Gene_PP.H4.abf = gene_h4)
    
    # Append results for the gene
    coloc_results <- bind_rows(coloc_results, snp_h4)
  }
  
  # Filter coloc results based on thresholds
  if (nrow(coloc_results) > 0 && "SNP.PP.H4" %in% names(coloc_results) && "Gene_PP.H4.abf" %in% names(coloc_results)) {
      coloc_results_filtered <- coloc_results %>%
        filter(SNP.PP.H4 > 0.8, Gene_PP.H4.abf > 0.8)
  
      # Merge with original data
      if (nrow(coloc_results_filtered) > 0) {
          merged_final <- merge(
            coloc_results_filtered, 
            merged_data_cleaned_no_na, 
            by.x = c("snp", "GeneSymbol"), 
            by.y = c("SNP", "symbol"), 
            all = FALSE
          )
  
          # Remove duplicate SNPs if specified
          if (remove_duplicate_snps) {
            merged_final_unique <- merged_final %>%
              group_by(snp) %>%
              slice_min(order_by = pvalue, n = 1) %>%
              ungroup()
          } else {
            merged_final_unique <- merged_final
          }

          # Calculate Zscore and new_Zscore
          merged_final_unique <- merged_final_unique %>%
            mutate(Zscore = beta / se_GTEX,
                   new_Zscore = ifelse(b * beta > 0, Zscore, -Zscore))

          # Define output path with disease-specific name
          filename <- basename(path)
          new_path <- file.path(output_dir, sub("\\.tsv$", paste0(".", disease_name, ".base_data.tsv"), filename))

          # Save to TSV with selected columns
          merged_final_unique %>%
            select(snp, A1, new_Zscore, Zscore, GeneSymbol) %>%
            write.table(file = new_path, sep = "\t", row.names = FALSE, quote = FALSE)
       } else {
          message(paste("No SNPs passed colocalization thresholds for", basename(path)))
          return(NULL)
        }
    } else {
      message(paste("No colocalization results available for", basename(path)))
      return(NULL)
      }
}

# Main loop to process each disease
for (i in seq_along(diseases)) {
  disease <- diseases[i]
  gwas_file <- gwas_files[i]
  
  message(paste("Processing", disease))
  
  # Read and process GWAS data
  gwas_data <- read.table(file.path(gwas_base_path, gwas_file),
                         header = TRUE, sep = '\t')
  
  # Filter rows where P < 1e-4
  gwas_filtered <- gwas_data %>% 
    filter(P < 1e-4) %>%
    mutate(maf = ifelse(freq > 0.5, 1 - freq, freq))
  
  # Get all .all.tsv files
  file_paths <- list.files(gtex_base_path, 
                          pattern = "\\.annotated.tsv$", 
                          full.names = TRUE)
  
  # Run coloc analysis for current disease
  lapply(file_paths, function(f) {
    run_coloc_analysis(f, 
                      gwas = gwas_filtered, 
                      disease_name = disease,
                      remove_duplicate_snps = TRUE)
  })
  
  message(paste("Completed", disease))
}

message("All analyses completed")