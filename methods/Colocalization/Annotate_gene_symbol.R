library(org.Hs.eg.db)
library(readr)
library(dplyr)



# Define base path for GTEx files
gtex_base_path <- "methods/Colocalization/GTEx/"

# Get all GTEx .all.tsv files
file_paths <- list.files(gtex_base_path, pattern = "\\.all.tsv$", full.names = TRUE)

# Function to annotate GTEx data
annotate_gtex <- function(file_path) {
  message(paste("Processing:", file_path))
  
  # Load GTEx data
  df <- read_delim(file_path, delim = "\t", col_types = cols())
  
  # Annotate gene symbols using ENSEMBL IDs
  df <- df %>%
    mutate(symbol = mapIds(org.Hs.eg.db,
                           keys = gene_id,
                           keytype = "ENSEMBL",
                           column = "SYMBOL",
                           multiVals = "first"))
  
  # Save the annotated file
  output_file <- sub("\\.all.tsv$", ".annotated.tsv", file_path)
  write_tsv(df, output_file)
  
  message(paste("Annotated file saved:", output_file))
}

# Apply annotation function to all GTEx files
lapply(file_paths, annotate_gtex)

message("All GTEx files have been annotated.")