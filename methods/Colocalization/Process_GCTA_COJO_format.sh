#!/bin/bash
outpath="methods/Colocalization/GWAS/"

# Read the list of files
while IFS=$'\t' read -r filename; do
    # Create output filename by adding _cojo suffix before the extension
    outfile="${outpath}${filename%.*}.GCTA-COJO.tsv"
    
    # Print header for the COJO format file
    echo -e "SNP\tA1\tA2\tfreq\tb\tse\tP\tN\ts" > "$outfile"
    
    # Process the file using awk
    # Skip the header line (NR>1) and reorder/rename columns as needed
    awk 'BEGIN {FS=OFS="\t"}
    NR>1 {
        print $4,     # SNP (SNPID)
              $6,     # A1 (Allele2)
              $5,     # A2 (Allele1)
              $8,     # freq (AF_Allele2)
              $11,    # b (BETA)
              $12,    # se (SE)
              $14,    # P (p.value)
              $10,    # N
              "NA"    # s (placeholder)
    }' "$filename" >> "$outfile"
    
    echo "Processed $filename -> $outfile"
done < list_sumstats_file.tsv
