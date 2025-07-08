diseases=(
    "Parkinson"
    "Osteoarthritis" 
    "BC"
    "CAD"
    "CKD"
    "Colorectal_cancer"
    "Gastric_cancer"
    "Hyperlipidemia"
    "Osteoporosis"
)

for disease in "${diseases[@]}"; do

echo "Processing $disease"

    disease="Osteoarthritis"
    bfile=""
    GTEx_path="methods/Colocalization/""$disease""/"
    result_path="./TRS_GTEx_results/""$disease""/"

    while read dataset_id study_label sample_group
    do
        plink --bfile "$bfile" \
        --out "$result_path""$dataset_id" \
        --score "$GTEx_path""$dataset_id".annotated."$disease".base_data.tsv 1 2 3 header
    done<list_GTEx.tsv

done