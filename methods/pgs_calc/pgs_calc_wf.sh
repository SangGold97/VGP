while read outdir trait_efo
do

nextflow run pgscatalog/pgsc_calc -r v2.0.0-beta \
	-profile docker \
	--input dataset_sheet_T2D.csv --target_build GRCh38 \
	--trait_efo "$trait_efo" \
	--efo_direct True \
	--genotypes_cache ./cache \
	--max_cpus 128 \
	--max_memory 750.GB \
	--min_overlap 0.01 \
	-resume \
	--outdir "$outdir"
done < list_trait.txt