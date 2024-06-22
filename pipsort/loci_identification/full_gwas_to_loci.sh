phen_name=$1
fname_gwas1=$2
fname_gwas2=$3
fsl1=5e-08
fsl2=1e-06
pop1=NOT_AFR
pop2=AFR


gsutil cp "${WORKSPACE_BUCKET}/gwas/$phen_name/$fname_gwas1" ./
gsutil cp "${WORKSPACE_BUCKET}/gwas/$phen_name/$fname_gwas2" ./

python find_peak_centers.py --infile $fname_gwas1 --outfile ${phen_name}_${pop1}_peaks.csv --pval_cutoff $fsl1 --pval_col p_value --pos_col pos --chr_col chrom
python find_peak_centers.py --infile $fname_gwas2 --outfile ${phen_name}_${pop2}_peaks.csv --pval_cutoff $fsl2 --pval_col p_value --pos_col pos --chr_col chrom

python add_rsid_col.py ${phen_name}_${pop1}_peaks.csv ,
python add_rsid_col.py ${phen_name}_${pop2}_peaks.csv ,


python construct_loci.py --s1_signals ${phen_name}_${pop1}_peaks.csv --s2_signals ${phen_name}_${pop2}_peaks.csv --delim , --outfile ${phen_name}_${pop1}_${pop2}_1MB_loci.txt --pval_col p_value --pos_col pos --chr_col chrom --id_col rsid --chr_name_prefix chr

python postprocess_loci.py ${phen_name}_${pop1}_${pop2}_1MB_loci.txt fixed_${phen_name}_loci.txt
mv fixed_${phen_name}_loci.txt ${phen_name}_${pop1}_${pop2}_1MB_loci.txt

gsutil cp ${phen_name}_${pop1}_${pop2}_1MB_loci.txt "${WORKSPACE_BUCKET}/pipsort/loci"



