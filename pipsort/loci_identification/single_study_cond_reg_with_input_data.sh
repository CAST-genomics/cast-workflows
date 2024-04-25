data=$1 #e.g. eur_data
gwas0=$2
phenocovar=$3
pvalthresh=$4
lead_snps=$5
max_num_reg=$6
logfile=cond_reg.log


nc=1 #num causal
plink_loc=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification
scripts=/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification

python $scripts/split_phenocovar.py $phenocovar f_phen f_covars

python $scripts/check_stopping_criteria_cond_reg.py --infile $gwas0 --pvalthresh $pvalthresh --pval_col p_value > stop_file 
stop=$( cat stop_file )
echo "STOP IS:"
echo $stop
rm stop_file
cntr=0
cp $gwas0 gwas${cntr}.phenotype.glm.linear
pval_col=p_value
while [ $stop -ne 0 ];
do
	#python script that identifies lead snp and LD friends
	python $scripts/get_lead_snp.py --infile gwas${cntr}.phenotype.glm.linear --pval_col $pval_col --rsid_col rsid > lead_snp_file
	pval_col=P #reset after evaluating first gwas
	lead_snp=$( cat lead_snp_file )
	echo $lead_snp
	echo $lead_snp >> $lead_snps
	rm lead_snp_file #remove after regressing out

	$plink_loc/plink2 --bfile $data --snps $lead_snp --recode A #output goes to plink2.raw
	python $scripts/add_covar.py --covars f_covars --genofile plink2.raw --snp ${lead_snp} #${lead_snp}_{some_capital_alphabet_letter} is the name of the col with genotypes from plink.raw

	#next round gwas
	cntr=$((cntr+1))
	$plink_loc/plink2 --glm hide-covar --bfile $data --pheno f_phen --covar f_covars --out gwas${cntr}


	python $scripts/check_stopping_criteria_cond_reg.py --infile gwas${cntr}.phenotype.glm.linear --pvalthresh $pvalthresh > stop_file 
	stop=$( cat stop_file )
	echo "STOP IS:"
	echo $stop
	rm stop_file
	if [[ $cntr -ge $max_num_reg ]]
	then
		echo $cntr
		echo $max_num_reg
		break
	fi
done


#rm remaining_snps.*
for (( i=0 ; i<=$cntr ; i++ ));
do
    rm gwas${i}.*
done
rm f_phen
rm f_covars
