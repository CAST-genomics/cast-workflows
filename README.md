# cast-workflows

Workflows used by the Center for Admixture Science and Technology

See individual folders for more info on each workflow:

* [targetTR](targetTR/README.md): targeted STR genotyping [[UKB launcher](targetTR/launch_ukb/README.md)] [[AoU launcher](targetTR/launch_aou/README.md)]
* [gwas](gwas/aou/README.md): Run SNP or TR GWAS in AoU.
* [subset_vcf](subset_vcf/README.md): subset AoU VCF files (prereq for our imputation/LAI jobs)
* [tr-imputation](tr-imputation/README.md): [Impute TRs from SNP data [UKB launcher](./tr-imputation/launch_ukb/README.md)] [[AoU launcher](./tr-imputation/launch_aou/README.md)]
* [local_ancestry](local_ancestry/README.md): Perform local ancestry inference with gnomix

Dockers that are used by WDL tasks can be found in the [docker](docker/README.md) folder.
