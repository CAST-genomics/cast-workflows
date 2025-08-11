import os
import pandas as pd
import numpy as np
import re
import scipy.stats as stats
from cyvcf2 import VCF
from collections import Counter

DEFAULTSAMPLECSV = "passing_samples_v7.1.csv"
SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.1.csv")

def GetPTCovarPath(phenotype):
    return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
        "phenotypes", "%s_phenocovar.csv"%phenotype)

def GetCohort(samplefile):
    cohort = "ALL" #if samplefile == DEFAULTSAMPLECSV
    if samplefile != DEFAULTSAMPLECSV:
        cohort = os.path.basename(samplefile[:-4]) #chop off .csv at the end
    return cohort

def CheckRegion(region):
    if region is None: return True
    return re.match(r"\w+:\d+-\d+", region) is not None

def Inverse_Quantile_Normalization(M):
    M = M.transpose()
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    Q = Q.transpose() 
    return Q

def NormalizeData(data, norm):
    # Add normalization quantile
    if norm == "quantile":
        data["phenotype"] = Inverse_Quantile_Normalization(data[["phenotype"]])
        return data

    # Add z-score normalization
    elif norm == "zscore":
        data["phenotype"]  = stats.zscore(data[["phenotype"]])
        return data

    else:
        ERROR("No normalization method specified")

def read_annotations(filename):
    annotations = []
    if filename:
        with open(filename, "r") as annotations_file:
            lines = annotations_file.readlines()
        for line in lines:
            if line.startswith("chrom"):
                # Header line
                continue
            chrom, start, gene = line.strip().split(",")[:3]
            annotations.append((str(chrom), str(start), str(gene)))
    return annotations

def get_alleles(record):
    ref_allele = record.REF
    alt_alleles = record.ALT
    info_fields = record.INFO
    ru_len = info_fields["PERIOD"]
    if ru_len is None:
        print("Error: Repeat Unit (RU) not provided in the info field")
        return
    #print("ru_len, ru, len ref allele", ru_len, ru, len(ref_allele))
    len_ref_allele = len(ref_allele)
    len_alt_alleles = [len(alt_allele) for alt_allele in alt_alleles]
    return len_ref_allele, len_alt_alleles, ru_len

def get_rc_from_allele(allele, ref_allele, alt_alleles, ru_len):
    if allele == 0:
        ref_rc = int(ref_allele/ru_len)
        #print("ref_rc ", ref_rc)
        return ref_rc
    else:
        rc = int(alt_alleles[allele-1]/ru_len)
        #print("call {} allele {}".format(allele, rc))
        return rc

def set_genotypes(data, annotations, cohort, samples, tr_vcf, outdir):
    print("Reading genotypes from the vcf file")
    # Extract the chrom from the tr_vcf file
    words = re.split("_|\.", tr_vcf)
    vcf_chrom = [word for word in words if word.startswith("chr")][0]
    samples_str = [str(sample) for sample in samples]
    vcf = VCF(tr_vcf, samples = samples_str)
    vcf_samples = vcf.samples

    relevant_annotations = []

    # Focus on a few loci of interest
    for chrom, start, gene in annotations:
        if chrom != vcf_chrom:
            continue

        relevant_annotations.append((chrom, start, gene))
        if not os.path.exists("genotypes/"):
            os.mkdir("genotypes")
        # Check if genotypes are extracted before
        df_load_path = "genotypes/genotypes_{}_{}_{}_{}.csv".format(gene, chrom, start, cohort)
        if not os.path.exists(df_load_path):
            # Load for all samples instead of cohort specific
            df_load_path = "genotypes/genotypes_{}_{}_{}_passing_samples_v7.1.csv".format(gene, chrom, start)
        if os.path.exists(df_load_path):
            ids_not_found = []
            data.loc[:, [gene]] = np.nan
            genotypes_df = pd.read_csv(df_load_path, index_col=0)
            for idx, row in data.iterrows():
                sample = idx
                if int(row["person_id"]) not in genotypes_df.index:
                    ids_not_found.append(int(row["person_id"]))
                else:
                    rc = float(genotypes_df.loc[int(row["person_id"])]["rc"])
                    data.at[idx, gene] = rc
            print("Number of samples where the genotype could not be loaded {}".format(len(ids_not_found)))
            print("For gene {} skipping {} samples".format(gene, len(ids_not_found)) +
                  "4 most common genotypes: ", genotypes_df["rc"].value_counts().head(4))
            continue

        # Compute genotypes from VCF file        
        all_alleles = []
        genotypes_df = pd.DataFrame(columns=["rc"])
        empty_calls, no_calls = 0, 0
        variant = list(vcf("{}:{}-{}".format(chrom, int(start)-1, int(start) + 1)))[0]
        print("Processing annotation {} {}:{}".format(gene, chrom, start))
        data.loc[:, [gene]] = np.nan
        samples_with_calls = set()
        print("Reading ref and alt alleles")
        ref_allele_len, alt_alleles_len, ru_len = get_alleles(variant)
        print("Reading calls")
        counter = 0
        for i, sample_call in enumerate(variant.genotypes):
            counter += 1
            sample = vcf_samples[i]
            if counter % 10000 == 0:
                print("reading call ", counter)
            rc_1 = get_rc_from_allele(int(sample_call[0]), ref_allele_len, alt_alleles_len, ru_len)
            rc_2 = get_rc_from_allele(int(sample_call[1]), ref_allele_len, alt_alleles_len, ru_len)
            rc = (rc_1 + rc_2) / 2.0
            all_alleles.extend([rc_1, rc_2])

            samples_with_calls.add(sample)
            data.loc[data["person_id"]==sample, gene] = rc
            genotypes_df.loc[sample] = rc
        print("Alleles count for the 4 most common alleles: ", Counter(all_alleles).most_common(4))
        if no_calls + empty_calls > 0:
            print("Skipping {} empty calls and {} no calls for {} on vcf".format(
                    empty_calls, no_calls, gene))
        genotypes_df.to_csv(df_load_path)
    return data, relevant_annotations

def print_stats(data, locus):
    print("For phenotype, mean {:.4f} and sd {:.4f}".format(
        data["phenotype"].mean(),
        data["phenotype"].std()))
    print(data["phenotype"].describe())
    print("For genotype, mean {:.4f} and sd {:.4f}".format(
        data[locus].mean(),
        data[locus].std()))
    print(data[locus].describe())
    min_height = data["phenotype"].min()
    max_height = data["phenotype"].max()
    samples_w_min_height = data[data["phenotype"]==min_height]
    print("samples_w_min_height: ", len(samples_w_min_height))
    samples_w_max_height = data[data["phenotype"]==max_height]
    print("samples_w_max_height: ", len(samples_w_max_height))
    
    print("median genotype at samples_w_min_height: ", samples_w_min_height[locus].median())
    print("median genotype at samples_w_max_height: ", samples_w_max_height[locus].median())

def load_snp_gwas_df(filename):
    snp_data = pd.read_csv(filename)
    print("Num snp variants loaded: ", len(snp_data))
    snp_data["POS"] = snp_data["POS"].astype(int)
    snp_data = snp_data.rename(columns={"POS": "pos",
                     "Pvalue_log10": "-log10pvalue",
                     "CHR": "chrom"})
    return snp_data

def load_gwas_tab(filename, max_nlp_val=100):
    data = pd.read_csv(filename, sep="\t")
    data["chrom"] = data["#CHROM"].apply(lambda x: "chr" + str(x))
    data["pos"] = data["POS"].astype(int)
    data["-log10pvalue"] = data["P"].apply(lambda x: -np.log10(float(x))).apply(lambda x: max_nlp_val if x == np.inf else x)
    
    return data