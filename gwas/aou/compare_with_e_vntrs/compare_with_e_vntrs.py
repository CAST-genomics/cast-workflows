import pandas as pd


def read_e_vntrs(filename):
    input_file = pd.ExcelFile(filename)
    data = pd.read_excel(input_file, "VNTRs")
    #data = data[["VNTR Coordinates (GRCh38)"]]
    # Exclude VNTRs that did not have an eQTL signal
    data = data.loc[data["eQTL"] == "Yes"]
    data["chrom_start"] = data["VNTR Coordinates (GRCh38)"].apply(lambda x: x.split("-")[0])
    data["chrom"] = data["chrom_start"].apply(lambda x: x.split(":")[0])
    data["start"] = data["chrom_start"].apply(lambda x: x.split(":")[1]).astype(int)
    return data

def read_gwas(filename, is_binary):
    if is_binary:
        columns = "phenotype,POS,ID,REF,ALT,PROVISIONAL_REF?,A1,OMITTED,A1_FREQ,FIRTH?,TEST,OBS_CT,OR,LOG(OR)_SE,Z_STAT,P,ERRCODE".split(",")
    else:
        columns = "phenotype,POS,ID,REF,ALT,PROVISIONAL_REF?,A1,OMITTED,A1_FREQ,TEST,OBS_CT,BETA,SE,T_STAT,P,ERRCODE".split(",")
    data = pd.read_csv(filename, sep="\t", header=0, names=columns, index_col=False)
    data["chrom_start"] = data["ID"].apply(lambda x: x.replace("_", ":"))
    data["chrom"] = data["chrom_start"].apply(lambda x: x.split(":")[0])
    data["start"] = data["chrom_start"].apply(lambda x: x.split(":")[1]).astype(int)
    return data

def main():
    e_vntrs = read_e_vntrs(filename="41467_2021_22206_MOESM3_ESM.xlsx")
    #gwas_sig_hits = read_gwas("results_v4_blood_significant_hits.tab", is_binary=False)
    #gwas_sig_hits = read_gwas("results_v3_eur_significant_hits.tab", is_binary=True)
    gwas_sig_hits = read_gwas("results_v2_all_significant_hits.tab", is_binary=True)
    print(f"Working with {len(e_vntrs)} e-vntrs and {len(gwas_sig_hits)} significant tr-gwas hits")
    gwas_e_vntrs = gwas_sig_hits.merge(e_vntrs, how="inner", on="chrom_start")
    print(gwas_e_vntrs[["chrom_start", "Nearest Gene", "phenotype", "P"]])

    # A one-off check for a Schizophrenia specific VNTR in the literature
    gwas_mir_137_dist = gwas_sig_hits.loc[gwas_sig_hits["chrom"] == "chr1"]["start"].apply(lambda x: abs(x - 98046174))
    #print("Distances of VNTRs with significant hit to the MIR137 VNTR: ", gwas_mir_137_dist.sort_values())

if __name__ == "__main__":
    main()
