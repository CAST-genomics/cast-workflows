import subprocess
import os
import sys

bucket = os.getenv('WORKSPACE_BUCKET')

phenotype = sys.argv[1]
full_gwas_cmd_file = sys.argv[2]

chr_regions_file = "chrom_regions_hg38.txt"
if (len(sys.argv) == 4):
    chr_regions_file = sys.argv[3] #for debugging on a smaller set of regions

with open(chr_regions_file, "r") as f:
    chr_regions = [line.strip() for line in f.readlines()]


f = open(full_gwas_cmd_file, 'r')
gwas_cmd = f.readline()
gwas_cmd = gwas_cmd.strip() #remove newline from cmd to concatenate region later
f.close()

for r in chr_regions:
    chrom = r.split(":")[0]
    full_cmd = gwas_cmd + " --region " + r
    subprocess.run(full_cmd, shell=True, check=True)
    files = os.listdir("./")
    filtered_files = [file for file in files if chrom in file]
    for f in filtered_files:
        copy_cmd = "gsutil cp " + f + " " + bucket+"/gwas/"+phenotype+"/"+f
        print(copy_cmd)
        subprocess.run(copy_cmd, shell=True, check=True)

