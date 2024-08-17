import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import pandas as pd
import subprocess
import sys

bucket = os.getenv('WORKSPACE_BUCKET')

fam = sys.argv[1]
fileprefix = fam[:-4]
kmeans_labels = sys.argv[2]
pca_components = sys.argv[3]
logfile = sys.argv[4]

labels = np.load(kmeans_labels)
pcs = np.load(pca_components)
pc1 = pcs[:,0]
pc2 = pcs[:,1]

fam = pd.read_csv(fam, sep="\t", header=None) #column labelled 1 is research ids
fam["label"] = labels
fam['pc1'] = pc1
fam['pc2'] = pc2


if not os.path.exists("/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/id_ancestry.tsv"):
    subprocess.run("gsutil cp "+bucket+"/pipsort/id_ancestry.tsv /home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/", shell=True, check=True)
ancestry = pd.read_csv("/home/jupyter/workspaces/impactofglobalandlocalancestryongenomewideassociationv7v6studies/cast-workflows/pipsort/loci_identification/id_ancestry.tsv", sep="\t")

fam_with_ancestry = fam.merge(ancestry, left_on=1, right_on="research_id", how="inner")
fam_with_ancestry = fam_with_ancestry[['research_id', 'label', 'ancestry_pred', 'ancestry_pred_other','pc1','pc2']]
#fam_with_ancestry.to_csv(region+"_cluster_labels.tsv", sep="\t", index=False)
assert(fam_with_ancestry.shape[0] == fam.shape[0])

log = file = open(logfile, 'w')
num_clusters = fam_with_ancestry['label'].nunique()
log.write("Ancestry pred counts per cluster\n")
for i in range(num_clusters):
    log.write("Cluster " + str(i) + " \n")
    subdf = fam_with_ancestry[fam_with_ancestry['label']==i]
    counts = subdf['ancestry_pred'].value_counts()
    csv_string = counts.to_csv(header=False, index_label='Value', path_or_buf=None)
    log.write(csv_string)
log.write("Ancestry pred other counts per cluster\n")
for i in range(num_clusters):
    log.write("Cluster " + str(i) + " \n")
    subdf = fam_with_ancestry[fam_with_ancestry['label']==i]
    counts = subdf['ancestry_pred_other'].value_counts()
    csv_string = counts.to_csv(header=False, index_label='Value', path_or_buf=None)
    log.write(csv_string)
log.close()


grouped = fam_with_ancestry.groupby(['label', 'ancestry_pred']).size().unstack().fillna(0)
grouped.plot(kind='bar', stacked=True)
plt.title('Histogram by Ancestry Prediction')
plt.xlabel('Cluster')
plt.ylabel('Count')
plt.legend(title='Ancestry')
plt.savefig(fileprefix+"_ancestry_pred.png")

plt.clf()

grouped = fam_with_ancestry.groupby(['label', 'ancestry_pred_other']).size().unstack().fillna(0)
grouped.plot(kind='bar', stacked=True)
plt.title('Histogram by Ancestry Prediction (with Other)')
plt.xlabel('Cluster')
plt.ylabel('Count')
plt.legend(title='Ancestry')
plt.savefig(fileprefix+"_ancestry_pred_other.png")
plt.clf()

#PCS


sns.scatterplot(x='pc1', y='pc2', hue='ancestry_pred', data=fam_with_ancestry, palette='Set1')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PC1 vs PC2 Colored by Ancestry')
plt.savefig(fileprefix+"_pcs_ancestry_pred.png")
plt.clf()
sns.scatterplot(x='pc1', y='pc2', hue='ancestry_pred_other', data=fam_with_ancestry, palette='Set1')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PC1 vs PC2 Colored by Ancestry')
plt.savefig(fileprefix+"_pcs_ancestry_pred_other.png")
plt.clf()
