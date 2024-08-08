import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from sklearn.preprocessing import normalize
from sklearn.cluster import MiniBatchKMeans
import subprocess
import time

region = "1_10583_1512464"

bucket = os.getenv('WORKSPACE_BUCKET')
saveloc = "/pipsort/preliminary"
'''
start = time.time()
m = np.load(region+".npy")
m = m.T

#m = m[:5000]

print(m.shape)

m = normalize(m, axis=0).astype(np.float16)
print(m.dtype)
end = time.time()
print("time to load and preprocess matrix = ", end - start)

start = time.time()

kmeans = MiniBatchKMeans(n_clusters=3,
                         random_state=0, n_init='auto').fit(m)

end = time.time()
print("fit time = ", end - start)

start = time.time()
labels = kmeans.predict(m)
end = time.time()
print("predict time = ", end - start)

fam = pd.read_csv(region+"_plink.fam", sep="\t", header=None) #column labelled 1 is research ids
fam["label"] = labels


if not os.path.exists("id_ancestry.tsv"):
    subprocess.run("gsutil cp "+bucket+"/pipsort/id_ancestry.tsv ./", shell=True, check=True)
ancestry = pd.read_csv("id_ancestry.tsv", sep="\t")

fam_with_ancestry = fam.merge(ancestry, left_on=1, right_on="research_id", how="inner")
fam_with_ancestry = fam_with_ancestry[['research_id', 'label', 'ancestry_pred', 'ancestry_pred_other']]
fam_with_ancestry.to_csv(region+"_cluster_labels.tsv", sep="\t", index=False)
assert(fam_with_ancestry.shape[0] == fam.shape[0])

grouped = fam_with_ancestry.groupby(['label', 'ancestry_pred']).size().unstack().fillna(0)
grouped.plot(kind='bar', stacked=True)
plt.title('Histogram by Ancestry Prediction')
plt.xlabel('Cluster')
plt.ylabel('Count')
plt.legend(title='Ancestry')
plt.savefig(region+"_ancestry_pred.png")

grouped = fam_with_ancestry.groupby(['label', 'ancestry_pred_other']).size().unstack().fillna(0)
grouped.plot(kind='bar', stacked=True)
plt.title('Histogram by Ancestry Prediction (with Other)')
plt.xlabel('Cluster')
plt.ylabel('Count')
plt.legend(title='Ancestry')
plt.savefig(region+"_ancestry_pred_other.png")
'''
print(bucket+saveloc)
exit(0)

subprocess.run("gsutil cp "+region+"_cluster_labels.tsv "+bucket+saveloc, shell=True, check=True)
subprocess.run("gsutil cp "+region+"_ancestry_pred.png "+bucket+saveloc, shell=True, check=True)
subprocess.run("gsutil cp "+region+"_ancestry_pred_other.png "+bucket+saveloc, shell=True, check=True)
