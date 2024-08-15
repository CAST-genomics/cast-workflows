import numpy as np
import os
import pandas as pd
from sklearn.preprocessing import normalize
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
import subprocess
import sys
import time

npy = sys.argv[1]
fileprefix = sys.argv[2]

start = time.time()
m = np.load(npy)
m = m.T #need to transpose because it will be numvariants x numsamples and we want the transpoe


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
np.save(fileprefix+"_kmeans_labels", labels)
end = time.time()
print("predict time = ", end - start)

#PCA
start = time.time()
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca.fit(m.T)
end = time.time()
print("pca time = ", end - start)
np.save(fileprefix+"_pcs", pca.components_.T)

