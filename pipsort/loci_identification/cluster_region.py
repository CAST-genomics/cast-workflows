import numpy as np
from sklearn.cluster import KMeans

m = np.load("first_region.npy")
m = m.T


print(m.shape)

kmeans = KMeans(n_clusters=3, random_state=0, n_init="auto").fit(m)
np.savetxt("first_region_3cluster_labels.txt", kmeans.labels_)


