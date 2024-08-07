import numpy as np
from sklearn.preprocessing import normalize
from sklearn.cluster import MiniBatchKMeans
import time

start = time.time()
m = np.load("1_10583_1512464.npy")
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

np.savetxt("first_region_3cluster_labels.txt", labels)


