import pandas as pd  
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
from sklearn.decomposition import PCA
import os

os.environ['MKL_NUM_THREADS'] = "4"
os.environ['OMP_NUM_THREADS'] = "4"

if __name__ == '__main__':

	# read in files and get rid of annoying end of line character from ccheminfolib
	arr_X = pd.read_csv('aso_combined.csv', index_col= 0, header = None)
	arr_X.dropna(axis = 1, inplace=True)
	# list of catalyst barcodes
	X_label = list(arr_X.index)

	# project to a 200 dimensional PCA space
	pca = PCA(n_components=200)
	arr_X = pca.fit_transform(arr_X) 

	print(sum(pca.explained_variance_ratio_))

	# open files to save clustering results
	nf = open('merck_box6cluster_PCA200_re.csv', 'w')
	nf2 = open('merck_box6cluster_PCA200_re_exemplars.csv', 'w')

	# do the KMeans
	kmeans = KMeans(n_clusters=6, random_state=0, n_init=1000).fit(arr_X)

	# idenitfy medoids from cluster exemplars
	closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, arr_X)

	# save cluster assignments and exemplars
	for i, j in enumerate(closest):
		nf2.write(X_label[j] +','+ str(i) + '\n')
	for i in range(len(kmeans.labels_.tolist())):
		nf.write(X_label[i] + ',' + str(kmeans.labels_.tolist()[i]) + ',' + '\n')
	nf.close()
	nf2.close()
