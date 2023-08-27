import numpy as np
import pandas as pd  
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import os

# from Andrew Zahrt

os.environ['MKL_NUM_THREADS'] = "4"
os.environ['OMP_NUM_THREADS'] = "4"

if __name__ == '__main__':

	# read in descriptors
	arr_X = pd.read_csv('desc_calc/aso_combined.csv', index_col= 0, header = None)
	arr_X.dropna(axis = 1, inplace=True)
	
	# do PCA projection
	pca = PCA(n_components=200)
	arr_X = pca.fit_transform(arr_X)

	print(sum(pca.explained_variance_ratio_))
	ev = sum(pca.explained_variance_ratio_)
	
	# calculate kmeans for k=1 to k=30
	distortions = []
	for k in range(1,30):
		print('started k=' + str(k))
		# calculate kmeans
		kmeanModel = KMeans(n_clusters=k, random_state=0, n_init=1000).fit(arr_X)
		kmeanModel.fit(arr_X)
		# calculate distortion, the sum of distances of ligand to closest cluster centroid for all ligands.
		distortions.append(np.sum(np.min(cdist(arr_X, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / arr_X.shape[0])
		print('did k=' + str(k))

	# make the plot
	plt.plot(range(1,30), distortions, 'bx-')
	plt.xlabel('k')
	plt.ylabel('Distortion')
	plt.title(f'Elbow Plot for PCA200 ASO Space (ev ={ev:.3f})')
	plt.savefig('elbow_plot_esp_PCA200.png')



