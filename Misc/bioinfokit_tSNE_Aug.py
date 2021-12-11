import time
import numpy as np
import pandas as pd

import matplotlib_venn
import adjustText

from bioinfokit.visuz import cluster

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import seaborn as sns


df = pd.read_table("D:/scRNASeq/GC/Stanford/5846_n1/matrix.mtx", header="infer")
#df = df.set_index(df.columns[0])
dft = df.T

pca_scores = PCA().fit_transform(dft)
df_pc = pd.DataFrame(pca_scores)
tsne_em = TSNE(n_components=2, perplexity=41.0, early_exaggeration=12, n_iter=1000, learning_rate=368, verbose=1).fit_transform(df_pc.loc[:,0:49])

cluster.tsneplot(score=tsne_em)
get_clusters = DBSCAN(eps=3, min_samples=10).fit_predict(tsne_em)

palette = sns.color_palette("viridis",len(set(get_clusters)))
x = cluster.tsneplot(score=tsne_em, colorlist=get_clusters, 
    colordot=palette, 
    legendpos='upper right', legendanchor=(1.15, 1))
