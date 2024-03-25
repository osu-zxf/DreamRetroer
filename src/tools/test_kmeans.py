from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def adaptive_clustering(data, max_clusters):
    # 初始化聚类数量和评估指标
    best_cluster_count = 2
    best_silhouette_score = -1

    # 迭代聚类数量
    for cluster_count in range(2, max_clusters + 1):
        # 进行K-means聚类
        kmeans = KMeans(n_clusters=cluster_count, random_state=0)
        labels = kmeans.fit_predict(data)

        # 计算轮廓系数
        silhouette_avg = silhouette_score(data, labels)

        # 判断是否得到更好的聚类结果
        if silhouette_avg > best_silhouette_score:
            best_silhouette_score = silhouette_avg
            best_cluster_count = cluster_count

    # 使用最佳聚类数量重新进行聚类
    best_kmeans = KMeans(n_clusters=best_cluster_count, random_state=0)
    best_labels = best_kmeans.fit_predict(data)

    return best_labels

# 示例用法
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
import torch

def load_molecules():
    routes = []
    route_file_name = '../dataset/retro190.txt'
    for line in open(route_file_name, "r"):
        routes.append(line)
    print('%d routes extracted from %s loaded' % (len(routes), route_file_name))
    return routes

def preprocess(X, fp_dim):
    # Compute fingerprint from mol to feature
    mol = Chem.MolFromSmiles(X)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=int(fp_dim),useChirality=True)
    onbits = list(fp.GetOnBits())
    arr = np.zeros(fp.GetNumBits())
    arr[onbits] = 1
    return arr

target_products = load_molecules()
target_products_fps = [preprocess(mol, 256) for mol in target_products]

max_clusters = 10  # 最大聚类数量

labels = adaptive_clustering(target_products_fps, max_clusters)
print(labels)