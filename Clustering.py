# /usr/bin/env python

#ライブラリの読み込み
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scrublet as scr
from statsmodels import robust
from statsmodels.robust.scale import mad_scale
import scnpy.external as sce
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
import anndata
import seaborn as sns
import sys
import os

# コマンドライン引数の処理
INPUT_DATA = sys.argv[1]
OUTPUT_DIR = sys.argv[2]
SAMPLE_COL = sys.argv[3]
os.makedirs(OUTPUT_DIR, exist_ok=True)

# データの読み込み
adata = sc.read_h5ad(INPUT_DATA)

#============================================
# 4, クラスタリング・バッチ処理のパラメータ最適化
#============================================

# ==================================
# 4.1, PCs数の最適化 (PCAの最適化)
# ==================================

print("\n[1] Optimizing number of PCs...")
n_pcs_list = [10, 15, 20, 25, 30, 40, 50]
pc_results = []

for n_pcs in n_pcs_list:
    print(f"  Testing {n_pcs} PCs...")
    
    sc.tl.pca(adata, n_comps=n_pcs)
    # バッチ処理（Harmony）
    sce.pp.harmony_integrate(adata, key=SAMPLE_COL, max_iter_harmony=20)
    
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, use_rep='X_pca_harmony')
    sc.tl.leiden(adata, resolution=0.6)
    
    if adata.n_obs > 5000:
        np.random.seed(42)
        sample_idx = np.random.choice(adata.n_obs, 5000, replace=False)
        X_sample = adata.obsm['X_pca_harmony'][sample_idx, :n_pcs]
        labels_sample = adata.obs['leiden'].iloc[sample_idx].astype(int).values
    else:
        X_sample = adata.obsm['X_pca_harmony'][:, :n_pcs]
        labels_sample = adata.obs['leiden'].astype(int).values
    
    sil_score = silhouette_score(X_sample, labels_sample)
    db_score = davies_bouldin_score(X_sample, labels_sample)
    ch_score = calinski_harabasz_score(X_sample, labels_sample)
    n_clusters = len(np.unique(labels_sample))
    
    pc_results.append({
        'n_pcs': n_pcs,
        'silhouette': sil_score,
        'davies_bouldin': db_score,
        'calinski_harabasz': ch_score,
        'n_clusters': n_clusters
    })
    
    print(f"    {n_pcs} PCs: Silhouette={sil_score:.4f}, DB={db_score:.4f}, Clusters={n_clusters}")

pc_df = pd.DataFrame(pc_results)
pc_df.to_csv(os.path.join(OUTPUT_DIR, 'pcs_results.csv'), index=False)

# PCs数の最適化の可視化 (Silhouette, Davies-Bouldin, Calinski-Harabasz)
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].plot(pc_df['n_pcs'], pc_df['silhouette'], 'o-', color='purple')
axes[0].axvline(pc_df.loc[pc_df['silhouette'].idxmax(), 'n_pcs'], color='red', linestyle='--', alpha=0.7)
axes[0].set_xlabel('Number of PCs')
axes[0].set_ylabel('Silhouette Score')
axes[0].set_title('Silhouette Score vs PCs')
axes[0].grid(True, alpha=0.3)

axes[1].plot(pc_df['n_pcs'], pc_df['davies_bouldin'], 'o-', color='green')
axes[1].axvline(pc_df.loc[pc_df['davies_bouldin'].idxmin(), 'n_pcs'], color='red', linestyle='--', alpha=0.7)
axes[1].set_xlabel('Number of PCs')
axes[1].set_ylabel('Davies-Bouldin Index')
axes[1].set_title('DB Index vs PCs (lower is better)')
axes[1].grid(True, alpha=0.3)

axes[2].plot(pc_df['n_pcs'], pc_df['calinski_harabasz'], 'o-', color='orange')
axes[2].axvline(pc_df.loc[pc_df['calinski_harabasz'].idxmax(), 'n_pcs'], color='red', linestyle='--', alpha=0.7)
axes[2].set_xlabel('Number of PCs')
axes[2].set_ylabel('Calinski-Harabasz Index')
axes[2].set_title('CH Index vs PCs (higher is better)')
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'pcs_optimization.png'), dpi=300)
plt.close()

best_n_pcs = int(pc_df.loc[pc_df['silhouette'].idxmax(), 'n_pcs'])
print(f"\nBest number of PCs: {best_n_pcs} (Silhouette: {pc_df['silhouette'].max():.4f})")

#============================================
# 4.2, k-neighborsの最適化 (Leidenの近傍数の最適化)
#============================================

print("\n[2] Optimizing k-neighbors...")
sc.tl.pca(adata, n_comps=best_n_pcs)

# バッチ処理（Harmony）
sce.pp.harmony_integrate(adata, key=SAMPLE_COL, max_iter_harmony=20)

k_values = [5, 10, 15, 20, 25, 30, 40, 50]
k_results = []

for k in k_values:
    print(f"  Testing k={k}...")
    
    sc.pp.neighbors(adata, n_neighbors=k, n_pcs=best_n_pcs, use_rep='X_pca_harmony')
    sc.tl.leiden(adata, resolution=0.6)
    
    if adata.n_obs > 5000:
        np.random.seed(42)
        sample_idx = np.random.choice(adata.n_obs, 5000, replace=False)
        X_sample = adata.obsm['X_pca_harmony'][sample_idx, :best_n_pcs]
        labels_sample = adata.obs['leiden'].iloc[sample_idx].astype(int).values
    else:
        X_sample = adata.obsm['X_pca_harmony'][:, :best_n_pcs]
        labels_sample = adata.obs['leiden'].astype(int).values
    
    sil_score = silhouette_score(X_sample, labels_sample)
    db_score = davies_bouldin_score(X_sample, labels_sample)
    ch_score = calinski_harabasz_score(X_sample, labels_sample)
    n_clusters = len(np.unique(labels_sample))
    
    k_results.append({
        'k': k,
        'silhouette': sil_score,
        'davies_bouldin': db_score,
        'calinski_harabasz': ch_score,
        'n_clusters': n_clusters
    })
    
    print(f"    k={k}: Silhouette={sil_score:.4f}, DB={db_score:.4f}, Clusters={n_clusters}")

k_df = pd.DataFrame(k_results)
k_df.to_csv(os.path.join(OUTPUT_DIR, 'k_neighbors_results.csv'), index=False)

# k-neighborsの最適化の可視化 (Silhouette, Davies-Bouldin, Calinski-Harabasz)
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].plot(k_df['k'], k_df['silhouette'], 'o-', color='orange')
axes[0].axvline(k_df.loc[k_df['silhouette'].idxmax(), 'k'], color='red', linestyle='--', alpha=0.7)
axes[0].set_xlabel('k-neighbors')
axes[0].set_ylabel('Silhouette Score')
axes[0].set_title('Silhouette Score vs k-neighbors')
axes[0].grid(True, alpha=0.3)

axes[1].plot(k_df['k'], k_df['davies_bouldin'], 'o-', color='green')
axes[1].axvline(k_df.loc[k_df['davies_bouldin'].idxmin(), 'k'], color='red', linestyle='--', alpha=0.7)
axes[1].set_xlabel('k-neighbors')
axes[1].set_ylabel('Davies-Bouldin Index')
axes[1].set_title('DB Index vs k-neighbors')
axes[1].grid(True, alpha=0.3)

axes[2].plot(k_df['k'], k_df['calinski_harabasz'], 'o-', color='blue')
axes[2].axvline(k_df.loc[k_df['calinski_harabasz'].idxmax(), 'k'], color='red', linestyle='--', alpha=0.7)
axes[2].set_xlabel('k-neighbors')
axes[2].set_ylabel('Calinski-Harabasz Index')
axes[2].set_title('CH Index vs k-neighbors')
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'k_neighbors_optimization.png'), dpi=300)
plt.close()

best_k = int(k_df.loc[k_df['silhouette'].idxmax(), 'k'])
print(f"\nBest k-neighbors: {best_k} (Silhouette: {k_df['silhouette'].max():.4f})")

#============================================
# 4.3, Resolution最適化　(Leidenの解像度の最適化)
#============================================

print("\n[3] Optimizing resolution...")
sc.pp.neighbors(adata, n_neighbors=best_k, n_pcs=best_n_pcs, use_rep='X_pca_harmony')

resolutions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
resolution_results = []

for res in resolutions:
    print(f"  Testing resolution={res:.2f}...")
    
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_res_{res}')
    
    if adata.n_obs > 5000:
        np.random.seed(42)
        sample_idx = np.random.choice(adata.n_obs, 5000, replace=False)
        X_sample = adata.obsm['X_pca_harmony'][sample_idx, :best_n_pcs]
        labels_sample = adata.obs[f'leiden_res_{res}'].iloc[sample_idx].astype(int).values
    else:
        X_sample = adata.obsm['X_pca_harmony'][:, :best_n_pcs]
        labels_sample = adata.obs[f'leiden_res_{res}'].astype(int).values
    
    sil_score = silhouette_score(X_sample, labels_sample)
    db_score = davies_bouldin_score(X_sample, labels_sample)
    ch_score = calinski_harabasz_score(X_sample, labels_sample)
    n_clusters = len(np.unique(labels_sample))
    
    resolution_results.append({
        'resolution': res,
        'silhouette': sil_score,
        'davies_bouldin': db_score,
        'calinski_harabasz': ch_score,
        'n_clusters': n_clusters
    })
    
    print(f"    Resolution {res:.2f}: Silhouette={sil_score:.4f}, DB={db_score:.4f}, Clusters={n_clusters}")

res_df = pd.DataFrame(resolution_results)
res_df.to_csv(os.path.join(OUTPUT_DIR, 'resolution_results.csv'), index=False)

# Resolution最適化の可視化 (Silhouette, Davies-Bouldin, Calinski-Harabasz)
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes[0, 0].plot(res_df['resolution'], res_df['silhouette'], 'o-', color='blue')
axes[0, 0].axvline(res_df.loc[res_df['silhouette'].idxmax(), 'resolution'], color='red', linestyle='--', alpha=0.7)
axes[0, 0].set_xlabel('Resolution')
axes[0, 0].set_ylabel('Silhouette Score')
axes[0, 0].set_title('Silhouette Score vs Resolution')
axes[0, 0].grid(True, alpha=0.3)

axes[0, 1].plot(res_df['resolution'], res_df['davies_bouldin'], 'o-', color='green')
axes[0, 1].axvline(res_df.loc[res_df['davies_bouldin'].idxmin(), 'resolution'], color='red', linestyle='--', alpha=0.7)
axes[0, 1].set_xlabel('Resolution')
axes[0, 1].set_ylabel('Davies-Bouldin Index')
axes[0, 1].set_title('DB Index vs Resolution')
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].plot(res_df['resolution'], res_df['calinski_harabasz'], 'o-', color='purple')
axes[1, 0].axvline(res_df.loc[res_df['calinski_harabasz'].idxmax(), 'resolution'], color='red', linestyle='--', alpha=0.7)
axes[1, 0].set_xlabel('Resolution')
axes[1, 0].set_ylabel('Calinski-Harabasz Index')
axes[1, 0].set_title('CH Index vs Resolution')
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].plot(res_df['resolution'], res_df['n_clusters'], 'o-', color='orange')
axes[1, 1].set_xlabel('Resolution')
axes[1, 1].set_ylabel('Number of Clusters')
axes[1, 1].set_title('Cluster Count vs Resolution')
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'resolution_optimization.png'), dpi=300)
plt.close()

best_resolution = res_df.loc[res_df['silhouette'].idxmax(), 'resolution']
print(f"\nBest resolution: {best_resolution:.2f} (Silhouette: {res_df['silhouette'].max():.4f})")

#============================================
# 4.4, パラメータ探索のサマリー作成
#============================================
print("\n=== Optimization Summary ===")
print(f"Best number of PCs: {best_n_pcs} (Silhouette: {pc_df['silhouette'].max():.4f})")
print(f"Best k-neighbors: {best_k} (Silhouette: {k_df['silhouette'].max():.4f})")
print(f"Best resolution: {best_resolution:.2f} (Silhouette: {res_df['silhouette'].max():.4f})")

summary_df = pd.DataFrame({
    'parameter': ['n_pcs', 'k_neighbors', 'resolution'],
    'optimal_value': [best_n_pcs, best_k, best_resolution],
    'max_silhouette': [pc_df['silhouette'].max(), 
                      k_df['silhouette'].max(), 
                      res_df['silhouette'].max()]
})
summary_df.to_csv(os.path.join(OUTPUT_DIR, 'optimization_summary.csv'), index=False)

#============================================
# 4.5, バッチ処理と最適化したパラメータで最終クラスタリング
#============================================
print("\n[4.5] Final clustering with optimal parameters...")

# PCA
sc.tl.pca(adata, n_comps=best_n_pcs)
# バッチ処理（Harmony）
sce.pp.harmony_integrate(adata, key=SAMPLE_COL, max_iter_harmony=20)
# 最適化したパラメータで近傍グラフ構築とクラスタリング
sc.pp.neighbors(adata, n_neighbors=best_k, n_pcs=best_n_pcs, use_rep='X_pca_harmony')
# 最適化した解像度でLeidenクラスタリング
sc.tl.leiden(adata, resolution=best_resolution)
# UMAPの計算
sc.tl.umap(adata)

# クラスタリング結果の保存
adata.write_h5ad(os.path.join(OUTPUT_DIR, 'adata_clustered.h5ad'))