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
SAMPLE_COL = sys.argv[3]  # サンプルIDを格納する列名
os.makedirs(OUTPUT_DIR, exist_ok=True)

# データの読み込み
adata = sc.read_h5ad(INPUT_DATA)

#============================================
# 5. クラスタリングの可視化
#============================================

# クラスタリング結果のUMAPプロット
sc.pl.umap(adata, color='leiden', save='_leiden_clusters.png')
# サンプルごとのUMAPプロット
sc.pl.umap(adata, color=SAMPLE_COL, save=f'_{SAMPLE_COL}.png')
# ミトコンドリア遺伝子の割合をUMAP上で可視化
sc.pl.umap(adata, color='pct_counts_mt', save='_pct_counts_mt.png')
