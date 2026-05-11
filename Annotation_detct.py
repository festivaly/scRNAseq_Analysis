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
import json

# コマンドライン引数の処理
INPUT_DATA = sys.argv[1]
OUTPUT_DIR = sys.argv[2]
SAMPLE_COL = sys.argv[3]
MARKER_FILE = sys.argv[4]
ANNOTATION_FILE = sys.argv[5]
os.makedirs(OUTPUT_DIR, exist_ok=True)

# データの読み込み
adata = sc.read_h5ad(INPUT_DATA)

#============================================
# 6.3, アノテーションのマッピング
#============================================

# クラスターごとにアノテーションをマッピング
with open(ANNOTATION_FILE, "r", encoding="utf-8") as annotation_handle:
	annotation_map = json.load(annotation_handle)

# クラスターIDを文字列に変換してアノテーションをマッピング
adata.obs['leiden'] = adata.obs['leiden'].astype(str)
adata.obs['cell_type'] = adata.obs['leiden'].map(annotation_map)
# UMAP上でアノテーションを可視化
sc.pl.umap(adata, color='cell_type', save='_annotation_umap.png')
# Cell typeごとのマーカー遺伝子の発現をドットプロットで確認
sc.pl.dotplot(adata, var_names=adata.var_names, groupby='cell_type', save='_annotation_dotplot.png')
# Cell typeごとのマーカー遺伝子の発現をヒートマップで確認
sc.pl.heatmap(adata, var_names=adata.var_names, groupby='cell_type', save='_annotation_heatmap.png')
# Cell typeごとのマーカー遺伝子の発現をバイオリンプロットで確認
sc.pl.violin(adata, var_names=adata.var_names, groupby='cell_type', save='_annotation_violin.png')

# データの保存
adata.write_h5ad(os.path.join(OUTPUT_DIR, 'annotated_data.h5ad'))