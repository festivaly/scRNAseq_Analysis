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
import ast

# コマンドライン引数の処理
INPUT_DATA = sys.argv[1]
OUTPUT_DIR = sys.argv[2]
MARKER_FILE = sys.argv[3]
os.makedirs(OUTPUT_DIR, exist_ok=True)

# データの読み込み
adata = sc.read_h5ad(INPUT_DATA)

#============================================
# 6.2, マーカーに基づくアノテーション
#============================================

# マーカー遺伝子のリストの読み込み
with open(MARKER_FILE, "r", encoding="utf-8") as marker_handle:
    markers_dict = ast.literal_eval(marker_handle.read())

all_genes = []
seen_genes = set()
for gene_list in markers_dict.values():
    for gene in gene_list:
        if gene not in seen_genes:
            seen_genes.add(gene)
            all_genes.append(gene)

# クラスターごとにマーカー遺伝子の発現をドットプロットで確認
sc.pl.dotplot(adata, var_names=all_genes, groupby='leiden', save='_marker_dotplot.png')
# クラスターごとにマーカー遺伝子の発現をヒートマップで確認
sc.pl.heatmap(adata, var_names=all_genes, groupby='leiden', save='_marker_heatmap.png')
# クラスターごとにマーカー遺伝子の発現をバイオリンプロットで確認
sc.pl.violin(adata, var_names=all_genes, groupby='leiden', save='_marker_violin.png')
#UMAP上でマーカー遺伝子の発現を可視化
for gene in all_genes:
    sc.pl.umap(adata, color=gene, save='_' + gene + '_umap.png')


# Module Scoreの計算
for cluster, gene_list in markers_dict.items():
    adata.obs[cluster + '_score'] = sc.tl.score_genes(
        adata,
        gene_list=gene_list,
        score_name=cluster + '_score',
        use_raw=False,
    )
# Module Scoreのドットプロット
score_cols = [col for col in adata.obs.columns if col.endswith('_score')]
sc.pl.dotplot(adata, var_names=score_cols, groupby='leiden',
                save='_module_score_dotplot.png')
# Module Scoreのヒートマップ
sc.pl.heatmap(adata, var_names=score_cols, groupby='leiden',
                save='_module_score_heatmap.png')
# Module Scoreのバイオリンプロット
sc.pl.violin(adata, var_names=score_cols, groupby='leiden',
                save='_module_score_violin.png')

# UMAP上でModule Scoreを可視化
for score in score_cols:
    sc.pl.umap(adata, color=score, save='_' + score + '_umap.png')

adata.write_h5ad(os.path.join(OUTPUT_DIR, 'adata_marker_annotated.h5ad'))