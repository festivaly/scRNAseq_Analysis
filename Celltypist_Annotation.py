# /# /usr/bin/env python

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
import celltypist
from celltypist import models

# コマンドライン引数の処理
INPUT_DATA = sys.argv[1]
OUTPUT_DIR = sys.argv[2]
MODEL_NAME = sys.argv[3]
os.makedirs(OUTPUT_DIR, exist_ok=True)

# データの読み込み
adata = sc.read_h5ad(INPUT_DATA)

#=============================================
# 6.1, CellTypistによる自動アノテーション
#=============================================

# CellTypistのモデルをダウンロードして読み込む
models.download_models(model=MODEL_NAME, force_update=True)
model = models.Model.load(model=MODEL_NAME)
predictions = celltypist.annotate(
    adata,
    model=model,
    majority_voting=True
)
adata = predictions.to_adata()

print("Celltypist_Process_finished")

adata.write_h5ad(os.path.join(OUTPUT_DIR, 'adata_annotated_celltypist.h5ad'))
