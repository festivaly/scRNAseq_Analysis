# /usr/bin/env python

#ライブラリの読み込み
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scrublet as scr
from statsmodels import robust
from statsmodels.robust.scale import mad_scale
import anndata
import seaborn as sns
import sys
import os

# コマンドライン引数の処理
INPUT_DATA = sys.argv[1]
OUTPUT_DIR = sys.argv[2]
SAMPLE_COL = sys.argv[3]
os.makedirs(OUTPUT_DIR, exist_ok=True)

#=============================================
# 2, データの読み込み
#=============================================

# データの読み込み
adata = sc.read_h5ad(INPUT_DATA)

#=============================================
# 3, 前処理 (正規化・log化・QC)
#=============================================

# QC指標: mt, ribo, hb の割合を計算
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS','RPL'))
adata.var['hb'] = adata.var_names.str.startswith(('HBA','HBB'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)

# QC指標 (総カウント数、遺伝子数) のlog変換
adata.obs['total_counts_log'] = np.log10(adata.obs['total_counts'])
adata.obs['n_genes_by_counts_log'] = np.log10(adata.obs['n_genes_by_counts'])

# サンプルを総カウント数で整列
df = adata.obs.copy()
df["total_counts_log"] = np.log10(df["total_counts"])
df["n_genes_by_counts_log"] = np.log10(df["n_genes_by_counts"])

order_counts = (
    df.groupby("sample")["total_counts_log"].median().sort_values().index.tolist()
)

# QC指標のバイオリンプロットをサンプルごとに作成
sc.pl.violin(
    adata,
    keys="total_counts_log",
    groupby="sample",
    stripplot=False,
    order=order_counts,
    show=False, 
    palette=["gray"]
)
plt.xticks(rotation=90)
plt.gcf().set_size_inches(30, 10)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "violin_total_counts_CON_raw.png"), dpi=300)
plt.close()

sc.pl.violin(
    adata,
    keys="n_genes_by_counts_log",
    groupby="sample",
    stripplot=False,
    order=order_counts,
    show=False,
    palette=["gray"]
)
plt.xticks(rotation=90)
plt.gcf().set_size_inches(30, 10)
plt.savefig(os.path.join(OUTPUT_DIR, "violin_n_genes_CON_raw.png"), dpi=300)
plt.close()

sc.pl.violin(
  adata,
  keys='pct_counts_mt',
  groupby='sample',
  stripplot=False,
  palette=["gray"],
  order=order_counts,
  show=False
)
plt.xticks(rotation=90)
plt.gcf().set_size_inches(30, 10)
plt.savefig(os.path.join(OUTPUT_DIR, "violin_pct_counts_mt_CON_raw.png"), dpi=300)
plt.close()

# 低発現遺伝子を除外
sc.pp.filter_genes(adata, min_cells=3)

# サンプル別にロバスト(MAD)なQCしきい値を適用
def _mad_limits(x, low_k=3.0, high_k=3.0):
    m = np.median(x)
    s = mad_scale(x, center=m)
    if not np.isfinite(s) or s == 0:
        low = np.percentile(x, 5)
        high = np.percentile(x, 95)
        return low, high
    return m - low_k*s, m + high_k*s

cap_mt = 20.0
qc_keep = pd.Series(True, index=adata.obs_names, name='qc_keep')
qc_summary = []
for smp, idx in adata.obs.groupby(SAMPLE_COL).groups.items():
    sub = adata.obs.loc[idx]
    low_g, high_g = _mad_limits(sub['n_genes_by_counts'], 3, 3)
    min_genes = max(200, low_g)
    max_genes = min(8000, high_g)
    low_c, high_c = _mad_limits(sub['total_counts'], 0, 4)
    max_counts = min(np.percentile(sub['total_counts'], 99.5), high_c)
    _, high_mt = _mad_limits(sub['pct_counts_mt'], 0, 3)
    max_mt = min(cap_mt, high_mt)
    keep = (
        (sub['n_genes_by_counts'] >= min_genes) &
        (sub['n_genes_by_counts'] <= max_genes) &
        (sub['total_counts'] <= max_counts) &
        (sub['pct_counts_mt'] <= max_mt)
    )
    qc_keep.loc[idx] = keep.values
    qc_summary.append({
        "sample": smp,
        "cells": int(len(sub)),
        "min_genes": float(min_genes),
        "max_genes": float(max_genes),
        "max_counts": float(max_counts),
        "max_mt(%)": float(max_mt),
        "keep_rate(%)": round(keep.mean()*100,1)
    })
qc_df = pd.DataFrame(qc_summary).round(2)
print(qc_df)
qc_df.to_csv(os.path.join(OUTPUT_DIR, "qc_thresholds_summary.csv"), index=False)
adata.obs['qc_keep'] = qc_keep
adata = adata[qc_keep].copy()

# サンプルごとにScrubletでダブレット検出
scrub_results = []

for sample in adata.obs[SAMPLE_COL].unique():
    adata_sample = adata[adata.obs[SAMPLE_COL] == sample].copy()

    counts_matrix = adata_sample.X
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    
    adata_sample.obs["doublet_score"] = doublet_scores
    adata_sample.obs["predicted_doublet"] = predicted_doublets
    adata_sample = adata_sample[adata_sample.obs["predicted_doublet"] == False].copy()
    
    if adata_sample.n_obs > 0:
        scrub_results.append(adata_sample)
        print(f"Sample {sample}: {adata_sample.n_obs} cells remain after doublet removal")
    else:
        print(f"Warning: Sample {sample} has 0 cells after doublet removal. Skipping.")

if scrub_results:
    sample_keys = [adata.obs[SAMPLE_COL].iloc[0] for adata in scrub_results]
    adata = anndata.concat(scrub_results, join="outer", label=SAMPLE_COL, keys=sample_keys)
else:
    raise ValueError("All samples were removed during doublet detection. Please check parameters.")

print("Doublet_Process_finished")

print("n_samples:", adata.obs[SAMPLE_COL].nunique())
print("min cells per sample:", adata.obs[SAMPLE_COL].value_counts().min())

# 生データを "counts" レイヤーに保存
adata.layers["counts"] = adata.X.copy()

# 正規化と対数変換
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 正規化前の生データを raw 属性に保存
adata.raw = adata.copy()

# HVG選択→スケーリング
sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor='seurat_v3', layer='counts')
sc.pl.highly_variable_genes(adata, show=False, save=os.path.join(OUTPUT_DIR, "highly_variable_genes.png"))

# 前処理後のデータを保存
adata.write_h5ad(os.path.join(OUTPUT_DIR, "adata_preprocessed.h5ad"))