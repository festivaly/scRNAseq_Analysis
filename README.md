# scRNA-seq解析パイプラインの構築

## パイプラインの流れ

### 1, アライメント・発現行列の取得
一般に生データはSRAに存在している
sam-toolsを用いてsraファイルをダウンロードする

**注意**

以下のパイプラインで可能なのは10x Genomicsのデータのみであるので、確認する必要がある。

ショートリード長が28bpであるものが基本的には対象である。またライブラリは確認する必要があり、v3以上が好ましい。ARCなどのATAC-seqがついているものも解析自体は可能であるが、少し工夫が必要であり、精度自体もあまりよくない

```bash
SRA_DIR=sra_dir
SRRID=sra_accesion_number
prefetch --output-directory $SRA_DIR --max-size 1000GB $SRRID
```

fastq-dump (fasterq-dump) でfastqファイルに加工する
```bash
SRA_DIR=sra_dir
SRRID=sra_accesion_number
FASTQ_DIR=fastq_dir
fastq-dump \
--split-files \
--gzip \
-O $FASTQ_DIR \
$SRA_DIR/${SRRID}.sra
```

CellRangerを用いて発現行列を取得する

このとき、リファレンスゲノムには十分注意を払う必要がある

そもそもリファレンスゲノムのビルドを行う必要があるが、このとき誤ったGTFファイルを用いてしまうと下流の解析を全てやり直す必要が出てくる

現在の一般的なヒトゲノムはhg38である

かなり古いヒトゲノムではあるもののそれだけ様々なツールの整備が進んでおり使いやすい

**注意**

ここでファイル名の変更が必要である

fastqファイルは2~4つ生成されるため、適切に名前を変更しないとCellRangerが読み込めない

ファイル数とファイルサイズから逐次判断する必要があるため、個別に対応しなければならない
```bash
SRA_DIR=sra_dir
FASTQ_DIR=fastq_dir
CellRanger=cell_ranger_path
source ${CellRanger}/sourceme.bash
CELL_RANGER_DIR=cell_ranger_ouutput_dir
mkdir -p $CELL_RANGER_DIR/${SRRID}
cd $CELL_RANGER_DIR/${SRRID}

${CellRanger}/cellranger count \
--id=$SRRID \
--fastqs=$FATSQ_DIR/$SRRID \
--sample=$SRRID \
--transcriptome=/path/to/reference_genome \
--localcores=16 \
--localmem=64 \ \
--create-bam=true \
```

**注意**

セクション1のみサンプルコマンドです

実際に実行する場合は調整が必要です

### 2, 読み込み

ここまではシェルの処理のみだが、ここからはpythonまたはRによる解析に移ることとなる

pythonの場合は、Scanpyが一般的である

Rの場合はSeuratが、一般的である

**Scanpyの利点**
処理が早いことが特徴である
大規模なデータセットになればなるほど、その差は歴然となっていく
そもそもRに比べ、pythonは圧倒的にはやい (RはC++が裏で動いているからと言われている)
発現行列はH5AD形式である

**Seuratの利点**
図が綺麗に出力される
また、古くからあるためライブラリやコードが豊富にある
図が綺麗なことから出力だけRで行うという人もいたが、最近互いにバージョンが進み、変換が複雑化しあまり見かけない
発現行列はRDS形式である

CellRangerの発現行列をH5ADまたはRDSにして読み込む

### 3, 前処理 (正規化・log化・QC)
ここでは品質の低い細胞・遺伝子の削減を行う

これより下流の解析は全て、正規化・log化を行った状態で行う

生の発現量が必要な場合は、レイヤーadata.layers['counts']などとして残しておくことが推奨される

一般に以下の閾値が知られている

発現遺伝子数が200未満の細胞の除外

発現細胞数が3未満の遺伝子の除外

ミトコンドリア由来の遺伝子が多すぎる細胞(5%以上) の除外

total countsの値が多すぎる細胞(2500以上) の除外

しかし、最近では動的な (分布に応じた) 閾値が用いられることが多い

例えば、ミトコンドリア割合が10%を超えていても、生物学的に意味のある細胞があることがわかってきている　(腫瘍細胞などで特に顕著)

また、最新のQCではリボソーム関連遺伝子やダブレットを除外することが多い

こののち、QCを可視化する

主にバイオリンプロットを用いて極端なサンプルがないかなどを観察する (脳の場合は二峰性になることが多い) 

サンプル除外なども視野に入れる

除外後、Highly-variable gene (HVG) をクリップする

HVGのみで解析を進める場合もあるが、できれば情報をロストしたくない場合などで残す

PCAやUMAPではHVG以外の遺伝子はノイズとなることが知られているが、現在のバージョンではデフォルトでHVGのみを対象として次元を削減してくれる

実行コマンドは以下
```bash
INPUT_DATA=/path/to/adata.h5ad
OUTPUT_DIR=/path/to/output_dir
SAMPLE_COL=sample_column_name
python Preprocessing.py $INPUT_DATA $OUTPUT_DIR $SAMPLE_COL
```

### 4, クラスタリング・バッチ処理

**クラスタリングについて**

次元削減は一般的にPCA (次元削減)→k-NNグラフ作成→Leiden (クラスタリング)→UMAP (次元圧縮)の順で行うことが多い

PCAとUMAPは共に次元削減手法ではあるが、PCAは計算を効率化しノイズを消し、UMAPはマッピングし可視化する手法 (すなわち2~3次元まで圧縮する)であり、意味合いが違うので使い分けた

一方、Leidenはクラスターに分類する手法であり、その他二つとはそもそもの概念が異なる

k-NNグラフは高次元空間の細胞同士をエッジで繋ぐ操作で、LeidenとUMAPの前処理としての意味を持つ

相当に高度な数学的処理が行われているが、その厳密な理解は運用に必要な部分ではない

かつてはクラスリングはLouvainや、次元圧縮t-SNEと呼ばれる手法がメインであったが、現在はUMAPが主流である

デフォルトの状態でもPCA、UMAPは共に出力される

しかし、パラメータが粗いため、状況を見ながら試して行くことが重要である
推奨される方法として、クラスタリング評価に有用ないくつかの指標 (Silhouette, Davies-Bouldin, Calinski-Harabasz) を用いて動的にパラメータを調整し、下流の結果を見ながらさらに微調整を加えて行く方法である。

すなわち、クラスタリング自体のクオリティーが善くても生物学的な妥当性に乏しい場合は (クラスターが粗すぎる、または細かすぎる)、以下を参考にパラメータを調整する (5,6を一旦実行したのち、検討する)


**バッチ処理について**

サンプルや実験手法ごとの偏りを是正していく

bbnkkやharmonyやcombatがpythonでは主流であるが、scVIなどGPUを活用する方法もある

注意しなければならないのは以下の3点

・どのタイミングで行うか

ツールによってPCA後に行うものとPCA前に行うものがある

・データはどのようなものが想定されているか

生データが必要なものや、正規化・log化したもので良いツールなど様々である

・データ自体が書き換わるか

ツールの仕様上、発現量そのものを上書きするものと細胞の座標だけ付記するものに大別される。下流の解析で発現量そのものが変えられることが必ずしも好ましくない場合は注意が必要。なお、前者はPCA前、後者はPCA後に行うものが多い


これらの3点については必ず使用前にドキュメントを参照する必要がある

今回はharmony-pyを用いる

harmon-pyは比較的標準的な方法で、R由来のコードで使用例が多い

また、PCA後にかけるため、発現量そのものが変わらない

scanpyのプラグインとして用いる場合、操作できるパラメータが少ないため、こだわる場合はharmony-py本体を呼び出すことが推奨されるが、ややコードが煩雑になる上、実際に必要とされるケースが少ない

実行コードは以下
```bash
INPUT_DATA=/path/to/adata_preprocessed.h5ad
OUTPUT_DIR=/path/to/output_dir
SAMPLE_COL=sample_column_name
python Clustering.py $INPUT_DATA $OUTPUT_DIR $SAMPLE_COL
```

### 5,可視化

PCAの寄与率グラフとUMAPを主に出力する

主に以下をマッピングする

・サンプル名

サンプルに異常な偏りがないか

・ミトコンドリア割合

極端にミトコンドリア割合が多いクラスターはないか

・クラスター

クラスターが多すぎないか、少なすぎないか

実行コードは以下
```bash
INPUT_DATA=/path/to/adata_clustered.h5ad
OUTPUT_DIR=/path/to/output_dir
SAMPLE_COL=sample_column_name
python Clustering_Visualization.py $INPUT_DATA $OUTPUT_DIR $SAMPLE_COL
```

### 6, アノテーション

アノテーションは様々な方法がある

主に自動アノテーションツールと手動アノテーションがある

自動アノテーションは有名なものにCelltypistやSingleRが存在しており、組織やCell Typeの粒度などから検討する

基本的に自動アノテーションツールのみでアノテーションが完結することは少なく、最終的には手動でアノテーションをつけることが推奨される

自動アノテーションツールは基本的に補助ツールとして用いることが好ましい
特にマーカー遺伝子のマッピングは必須である
データの性質や結果に併せて手法を組み合わせて行く必要がある

例えばT細胞を更に細かく分類したい場合、PBMCを一度クラスタリングするだけでアノテーションを行うことはほぼ不可能なので、T細胞を更にサブクラスタリングしてアノテーションを行うことが多い

Celltypistの実行コードは以下
```bash
INPUT_DATA=/path/to/adata_clustered.h5ad
OUTPUT_DIR=/path/to/output_dir
MODEL_NAME=model_name
python Celltypist_Annotation.py $INPUT_DATA $OUTPUT_DIR $MODEL_NAME
```

またマーカー遺伝子ファイルMarker_List.jsonを編集して以下のコマンドを実行すると、細胞ごとのCelltype Module Scoreがマッピングできる
```bash
INPUT_DATA=/path/to/adata_annotated_celltypist.h5ad
OUTPUT_DIR=/path/to/output_dir
MARKER_FILE=/path/to/marker_list.txt
python Marker_Annotation.py $INPUT_DATA $OUTPUT_DIR $MARKER_FILE
```

更にこのデータに基づき、Annotation.jsonでクラスターごとにCelltypeを割り当て以下のコマンドを実行すると、CelltypeをマッピングしたUMAPが作成できる
```bash
INPUT_DATA=/path/to/adata_marker_annotated.h5ad
OUTPUT_DIR=/path/to/output_dir
SAMPLE_COL=sample_column_name
MARKER_FILE=/path/to/marker_list.txt
ANNOTATION_FILE=/path/to/annotation_list.json
python Annotation_detct.py $INPUT_DATA $OUTPUT_DIR $SAMPLE_COL $MARKER_FILE $ANNOTATION_FILE
```

### 7, 下流解析及び留意事項

下流では様々な解析を行うことができる

主な手法は以下

**DEG解析**
二群の遺伝子発現を比較し、高差次発現遺伝子を同定することができる

**TF activity解析**
転写活性を同定することができ、二群間での比較も可能

**Pseudotime解析**
擬似時間をマッピングすることができる (例えば腫瘍の進展や細胞の分化経路を辿ることができる)

**Cell-Cell Communication (CCC)**
細胞間でのコミュニケーションを推定することができる

**Pathway解析 (及びGSE解析)**
DEGの結果を用いて、EnrichされているPathwayを同定したり、PathwayでEnrichされている遺伝子を同定することができる

**SNV解析**
一般的にはWGSデータが必要であるが、GATKなどのツールを活用することで、somaticSNVやgermlineSNVを同定することが可能である

**CNV解析**
コピー数をマッピングすることで増幅または欠失している染色体や領域を特定することができる

**注意点**
下流解析に応じて必要とされるデータが異なる

例えば、SNV解析ではCellRangerで出力されるBamファイルが要求される

また例えばこれらの解析をアノテーションに使用しなければならないケースもあるため (例えば腫瘍細胞の解析において、非常に遺伝子プロファイルが似ているCelltypeを同定する際にCNV解析の結果を使うことなどがある)、適切なタイミングで解析の順序を挿入したり、除外したりしてアノテーションを進める必要がある

また複雑なケースでは、データのコンディションに併せてアノテーションをかけたデータ同士を統合することなどがある (このとき、scANVIを使う)

## 使い方
**pythonコードの使い方**

以下に実行コマンドを整理して列挙する

/path/to/hogehogeとなっていたり、大文字=hogehogeとなっているものは適宜環境を書き換えることで、実行できる

3) 前処理
```bash
INPUT_DATA=/path/to/adata.h5ad
OUTPUT_DIR=/path/to/output_dir
SAMPLE_COL=sample_column_name
python Preprocessing.py $INPUT_DATA $OUTPUT_DIR $SAMPLE_COL
```

4) クラスタリング・バッチ処理
```bash
INPUT_DATA=/path/to/adata_preprocessed.h5ad
OUTPUT_DIR=/path/to/output_dir
SAMPLE_COL=sample_column_name
python Clustering.py $INPUT_DATA $OUTPUT_DIR $SAMPLE_COL
```

5) 可視化
```bash
INPUT_DATA=/path/to/adata_clustered.h5ad
OUTPUT_DIR=/path/to/output_dir
SAMPLE_COL=sample_column_name
python Clustering_Visualization.py $INPUT_DATA $OUTPUT_DIR $SAMPLE_COL
```

6.1) Celltypist
```bash
INPUT_DATA=/path/to/adata_clustered.h5ad
OUTPUT_DIR=/path/to/output_dir
MODEL_NAME=model_name
python Celltypist_Annotation.py $INPUT_DATA $OUTPUT_DIR $MODEL_NAME
```

6.2) マーカー遺伝子のスコアリング
```bash
INPUT_DATA=/path/to/adata_annotated_celltypist.h5ad
OUTPUT_DIR=/path/to/output_dir
MARKER_FILE=marker_list
python Marker_Annotation.py $INPUT_DATA $OUTPUT_DIR $MARKER_FILE
```

6.3) クラスターアノテーションのマッピング
```bash
INPUT_DATA=/path/to/adata_marker_annotated.h5ad
OUTPUT_DIR=/path/to/output_dir
SAMPLE_COL=sample_column_name
MARKER_FILE=marker_list
ANNOTATION_FILE=annotation_list.json
python Annotation_detct.py $INPUT_DATA $OUTPUT_DIR $SAMPLE_COL $MARKER_FILE $ANNOTATION_FILE
```

**jupyter labでの実行方法**

scRNAseq_pipeline_Run.ipynbで全く同じコードが実行できる

逐次的に結果を確認したい場合に使って欲しい

Jupyterで実行する場合は環境構築ガイド (後述) を参考にjupyter labを起動し、セル2でパスを設定し実行して欲しい

ローカル環境なので当然パスも同じ

GoogleColabでも実行は可能だが細かい仕様が異なるため、scRNAseq_pipeline_Run_For_GoogleColab.ipynbを使用する

パスの設定などは同様にセル2で行える

但し、データが大規模の場合メモリがフローして解析ができない可能性がある

また、ファイルやフォルダを作成してアップロード・指定する必要がある

## 前提知識
**Gitをクローンする方法**

GitはGitアカウントがなくてもクローンすることができる

以下をターミナルで実行し、2.x.xなどと出てくれば準備は完了している

「インストールしますか？」と問われた場合は許可する
```bash
git --version
```

次に、以下を実行するとパスを問われ、指定した場所に本レポジトリがクローンされ、コードが自由に使用できる
```bash
git clone htpss
```

**環境構築ガイド**

本パイプラインはpython環境を必要とする

venvというツールでビルドすることも十分に可能だが、一般に環境が汚染されやすいvenvではなくconda環境が使用されることが多い

conda環境にもAnacondaとMinicondaが存在するが、よりサイズが小さく環境を作りやすいMinicondaが推奨される (Anacondaの方が初心者向きではある)

古くはcondaコマンドによるビルドが好まれたが、更に依存関係 (ライブラリ間の関係) の解決が非常に早いことで知られるmambaコマンドを使用する

Minicondaに代わり、商用利用が可能でありmambaをネイティブとするMiniforgeが最近登場した

理由がなければMinicondaかMiniforgeが推奨される

以下はMiniforgeのビルド方法である

まずターミナルを開き、ビルドスクリプトをダウンロードする
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
```

次にビルドを実行する
```bash
bash Miniforge3-MacOSX-arm64.sh -b -p ~/miniforge3
```

基本的にはオートになるコマンドではあるが、何か聞かれた場合はEnterで大丈夫であるが、conda init (またはmamba init) の許可を求められた場合はnoと答えることをお勧めする

理由として最もクリーンに環境を構築できる方法だからである

自分である程度、自由にpythonを動かせるようになってから、汚い環境の整理に時間を使わなければならなくなる可能性が高い

スクリプト実行後、ダウンロードしたファイルは削除可能である
ターミナル起動後毎回以下のコマンドを実行するとmamba環境が使えるようになる
```bash
. ~/miniforge3/etc/profile.d/mamba.sh
```

mamba環境起動後、以下のコマンドを実行するとパイプライン用の環境を構築・起動できる
```bash
mamba env create -f environment.yml
mamba activate scRNAseq
```

2回目以降は以下のコマンドのみで構わない
```bash
. ~/miniforge3/etc/profile.d/mamba.sh
mamba activate scRNAseq
```

起動した環境にはjupyter labが入っているため、以下のコマンドでjupyter labがブラウザ上で起動する (ローカルでの処理のため、オフライン環境でも動作する)

そのまま、ipynbファイルが実行可能である
```bash
jupyter lab
```

使用終了後はターミナルに戻ってCtrl+Cでシャットダウンできる


