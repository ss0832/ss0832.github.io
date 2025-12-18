---
title: 【計算化学】MultiOptPy-v1.20.3の機能追加内容
published: 2025-12-18
tags: [MultiOptPy, python]
category: Computational Chemistry
draft: false
---
最終更新：2025-12-18


### 1. モード追跡型 RS-I-RFO (MF-RS-I-RFO)


概要

遷移状態探索（Saddle Point Optimization）において、特定のへシアンのモードを追跡しながら最適化を進めるアルゴリズムです。従来のRS-I-RFO法をベースに、ステップ間でのモードの入れ替わりを防ぐための追跡ロジックが組み込まれています。

実装詳細モード追跡戦略: 

ModeFollowing クラスを使用し、先行ステップのへシアン固有ベクトルとのオーバーラップを計算します。

質量重み付けオーバーラップ (MWO): 原子リストが提供された場合、デカルト座標ではなく質量重み付け座標系での投影を行い、物理的に妥当な追跡を実現します。

適応型更新 (EMA): 指数移動平均（EMA）を用いて参照ベクトルを動的に更新することで、モードの回転に追従します。

勾配バイアス: 現在の力（勾配）の方向へのオーバーラップをスコアに加算し、エネルギー的に重要なモードの選択を優先できます。

コマンドライン引数 (-opt)-opt 引数の文字列内でコロン : を使用して詳細なパラメータを指定します。

基本形式: -opt mf_rsirfo:<target_index>:<ema_val>:<grad_val>

パラメータ:

target_index: 追跡対象とする固有値のインデックス（0から開始）。

ema<val>: 更新レート（0.0〜1.0）。1.0で完全置換（適応型）、0.0で固定。

grad<val>: 勾配方向への重み付け。


使用例:
```bash
# 質量荷重モード1を追跡し、適応更新レート0.5、勾配バイアス0.3で遷移状態構造最適化を実行

python optmain.py input.xyz -opt mwmf_rsirfo_fsb:1:ema0.5:grad0.3 -fc 5 -order 1 -freq -tcc
```


### 2. 拘束条件付き RS-I-RFO (C-RS-I-RFO)

概要

幾何学的な拘束条件（結合距離、角度など）を課した状態で、RS-I-RFOによる鞍点探索を行う手法です。

実装詳細部分

空間投影: 特異値分解（SVD）を用いて、課された拘束条件と直交する「ヌル空間基底」を構築します。

制約付き最適化: 全空間の勾配とへシアンをこの部分空間に投影し、その中でRFOステップを計算します。これにより、拘束を維持したまま鞍点への移動ベクトルを算出します。

SHAKE様補正: 数値誤差による拘束のズレを補正するための幾何学的調整ロジックが含まれています。

コマンドライン引数 (-opt, -pc) -opt で crsirfo 系のメソッドを指定し、-pc で拘束条件を定義します。

使用例:
```Bash
# 原子1と2の距離、および原子2-3-4の角度を固定したまま構造最適化
python optmain.py input.xyz -opt crsirfo_block_fsb -pc bond 1,2 bend 2,3,4

```
3. 複数PES情報の統合とモデル関数最適化 (BITSS等)


概要


複数の電子状態（PES）のエネルギーや勾配を組み合わせた「有効ポテンシャル」上で構造最適化を行う機能です。

特にv1.20.3では、Binary-Image Transition State Search (BITSS) 法が実装されています。

実装詳細

独立した計算機インスタンス: ModelFunctionHandler が State1 と State2 用に独立したディレクトリと計算機を生成し、状態間の干渉を防ぎます。

BITSSの実装:

6N次元拡張: 2つの構造（Image 1 & 2）を連結し、$6N \times 6N$ のへシアンおよび $2N$ 原子系として扱います。

二階微分 (Hessian): 距離およびエネルギー等価拘束項の二階微分を数学的に導出し、calc_hess メソッドとして実装しています。

コマンドライン引数 (-mf) -mf 引数を使用して、モデル関数の種類と付随するパラメータ（BITSSの場合は参照構造のファイルパス）を指定します。

使用例:
```Bash
# BITSS法を用いて、target.xyzを目標構造として二点間をつなぐ経路の探索
python optmain.py start.xyz -mf bitss target.xyz

# Seam モデル関数を使った最小エネルギー交差点(Minimum energy seam of crossing)（電荷0, 多重度1と3の状態間）
python optmain.py input.xyz -mf seam 0 3 -elec 0 -spin 3
```

備忘録：内部的なデータ構造の変更

OptimizationState: BITSSモード時には element_list や geometry が自動的に2倍（2N）へ拡張されます。

Hessian統合: Model_hess (PES由来) と bias_hessian (拘束・ペナルティ項由来) が個別に管理され、オプティマイザーに渡される直前に合算されます。

### 参考
- _J. Chem. Phys._ 157, 124107 (2022)
- _J. Am. Chem. Soc._ 2015, 137, 10, 3433–3445

