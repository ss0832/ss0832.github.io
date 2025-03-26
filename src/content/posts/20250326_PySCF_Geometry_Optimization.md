---
title: 【PySCF】PySCFを用いて構造最適化する方法
published: 2025-03-26
description: ""
tags: [Pyscf]
category: Computational Chemistry
draft: false
---
最終更新：2025-03-26

## 概要

PySCFはPythonで書かれたオープンソースの第一原理電子構造計算ソフトウェアである。
本稿では、最新版（PySCF 2.3.0）に準拠した構造最適化の手順を、Hartree-Fock (HF) およびDFTの両方の例とともに簡潔に示す。
デフォルトの設定では、構造最適化に使うsolverはgeomeTRIC(https://github.com/leeping/geomeTRIC)である。

## 1. Hartree-Fock法による構造最適化

```python 
# 1. 必要なモジュールをインポート
from pyscf import gto, scf
from pyscf.geomopt.geometric_solver import optimize
# geometric_solverからoptimize関数を呼び出す

# 2. 分子情報を定義
mol = gto.Mole()
mol.atom = '''
O  0.000000   0.000000   0.000000
H  0.757000   0.586000   0.000000
H -0.757000   0.586000   0.000000
'''
# たとえば水分子の座標を定義
mol.charge = 0             # 電荷（中性なら0）
mol.spin = 0               # スピン多重度1 (singlet) のため spin=0
mol.basis = '6-31G*'       # 基底関数を指定
mol.build()

# 3. SCF計算（Restricted Hartree-Fock）の設定
mf = scf.RHF(mol)
mf.max_cycle = 100         # SCF計算の反復回数を設定
mf.conv_tol = 1e-9         # SCF計算のエネルギー収束閾値
mf.verbose = 4             # 計算ログの出力レベル(大きいほど詳細)
mf.diis = True             # 収束を加速するDIISアルゴリズムを使用

# 4. 構造最適化を実行
# optimize関数にSCFインスタンスを渡して呼び出すと、分子の座標を最適化したうえで
# 収束した最終座標を返す
optimized_coords = optimize(mf)
print("Optimized Coordinates (Angstrom) with HF:", optimized_coords)
```

行ごとの意味:  
1) pyscf.geomopt.geometric_solverからoptimize関数をインポートしている。  
2) Moleオブジェクトを生成し、原子座標・基底関数などを定義している。  
3) RHFを用いてSCF計算のセットアップを行い、反復回数や収束条件などを指定している。  
4) optimize関数を呼び出して構造最適化を実行し、最終的な座標を取得している。  

### SCFのアルゴリズムを変更する場合

- mf.diis = False とするとDIISを無効化する。  

## 2. DFTによる構造最適化

```python 
from pyscf import gto, dft
from pyscf.geomopt.geometric_solver import optimize

# 1. 分子情報を定義
mol = gto.Mole()
mol.atom = '''
O  0.000000   0.000000   0.000000
H  0.757000   0.586000   0.000000
H -0.757000   0.586000   0.000000
'''
# 例: トリプレット系の場合 spin=2 (n-1) にする点に注意
mol.charge = 0
mol.spin = 0
mol.basis = '6-31G*'
mol.build()

# 2. DFTの設定 (Restricted Kohn-Sham)
mf_dft = dft.RKS(mol)
mf_dft.xc = 'B3LYP'       # XC汎関数を指定 (例: B3LYP)
mf_dft.grids.level = 3    # DFT計算用グリッドの細かさ(0～9程度)
mf_dft.max_cycle = 100    # SCF計算の反復回数
mf_dft.conv_tol = 1e-8    # SCFエネルギー収束閾値
mf_dft.verbose = 4

# 3. 構造最適化の実行
opt_coords_dft = optimize(mf_dft)
print("Optimized Coordinates (Angstrom) with DFT:", opt_coords_dft)
```

行ごとの意味:  
1) SCFの場合と同様、分子情報を定義する。  
2) `dft.RKS`を立ち上げ、汎関数(xc)や計算精度、グリッド細分化レベルなどを指定する。  
3) HFと同様にoptimize関数を呼び出すことで座標最適化を実行し、最終的な原子座標を取得する。  

### DFTのグリッド設定

- `mf_dft.grids.level`の数字を上げるほど細かいグリッドを使い、計算精度が向上する一方で計算コストも増大する。  
- `mf_dft.grids.atom_grid = (75,302)` のようにより詳細なグリッドを個別に設定することも可能。  

mf_dft.grids.atom_gridでタプル (R, A) の形でグリッドの細かさを指定する場合、最初の要素 R は“各原子に対する半径方向の格子点数（radial grid points）”を、2番目の要素 A は“角度方向の格子点数（angular grid points）”を表す。たとえば (75, 302) と指定した場合、原子ごとに半径方向に 75 個の積分点、角度方向に 302 個の積分点を設定することになる。

## 注意事項

1. PySCF バージョン2.2.0 以降では、構造最適化関数は「pyscf.geomopt.geometric_solver」から`optimize`を呼び出す。  
2. スピン多重度がn重項の場合は、コード上で「mol.spin = n - 1」と設定しなければならない。
(2.の補足)　ラジカル等、開核系の計算でスピン多重度が2以上になる場合や、一重項ラジカル等を取り扱う場合は、`scf.RHF`から`scf.UHF`（DFTの場合は`scf.RKS`から`scf.UKS`）に書き換えたうえで構造最適化を行う。 
3. SCF計算での収束がうまくいかない場合は、反復回数 (`max_cycle`) やアルゴリズム (`DIIS`, `level_shift`など) を適宜調整する。  
4. DFT計算では、グリッドの細かさ (`grids.level` など) を調整することで精度と計算時間のトレードオフを制御できる。
(4.の補足)　M06等のmeta-GGA汎関数はHessianの計算時にグリッドが十分細かくとらないと、不正確な結果が出力されるので注意。（少なくともデフォルトの`mf_dft.grids.level`よりも細かくとる必要がある）

