---
title: 【Psi4】電子状態計算の手法の指定方法
published: 2025-04-01
description: "量子化学計算ソフトウェアPsi4を用いた電子状態計算の指定方法の解説"
tags: [Psi4]
category: Computational Chemistry
draft: false
---
最終更新：2025-04-01

## 概要

Psi4は、高性能な量子化学計算ソフトウェアであり、Pythonインターフェースを持つオープンソースのプログラムである。本記事では、Psi4を用いた様々な電子状態計算法の実行方法について解説する。Hartree-Fock法（HF）、密度汎関数理論（DFT）、MP2、結合クラスター法（CC）などの異なる理論レベルでの計算方法を紹介する。


## 使用したバージョン
```
1.9.1
```

## Psi4を用いた電子状態計算の方法

全般的に、`(計算手法)/(基底関数)`で計算手法の指定が可能である。

### Hartree-Fock法（HF法）

HF法は電子状態計算の基礎となる手法であり、多電子系を解くための基本的な近似法である。

```python
import psi4
import numpy as np

# 出力ファイルとメモリ設定
psi4.set_output_file('hf_calc.out')
psi4.set_memory('2 GB')

# 分子の定義
h2o = psi4.geometry("""
    O  0.000000  0.000000  0.000000
    H  0.757000  0.586000  0.000000
    H -0.757000  0.586000  0.000000
""")

# HF計算の実行（基底関数はcc-pVDZ）
energy_hf = psi4.energy('scf/cc-pvdz')  # scfはHF法と同義
print(f'HF/cc-pVDZ Energy: {energy_hf:.10f} Hartree')

# 異なる基底関数での計算
energy_hf_631g = psi4.energy('scf/6-31g')
print(f'HF/6-31G Energy: {energy_hf_631g:.10f} Hartree')

# 計算結果の詳細出力
psi4.core.print_variables()
```

### 密度汎関数理論（DFT）

DFTは電子相関を効率的に取り入れる手法であり、計算コストと精度のバランスが良い。

```python
# B3LYP汎関数での計算
energy_b3lyp = psi4.energy('b3lyp/cc-pvdz')
print(f'B3LYP/cc-pVDZ Energy: {energy_b3lyp:.10f} Hartree')

# 分散力補正を加えたwB97X-D汎関数
energy_wb97xd = psi4.energy('wb97x-d/cc-pvdz')
print(f'wB97X-D/cc-pVDZ Energy: {energy_wb97xd:.10f} Hartree')

# 異なる分子での計算（ベンゼン）
benzene = psi4.geometry("""
    0 1
    C        0.000000    1.396792    0.000000
    C        1.209612    0.698396    0.000000
    C        1.209612   -0.698396    0.000000
    C        0.000000   -1.396792    0.000000
    C       -1.209612   -0.698396    0.000000
    C       -1.209612    0.698396    0.000000
    H        0.000000    2.484216    0.000000
    H        2.151390    1.242108    0.000000
    H        2.151390   -1.242108    0.000000
    H        0.000000   -2.484216    0.000000
    H       -2.151390   -1.242108    0.000000
    H       -2.151390    1.242108    0.000000
""")

# DFT計算と解析
psi4.set_options({'PRINT': 2})  # 詳細な出力を有効化
energy_benzene = psi4.energy('b3lyp/def2-svp')
print(f'ベンゼンのB3LYP/def2-SVP Energy: {energy_benzene:.10f} Hartree')
```

## 電子相関を考慮した高度な計算法

### MP2（二次のMøller-Plesset摂動論）

MP2法はHF法に電子相関を摂動論的に加えた方法である。

```python
# MP2計算
energy_mp2 = psi4.energy('mp2/cc-pvdz')
print(f'MP2/cc-pVDZ Energy: {energy_mp2:.10f} Hartree')

# MP2とHFのエネルギー差（相関エネルギー）
variables = psi4.core.variables()
correlation_energy = variables["MP2 CORRELATION ENERGY"]
print(f'MP2 相関エネルギー: {correlation_energy:.10f} Hartree')

# RI近似を使用した高速MP2計算
psi4.set_options({'MP2_TYPE': 'DF'})  # 密度fitting近似
energy_df_mp2 = psi4.energy('mp2/cc-pvdz')
print(f'DF-MP2/cc-pVDZ Energy: {energy_df_mp2:.10f} Hartree')
```

### 結合クラスター法（CC）

CCSDは高精度な電子相関法である。

```python
# 小さな分子でCCSD計算（計算コストが高い）
h2 = psi4.geometry("""
    0 1
    H
    H 1 0.74
""")

# CCSD計算
energy_ccsd = psi4.energy('ccsd/cc-pvdz', molecule=h2)
print(f'CCSD/cc-pVDZ Energy: {energy_ccsd:.10f} Hartree')

# CCSD(T)計算（"金の標準"とも呼ばれる高精度計算）
energy_ccsdt = psi4.energy('ccsd(t)/cc-pvdz', molecule=h2)
print(f'CCSD(T)/cc-pVDZ Energy: {energy_ccsdt:.10f} Hartree')
```

### 開殻系の計算

開殻系では非制限法を用いる。計算手法の前に「U」か「u」をつけるだけで使用可能である。例えばHartree-Forck法では`uhf`と入力すればよい。

```python
# メチルラジカル（開殻系）
methyl_radical = psi4.geometry("""
    0 2  # 電荷0、二重項（スピン多重度2）
    C  0.000000  0.000000  0.000000
    H  0.000000  1.078000  0.000000
    H  0.933000 -0.539000  0.000000
    H -0.933000 -0.539000  0.000000
""")

# UHF計算（Unrestricted HF）
psi4.set_options({'REFERENCE': 'UHF'})
energy_uhf, wfn = psi4.energy('scf/cc-pvdz', molecule=methyl_radical, return_wfn=True)
print(f'UHF/cc-pVDZ Energy: {energy_uhf:.10f} Hartree')


```


### 参考
- Psi4公式サイト(https://psicode.org/)

