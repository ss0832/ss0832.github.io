---
title:  【Psi4】RCCSD法とUCCSD法による高精度電子相関計算
published: 2025-04-07
description: "Psi4を用いた制限・非制限CCSD法およびCCSD(T)法の実装方法と各種オプションの解説"
tags: [psi4]
category: Computational Chemistry
draft: false
---
最終更新：2025-04-07

## 概要

結合クラスター法（Coupled Cluster Method）は、波動関数理論に基づく高精度電子相関計算手法である。特に単一・二電子励起を含むCCSD法と摂動論的三電子励起を加えたCCSD(T)法は、化学的精度で電子相関を考慮できる優れた手法として知られている。本記事では、閉殻系に適用する制限CCSD（RCCSD）と開殻系に適用する非制限CCSD（UCCSD）、およびそれぞれの摂動論的三電子励起を含む方法（RCCSD(T)およびUCCSD(T)）について、Psi4での実装方法を解説する。

## 基本的な使い方

### モジュールのインポート

```python
import psi4
import numpy as np
```

### RCCSD計算の基本例（閉殻系）

```python
# 水分子の定義（閉殻系）
h2o = psi4.geometry("""
0 1
O
H 1 0.96
H 1 0.96 2 104.5
""")

# RHF参照のCCSD計算（RCCSD）
psi4.set_options({'basis': 'cc-pvdz', 'reference': 'rhf'})
rccsd_energy = psi4.energy('ccsd')

print(f"RCCSD total energy: {rccsd_energy}")
```

### RCCSD(T)計算（閉殻系）

```python
# RCCSD(T)計算
rccsd_t_energy = psi4.energy('ccsd(t)')
print(f"RCCSD(T) total energy: {rccsd_t_energy}")
```

### UCCSD計算の基本例（開殻系）

```python
# メチルラジカルの定義（開殻系）
methyl = psi4.geometry("""
0 2
C  0.0000  0.0000  0.0000
H  0.0000  1.0789  0.0000
H  0.9343 -0.5394  0.0000
H -0.9343 -0.5394  0.0000
""")

# UHF参照のCCSD計算（UCCSD）
psi4.set_options({'basis': 'cc-pvdz', 'reference': 'uhf'})
uccsd_energy = psi4.energy('ccsd', molecule=methyl)

print(f"UCCSD total energy: {uccsd_energy}")
```

### UCCSD(T)計算（開殻系）

```python
# UCCSD(T)計算
uccsd_t_energy = psi4.energy('ccsd(t)', molecule=methyl)
print(f"UCCSD(T) total energy: {uccsd_t_energy}")
```

## 制限法と非制限法の選択

### 制限CCSD（RCCSD）

閉殻系の分子（一重項状態など）に適している。スピン対称性が保たれる。

```python
psi4.set_options({'reference': 'rhf'})
```

### 非制限CCSD（UCCSD）

ラジカル種や結合解離など、開殻系や強い静的相関がある系に適している。

```python
psi4.set_options({'reference': 'uhf'})
```

## 主要なオプション

### 凍結軌道の設定

内殻電子を相関計算から除外することで計算量を削減する。

```python
# 内殻軌道を自動的に凍結（デフォルト）
psi4.set_options({'freeze_core': 'true'})

# すべての軌道を相関計算に含める
psi4.set_options({'freeze_core': 'false'})

# 凍結軌道を明示的に指定
psi4.set_options({
    'freeze_core': 'false',
    'num_frozen_docc': [2]  # 2つの軌道を凍結
})
```



## ハンズオン：RCCSDとUCCSD計算

以下に、閉殻系（水分子）と開殻系（酸素原子）に対する様々なCCSD計算の例を示す。

```python
import psi4
import numpy as np

# 出力レベルの設定
psi4.core.set_output_file('ccsd_example.out', False)
psi4.set_memory('4 GB')

# 閉殻分子の例：水分子
h2o = psi4.geometry("""
0 1
O
H 1 0.96
H 1 0.96 2 104.5
""")

# 開殻分子の例：三重項酸素分子
o2 = psi4.geometry("""
0 3
O
O 1 1.208
""")

print("\n===== 閉殻系の計算（水分子） =====")

# 基底関数設定
psi4.set_options({'basis': 'aug-cc-pvdz'})

# 1. RCCSD計算（凍結軌道あり）
psi4.set_options({
    'reference': 'rhf',
    'freeze_core': 'true'
})
rccsd_e = psi4.energy('ccsd', molecule=h2o)
rccsd_corr = psi4.core.variable("CCSD CORRELATION ENERGY")
print("\n1. RCCSD計算（凍結軌道あり）")
print(f"RCCSD Total Energy: {rccsd_e}")
print(f"RCCSD Correlation Energy: {rccsd_corr}")
psi4.core.clean()
# 2. RCCSD計算（凍結軌道なし）
psi4.set_options({
    'reference': 'rhf',
    'freeze_core': 'false'
})
rccsd_full_e = psi4.energy('ccsd', molecule=h2o)
rccsd_full_corr = psi4.core.variable("CCSD CORRELATION ENERGY")
print("\n2. RCCSD計算（凍結軌道なし）")
print(f"RCCSD Total Energy: {rccsd_full_e}")
print(f"RCCSD Correlation Energy: {rccsd_full_corr}")

# 3. RCCSD(T)計算
rccsd_t_e = psi4.energy('ccsd(t)', molecule=h2o)
rccsd_t_corr = psi4.core.variable("CCSD(T) CORRELATION ENERGY")
t_contrib = psi4.core.variable("(T) CORRECTION ENERGY")
print("\n3. RCCSD(T)計算")
print(f"RCCSD(T) Total Energy: {rccsd_t_e}")
print(f"RCCSD(T) Correlation Energy: {rccsd_t_corr}")
print(f"(T) Contribution: {t_contrib}")


# 4. 厳密な収束条件でのRCCSD計算
psi4.set_options({
    'reference': 'rhf',
    'freeze_core': 'true',
    'e_convergence': 1e-9,
    'r_convergence': 1e-8,
    'maxiter': 150
})
rccsd_tight_e = psi4.energy('ccsd', molecule=h2o)
print("\n5. 厳密な収束条件でのRCCSD計算")
print(f"RCCSD (tight) Total Energy: {rccsd_tight_e}")
psi4.core.clean()
print("\n\n===== 開殻系の計算（三重項酸素分子） =====")

# 5. UCCSD計算
psi4.set_options({
    'reference': 'uhf',
    'freeze_core': 'true',
    'e_convergence': 1e-8,
    'r_convergence': 1e-7
})
uccsd_e = psi4.energy('ccsd', molecule=o2)
uccsd_corr = psi4.core.variable("CCSD CORRELATION ENERGY")
print("\n6. UCCSD計算")
print(f"UCCSD Total Energy: {uccsd_e}")
print(f"UCCSD Correlation Energy: {uccsd_corr}")
# 6. UCCSD(T)計算
uccsd_t_e = psi4.energy('ccsd(t)', molecule=o2)
uccsd_t_corr = psi4.core.variable("CCSD(T) CORRELATION ENERGY")
ut_contrib = psi4.core.variable("(T) CORRECTION ENERGY")
print("\n7. UCCSD(T)計算")
print(f"UCCSD(T) Total Energy: {uccsd_t_e}")
print(f"UCCSD(T) Correlation Energy: {uccsd_t_corr}")
print(f"(T) Contribution: {ut_contrib}")

```

## 参考サイト

- [Psi4マニュアル: CC計算](https://psicode.org/psi4manual/master/methods.html#cc)
- [GitHub: Psi4プロジェクト](https://github.com/psi4/psi4)
