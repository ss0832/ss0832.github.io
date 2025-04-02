---
title:  【PySCF】Møller-Plesset法を用いた電子相関計算の実装
published: 2025-04-02
description: "PySCFを用いたMøller-Plesset法の計算の使用方法と各種オプションの解説"
tags: [pyscf]
category: Computational Chemistry
draft: false
---
最終更新：2025-04-02

## 概要

Møller-Plesset摂動論（MPn法）は、Hartree-Fock法で得られたハミルトニアンと厳密なハミルトニアンの差を摂動として取り扱い、電子相関効果を段階的に取り込む手法である。PySCF ver.2.8.0では、2次の摂動計算のみが公式に実装されている。本記事では、PySCFにおけるMPn法の基本的な使い方と様々なオプションについて解説する。

## 使用したPySCFのバージョン
```
2.8.0
```

## 基本的な使い方

### モジュールのインポート

```python
from pyscf import gto, scf, mp
```

### MP2計算の基本例

```python
# 分子の定義
mol = gto.M(
    atom = '''
    O  0.000000  0.000000  0.000000
    H  0.000000  0.757160  0.586260
    H  0.000000 -0.757160  0.586260
    ''',
    basis = 'cc-pvdz'
)

# Hartree-Fock計算
mf = scf.RHF(mol).run()

# MP2計算の実行
mp2 = mp.MP2(mf)
e_corr, t2 = mp2.kernel()
mp2_total_energy = mf.e_tot + e_corr

print(f"HF energy: {mf.e_tot:.8f}")
print(f"MP2 correlation energy: {e_corr:.8f}")
print(f"MP2 total energy: {mp2_total_energy:.8f}")
```


## 様々な参照状態に基づくMPn計算

PySCFでは、閉殻系（Restricted、R）、開殻系（Unrestricted、U）のHartree-Fock波動関数に基づくMPn計算が可能である。

### 非制限MPn計算（UMPn）

```python
# 非制限HF計算
umf = scf.UHF(mol).run()

# UMP2計算
ump2 = mp.UMP2(umf)
ump2_corr = ump2.kernel()[0]
ump2_total = umf.e_tot + ump2_corr

print(f"UHF energy: {umf.e_tot:.8f}")
print(f"UMP2 correlation energy: {ump2_corr:.8f}")
print(f"UMP2 total energy: {ump2_total:.8f}")
```


## 主要なオプション

### フローズンコア近似

内殻電子を相関計算から除外することで計算負荷を軽減できる。

```python
# フローズンコア近似を使用したMP2計算
mp2_fc = mp.MP2(mf)
mp2_fc.frozen = 1  # 酸素の1s軌道をフリーズ
# または自動設定
# mp2_fc.frozen = 'auto'
mp2_fc_corr = mp2_fc.kernel()[0]
mp2_fc_total = mf.e_tot + mp2_fc_corr

print(f"MP2 (frozen core) correlation energy: {mp2_fc_corr:.8f}")
print(f"MP2 (frozen core) total energy: {mp2_fc_total:.8f}")
```

### 直接法（ディスク使用量の削減）

中間データをディスクに保存せず、必要に応じて再計算する方法。

```python
# 直接MP2法
mp2_direct = mp.MP2(mf)
mp2_direct.direct = True
mp2_direct_corr = mp2_direct.kernel()[0]

print(f"Direct MP2 correlation energy: {mp2_direct_corr:.8f}")
```

### メモリ使用量の制限

```python
# メモリ使用量の制限
mp2_mem = mp.MP2(mf)
mp2_mem.max_memory = 2000  # 単位: MB
mp2_mem_corr = mp2_mem.kernel()[0]

print(f"MP2 correlation energy (memory limited): {mp2_mem_corr:.8f}")
```




## 解析的勾配計算

MP2では解析的勾配が実装されており、構造最適化に利用できる。

```python
# MP2勾配計算
mp2_grad = mp.MP2(mf)
mp2_grad.kernel()
gradients = mp2_grad.nuc_grad_method().kernel()

print("MP2 gradients:")
print(gradients)
```

## ハンズオン：様々なMPn計算の例

以下は、水分子に対する様々なMPn計算を実行し、結果を解析する包括的な例である。

```python
import numpy as np
from pyscf import gto, scf, mp

# 水分子の定義
mol = gto.M(
    atom = '''
    O  0.000000  0.000000  0.000000
    H  0.000000  0.757160  0.586260
    H  0.000000 -0.757160  0.586260
    ''',
    basis = 'cc-pvdz',
    verbose = 3
)

# RHF計算
mf = scf.RHF(mol)
mf.kernel()
print(f"RHF energy: {mf.e_tot:.8f} a.u.")

# MP2計算
mp2 = mp.MP2(mf)
mp2_corr = mp2.kernel()[0]
mp2_total = mf.e_tot + mp2_corr
print(f"\n1. 基本MP2計算")
print(f"MP2 correlation energy: {mp2_corr:.8f} a.u.")
print(f"MP2 total energy: {mp2_total:.8f} a.u.")


# フローズンコアMP2計算
mp2_fc = mp.MP2(mf)
mp2_fc.frozen = 1  # 酸素の1s軌道をフリーズ
mp2_fc_corr = mp2_fc.kernel()[0]
mp2_fc_total = mf.e_tot + mp2_fc_corr
print(f"\n2. フローズンコアMP2計算")
print(f"MP2 (frozen core) correlation energy: {mp2_fc_corr:.8f} a.u.")
print(f"MP2 (frozen core) total energy: {mp2_fc_total:.8f} a.u.")

# 直接MP2計算
mp2_direct = mp.MP2(mf)
mp2_direct.direct = True
mp2_direct_corr = mp2_direct.kernel()[0]
mp2_direct_total = mf.e_tot + mp2_direct_corr
print(f"\n3. 直接MP2計算")
print(f"Direct MP2 correlation energy: {mp2_direct_corr:.8f} a.u.")
print(f"Direct MP2 total energy: {mp2_direct_total:.8f} a.u.")


# MP2勾配計算
mp2_grad = mp.MP2(mf)
mp2_grad.kernel()
gradients = mp2_grad.nuc_grad_method().kernel()
print(f"\n4. MP2勾配計算")
print("MP2 gradients (a.u.):")


# スピン開殻系（三重項酸素分子）のUMP2計算
o2_mol = gto.M(
    atom = '''
    O  0.000000  0.000000  0.000000
    O  0.000000  0.000000  1.207
    ''',
    basis = 'cc-pvdz',
    spin = 2,  # 三重項状態
    verbose = 3
)

# UHF計算
umf = scf.UHF(o2_mol)
umf.kernel()

# UMP2計算
ump2 = mp.UMP2(umf)
ump2_corr = ump2.kernel()[0]
ump2_total = umf.e_tot + ump2_corr
print(f"\n5. 三重項酸素分子のUMP2計算")
print(f"UHF energy: {umf.e_tot:.8f} a.u.")
print(f"UMP2 correlation energy: {ump2_corr:.8f} a.u.")
print(f"UMP2 total energy: {ump2_total:.8f} a.u.")
```


## 参考サイト

- [PySCF公式ドキュメント: MP2モジュール](https://pyscf.org/user/mp.html)
- [PySCF GitHub: MPnコード](https://github.com/pyscf/pyscf/tree/master/pyscf/mp)
- [PySCF API Documentation](https://pyscf.org/pyscf_api_docs/pyscf.mp.html)
- [PySCF Examples Repository](https://github.com/pyscf/examples)