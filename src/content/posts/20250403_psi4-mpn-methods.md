---
title:  【Psi4】MPn法による相関エネルギー計算の使用方法
published: 2025-04-03
description: "Psi4を用いたMP2-MP4摂動論計算の使用方法の解説"
tags: [psi4]
category: Computational Chemistry
draft: false
---
最終更新：2025-04-03

## 概要

Møller-Plesset摂動論（MPn法）は、Hartree-Fock法で得られた波動関数を基に電子相関エネルギーを計算する手法である。Psi4では、MP2からMP4までの計算が可能で、分子系の高精度エネルギー計算に広く用いられている。本記事では、Psi4におけるMPn法の実装方法と主要なオプションについて解説する。

## 使用したバージョン
```
1.9.1
```

## 基本的な使い方

### モジュールのインポート

```python
import psi4
import numpy as np
```

### MP2計算の基本例

```python
# 分子の定義
mol = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

# 計算実行
psi4.set_options({'basis': 'cc-pvdz'})
mp2_energy = psi4.energy('mp2')

print(f"MP2 total energy: {mp2_energy}")

```

### MP3、MP4の実行

```python
# MP3計算
mp3_energy = psi4.energy('mp3')
print(f"MP3 total energy: {mp3_energy}")

# MP4計算
mp4_energy = psi4.energy('mp4')
print(f"MP4 total energy: {mp4_energy}")
```



## 主要なオプション

### 凍結軌道の設定

内殻電子を電子相関の計算から除外することで計算量を削減する。

```python
# デフォルト（元素ごとに自動決定）
psi4.set_options({'freeze_core': 'true'})

# 凍結軌道を明示的に指定
psi4.set_options({
    'freeze_core': 'false',
    'num_frozen_docc': [2]  # 2つの軌道を凍結
})
```

### ダイレクト計算

ディスク使用量を削減し、メモリ上で計算を完結させる。

```python
psi4.set_options({'mp2_type': 'direct'})
```



## ハンズオン：MPn計算例

以下に、水分子に対する様々なMPn計算の例を示す。

```python
import psi4
import numpy as np

# 出力レベルの設定
psi4.core.set_output_file('mpn_example.out', False)
psi4.set_options({'basis': 'aug-cc-pvdz'})
psi4.set_memory('4 GB')
# 水分子の定義
mol = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

# SCF計算
psi4.set_options({'reference': 'rhf'})
scf_e = psi4.energy('scf')

print(f"\nSCF Energy: {scf_e}")

# MP2計算（凍結軌道あり、デフォルト）
psi4.set_options({'freeze_core': 'true'})
mp2_e = psi4.energy('mp2')
mp2_corr = psi4.core.variable("MP2 CORRELATION ENERGY")
print("\n1. MP2計算（凍結軌道あり）")
print(f"MP2 Total Energy: {mp2_e}")
print(f"MP2 Correlation Energy: {mp2_corr}")
mp2_grad = psi4.gradient('mp2')
print("MP2 Gradient :\n", np.array(mp2_grad))
mp2_hess = psi4.hessian('mp2')
print("MP2 Hessian :\n", np.array(mp2_hess))

# MP2計算（凍結軌道なし）
psi4.set_options({'freeze_core': 'false'})
mp2_full_e = psi4.energy('mp2')
mp2_full_corr = psi4.core.variable("MP2 CORRELATION ENERGY")
print("\n2. MP2計算（凍結軌道なし）")
print(f"MP2 Total Energy: {mp2_full_e}")
print(f"MP2 Correlation Energy: {mp2_full_corr}")
mp2_grad = psi4.gradient('mp2')
print("MP2 Gradient :\n", np.array(mp2_grad))
mp2_hess = psi4.hessian('mp2')

# MP3計算
mp3_e = psi4.energy('mp3')
mp3_corr = psi4.core.variable("MP3 CORRELATION ENERGY")
print("\n3. MP3計算")
print(f"MP3 Total Energy: {mp3_e}")
print(f"MP3 Correlation Energy: {mp3_corr}")
mp3_grad = psi4.gradient('mp3')
print("MP3 Gradient :\n", np.array(mp3_grad))
mp3_hess = psi4.hessian('mp3')
print("MP3 Hessian :\n", np.array(mp3_hess))

# MP4計算
mp4_e = psi4.energy('mp4')
mp4_corr = psi4.core.variable("MP4 CORRELATION ENERGY")
mp4_sdq = psi4.core.variable("MP4(SDQ) CORRELATION ENERGY")
print("\n4. MP4計算")
print(f"MP4 Total Energy: {mp4_e}")
print(f"MP4 Correlation Energy: {mp4_corr}")
print(f"MP4(SDQ) Correlation Energy: {mp4_sdq}")
mp4_grad = psi4.gradient('mp4')
print("MP4 Gradient :\n", np.array(mp4_grad))
mp4_hess = psi4.hessian('mp4')
print("MP4 Hessian :\n", np.array(mp4_hess))

# メソッド比較
print("\n5. 水分子のエネルギー比較")
methods = ["SCF", "MP2", "MP3", "MP4"]
energies = [scf_e, mp2_e, mp3_e, mp4_e]
for method, energy in zip(methods, energies):
    print(f"{method}: {energy}")
```


## 参考サイト

- [Psi4マニュアル: MPn計算](https://psicode.org/psi4manual/master/methods.html#mp2)
- [Psi4マニュアル: 入門チュートリアル](https://psicode.org/psi4manual/master/tutorial.html)
- [GitHub: Psi4プロジェクト](https://github.com/psi4/psi4)
- [Psi4マニュアル: コア関連API](https://psicode.org/psi4manual/master/api/psi4.core.html)