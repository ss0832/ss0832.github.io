---
title: 【PySCF】PCM法を用いた溶媒効果の考慮方法
published: 2025-03-31
description: "PySCFで分子の溶媒効果をPCM法によって計算する方法について解説する"
tags: [pyscf, PCM]
category: Computational Chemistry
draft: false
---
最終更新：2025-03-31

## 概要

PySCFは、量子化学計算のためのPythonベースのオープンソースソフトウェアである。バージョン2.8.0では、溶媒効果を考慮するための様々な方法が実装されており、その中でもPCM（Polarizable Continuum Model）法は広く使われている手法である。PCM法は分子を溶媒中の空洞内に配置し、溶媒を連続的な誘電体として扱うことで、溶媒効果を効率的に計算する方法である。本稿では、PySCFでPCM法を用いて溶媒効果を考慮する方法について解説する。

(注意)PCM法のPySCFの実装が実験的であるため、ここで示した方法が今後使用不可になる可能性がある。

## 必要なモジュール

PySCFのPCMモジュールを使用するには、以下のライブラリが必要である：
- pyscf
- geometric（構造最適化時に利用）

## 基本的な使用方法

PySCFでPCM法を利用する基本的な流れは以下の通りである：

1. PCMオブジェクトを作成する
2. SCF計算に溶媒効果を組み込む
3. 必要な計算を実行する

以下は、水分子における基本的なPCM-HF計算の例である：

```python
import numpy as np
from pyscf import gto, scf, solvent

# 分子を定義
mol = gto.M(
    atom='''
    O  0.0  0.0  0.0
    H  0.0  0.0  0.9
    H  0.0  0.8  -0.3
    ''',
    basis='6-31g'
)

# PCMオブジェクトを作成（水を溶媒として設定）
pcm = solvent.PCM(mol)
pcm.eps = 78.3553  # 水の誘電率
pcm.method = 'C-PCM'  # C-PCM法を使用

# PCM法を用いたHF計算
mf = scf.RHF(mol).PCM(pcm)
e_pcm = mf.kernel()

# 結果の出力
print(f'PCM-HF energy: {e_pcm:.8f}')
```

## PCMパラメーターの設定

PySCFでPCMを利用する際に設定できる主要なパラメーターは以下の通りである：

1. `eps`: 溶媒の誘電率（例：水は78.3553、アセトニトリルは35.688）
2. `method`: PCMの種類（'C-PCM'、'IEF-PCM'、'SS(V)PE'など）
3. `lebedev_order`: Lebedev格子の次数（精度と計算コストのバランス）
4. `cavity_radii`: 空洞面を構築する際の原子半径

以下は、異なる溶媒（エタノール）でのPCM-DFT計算の例である：

```python
import numpy as np
from pyscf import gto, scf, dft, solvent

# 分子を定義
mol = gto.M(
    atom='''
    C  -0.755  0.022  -0.046
    C  0.764  -0.001  0.036
    O  1.327  1.303  -0.029
    H  -1.131  -0.961  0.287
    H  -1.095  0.262  -1.056
    H  -1.175  0.773  0.634
    H  1.142  -0.549  0.917
    H  1.154  -0.504  -0.856
    H  2.294  1.263  0.025
    ''',
    basis='def2-svp'
)

# PCMオブジェクトを作成（エタノールを溶媒として設定）
pcm = solvent.PCM(mol)
pcm.eps = 24.852  # エタノールの誘電率
pcm.method = 'IEF-PCM'  # IEF-PCM法
pcm.lebedev_order = 29  # Lebedev格子の次数を上げて精度向上
pcm.cavity_radii = 'bondi'  # Bondi半径セットを使用

# PCM法を用いたDFT計算（B3LYP汎関数）
mf = dft.RKS(mol).PCM(pcm)
mf.xc = 'B3LYP'
e_pcm = mf.kernel()

print(f'PCM-B3LYP energy: {e_pcm:.8f}')
```


## 溶媒効果を考慮した構造最適化

PySCFとgeometricライブラリを組み合わせることで、PCMを考慮した構造最適化も可能である：

```python
from pyscf.geomopt.geometric_solver import optimize

# PCMを考慮したDFT計算オブジェクトを作成
mf = dft.RKS(mol).PCM(pcm)
mf.xc = 'B3LYP'

# PCMを考慮した構造最適化
mol_opt = optimize(mf)

print("最適化された構造(Bohr単位):")
for i in range(mol_opt.natm):
    atom = mol_opt.atom_symbol(i)
    coord = mol_opt.atom_coords()[i]
    print(f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")
```

## ハンズオン：PCM法の様々なオプションを試す

以下にPCM法での様々なオプションを含むハンズオン例をまとめて示す：

```python
import numpy as np
from pyscf import gto, scf, dft, solvent, mp, cc
from pyscf.geomopt.geometric_solver import optimize

# 分子を定義（アセトンを例として）
mol = gto.M(
    atom='''
    C  0.000  0.000  0.000
    O  0.000  0.000  1.220
    C -1.290  0.000 -0.800
    C  1.290  0.000 -0.800
    H -2.180  0.000 -0.167
    H -1.320  0.890 -1.433
    H -1.320 -0.890 -1.433
    H  2.180  0.000 -0.167
    H  1.320  0.890 -1.433
    H  1.320 -0.890 -1.433
    ''',
    basis='cc-pVDZ',
    verbose=4
)

# 1. 様々な溶媒でのPCM計算

# 水中での計算
pcm_water = solvent.PCM(mol)
pcm_water.eps = 78.3553
pcm_water.method = 'C-PCM'
mf_water = scf.RHF(mol).PCM(pcm_water)
e_water = mf_water.kernel()
print(f'水中でのPCM-HF energy: {e_water:.8f}')

# アセトニトリル中での計算
pcm_acn = solvent.PCM(mol)
pcm_acn.eps = 35.688
pcm_acn.method = 'C-PCM'
mf_acn = scf.RHF(mol).PCM(pcm_acn)
e_acn = mf_acn.kernel()
print(f'アセトニトリル中でのPCM-HF energy: {e_acn:.8f}')

# 2. 異なるPCM法の比較

# IEF-PCM法
pcm_ief = solvent.PCM(mol)
pcm_ief.eps = 78.3553
pcm_ief.method = 'IEF-PCM'
mf_ief = scf.RHF(mol).PCM(pcm_ief)
e_ief = mf_ief.kernel()
print(f'IEF-PCM法によるHF energy: {e_ief:.8f}')

# SS(V)PE法
pcm_ssvpe = solvent.PCM(mol)
pcm_ssvpe.eps = 78.3553
pcm_ssvpe.method = 'SSVPE'
mf_ssvpe = scf.RHF(mol).PCM(pcm_ssvpe)
e_ssvpe = mf_ssvpe.kernel()
print(f'SS(V)PE法によるHF energy: {e_ssvpe:.8f}')

# 3. Lebedev格子の次数による精度比較

# 低次（高速だがやや粗い）
pcm_low = solvent.PCM(mol)
pcm_low.eps = 78.3553
pcm_low.lebedev_order = 17
mf_low = scf.RHF(mol).PCM(pcm_low)
e_low = mf_low.kernel()
print(f'Lebedev次数17でのPCM-HF energy: {e_low:.8f}')

# 高次（精度良いが遅い）
pcm_high = solvent.PCM(mol)
pcm_high.eps = 78.3553
pcm_high.lebedev_order = 41
mf_high = scf.RHF(mol).PCM(pcm_high)
e_high = mf_high.kernel()
print(f'Lebedev次数41でのPCM-HF energy: {e_high:.8f}')

# 4. 異なる空洞半径の設定

# Bondi半径
pcm_bondi = solvent.PCM(mol)
pcm_bondi.eps = 78.3553
pcm_bondi.cavity_radii = 'bondi'
mf_bondi = scf.RHF(mol).PCM(pcm_bondi)
e_bondi = mf_bondi.kernel()
print(f'Bondi半径を使用したPCM-HF energy: {e_bondi:.8f}')

# UFF半径
pcm_uff = solvent.PCM(mol)
pcm_uff.eps = 78.3553
pcm_uff.cavity_radii = 'uff'
mf_uff = scf.RHF(mol).PCM(pcm_uff)
e_uff = mf_uff.kernel()
print(f'UFF半径を使用したPCM-HF energy: {e_uff:.8f}')

# 5. DFT汎関数との組み合わせ

# B3LYP/PCM
pcm = solvent.PCM(mol)
pcm.eps = 78.3553
mf_b3lyp = dft.RKS(mol).PCM(pcm)
mf_b3lyp.xc = 'B3LYP'
e_b3lyp = mf_b3lyp.kernel()
print(f'PCM-B3LYP energy: {e_b3lyp:.8f}')

# ωB97X-D/PCM
mf_wb97 = dft.RKS(mol).PCM(pcm)
mf_wb97.xc = 'wb97xd'
e_wb97 = mf_wb97.kernel()
print(f'PCM-ωB97X-D energy: {e_wb97:.8f}')

# 6. 後処理計算

# PCM-MP2
mp2 = mp.MP2(mf_water)
e_mp2_corr, _ = mp2.kernel()
print(f'PCM-MP2 correlation energy: {e_mp2_corr:.8f}')
print(f'PCM-MP2 total energy: {mp2.e_tot:.8f}')

# PCM-CCSD
cc_obj = cc.CCSD(mf_water)
cc_obj.kernel()
print(f'PCM-CCSD correlation energy: {cc_obj.e_corr:.8f}')
print(f'PCM-CCSD total energy: {cc_obj.e_tot:.8f}')

# 7. 構造最適化（小さな分子で実行）
small_mol = gto.M(
    atom='''
    O  0.0  0.0  0.0
    H  0.0  0.0  0.9
    H  0.0  0.8  -0.3
    ''',
    basis='6-31g'
)

pcm_opt = solvent.PCM(small_mol)
pcm_opt.eps = 78.3553
mf_opt = dft.RKS(small_mol).PCM(pcm_opt)
mf_opt.xc = 'B3LYP'

# 構造最適化を実行
mol_opt = optimize(mf_opt)
print("水中での最適化構造(Bohr単位):")
for i in range(mol_opt.natm):
    atom = mol_opt.atom_symbol(i)
    coord = mol_opt.atom_coords()[i]
    print(f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")

# 溶媒中と気相中でのエネルギー差（溶媒和自由エネルギー）
mf_gas = scf.RHF(mol)
e_gas = mf_gas.kernel()
solvation_energy = e_water - e_gas
print(f'溶媒和エネルギー: {solvation_energy:.8f} Hartree')
print(f'溶媒和エネルギー: {solvation_energy*627.5095:.4f} kcal/mol')
```

上記のコード例を実行することで、様々な溶媒条件下での計算結果を比較が可能である。パラメータを変更して検討することで、PCM法のパラメータを選択するための参考になるだろう。PySCFに実装されているPCM法はHF法だけでなく、DFT、MP2、CCSDなど様々な電子状態計算法と組み合わせて利用可能である。

