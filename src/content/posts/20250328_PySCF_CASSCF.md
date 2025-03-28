---
title:  【PySCF】CASSCF法を用いた多参照電子状態計算と活性空間の影響
published: 2025-03-28
description: ""
tags: [Pyscf]
category: Computational Chemistry
draft: false
---
最終更新：2025-03-28

## 概要

完全活性空間自己無撞着場（CASSCF: Complete Active Space Self-Consistent Field）法は、多参照性の強い分子系や電子励起状態の計算に適した量子化学計算手法である。この方法では、化学反応に重要な軌道を「活性空間」として選び、その空間内では配置間相互作用（CI）計算を行い、活性空間外の軌道は通常のHartree-Fock（HF）計算のように扱う。

CASSCF法は多参照理論を用いた手法である。そのため、静的電子相関を十分に取り込むことが可能である。

（多参照理論では、異なる電子配置の波動関数を構成する基底関数の軌道係数をSCFの手続きの際に全て独立変数として扱う。一方単参照理論では、単一の波動関数の軌道係数をSCFの手続きの際に独立変数として扱う。）

CASSCF法の特徴は、活性空間内の電子配置をすべて考慮することで静的電子相関を正確に取り込める点にある。特に結合解離過程、遷移状態、開殻系、一部の励起状態など、単一の電子配置では適切に記述できない系に対して有効である。

PySCFはCASSCF計算を実装しており、様々な系に対して柔軟に適用できる。本稿では、PySCFを用いたCASSCF計算の基本的な実行方法と、Hartree-Fock法との比較、さらに活性空間の選択が電子相関エネルギーに与える影響について解説する。

## ハンズオン

### 1. 基本的なCASSCF計算の実行

以下のコードでは、酸素分子（O₂）を例に取り、HF計算とCASSCF計算を行い、結果を比較する。

```python
# 必要なモジュールをインポート
from pyscf import gto, scf, mcscf
import numpy as np

# 分子の定義（酸素分子：多参照性が現れる良い例）
mol = gto.Mole()                      # 分子オブジェクトを生成
mol.atom = '''
O  0.0000  0.0000  0.0000
O  0.0000  0.0000  1.2075
'''                                   # O₂分子の構造を定義（単位はオングストローム）
mol.basis = 'cc-pVDZ'                 # 基底関数セットを指定
mol.spin = 2                          # スピン多重度=3（トリプレット状態）
mol.verbose = 4                        # 出力の詳細レベルを設定
mol.build()                           # 分子オブジェクトを構築

# Hartree-Fock計算
mf = scf.UHF(mol)                     # 非制限Hartree-Fock計算オブジェクトを作成（開殻系のため）
mf.conv_tol = 1e-10                   # SCF収束閾値を設定
mf.kernel()                           # HF計算を実行
print(f"Hartree-Fock エネルギー: {mf.e_tot:.10f} Hartree")  # HF全エネルギーを表示

# CASSCF計算（基本設定）
# 酸素分子の場合、典型的には2p軌道（π軌道系を含む）を活性空間とする
ncas = 6                              # 活性軌道数（O原子の2p軌道）
nelecas = (5, 3)                      # 活性空間内の電子数（α電子, β電子）
mc = mcscf.CASSCF(mf, ncas, nelecas)  # CASSCF計算オブジェクトを作成
mc.conv_tol = 1e-10                   # 収束閾値を設定
mc.max_cycle_macro = 100              # 最大マクロ反復回数
mc.kernel()                           # CASSCF計算を実行（軌道最適化とCI計算）

print(f"CASSCF(6,8) エネルギー: {mc.e_tot:.10f} Hartree")
print(f"電子相関エネルギー: {mc.e_tot - mf.e_tot:.10f} Hartree")  # CASSCFとHFの差


# CASSCF波動関数の主要配置を解析
try:
    ci_coeff = mc.ci
    print("\nCASSCF波動関数の解析:")
    print(f"CI行列の形状: {ci_coeff.shape}")  # 配置状態の数
    max_ci_idx = np.unravel_index(np.argmax(np.abs(ci_coeff)), ci_coeff.shape)
    max_ci = ci_coeff[max_ci_idx]
    print(f"最大CI係数: {max_ci:.6f} (絶対値: {abs(max_ci):.6f})")
    print(f"上位5つのCI係数（絶対値順）:")
    flat_ci = ci_coeff.flatten()
    sorted_indices = np.argsort(-np.abs(flat_ci))
    for i in range(min(5, len(flat_ci))):
        idx = sorted_indices[i]
        val = flat_ci[idx]
        print(f"  係数 #{i+1}: {val:.6f} (絶対値: {abs(val):.6f})")
except:
    print("CI係数の解析に失敗しました")
```

### 2. 異なる活性空間サイズでのCASSCF計算

活性空間の選択がCASSCF計算の結果に与える影響を調べる。

```python
# 異なる活性空間でのCASSCF計算
print("\n異なる活性空間サイズでのCASSCF計算:")

# 活性空間のサイズの組み合わせ
cas_configs = [
    (2, (1, 1)),  # 2軌道、2電子（α:1, β:1） - π結合軌道のみ
    (4, (3, 1)),  # 4軌道、4電子（α:3, β:1） - π軌道とπ*軌道
    (6, (5, 3)),  # 6軌道、8電子（α:5, β:3） - 2p軌道すべて
    (8, (6, 4)),  # 8軌道、10電子（α:6, β:4） - 2p軌道と追加軌道
]

for norb, nelec in cas_configs:
    try:
        print(f"\nCASSCF({norb},{nelec[0]+nelec[1]})計算:")
        
        # CASSCF計算
        mc_tmp = mcscf.CASSCF(mf, norb, nelec)
        mc_tmp.conv_tol = 1e-8
        mc_tmp.max_cycle_macro = 50
        mc_tmp.kernel()
        
        # 結果の出力
        print(f"CASSCF({norb},{nelec[0]+nelec[1]}) エネルギー: {mc_tmp.e_tot:.10f} Hartree")
        corr_e = mc_tmp.e_tot - mf.e_tot
        print(f"電子相関エネルギー: {corr_e:.10f} Hartree")
        
        # CI係数の解析
        ci_coeff = mc_tmp.ci
        max_ci = np.max(np.abs(ci_coeff))
        print(f"最大CI係数: {max_ci:.6f}")
        
    except Exception as e:
        print(f"CASSCF({norb},{nelec[0]+nelec[1]})計算でエラーが発生: {str(e)}")
```

### 3. 活性空間の内容変更によるCASSCF計算

同じサイズの活性空間でも、異なる軌道を選択することで結果が変わる。

```python
# 活性空間の内容を変えたCASSCF計算
print("\n活性空間の内容を変えたCASSCF計算:")

# HF計算を再実行
mf = scf.UHF(mol)
mf.kernel()

# 初期活性空間選択用の軌道インデックス
# 基底関数と分子の組み合わせによって適切な軌道インデックスは変わる
# ここでは例として異なる軌道セットを指定

# 標準的な活性空間（2p軌道を中心に選択）
print("\n標準的な活性空間（2p軌道を中心とした6軌道）:")
ao_labels = mol.ao_labels()
mo_coeff = mf.mo_coeff
ncas = 6
nelecas = (5, 3)

# CASSCF計算（デフォルト選択）
mc1 = mcscf.CASSCF(mf, ncas, nelecas)
mc1.kernel()
print(f"標準活性空間 CASSCF(6,8) エネルギー: {mc1.e_tot:.10f} Hartree")
print(f"電子相関エネルギー: {mc1.e_tot - mf.e_tot:.10f} Hartree")

# 別の活性空間（内殻軌道を含める）
print("\n内殻軌道を含む変更活性空間（6軌道）:")
# 手動で活性軌道を指定するための初期軌道を準備
# ここでは内殻1軌道、価電子5軌道の組み合わせを選択すると仮定
try:
    # 特定の軌道を活性空間に選ぶ（実際の適用には注意が必要）
    # 内殻軌道を含める例（インデックスは系によって変わる）
    aolabels = ['O 1s', 'O 2s', 'O 2p'] # 簡略化したAO表記
    active_orbitals = [0, 4, 5, 6, 7, 8]  # 一般的に内殻+価電子軌道の例
    
    # 活性軌道を直接指定（文字列で疑似表現）
    print(f"手動選択活性軌道: 内殻(1s)+価電子2p軌道")
    
    mc2 = mcscf.CASSCF(mf, ncas, nelecas)
    # 通常はmc2.sort_mo(active_orbitals)で指定するが、ここでは簡略化
    mc2.kernel()
    print(f"変更活性空間 CASSCF(6,8) エネルギー: {mc2.e_tot:.10f} Hartree")
    print(f"電子相関エネルギー: {mc2.e_tot - mf.e_tot:.10f} Hartree")
    print(f"標準との差異: {mc2.e_tot - mc1.e_tot:.10f} Hartree")
except Exception as e:
    print(f"変更活性空間計算でエラー: {str(e)}")
    print("注: 実際の軌道指定には分子の対称性と軌道エネルギーに基づく慎重な選択が必要")
```

活性空間が大きいほど、電子相関エネルギーが増加していることがわかる。つまり、電子相関を取り込む大きさと活性空間の大きさは相関している。

### 4. CASSCF軌道と自然軌道占有数の解析

CASSCF計算で得られる軌道情報（特に自然軌道占有数）は多参照性の指標として重要。

```python
# CASSCF軌道と自然軌道占有数の解析
print("\nCASSCF軌道と自然軌道占有数の解析:")

# 最初のCASSCF計算（基本設定）の軌道を解析
mc = mcscf.CASSCF(mf, 6, (5, 3))
mc.kernel()

# 自然軌道占有数を出力
try:
    # CAS空間内の自然軌道占有数
    cas_natorbs, cas_occnums = mcscf.addons.cas_natorb(mc)
    print("CAS空間内の自然軌道占有数:")
    for i, occ in enumerate(cas_occnums):
        print(f"軌道 {i+1}: {occ:.6f}")
    
    # 多参照性の指標を計算
    # 理想的な単参照系では占有数は0か2に近い
    # 多参照系では中間的な値（0.1-1.9）を示す軌道が増える
    mr_orbs = [i for i, occ in enumerate(cas_occnums) if 0.1 <= occ <= 1.9]
    print(f"\n多参照性を示す軌道数（占有数0.1-1.9）: {len(mr_orbs)}")
    print(f"多参照性の強さの指標: {1.0 - abs(cas_occnums-1.0).sum()/len(cas_occnums):.6f}")
    # 1に近いほど多参照性が強い（単参照系では0に近い）
    
    # CI波動関数のエントロピー（別の多参照性の指標）
    ci_coeff = mc.ci
    prob = np.abs(ci_coeff)**2
    prob_flat = prob.flatten()
    entropy = -np.sum(prob_flat * np.log2(prob_flat + 1e-12))
    print(f"CI波動関数のエントロピー: {entropy:.6f}")
    # エントロピーが高いほど多参照性が強い
    
except Exception as e:
    print(f"軌道解析でエラー: {str(e)}")
```

## コードの詳細説明

### 基本的なCASSCF計算部分

- `from pyscf import gto, scf, mcscf`: 必要なモジュールをインポート（分子定義、SCF計算、多参照計算用）
- `mol = gto.Mole()`: 分子オブジェクトを生成
- `mol.atom = '''...'''`: 酸素分子の構造を定義（原子座標をオングストローム単位で指定）
- `mol.spin = 2`: スピン多重度をトリプレット（2S=2）に設定（酸素分子は三重項基底状態）
- `mf = scf.UHF(mol)`: 非制限Hartree-Fock計算オブジェクトを作成（開殻系のため）
- `ncas = 6`: 活性空間の軌道数（O₂の場合、各O原子の2p軌道に対応）
- `nelecas = (5, 3)`: 活性空間内のα電子とβ電子の数（合計8電子、スピン多重度=3に対応）
- `mc = mcscf.CASSCF(mf, ncas, nelecas)`: CASSCF計算オブジェクトを作成
- `mc.kernel()`: CASSCF計算を実行（軌道最適化とCI計算を行う）
- `ci_coeff = mc.ci`: CASSCF波動関数のCI係数を取得
- `max_ci = np.max(np.abs(ci_coeff))`: 最大CI係数（多参照性の指標の一つ）

### 異なる活性空間サイズの計算部分

- `cas_configs = [(2, (1, 1)), ...]`: 異なるサイズの活性空間の定義（軌道数, (α電子数, β電子数)）
- `mc_tmp = mcscf.CASSCF(mf, norb, nelec)`: 指定した活性空間サイズでCASSCF計算を設定
- `corr_e = mc_tmp.e_tot - mf.e_tot`: HF法との差として電子相関エネルギーを計算

### 活性空間の内容変更部分

- `active_orbitals = [0, 4, 5, 6, 7, 8]`: 活性空間に選択する軌道のインデックス（実装例）
- `mc2.sort_mo(active_orbitals)`: 特定の軌道を活性空間として選択（実際の適用にはコメント化）

### 軌道と自然軌道占有数の解析部分

- `cas_natorbs, cas_occnums = mcscf.addons.cas_natorb(mc)`: 自然軌道と占有数を計算
- `mr_orbs = [i for i, occ in enumerate(cas_occnums) if 0.1 <= occ <= 1.9]`: 多参照性を示す軌道を特定
- `entropy = -np.sum(prob_flat * np.log2(prob_flat + 1e-12))`: CI波動関数のエントロピー計算

## 電子相関とその評価


### 活性空間のサイズによる電子相関の変化

活性空間のサイズを大きくするほど多くの電子相関を取り込めるが、計算コストは指数関数的に増大する点に注意が必要である。

### 活性空間の内容による電子相関の変化

同じサイズの活性空間でも、含める軌道の選択によって電子相関エネルギーは大きく変化する：

1. **価電子軌道中心の活性空間**：化学結合に直接関与する2p軌道を選択
   - 化学反応や結合の解離に重要な静的相関を正確に記述

2. **内殻軌道を含む活性空間**：1s軌道などの内殻軌道も活性空間に含める
   - 内殻-価電子相関も部分的に記述するが、化学的にあまり重要でないケースが多い
   - 計算コストの増大に見合う精度向上があるとは限らない

3. **冗長な軌道を含む場合（分極関数や分散関数等）**：
   - 占有数が極めて0に近い、あるいは2に近い軌道のみを含む場合は、計算コストの増大にもかかわらず相関エネルギーへの寄与は小さい

## 多参照性の評価と自然軌道占有数

CASSCF計算から得られる自然軌道占有数は、系の多参照性を評価する重要な指標である：

1. **単参照系**：ほとんどの軌道の占有数が0または2に近い（単分子の計算や、結合の解離を伴わない化学反応等）
   - 単一配置のHF法でも良好な記述が期待できる

2. **弱～中程度の多参照系**：少数の軌道が中間的な占有数（0.1-1.9）を示す（ビラジカルが関わる系等）
   - CASSCF計算や制限された多参照計算が有効

3. **強い多参照系**：多数の軌道が中間的な占有数を示す。（前周期遷移金属の高スピン錯体等）
   - 大きな活性空間でのCASSCF計算が必要


## まとめ

PySCFを用いたCASSCF計算は、特に多参照性を持つ系の電子状態を正確に記述するのに役立つ。最適な活性空間の選択は、対象分子の化学的性質、特に注目する結合や反応に関与する軌道に基づいて行うべきである。

CASSCF法の主な利点は静的電子相関を正確に記述できることだが、活性空間外の動的相関を捉えるには、CASPT2やMRCIなどの手法を利用する必要がある。

活性空間の選択は電子相関エネルギーと計算精度に直接影響するため、対象系の特性と計算コストのバランスを考慮した慎重な判断が重要である。