---
title:  【PySCF】Full CI法による高精度電子状態計算と電子相関評価
published: 2025-03-27
description: ""
tags: [Pyscf]
category: Computational Chemistry
draft: false
---
最終更新：2025-03-27

## 概要

Full CI（完全配置間相互作用）法は、与えられた基底関数セット内で厳密解を求めることができる量子化学計算手法である。postHF法の１つである。

この方法では、可能なすべての電子配置を考慮した波動関数を構築するため、電子相関を取り込むことができる。

電子相関とは、Hartree-Fock法に使用した近似が原因で生じた誤差のことを指す。

2つの分類が存在する。

電子間の衝突によるもの（HF法は平均場近似により電子間衝突を考慮不可）を動的電子相関と呼ぶ。

分子軌道の擬縮退によるもの(HF法では全分子軌道に対して１つの電子配置しか考慮していない)を静的電子相関と呼ぶ。

Full CI法では主に動的電子相関を考慮することが可能である。（静的電子相関も考慮したい場合は、多参照理論(MR)を用いた計算手法が必要である。）

ただし、計算量が電子数と基底関数数に対して指数関数的に増大するという欠点がある。そのため、実用的には小さな分子系や小さな基底関数セットでの使用に限られる。

PySCFライブラリは、Full CI計算を実行するための機能を提供している。
本記事では、PySCFを使用してFull CI計算を実行し、Hartree-Fock（HF）法との比較を通じて電子相関の大きさを評価する方法について解説する。


## ハンズオン

### 1. 基本的なFull CI計算

以下のコードでは、水分子を例にとり、Hartree-Fock計算とFull CI計算を順番に実行して結果を比較する。

```python
# 必要なモジュールをインポート
from pyscf import gto, scf, fci
import numpy as np
import time
# 分子の定義（水分子の例）
mol = gto.Mole()                      # 分子オブジェクトを生成
mol.atom = '''
O  0.000000  0.000000  0.000000
H  0.000000  0.757000  0.586000
H  0.000000 -0.757000  0.586000
'''                                   # 水分子の構造を定義
mol.basis = 'STO-3G'                  # 小さな基底関数セットを指定（Full CIのため）
mol.verbose = 4                        # 詳細な出力レベルを設定
mol.build()                           # 分子オブジェクトを構築

# Hartree-Fock計算の実行
mf = scf.RHF(mol)                     # 制限付きHartree-Fockオブジェクトを作成
mf.conv_tol = 1e-10                   # SCF収束閾値を設定
mf.kernel()                           # SCF計算を実行
print(f"Hartree-Fock エネルギー: {mf.e_tot:.10f} Hartree")  # HF全エネルギーを出力

# Full CI計算の実行
cisolver = fci.FCI(mf)                # Full CIソルバーをHF結果に基づいて初期化
cisolver.max_cycle = 200              # 最大反復回数を設定
cisolver.conv_tol = 1e-10             # CI収束閾値を設定
e_fci, fcivec = cisolver.kernel()     # Full CI計算を実行し、エネルギーと波動関数ベクトルを取得
print(f"Full CI エネルギー: {e_fci:.10f} Hartree")

# 電子相関エネルギーの計算と比較
corr_energy = e_fci - mf.e_tot        # 電子相関エネルギー = FCI - HF
print(f"電子相関エネルギー: {corr_energy:.10f} Hartree")
print(f"電子相関エネルギー: {corr_energy*627.5095:.6f} kcal/mol")  # Hartreeからkcal/molに変換

# 相対的な電子相関の寄与
rel_corr = abs(corr_energy / mf.e_tot) * 100
print(f"全エネルギーに対する電子相関の寄与: {rel_corr:.6f}%")
```

### 2. 配置状態関数（CSF）の解析

Full CI計算では、さまざまな電子配置の寄与を分析することができる。以下のコードでは、主要な配置の寄与を抽出して表示する。

```python
### 前のコードの続き
# 主要な配置状態関数の寄与を分析
print("\n主要な配置状態関数の寄与:")

# Full CI波動関数の解析
norb = mf.mo_coeff.shape[1]           # 分子軌道の数を取得
nelec = mol.nelec                      # α電子とβ電子の数を取得
na, nb = nelec                         # α電子とβ電子の数を分離

# 主要な配置状態を表示する関数
def print_configurations(fcivec, norb, nelec, threshold=0.005):
    from pyscf import fci
    na, nb = nelec
    # 波動関数の形状を(αの配置数, βの配置数)に変換
    fcivec_reshape = fcivec.reshape(fci.cistring.num_strings(norb, na), 
                                   fci.cistring.num_strings(norb, nb))
    # 各配置の寄与を計算
    for ia, stringa in enumerate(fci.cistring.gen_strings4orblist(range(norb), na)):
        for ib, stringb in enumerate(fci.cistring.gen_strings4orblist(range(norb), nb)):
            coeff = fcivec_reshape[ia, ib]
            # 閾値以上の係数を持つ配置を表示
            if abs(coeff) >= threshold:
                occ_a = bin(stringa)[2:].zfill(norb)
                occ_b = bin(stringb)[2:].zfill(norb)
                print(f"配置: α={occ_a}, β={occ_b}, 係数={coeff:.6f}, 寄与={coeff**2:.6f}")

# 0.05以上の寄与を持つ配置を表示
try:
    print_configurations(fcivec, norb, nelec, threshold=0.005)
except:
    print("配置状態関数の解析はより複雑な実装を要します")

# 第一励起状態の計算（必要に応じて）
try:
    e_ex, c_ex = cisolver.kernel(nroots=2)[0:2]  # 基底状態と第一励起状態を計算
    print(f"\n第一励起エネルギー: {e_ex[1]:.10f} Hartree")
    print(f"励起エネルギー: {(e_ex[1]-e_ex[0])*27.2114:.6f} eV")  # Hartreeから電子ボルトに変換
except:
    print("\n励起状態の計算はより複雑な設定を要します")
```

### 3. 異なる基底関数セットでの比較

基底関数セットの違いによる電子相関の変化を調査する。

```python
### 前のコードの続き
# 異なる基底関数セットを用いた比較
basis_sets = ['STO-3G', 'STO-6G', '3-21G', '6-31G']  # 比較する基底関数セット

print("\n基底関数セットによる電子相関エネルギーの比較:")
print("基底関数    HF Energy       FCI Energy      相関エネルギー  　Full CI計算時間")
print("---------------------------------------------------------------------------")

for basis in basis_sets:
    try:

        # 新しい分子オブジェクトを作成
        mol_tmp = gto.Mole()
        mol_tmp.atom = mol.atom
        mol_tmp.basis = basis      # 基底関数セットを変更
        mol_tmp.build()
        
        # HF計算
        mf_tmp = scf.RHF(mol_tmp)
        mf_tmp.kernel()
        
        # Full CI計算（小さな分子でのみ実行可能）
        cisolver_tmp = fci.FCI(mf_tmp)
        full_ci_start_time = time.time()
        e_fci_tmp = cisolver_tmp.kernel()[0]
        full_ci_end_time = time.time()
        # 相関エネルギーの計算
        corr_tmp = e_fci_tmp - mf_tmp.e_tot
        elapsed_time = full_ci_end_time - full_ci_start_time
        # 結果の表示
        print(f"{basis:8s}  {mf_tmp.e_tot:12.8f}  {e_fci_tmp:12.8f}  {corr_tmp:12.8f}  {elapsed_time:12.8f}")
    except Exception as e:
        print(f"{basis:8s}  計算エラー: 基底関数が大きすぎるかメモリ不足です - {str(e)}")
```

## コードの詳細説明

### 基本的なFull CI計算部分

- `from pyscf import gto, scf, fci`: 分子定義、SCF計算、Full CI計算のモジュールをインポート
- `mol = gto.Mole()`: 分子オブジェクトを生成
- `mol.atom = '''...'''`: 水分子の幾何構造を定義（O原子が原点、H原子が適切な位置）
- `mol.basis = 'STO-3G'`: 小さな基底関数セットを選択（Full CIの計算コスト削減のため）
- `mf = scf.RHF(mol)`: 制限付きHartree-Fock計算オブジェクトの作成
- `mf.kernel()`: Hartree-Fock SCF計算の実行
- `cisolver = fci.FCI(mf)`: HF結果に基づいたFull CIソルバーの初期化
- `e_fci, fcivec = cisolver.kernel()`: Full CI計算の実行とエネルギー・波動関数の取得
- `corr_energy = e_fci - mf.e_tot`: 電子相関エネルギーを計算（FCI - HF）
- `rel_corr = abs(corr_energy / mf.e_tot) * 100`: 全エネルギーに占める相関エネルギーの割合を百分率で計算

### 配置状態関数の解析部分

- `norb = mf.mo_coeff.shape[1]`: 分子軌道の総数を取得
- `nelec = mol.nelec`: 系内の電子数（α電子とβ電子）を取得
- `print_configurations(fcivec, norb, nelec, threshold=0.05)`: 寄与の大きい配置状態を表示（係数の絶対値が0.05以上のもの）
- `fcivec_reshape = fcivec.reshape(...)`: CI波動関数ベクトルを適切な形状に変換
- `for ia, stringa in enumerate(...)`: α電子の配置を列挙
- `for ib, stringb in enumerate(...)`: β電子の配置を列挙
- `coeff = fcivec_reshape[ia, ib]`: 特定の配置の係数を取得
- `coeff**2`: 配置の寄与（確率）を計算
- `e_ex, c_ex = cisolver.kernel(nroots=2)[0:2]`: 基底状態と第一励起状態を計算
- `(e_ex[1]-e_ex[0])*27.2114`: 励起エネルギーをHartreeから電子ボルト（eV）に変換

### 異なる基底関数セットでの比較部分

- `basis_sets = ['STO-3G', 'STO-6G', '3-21G', '6-31G']`: 比較する基底関数セットのリスト
- `mol_tmp = gto.Mole()`: 新しい分子オブジェクトを生成
- `mol_tmp.basis = basis`: 基底関数セットを変更
- `cisolver_tmp = fci.FCI(mf_tmp)`: 新しい基底関数セットでのFull CIソルバーを作成
- `corr_tmp = e_fci_tmp - mf_tmp.e_tot`: 新しい基底関数セットでの相関エネルギーを計算
- `try/except`: 大きな基底関数セットでのFull CI計算が失敗した場合のエラーハンドリング

## 電子相関とその評価

Full CI法は、与えられた基底関数セット内で可能なすべての電子配置を考慮するため、電子相関を完全に取り込むことができる。HF法と比較することで、以下のような電子相関の特性が評価できる：

1. **電子相関エネルギー**：HF法と比べてFull CI法ではエネルギーがどれだけ低下するかを示す。通常、負の値となる。

2. **相関エネルギーの相対的寄与**：全エネルギーに対する相関エネルギーの割合。これが大きいと電子相関を考慮出来ていることとなる。

3. **基底関数依存性**：基底関数セットが大きくなるにつれて、相関エネルギーの絶対値も大きくなる傾向がある。これは、より多くの仮想軌道が利用可能になり、より多くの励起配置が考慮されるためである。


## まとめ

PySCFを使用したFull CI計算により、与えられた基底関数セットにおける厳密な電子状態を求めることができる。これにより、HF法では記述できない電子相関効果を完全に取り込むことが可能になる。ただし、計算コストが非常に高いため、実際の応用は小さな分子系に限られる。

電子相関エネルギーの定量的な評価は、量子化学計算の精度を理解するうえで重要であり、Full CI計算はその参照点として価値がある。より大きな系では、CCSD(T)や多参照法など、計算効率と精度のバランスがとれた他の高精度手法が代替として使用される。


### 参考
- https://www.jstage.jst.go.jp/article/oubutsu/86/8/86_720/_pdf
- https://www2.itc.nagoya-u.ac.jp/pub/pdf/pdf/vol07_03/333_350kouza.pdf
