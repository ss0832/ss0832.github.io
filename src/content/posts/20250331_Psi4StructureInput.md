---
title: 【Psi4】分子構造の入力方法と活用テクニック
published: 2025-03-31
description: "量子化学計算ソフトウェアPsi4における分子構造の入力方法を解説する"
tags: [Psi4]
category: Computational Chemistry
draft: false
---
最終更新：2025-03-31

## 概要

Psi4は、高性能な量子化学計算ソフトウェアであり、分子構造の指定方法として複数の便利な選択肢を提供している。本記事では、Psi4で計算を行うために必要な分子構造の入力方法について、基本から応用まで解説する。Psi4では、Z-matrix形式、直交座標、外部のデータベースからの読み込みなど、様々な形式で分子構造を指定できる。適切な入力方法を選ぶことで、計算の効率化や複雑な系の設定が容易になる。


## Psi4のバージョン

1.9.1

## 基本的な分子構造の定義

Psi4では`molecule`ブロックを使って分子構造を定義する。基本的な構文は以下のようになる：

```python
molecule mol_name {
    [分子の電荷] [スピン多重度]
    [分子構造データ]
}
```

### 1. カルテシアン座標による入力

最も一般的な入力方法は、原子記号とXYZ座標（単位はÅ）を指定する方法である：

```python
molecule water {
    0 1
    O  0.000000  0.000000  0.000000
    H  0.757000  0.586000  0.000000
    H -0.757000  0.586000  0.000000
}
```

### 2. Z-matrix形式による入力

内部座標を使用するZ-matrix形式も広く使われている：

```python
molecule methane {
    0 1
    C
    H 1 1.089
    H 1 1.089 2 109.5
    H 1 1.089 2 109.5 3 120.0
    H 1 1.089 2 109.5 3 -120.0
}
```

ここで、最初の数字は結合する原子番号、二番目は結合距離（Å）、三番目は結合角（度）、四番目は二面角（度）を表す。


## 高度な構造指定オプション

### 1. 単位の指定

デフォルトではÅとdegree（度）が使用されるが、明示的に単位を指定することも可能である：

```python
molecule h2o_bohr {
    0 1
    O  0.000000  0.000000  0.000000
    H  1.430000  1.107000  0.000000
    H -1.430000  1.107000  0.000000
    units bohr
}
```

使用可能な単位：
- `angstrom` (または `ang`): オングストローム（デフォルト）
- `bohr` (または `au`): ボーア単位
- `nm`: ナノメートル
- `pm`: ピコメートル

### 2. 対称性の指定

分子の対称性を明示的に指定したり、無視したりすることができる：

```python
molecule benzene {
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
    symmetry d6h
}
```

または対称性を無視する場合：

```python
molecule no_symmetry {
    0 1
    O  0.000000  0.000000  0.000000
    H  0.757000  0.586000  0.000000
    H -0.757000  0.586000  0.000000
    symmetry c1
}
```

### 3. 変数を使ったパラメトリック指定

Z-matrix内で変数を使用することで、構造パラメータを柔軟に設定できる：

```python
molecule h2o_param {
    0 1
    O
    H 1 R
    H 1 R 2 A
    
    R = 0.96
    A = 104.5
}
```

## 様々な構造指定方法


### 1. 化学データベースからの取得

PubChemなどのオンラインデータベースから直接分子構造を取得できる：

```python
molecule caffeine {
    pubchem:caffeine
}
```

### 2. フラグメント指定

分子をフラグメントに分けて指定することもできる：

```python
molecule water_dimer {
    0 1
    --
    0 1
    O  -1.551007  -0.114520   0.000000
    H  -1.934259   0.762503   0.000000
    H  -0.599677   0.040712   0.000000
    --
    0 1
    O   1.350625   0.111469   0.000000
    H   1.680398  -0.373741  -0.758561
    H   1.680398  -0.373741   0.758561
}
```

## ハンズオン：すべてのオプションを試す

以下にPsi4での分子構造入力に関するすべてのオプションを網羅したハンズオン例を示す：

```python
import psi4
import numpy as np

# Psi4のバージョン確認
print(psi4.__version__)  # 1.9.1

# 1. 基本的なカルテシアン座標による入力
psi4.set_memory('2 GB')
psi4.set_output_file('basic_cartesian.out')

h2o_cart = psi4.geometry("""
    0 1
    O  0.000000  0.000000  0.000000
    H  0.757000  0.586000  0.000000
    H -0.757000  0.586000  0.000000
""")

energy_cart = psi4.energy('scf/cc-pvdz', molecule=h2o_cart) # scfと指定するとHartree-Fock法で計算を実行する。
print(f"H2O カルテシアン座標でのSCFエネルギー: {energy_cart:.8f} Hartree")

# 2. Z-matrix形式による入力
psi4.set_output_file('zmat_format.out')

h2o_zmat = psi4.geometry("""
    0 1
    O
    H 1 R1
    H 1 R2 2 A1
    
    R1 = 0.96
    R2 = 0.96
    A1 = 104.5
""")

energy_zmat = psi4.energy('scf/cc-pvdz', molecule=h2o_zmat)
print(f"H2O Z-matrix形式でのSCFエネルギー: {energy_zmat:.8f} Hartree")

# 3. 単位の指定
psi4.set_output_file('bohr_units.out')

h2o_bohr = psi4.geometry("""
    0 1
    O  0.000000  0.000000  0.000000
    H  1.430000  1.107000  0.000000
    H -1.430000  1.107000  0.000000
    units bohr
""")

energy_bohr = psi4.energy('scf/cc-pvdz', molecule=h2o_bohr)
print(f"H2O ボーア単位でのSCFエネルギー: {energy_bohr:.8f} Hartree")

# 4. 対称性の指定
psi4.set_output_file('symmetry.out')

h2o_sym = psi4.geometry("""
    0 1
    O  0.000000  0.000000  0.000000
    H  0.757000  0.586000  0.000000
    H -0.757000  0.586000  0.000000
    symmetry c2v
""")

energy_sym = psi4.energy('scf/cc-pvdz', molecule=h2o_sym)
print(f"H2O 対称性C2vを指定したSCFエネルギー: {energy_sym:.8f} Hartree")

# 5. 対称性を無視
psi4.set_output_file('no_symmetry.out')

h2o_no_sym = psi4.geometry("""
    0 1
    O  0.000000  0.000000  0.000000
    H  0.757000  0.586000  0.000000
    H -0.757000  0.586000  0.000000
    symmetry c1
""")

energy_no_sym = psi4.energy('scf/cc-pvdz', molecule=h2o_no_sym)
print(f"H2O 対称性を無視したSCFエネルギー: {energy_no_sym:.8f} Hartree")


# 6. 外部データベースからの取得
psi4.set_output_file('pubchem.out')

try:
    # Pythonで実行する場合の正しい構文
    aspirin_db = psi4.geometry("pubchem:aspirin")
    
    energy_aspirin = psi4.energy('scf/sto-3g', molecule=aspirin_db)
    print(f"アスピリン PubChemから取得したSCFエネルギー: {energy_aspirin:.8f} Hartree")
except:
    print("PubChem接続エラー: インターネット接続を確認してください")


# 7. フラグメント指定
psi4.set_output_file('water_dimer.out')

water_dimer = psi4.geometry("""
    0 1
    --
    0 1  # 第1フラグメントの電荷とスピン多重度
    O  -1.551007  -0.114520   0.000000
    H  -1.934259   0.762503   0.000000
    H  -0.599677   0.040712   0.000000
    --
    0 1  # 第2フラグメントの電荷とスピン多重度
    O   1.350625   0.111469   0.000000
    H   1.680398  -0.373741  -0.758561
    H   1.680398  -0.373741   0.758561
    no_com         # 重心を原点に移動しない
    no_reorient    # 分子を再配向しない
""")

energy_dimer = psi4.energy('scf/cc-pvdz', molecule=water_dimer)
print(f"水二量体のSCFエネルギー: {energy_dimer:.8f} Hartree")

# 8. 構造情報の出力
psi4.set_output_file('mol_info.out')

# 分子情報の取得
print("\n分子情報:")
print(f"原子数: {h2o_cart.natom()}")
print(f"電荷: {h2o_cart.molecular_charge()}")
print(f"スピン多重度: {h2o_cart.multiplicity()}")
print(f"対称性: {h2o_cart.schoenflies_symbol()}")

# 構造の出力
print("\n原子座標:")
for i in range(h2o_cart.natom()):
    symbol = h2o_cart.symbol(i)
    x, y, z = h2o_cart.x(i), h2o_cart.y(i), h2o_cart.z(i)
    print(f"{symbol} {x:.6f} {y:.6f} {z:.6f}")
    
# 9. 特殊なケース: アニオン・カチオン
psi4.set_output_file('ion.out')

water_anion = psi4.geometry("""
    -1 1  # 電荷: -1
    O  0.000000  0.000000  0.000000
    H  0.757000  0.586000  0.000000
""")

energy_anion = psi4.energy('scf/cc-pvdz', molecule=water_anion)
print(f"HO- のSCFエネルギー: {energy_anion:.8f} Hartree")

# 10. Python変数を使った高度なパラメータ制御
r_oh = 0.96
angle = 104.5

# 変数を使って構造を定義
h2o_python_var = psi4.geometry(f"""
0 1
O  
H  1 {r_oh} 
H  1 {r_oh*np.cos(np.radians(angle))} 2 {r_oh*np.sin(np.radians(angle))} 

""")

psi4.set_output_file('python_var.out')
energy_python_var = psi4.energy('scf/cc-pvdz', molecule=h2o_python_var)
print(f"Python変数を使った構造定義のSCFエネルギー: {energy_python_var:.8f} Hartree")

# 11. ゴーストアトム（バーチャルサイト）の指定
psi4.set_output_file('ghost_atom.out')

h2o_ghost = psi4.geometry("""
    0 1
    O  0.000000  0.000000  0.000000
    H  0.757000  0.586000  0.000000
    H -0.757000  0.586000  0.000000
    @He 0.0 0.0 1.0  # ヘリウムのゴーストアトム
""")

energy_ghost = psi4.energy('scf/cc-pvdz', molecule=h2o_ghost)
print(f"ゴーストアトムを含む計算のSCFエネルギー: {energy_ghost:.8f} Hartree")
```

## 参考サイト

1. Psi4公式ドキュメント：https://psicode.org/psi4manual/1.9.1/index.html
2. Psi4 GitHub リポジトリ：https://github.com/psi4/psi4
3. Psi4 分子入力ガイド：https://psicode.org/psi4manual/1.9.1/psithon_input.html#psithon-molecule
4. Psi4教育ワークショップ：http://www.psicode.org/workshops.html
5. MolSSI（The Molecular Sciences Software Institute）：https://molssi.org/software/psi4/
6. オンラインZ-matrix生成ツール：https://www.webmo.net/