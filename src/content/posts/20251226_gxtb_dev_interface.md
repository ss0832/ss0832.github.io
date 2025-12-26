---
title: 【計算化学】 正式リリース前のg-xTBをASEで利用するための非公式ラッパー「pygxtb」
published: 2025-12-26
tags: [pygxtb, python, ase, g-xTB]
category: Computational Chemistry
draft: false
---
最終更新：2025-12-26


### 概要


`pygxtb`の開発段階のコミット時点(https://github.com/grimme-lab/g-xtb, Commit: `df43e00`)の半経験的量子化学計算手法であるg-xTB (General-Purpose Extended Tight-Binding) を、Pythonスクリプト内から容易に利用するための非公式インターフェースライブラリです。

レポジトリのURL：　https://github.com/ss0832/pygxtb/


ASE (Atomic Simulation Environment) のインターフェースとして実装されているため、既存のASEワークフロー（構造最適化、MD等）にスムーズに組み込むことが可能です。

Grimme Labが提供する g-xTB の開発段階のバイナリおよびパラメータファイルがパッケージに含まれており、パッケージ内部で自動的にパス解決を行います。これにより、別途環境変数を設定する手間が省かれています。

（注意）本ライブラリがラップしている `g-xTB` は、2025年12月時点で開発段階（プレリリース）のものです。本ツールは利便性のためにバイナリを同梱していますが、将来的に公式リリースされるバージョンとは仕様やパラメータが変更される可能性があります。

### g-xTBについて：Full DFTに迫る次世代半経験的手法

本ライブラリがインターフェースを提供する g-xTB は、Bonn大学のGrimme研究室によって開発された最新の半経験的量子化学計算手法です。

2025年に公開されたプレプリント によると、本手法は従来のGFNn-xTBファミリーの後継として、「半経験的手法（SQM）の速度」を維持しつつ「DFT（密度汎関数法）の精度」に近づくことを目指して設計されました。

論文で報告されている主な特徴は以下の通りです。

1. 設計思想とターゲット

g-xTBは、領域分割ハイブリッド汎関数である $\omega$B97M-V/def2-TZVPPD レベルの精度を、Tight-Binding（TB）法の計算コストで実現することを目標としています。これを実現するために、以下の新しい物理モデルが導入されています。

分子内原子適応型AO基底 (q-vSZP): 原子の電荷や配位数に応じて基底関数が伸縮し、イオンや分極した系を柔軟に記述します。

領域分割近似Fock交換: 従来のTB法では難しかった反応障壁やバンドギャップの精度を向上させるため、近似的な非局所Fock交換項を取り入れています。


原子補正ポテンシャル (ACPs): 基底関数の不足を補うための経験的な補正項です。

2. ベンチマーク性能

論文内では、約32,000点のデータポイントを用いた大規模なベンチマークが行われています。

総合精度 (GMTKN55): 主要なベンチマークセットであるGMTKN55において、g-xTBは WTMAD-2（重み付き平均絶対偏差）で 9.3 kcal/mol を記録しました。これは、前身であるGFN2-xTB (25.0 kcal/mol) から劇的な精度向上であり、低コストなDFT手法（B97-3cなど）に匹敵する数値と報告されています。

遷移金属と反応障壁: GFN2-xTBの弱点とされていた反応障壁の高さ（過小評価される傾向があった）や、遷移金属錯体の構造最適化においても、DFTに近い精度が達成されているとされています。例えば、反応障壁のベンチマーク（BH9）では、GFN2-xTBと比較して平均誤差が半分以下に低減しています。

3. 計算速度

計算コストに関しては、GFN2-xTBと比較して約30%〜50%程度のオーバーヘッド（速度低下）が生じるとされていますが、依然としてDFTと比較すると数桁（1000倍以上）高速であり、ハイスループットスクリーニング等に適した特性を維持しています。


`pygxtb` は、この正式リリース前の新しい計算エンジンを検証のためにPythonエコシステムへ手軽に接続することを目的としています。


### インストール

Pythonのパッケージ管理システムである pip を介してインストール可能です。

```Bash

pip install pygxtb
```

### 使用例

以下は、水分子（H2O）を生成し、pygxtb を計算機として割り当ててポテンシャルエネルギー、力、双極子モーメントを計算するスクリプト例です。

```Python

from ase.build import molecule
from pygxtb import PygxTB

# 水分子の作成
atoms = molecule("H2O")

# Calculatorの割り当て
# パッケージ内のバイナリ/パラメータを自動的に検出します
calc = PygxTB()
atoms.calc = calc

# ポテンシャルエネルギーの計算
energy = atoms.get_potential_energy()
print(f"Energy: {energy:.6f} eV")

# 力の計算
forces = atoms.get_forces()
print("Forces (eV/Ang):")
print(forces)

# 双極子モーメントの取得
dipole = atoms.get_dipole_moment()
print(f"Dipole: {dipole}")
```

### 使用上の注意

（既知の問題）実行パスの長さ制限について 実行ディレクトリの絶対パスが長すぎると、g-xTBのバイナリがパスを正しく認識できずエラーが発生する場合があります。 これは、g-xTB内部（Fortranコード）におけるパス文字列の長さ制限に起因すると推測されます。

回避策： エラーを避けるため、可能な限り**浅い階層（絶対パスの文字数が短いディレクトリ）**でプログラムを実行することを推奨します。

### ライセンスと依存関係

- ライセンス: GNU General Public License v3 (GPLv3)

- 依存関係: ASE (Atomic Simulation Environment) が必須となります。

### 参考
- Preprint: g-xTB: A General-Purpose Extended Tight-Binding Electronic Structure Method For the Elements H to Lr (Z=1–103) https://doi.org/10.26434/chemrxiv-2025-bjxvt

- Source Repository: https://github.com/grimme-lab/g-xtb  (Commit: df43e00)

- Reference Implementation: https://github.com/stefanbringuier/g-xtb-ase/tree/main?tab=readme-ov-file

