---
title:  【計算化学】ASEでAFIR法を使用する
published: 2025-03-04
description: "ASEを使った人工力誘起反応法の実装について解説"
tags: [AFIR, ASE, Python, Quantum Chemistry, Reaction Coordinate]
category: Computational Chemistry
draft: false
---
最終更新：2025-03-06


ASE_AFIR（Atomic Simulation Environment - Artificial Force Induced Reaction）と呼ばれるパッケージについて説明する。これは化学反応経路の自動探索を可能にする人工力誘起反応（AFIR）法をPythonの原子シミュレーション環境（ASE）に実装した自作のパッケージである。このパッケージでは、量子化学計算と分子動力学シミュレーションを組み合わせた反応経路探索機能を提供している。


Repositoryはこちら：https://github.com/ss0832/ASE_AFIR

## AFIR法とは

AFIR法は北海道大学の前田理教授らによって開発された計算化学手法である。この手法により、分子間に人工的な力を加えて適切に構造最適化すると、反応系から生成系までの情報を含むtrajectoryが得られる。これをAFIR経路と呼ぶ。

このAFIR経路に対して、必要に応じてNEB(Nudged Elastic Band)法などの経路最適化手法を用いて経路を緩和する。

最終的に得られる経路のエネルギー極大値をもつ分子構造から、良質な近似遷移状態構造を得ることが可能である。これを用いて遷移状態及びIRCを求める。これにより、化学反応の遷移状態や反応経路を効率的に探索可能である。

緩和スキャンといった従来の方法では見つけることが難しい反応経路も、人工力を導入することで、発見が容易になるという利点がある。

また、AFIR関数のガンマの値を調節することで、超えうる活性化障壁の高さを制限して、反応経路の探索が可能である。

例えば、（化学反応ではないが）ガンマの値を50kJ/mol程度にすることで、立体配座探索を行うことが可能である。

具体的には、反応に関わる分子フラグメント間に引力または斥力を加えることで、エネルギー障壁を越え、近似反応経路を効率よく探索する。この方法は特に、遷移状態の予測が難しい複雑な反応系で威力を発揮する。

## ASE_AFIRの特徴

- ASEの最適化アルゴリズムとの統合
- 様々な量子化学計算ソフト（Psi4、ORCA、NWChemなど）との互換性
- AFIR法のパラメータ（力の強さ、原子選択など）の制御
- 引力型・斥力型両方の人工力のサポート


## インストールと必要環境
ASE_AFIRはPython環境で動作し、主要な依存パッケージはASE、NumPy、PyTorch、Matplotlibである。インストールは以下のコマンドで行う。：

```bash
git clone https://github.com/ss0832/ASE_AFIR.git
cd ASE_AFIR
pip install -r requirements.txt
```
Psi4などの量子化学計算ソフトウェアは別途インストールが必要です：

```bash
conda install -c psi4 psi4
```

## 実装の技術的詳細

ASE_AFIRはPyTorchを用いた自動微分機能を活用し、AFIRポテンシャルの効率的な計算を実現している。これにより、量子化学計算との組み合わせによる反応経路探索が可能である。コードは完全にPythonで実装されており、モジュール構造によって拡張性も確保されている。

AFIRのエネルギー項は次の式で表される：

$$E_{AFIR} = \alpha \frac{\sum_{i \in A}\sum_{j \in B} \omega_{ij}r_{ij}}{\sum_{i \in A}\sum_{j \in B} \omega_{ij}}$$

ここで、$\alpha$はAFIR力の強さを制御するパラメータ、$\omega_{ij}$は原子間の距離に依存する重み、$r_{ij}$は原子間距離である。ASE_AFIRでは、このエネルギー項を量子化学計算のエネルギーに加えることで、人工力を導入している。

## 応用例：Claisen転移反応

反応物にAFIR関数を加えることでClaisen転移の近似反応経路を得られる。

```python
# Claisen転位のASE_AFIR実装例
import numpy as np
from ase import Atoms
from ase.calculators.psi4 import Psi4
from ase_afir_calc import ASE_AFIR

# 分子構造の定義
symbols = ['C', 'O', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']
positions = np.array([
        [-0.03989051,   -1.44876143,   -4.10172098],  # 1: C
        [-0.24675417,   -0.16352192,   -4.69355444],  # 2: O
        [ 0.99900999,    0.34386140,   -5.17887270],  # 3: C
        [ 1.99215654,    0.47134851,   -4.00882699],  # 4: C
        [ 3.32183313,    0.33200893,   -4.23043860],  # 5: C
        [-0.83456106,   -1.86233332,   -3.08485125],  # 6: C
        [ 0.74232880,   -2.08390788,   -4.46175144],  # 7: H
        [-0.67977551,   -2.82401638,   -2.64201141],  # 8: H
        [-1.61677881,   -1.22718586,   -2.72481917],  # 9: H
        [ 1.39589792,   -0.32674955,   -5.91212635],  # 10: H
        [ 0.84422389,    1.30554411,   -5.62171311],  # 11: H
        [ 1.63235197,    0.66994296,   -3.02089940],  # 12: H
        [ 3.68163775,    0.13341526,   -5.21836632],  # 13: H
        [ 4.01187671,    0.42059024,   -3.41748521]   # 14: H
])
atoms = Atoms(symbols=symbols, positions=positions)

# 量子化学計算の設定
qm_calc = Psi4(method='b3lyp', basis='6-31g')

# AFIR計算の設定 - C=C結合形成のために人工力を適用
afir_calc = ASE_AFIR(
    calculator=qm_calc,
    AFIR_gamma=150.0,  # kJ/mol 
    AFIR_Fragm_1=[5, 7, 8],  # ビニル基のC原子とH原子(原子を0からラベリングしている, ASEのフォーマットに合わせている)
    AFIR_Fragm_2=[4, 12, 13],  # アリル基のC原子とH原子(原子を0からラベリングしている)
    p=6.0
)

# 計算機の割り当て
atoms.calc = afir_calc

# 構造最適化の実行
from ase.optimize import BFGS
opt = BFGS(atoms, trajectory='claisen_reaction.traj')
opt.run(fmax=0.02, steps=200)

# 反応軌跡の解析
from analyze_trajectory import analyze_trajectory
analyze_trajectory()

```
上記のコードでは、アリルビニルエーテルのClaisen転位反応を、AFIRを用いてシミュレーションしている。ビニル基（6: C）とアリル基（5: C）の間に人工引力を加えることで、C-C結合形成と同時に進行するC-O結合切断を含む反応経路を自動的に探索できる（近似反応経路が得られる）。

## 反応軌跡の解析
反応軌跡の解析には、専用のスクリプトが用意されている：

```Python
# analyze_trajectory.pyの主な機能
# 1. 重要な結合距離（C-C形成、C-O切断など）の追跡
# 2. エネルギープロファイルのプロット
# 3. 反応過程の主要構造の3D表示
# 4. 全軌跡のXYZファイル出力
```
解析の結果、結合距離の変化や反応エネルギープロファイルのグラフなどが得られる。

得られたAFIR経路のエネルギー極大値が遷移状態構造を求めるための候補となる。さらに、NEB(Nudged Elastic Band)法もしくはLUP(Locally Plains Update)法により、AFIR経路を緩和することで、より正確な遷移状態構造に近い候補を見つけることが可能である。

以下にNEB法による経路の緩和を行うプログラムを示す。
```python
"""
XYZファイルから選択されたフレームに対してNEB計算を実行するスクリプト
電子状態計算にはASEとPsi4を使用します。
NEBの全最適化iterationのエネルギープロファイルを出力します。
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import BFGS
from ase.calculators.psi4 import Psi4


# パラメータ
XYZ_FILE = 'claisen_reaction_animation.xyz'
NUM_IMAGES = 30  # NEBに使用する画像数
MAX_NEB_ITERATIONS = 10
OUTPUT_DIR = 'neb_results'
TS_DIR = 'ts_optimization'
PROFILES_DIR = f'{OUTPUT_DIR}/energy_profiles'  # エネルギープロファイル保存ディレクトリ
eV2kcalmol = 23.06031

# 出力ディレクトリの作成（存在しない場合）
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TS_DIR, exist_ok=True)
os.makedirs(PROFILES_DIR, exist_ok=True)

# ステップ1: XYZファイルからすべての構造を読み込む
print("Reading structures from XYZ file...")
all_structures = read(XYZ_FILE, index=':')
total_structures = len(all_structures)
print(f"Found {total_structures} structures in the file.")

# ステップ2: NEB経路用にNUM_IMAGES個の等間隔構造を選択
if total_structures <= NUM_IMAGES:
    images = all_structures
    print("Using all available structures as there are fewer than requested.")
else:
    # 等間隔で構造を選択するためのインデックスを計算
    indices = np.round(np.linspace(0, total_structures - 1, NUM_IMAGES)).astype(int)
    images = [all_structures[i] for i in indices]
    print(f"Selected {NUM_IMAGES} images with indices: {indices}")

# ステップ3: B3LYP/6-31GでのPsi4計算機設定
psi4_template = {
    'basis': '6-31G',  # 基底関数
    'method': 'B3LYP',     # 交換相関汎関数
    'reference': 'RHF',  # 参照波動関数
    'scf_type': 'pk',  # SCF収束のためのアルゴリズムの指定
    'maxiter': 100,    # SCF最大iteration回数
    'print': 1,        # 出力レベル
    'mem': '2GB',       # メモリ割り当て
    'num_threads': 8,       # 使用するCPUスレッド数
    'psi4_options': {
        'nthreads': 8, }
}

# ステップ4: NEB計算のセットアップ
print("Setting up NEB calculation...")
for i, image in enumerate(images):
    calc = Psi4(**psi4_template)
    image.calc = calc
    # 初期構造を保存
    write(f'{OUTPUT_DIR}/initial_image_{i:02d}.xyz', image)

# NEBインスタンスの作成
neb = NEB(images, k=0.0)  # 標準NEB
### AFIR経路を内挿補間する必要はない

# 各iterationでのエネルギープロファイル保存用関数
iteration_count = [0]  # 関数内で変更できるようリストを使用

def save_energy_profile():
    """各iterationでエネルギープロファイルを保存"""
    current_iter = iteration_count[0]
    
    # エネルギー計算
    energies = [image.get_potential_energy() for image in images]
    min_energy = min(energies)
    rel_energies = [(e - min_energy) * eV2kcalmol for e in energies]
    
    # データをファイルに保存
    with open(f'{PROFILES_DIR}/energy_profile_iter_{current_iter:03d}.dat', 'w') as f:
        f.write("# Image Index, Energy (eV), Relative Energy (kcal/mol)\n")
        for i, (e, rel_e) in enumerate(zip(energies, rel_energies)):
            f.write(f"{i} {e} {rel_e}\n")
    
    # プロットの作成
    plt.figure(figsize=(10, 6))
    plt.plot(range(len(images)), rel_energies, 'o-')
    plt.xlabel('Image Index')
    plt.ylabel('Energy (kcal/mol)')
    plt.title(f'NEB Energy Profile - Iteration {current_iter}')
    plt.grid(True)
    plt.savefig(f'{PROFILES_DIR}/energy_profile_iter_{current_iter:03d}.png')
    plt.close()
    
    # 現在の経路のコピーを保存
    write(f'{PROFILES_DIR}/neb_path_iter_{current_iter:03d}.xyz', images)
    
    # iterationカウンターをインクリメント
    iteration_count[0] += 1
    
    # 現在の最大エネルギーポイントを出力
    max_idx = np.argmax(rel_energies)
    print(f"Iteration {current_iter}: Max energy at image {max_idx} = {rel_energies[max_idx]:.2f} kcal/mol")

# 初期エネルギープロファイルを保存（iteration0）
save_energy_profile()

# ステップ5: NEB最適化の実行
print(f"Starting NEB optimization for {MAX_NEB_ITERATIONS} iterations...")
optimizer = BFGS(neb, trajectory=f'{OUTPUT_DIR}/neb_trajectory.traj')

# エネルギープロファイル保存関数をオプティマイザに接続
optimizer.attach(save_energy_profile)

# オプティマイザの実行
optimizer.run(fmax=0.05, steps=MAX_NEB_ITERATIONS)  # fmaxは力の許容値

# ステップ6: 最終結果の分析
print("NEB optimization completed. Creating final analysis...")

# 最終NEB経路を保存
write(f'{OUTPUT_DIR}/final_neb_path.xyz', images)

# 分析用の最終エネルギーを計算
final_energies = [image.get_potential_energy() for image in images]
min_energy = min(final_energies)
final_rel_energies = [(e - min_energy) * eV2kcalmol for e in final_energies]

# 視覚性を向上させた最終エネルギープロファイルプロットの作成
plt.figure(figsize=(10, 6))
plt.plot(range(len(images)), final_rel_energies, 'o-', linewidth=2, color='blue')
plt.xlabel('Image Index', fontsize=12)
plt.ylabel('Energy (kcal/mol)', fontsize=12)
plt.title('Final NEB Energy Profile', fontsize=14)
plt.grid(True)

# 反応物、遷移状態、生成物のマーカーを追加
max_idx = np.argmax(final_rel_energies)
plt.scatter([0, max_idx, len(images)-1], 
            [final_rel_energies[0], final_rel_energies[max_idx], final_rel_energies[-1]], 
            color=['green', 'red', 'purple'], 
            s=100, zorder=5)
plt.annotate('Reactant', (0, final_rel_energies[0]), textcoords="offset points", 
             xytext=(0,10), ha='center', fontsize=10)
plt.annotate('TS', (max_idx, final_rel_energies[max_idx]), textcoords="offset points", 
             xytext=(0,10), ha='center', fontsize=10)
plt.annotate('Product', (len(images)-1, final_rel_energies[-1]), textcoords="offset points", 
             xytext=(0,10), ha='center', fontsize=10)

plt.savefig(f'{OUTPUT_DIR}/final_energy_profile.png', dpi=300, bbox_inches='tight')
print(f"Final energy profile saved to {OUTPUT_DIR}/final_energy_profile.png")

# 最終エネルギーデータをファイルに保存
with open(f'{OUTPUT_DIR}/final_energy_profile.dat', 'w') as f:
    f.write("# Image Index, Energy (eV), Relative Energy (kcal/mol)\n")
    for i, (e, rel_e) in enumerate(zip(final_energies, final_rel_energies)):
        f.write(f"{i} {e} {rel_e}\n")

# ステップ7: エネルギープロファイル上のエネルギー最大値の特定
max_energy_idx = np.argmax(final_rel_energies)
max_energy = final_rel_energies[max_energy_idx]
print(f"Identified global maximum at image {max_energy_idx} with relative energy {max_energy:.2f} kcal/mol")

# 遷移状態最適化用に画像のコピーを作成
ts_guess = images[max_energy_idx].copy()
ts_guess.calc = Psi4(**psi4_template)

# 近似遷移状態構造を保存
write(f'{TS_DIR}/ts_guess.xyz', ts_guess)

```
これにより、AFIR経路から得られる遷移状態構造の候補をより正確な遷移状態構造に近づけることが可能である。

Repositoryはこちら：https://github.com/ss0832/ASE_AFIR



## 参考文献
- S. Maeda, K. Ohno, K. Morokuma, "Systematic exploration of the mechanism of chemical reactions: the global reaction route mapping (GRRM) strategy using the ADDF and AFIR methods", Phys. Chem. Chem. Phys., 2013, 15, 3683-3701.
- S. Maeda, K. Morokuma, "Finding Reaction Pathways of Type A + B → X: Toward Systematic Prediction of Reaction Mechanisms", J. Chem. Theory Comput., 2011, 7, 2335-2345.
- A. H. Larsen et al., "The atomic simulation environment—a Python library for working with atoms", J. Phys.: Condens. Matter, 2017, 29, 273002.
- R. M. Parrish et al., "Psi4 1.1: An Open-Source Electronic Structure Program Emphasizing Automation, Advanced Libraries, and Interoperability", J. Chem. Theory Comput., 2017, 13, 3185-3197.

### 追記・修正

2025/3/5 プログラムの実装ミスを修正。analyze_trajectory()で得られるエネルギーに人工力ポテンシャルが含まれるため、人工力なしのポテンシャルを取り出せるように実装しなおす予定である。
2025/3/6 NEB法を実行するプログラムの追加
