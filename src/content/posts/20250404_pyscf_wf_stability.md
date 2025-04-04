---
title:  【PySCF】波動関数の安定性解析
published: 2025-04-04
description: "PySCFでの波動関数安定性の解析方法"
tags: [pyscf]
category: Computational Chemistry
draft: false
---
最終更新：2025-04-04

# 概要
PySCFにおいて、UHFやUKSで得られた波動関数が安定解であるかどうかをチェックする方法を示す。内部安定性と外部安定性の両方を評価し、必要に応じて新たな軌道を求めることが可能である。

## 波動関数の安定性解析の基本
PySCFでは、以下を満たすことで波動関数の安定性を評価できる。  
- 内部安定性 (internal=True)  
- 外部安定性 (external=True)  

解析の結果、不安定な場合には軌道が再計算される。

## コード例
下記に、PySCFでUHF計算を行った後に、内部安定性と外部安定性を解析する一連の流れを示す。SCF計算→安定性解析→再計算の順で進める。

```python
import numpy as np
from pyscf import gto, scf

# 分子の定義
mol = gto.M(
    atom="O 0 0 0; O 0 0 1.2",
    basis="cc-pvtz",
    spin=2,
    charge=0,
    verbose=4
)

# UHF計算の実行
mf = scf.UHF(mol)
mf.kernel()

# 安定性解析（内部+外部）
# return_status=True で安定性フラグを取得し、内部・外部それぞれを判定する
internal_orb, external_orb, stable_i, stable_e = mf.stability(
    internal=True, 
    external=True,
    return_status=True
)

print(f"内部安定性: {'安定' if stable_i else '不安定'}")
print(f"外部安定性: {'安定' if stable_e else '不安定'}")

# 内部安定性が不十分なら、内部安定化で得られた軌道で再計算
if not stable_i:
    print("内部安定性を改善する...")
    dm = mf.make_rdm1(internal_orb, mf.mo_occ)  # 新しい密度行列を生成
    mf_new = mf.copy()  # SCFオブジェクトをコピー
    mf_new.kernel(dm)  # 新しい密度行列から再スタート
    print(f"改善後のエネルギー: {mf_new.e_tot}")

# 外部安定性が不十分なら、外部安定化で得られた軌道で再計算
elif not stable_e:
    print("外部安定性を改善する...")
    dm = mf.make_rdm1(external_orb, mf.mo_occ)  # 新しい密度行列を生成
    mf_new = mf.copy()  # SCFオブジェクトをコピー
    mf_new.kernel(dm)  # 新しい密度行列から再スタート
    print(f"改善後のエネルギー: {mf_new.e_tot}")

# 内部と外部の両方が不安定の場合
else:
    print("内部安定性を改善する...")
    dm = mf.make_rdm1(internal_orb, mf.mo_occ)  # 新しい密度行列を生成
    mf_new = mf.copy()  # SCFオブジェクトをコピー
    mf_new.kernel(dm)  # 新しい密度行列から再スタート
    print(f"改善後のエネルギー: {mf_new.e_tot}")  
    _, external_orb, _, stable_e = mf_new.stability(
    internal=False, 
    external=True,
    return_status=True)
    
    print("外部安定性を改善する...")
    dm = mf.make_rdm1(external_orb, mf_new.mo_occ)  # 新しい密度行列を生成
    mf_new = mf.copy()  # SCFオブジェクトをコピー
    mf_new.kernel(dm)  # 新しい密度行列から再スタート
    print(f"改善後のエネルギー: {mf_new.e_tot}")

# 再チェック（オプション）:
_, _, stable_i_2, stable_e_2 = mf.stability(
    internal=True, 
    external=True,
    return_status=True
)
print(f"再チェック後の内部安定性: {'安定' if stable_i_2 else '不安定'}")
print(f"再チェック後の外部安定性: {'安定' if stable_e_2 else '不安定'}")
```

## ケーススタディ

1. 内部安定性のみチェック（`internal=True, external=False`）  
   ```python
   internal_orb, _, stable_i, _ = mf.stability(internal=True, external=False, return_status=True)
   print("内部安定性:", "安定" if stable_i else "不安定")
   ```

2. 外部安定性のみチェック（`internal=False, external=True`）  
   ```python
   _, external_orb, _, stable_e = mf.stability(internal=False, external=True, return_status=True)
   print("外部安定性:", "安定" if stable_e else "不安定")
   ```

3. 内部と外部を同時にチェック（`internal=True, external=True`）  
   ```python
   internal_orb, external_orb, stable_i, stable_e = mf.stability(
       internal=True, 
       external=True,
       return_status=True
   )
   print("内部:", "安定" if stable_i else "不安定", "外部:", "安定" if stable_e else "不安定")
   ```


## 参考サイト
- [PySCF公式リポジトリ](https://github.com/pyscf/pyscf)
- [PySCF公式ドキュメント](https://pyscf.org)

