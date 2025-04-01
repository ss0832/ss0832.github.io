---
title:  【Psi4】利用可能なDFT汎関数の一覧
published: 2025-04-01
description: "Psi4で利用できるDFT汎関数の種類と特徴の解説"
tags: [Psi4]
category: Computational Chemistry
draft: false
---
最終更新：2025-04-01

## 概要

Psi4は様々な密度汎関数理論（DFT）汎関数を実装している。これらの汎関数は交換汎関数と相関汎関数の組み合わせによって構成され、異なる分子系や物性の計算に対して異なる精度と計算コストを持っている。本記事ではPsi4で利用可能なDFT汎関数の種類を分類して紹介する。

## LDA汎関数（局所密度近似）

最も単純なDFT汎関数であり、一様電子ガスに基づく近似を用いる。

- `SVWN`: Slater交換 + VWN5相関
- `SPW92`: Slater交換 + PW92相関

## GGA汎関数（一般化勾配近似）

電子密度に加えてその勾配も考慮する汎関数。

- `PBE`: PBE交換 + PBE相関
- `BLYP`: B88交換 + LYP相関
- `BP86`: B88交換 + P86相関
- `PW91`: PW91交換 + PW91相関
- `B97-D`: Grimmeの分散補正付きB97
- `SOGGA`: 第二次勾配近似
- `SOGGA11`: SOGGA汎関数の改良版
- `GAM`: GAM交換相関汎関数
- `OP`: OPカップリング補正
- `N12`: ミネソタN12汎関数
- `FT97`: FT97交換相関汎関数

## メタGGA汎関数

電子密度、その勾配に加えて、運動エネルギー密度も考慮する汎関数。

- `M06-L`: ミネソタM06-L汎関数
- `M11-L`: ミネソタM11-L汎関数
- `MN12-L`: ミネソタMN12-L汎関数
- `M08-HX`: ミネソタM08-HX汎関数
- `TPSS`: TPSS交換相関汎関数
- `revTPSS`: 改訂版TPSS
- `B97M-V`: B97ベースメタGGA + VV10非局所相関
- `B97M-D3BJ`: B97ベースメタGGA + D3(BJ)分散補正

## ハイブリッドGGA汎関数

厳密な交換（HF交換）をある割合で混合したGGA汎関数。

- `PBE0`: PBE汎関数に25%のHF交換を混合
- `B3LYP`: 最も広く使われるハイブリッド汎関数（20%のHF交換を含む）
- `B3LYP5`: 別バージョンのB3LYP（VWN5相関を使用）
- `B97`: Beckeの1997年ハイブリッド汎関数
- `SOGGA11-X`: SOGGA11のハイブリッド版
- `PW6B95`: PW6交換 + B95相関のハイブリッド
- `B97-1`: B97の改良版
- `B97-2`: さらに改良されたB97
- `B97-K`: 運動論特性向けに最適化されたB97
- `B98`: B97の改良版
- `HSE-HJS`: Heyd-Scuseria-Ernzerhofの交換相関汎関数
- `wB97`: ロングレンジ補正付きB97
- `wB97X`: HF交換成分を増やした長距離補正B97
- `wB97X-D`: 分散補正付きwB97X
- `wB97X-V`: VV10非局所相関付きwB97X
- `wB97M-V`: メタGGA成分を持つwB97汎関数

## ハイブリッドメタGGA汎関数

HF交換とメタGGA汎関数を組み合わせたもの。

- `M05`: ミネソタM05汎関数
- `M05-2X`: 2倍のHF交換を含むM05
- `M06`: ミネソタM06汎関数
- `M06-2X`: 2倍のHF交換を含むM06
- `M06-HF`: 100%のHF交換を含むM06
- `M08-SO`: ミネソタM08-SO汎関数
- `M11`: ミネソタM11汎関数（レンジセパレート型）
- `MN15`: ミネソタMN15汎関数
- `TPSSh`: TPSS汎関数に10%のHF交換を混合
- `revTPSSh`: 改訂版TPSSのハイブリッド
- `PW6B95`: Zhao-Truhlar汎関数
- `PWB6K`: PW6B95の修正版

## 二重ハイブリッド汎関数

HF交換に加えてMP2型の相関も考慮する高精度汎関数。

- `B2PLYP`: B88交換とLYP相関の二重ハイブリッド
- `DSD-PBEP86`: P86相関を使用したスピン成分依存二重ハイブリッド
- `PBE0-DH`: PBE0ベースの二重ハイブリッド
- `PBE0-2`: PBE0の別バージョン二重ハイブリッド

## 分散補正

多くの汎関数はファンデルワールス分散力を正確に記述できないため、以下の分散補正を追加できる。

- `-D`: Grimmeの最初の分散補正（D1）
- `-D2`: Grimmeの第二世代分散補正
- `-D3`: Grimmeの第三世代分散補正
- `-D3(BJ)`: Becke-Johnson減衰関数を用いたD3補正
- `-D3M`: より最適化されたD3パラメータ
- `-D3M(BJ)`: Becke-Johnson減衰関数を用いた修正D3

## 非局所汎関数

非局所相関を含む汎関数。

- `-NL`: 非局所相関補正
- `-VV10`: Vydrov-Van Voorhisの非局所相関

## 使用例

Psi4でDFT計算を実行するには以下のように汎関数名を指定する：

```python
# B3LYP汎関数を使用した計算の例
import psi4

mol = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

# B3LYP/cc-pVDZレベルでのエネルギー計算
energy = psi4.energy('b3lyp/cc-pVDZ')
print(f"B3LYP/cc-pVDZ Energy: {energy} Hartree")

# PBE0-D3汎関数を使用した場合
energy = psi4.energy('pbe0-d3/aug-cc-pVTZ')
print(f"PBE0-D3/aug-cc-pVTZ Energy: {energy} Hartree")
```

## 参考サイト

- [Psi4マニュアル: DFT By Functional](https://psicode.org/psi4manual/master/dft_byfunctional.html)
- [Psi4マニュアル: DFT](https://psicode.org/psi4manual/master/dft.html)
- [Psi4公式GitHub](https://github.com/psi4/psi4)