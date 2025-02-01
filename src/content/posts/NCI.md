---
title: Non-Covalent Interaction(NCI)の解析方法
published: 2025-02-01
tags: [NCIplot]
category: Computational Chemistry
draft: false
---
最終更新：2025-02-01

※2025年1月現在再現可能な方法である。これを備忘録として残す。

## 1, 概要

分子構造の電子状態計算の結果得た分子軌道から、非共有結合相互作用が存在する箇所を可視化する手順を解説する。

## 2, わかること

・電子状態計算の結果得られた系の分子間相互作用が生じている場所がわかる。

・分子間相互作用が生じていると考えられる箇所で、反発しているのか、引き合っているのかが分かる。


## 3, 解析方法の手順

NCIの計算にMultiwfnを用い、可視化には、VMDを使用する。

初めにNCIの計算の手順について解説する。

### MultiwfnのLinuxへの導入方法

a, curlコマンドを実行してソフトウェアをダウンロードする。
```
curl -O http://sobereva.com/multiwfn/misc/Multiwfn_3.8_dev_bin_Linux.zip 
```

（リンク切れであればhttp://sobereva.com/multiwfn/ にブラウザ等でアクセスしてDevelopment version: 3.8(dev)のBinary packageのLinux 64bitをダウンロードする） 

b, unzip等で適当なディレクトリで解凍する。

c, Multiwfn_3.7_bin_Linuxへ移動し、Multiwfnのバイナリファイルが存在することを確認する。

d, Multiwfnを実行できれば導入成功。
（Windowsでも導入できるが手順は似ているので省略する。）

### Multiwfn version: 3.8(dev)を使った手順

a, NCIを計算したい電子状態構造の情報を含むGaussianのchkファイルを用意する。

b, formchkコマンドなどでchkファイル(xxx_MO.rrm)をfchkファイルに変換する。

c, Multiwfnが存在するディレクトリに移動し、以下のコマンドを実行する。
```
./Multiwfn xxx.fchk
```
(成功すると************ Main function menu ************が見える。)

d, 20 Visual study of weak interactionを選択する。（20を入力してEnterを押す。成功すると、============ Visual study of weak interaction ============が見える。）

e, 1 NCI analysis (also known as RDG analysis. JACS, 132, 6498)を選択する。(NCIplotと同じ手法)

f, Please select a method to set up gridが表示されるので、適当に選択する。

g, 処理が終わるまでしばらく待ち、 
3 Output cube files to func1.cub and func2.cub in current folder
を選択する。(func1.cubとfunc2.cubが生成されていれば成功)

（こちらの動画で(https://www.youtube.com/watch?v=Z9xN1-o7OUY)使い方が説明されている。）

NCIを計算した後、VMDを用いて可視化する。
### VMDの導入方法(windows11)

a, 次のサイトにアクセスする。https://www.ks.uiuc.edu/Research/vmd/

b, VMD ver.1.9.3をダウンロードする。

c, https://github.com/stecue/gMultiwfn/blob/master/examples/RDGfill.vmd
から、RDGfill.vmdをダウンロードする。

d, VMD.exeが存在するディレクトリにRDGfill.vmdを置く。

### 可視化方法(windows11)

a, VMD.exeが存在するディレクトリに先ほど生成したfunc1.cubとfunc2.cubを置く。

b, VMD.exeを実行する。 

c, コマンドプロンプトに”source RDGfill.vmd”と入力する。

d, VMD x.x.x OpenGL Displayが表示されれば成功。

※注目したい相互作用が見えにくいときは、VMD mainのGraphics > Presentationsをクリックし、Isovalueの値を調節することで、見やすくできる。これは、可視化する電子密度の下限の閾値を調節している。値が大きいほど電子密度が小さい箇所も可視化する。

### NCIが可視化されたときに表示されている色についてのメモ：

赤＝分子構造間で電子反発がある。

緑＝vDW相互作用等、弱い相互作用が存在する。

青＝水素結合等の強い電子誘引効果が存在する。


#### 参考:
- _J. Am. Chem. Soc._ 2010, 132 (18), 6498–6506. DOI:10.1021/ja100936w. (NCIplotの論文1)
- _J. Chem. Inf. Model._ 2020, 60 (1), 6–10. DOI:10.1021/acs.jcim.9b00950. (NCIplotの論文2)
- _The Journal of Chemical Physics_ 2024, 161 (8), 082503. DOI: 10.1063/5.0216272. (Multiwfnの最新の論文)
- http://sobereva.com/multiwfn/misc/Multiwfn_3.8_dev.pdf (Multiwfnのマニュアルのダウンロードリンク)
- https://computational-chemistry.com/top/blog/2017/04/04/non-covalent-interaction/ (日本語のサイト。参考文献で元論文をたどれる。)
- _Adv. Synth. Catal._ 2022, 364 (14), 2333–2339. DOI:10.1002/adsc.202200327.（活用例1）
- _J. Am. Chem. Soc._ 2022, 144 (2), 798–806. DOI:10.1021/jacs.1c09889. （活用例2）


