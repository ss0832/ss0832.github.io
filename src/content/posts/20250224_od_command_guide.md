---
title:  【Unix】odの解説
published: 2025-02-24
description: ""
tags: [od]
category: Unix
draft: false
---
最終更新：2025-02-24

# 概要

`od`コマンドは、バイナリファイルや任意のファイルの内容を、8進数、16進数、10進数などの形式で表示するためのツールである。ファイルの内部構造を詳しく調査する際に有用である。

使用方法は以下の通りである。
```
od [オプション] [ファイル]
```
# オプション一覧
  
- `-c`      :可読な文字を表示
- `-b`      :8進数表示
- `-d`, `-x`: 2バイト単位の10進数および16進数表示 
- `-t`      : 細かい出力形式の指定
- `-j`, `-N`: 出力範囲を制限



# オプションとハンズオン

以下、主要なオプションについて解説する。

## 1. -c オプション（文字表示）

・機能: ファイルの内容を、対応する文字（印字可能な文字およびエスケープシーケンス）として表示する。  
・使用例:
```bash
od -c sample.txt
```
・出力例:
```
0000000   H   e   l   l   o       W   o   r   l   d  \n
```

## 2. -b オプション（8進数バイト表示）

・機能: 各バイトを8進数形式で表示する。  
・使用例:
```bash
od -b sample.txt
```
・出力例:
```
0000000 110 145 154 154 157 040 127 157 162 154 144 012
```

## 3. -d オプション（10進数表示：符号なし2バイト単位）

・機能: 2バイト単位で符号なし10進数として表示する。  
・使用例:
```bash
od -d sample.txt
```
・出力例:
```
0000000 18533 22316 277
```

## 4. -x オプション（16進数表示：2バイト単位）

・機能: 2バイト単位で16進数として表示する。  
・使用例:
```bash
od -x sample.txt
```
・出力例:
```
0000000 4865 6c6c 6f20 576f 726c 640a
```

## 5. -t オプション（出力形式の指定）

・機能: `-t`オプションにより、出力形式を柔軟に指定できる。例えば、1バイトごとの16進数表示は「x1」と指定する。  
・使用例:
```bash
od -t x1 sample.txt
```
・出力例:
```
0000000 48 65 6c 6c 6f 20 57 6f 72 6c 64 0a
```

## 6. -j オプション（スキップ開始位置の指定）

・機能: 指定したバイト数分をスキップして出力を開始する。  
・使用例:
```bash
od -j 5 -c sample.txt
```
・出力例:
```
0000005   o       W   o   r   l   d  \n
```

## 7. -N オプション（最大出力バイト数の指定）

・機能: 最大で指定したバイト数のみ出力する。  
・使用例:
```bash
od -N 10 -c sample.txt
```
・出力例:
```
0000000   H   e   l   l   o       W   o   r
```
