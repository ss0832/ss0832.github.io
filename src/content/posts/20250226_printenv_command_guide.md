---
title:  【Unix】printenvの解説
published: 2025-02-26
description: ""
tags: [printenv]
category: Unix
draft: false
---
最終更新：2025-02-26

# 概要

`printenv`コマンドは、現在のシェル環境変数を表示するためのコマンドである。引数が指定されなければ、すべての環境変数を一覧表示し、特定の環境変数名を引数に与えると、その環境変数の値のみを表示する。

## コマンドの使用方法

```
printenv [オプション] [環境変数名...]
```

## オプション一覧と詳細

`printenv`には多くのオプションがない場合が多いが、以下は代表的な例である。  
なお、環境によりオプションの有無が異なる場合がある。

- (オプション未指定)  
  すべての環境変数とその値を一覧表示する。  
- variable  
  指定された変数の値のみを表示する。複数指定することもできる。  
- `-0` (ゼロ終端)  
  出力を改行ではなくヌル文字（\0）で区切る。スクリプトやパイプで処理する際に使用する場合がある。  
- `--help`  
  簡単なヘルプメッセージを表示する。  
- `--version`  
  バージョン情報を表示する。

---

## ハンズオンと出力例

### 1. すべての環境変数を表示
```bash
printenv
```
出力例:
```
PATH=/usr/local/bin:/usr/bin:/bin
HOME=/home/user
SHELL=/bin/bash
LANG=C.UTF-8
...
```
（環境によって異なる。）

### 2. 特定の環境変数（例：HOME）を表示
```bash
printenv HOME
```
出力例:
```
/home/user
```

### 3. 複数の環境変数を指定
```bash
printenv PATH SHELL
```
出力例:
```
/usr/local/bin:/usr/bin:/bin
/bin/bash
```

### 4. `-0`オプション（ゼロ終端で変数を表示）
```bash
printenv -0
```
出力例:
```
PATH=/usr/local/bin:/usr/bin:/bin\0HOME=/home/user\0SHELL=/bin/bash\0...
```
改行の代わりにヌル文字 `\0` が挿入される。

### 5. `--help`オプション
```bash
printenv --help
```
出力例（環境依存）:
```
Usage: printenv [OPTION]... [NAME]...
Print the values of the specified environment VARIABLE(s).  If none specified, print
...
```

### 6. `--version`オプション
```bash
printenv --version
```
出力例（環境依存）:
```
printenv (GNU coreutils) 8.30
Copyright (C) ...
```
