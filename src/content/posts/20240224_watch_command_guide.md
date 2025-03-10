---
title:  【Unix】watchの解説
published: 2025-02-24
description: ""
tags: [watch]
category: Unix
draft: false
---
最終更新：2025-02-24

# 概要

`watch`コマンドは、指定したコマンドを一定間隔で繰り返し実行し、その結果を常に表示し続けるためのコマンドである。システムリソースやファイルの変化をリアルタイムで監視したい場合に有用である。

使用方法は以下の通りである。  
```
watch [オプション] [コマンド]
```


## オプション一覧とハンズオン

### 1. -n (更新間隔の指定)

- 機能: コマンドを再実行する間隔を秒単位で指定する。  
- 使用例:
```bash
watch -n 5 ls -l
```
上記では、`ls -l`コマンドが5秒おきに実行され、一覧が更新される。  
- 出力例:
  (5秒ごとに下記のようなディレクトリ情報が更新される)
  ```
  total 12
  -rw-r--r--  1 user group  104 Feb 24 12:34 file1.txt
  -rw-r--r--  1 user group  256 Feb 24 13:00 file2.log
  ...
  ```

### 2. -d (差分表示モード)

- 機能: 直前の結果との違い（差分）をハイライト表示する。  
- 使用例:
```bash
watch -d df -h
```
上記では、ディスク使用量の変化部分がハイライトされる。  
- 出力例:  
  ```
  Filesystem   Size  Used Avail Use% Mounted on
  /dev/sda1     20G   15G    5G  75% /
  ...
  ```
  (変化があれば数字の部分がハイライトされる)

### 3. -t (タイトル行を非表示)

- 機能: 画面上部に表示される更新間隔の情報や時刻の表示を隠す。  
- 使用例:
```bash
watch -t uname -a
```
上記では、`uname -a`の結果だけが表示され、画面上部のヘッダが非表示になる。  
- 出力例:
  ```
  Linux hostname 5.10.0-13-amd64 #1 SMP Debian 5.10.106-1 x86_64 GNU/Linux
  ```

### 4. -x (コマンドと引数を1つの文字列として扱う)

- 機能: watchの引数として、スペースを含むコマンドを1つの引数として扱う際に使用する。  
- 使用例:
```bash
watch -x echo "Hello World"
```
（環境によっては効果が分かりにくい場合がある。通常はコマンドの引数にスペースがある場合、クォートで対応することが多い）

- 出力例:
  ```
  Hello World
  ```

