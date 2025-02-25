---
title:  【Unix】timeの解説
published: 2025-02-25
description: ""
tags: [time]
category: Unix
draft: false
---
最終更新：2025-02-25

# 概要

`time`コマンドは、指定したコマンドの実行に要した実行時間やCPU使用率を計測するためのコマンドである。プログラムのパフォーマンスや処理効率を確認したい場合に有用である。

使用方法は以下の通りである。  
```
time [オプション] コマンド
```
## オプション一覧
  
- `-p`: POSIX形式（秒単位）で結果を簡潔に表示できる。  
- `-v`: 詳細な統計情報を得ることができる。(GNU版, 内訳は実行して確かめてみるとよい。)  

## オプションとハンズオン

### (オプションなし)

・機能: デフォルトでは、実行したコマンドの「実行に要した実時間（real）」「ユーザモードCPU時間（user）」「カーネルモードCPU時間（sys）」を表示する。  
・使用例:  
```bash
time ls
```
・出力例（環境により異なる）:  
```
file1.txt  file2.txt
real    0m0.003s
user    0m0.001s
sys     0m0.002s
```

### -p オプション (POSIX形式での出力)

・機能: POSIX形式の簡潔な出力を行う。`real`, `user`, `sys`の数値を秒数のみで表示する。  
・使用例:  
```bash
time -p sleep 2
```
・出力例:  
```
real 2.00
user 0.00
sys  0.00
```
POSIX形式では、経過時間・ユーザCPU時間・システムCPU時間がそれぞれ「real」「user」「sys」のキーで秒数表示されることが特徴である。

### GNU timeコマンド拡張例

実装によっては、Linuxなどで以下のような拡張オプションが利用できる。

#### -v オプション (詳細統計の表示)

・機能: コマンド実行時のページフォルト数やCPU占有率など、詳細な情報を表示する。`GNU time`で有効な場合が多い。  
・使用例:  
```bash
/usr/bin/time -v ls
```
・出力例:  
```
Command being timed: "ls"
User time (seconds): 0.00
System time (seconds): 0.00
Percent of CPU this job got: 25%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.00
...
```

