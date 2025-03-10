---
title:  【Unix】rloginの解説
published: 2025-02-27
description: ""
tags: [rlogin]
category: Unix
draft: false
---
最終更新：2025-02-27

# 概要

`rlogin`（`rsh`）コマンドは、リモートホストへのログインを行うためのコマンドである。ネットワークを介して遠隔地のシステムにログインし、端末操作を可能にする。
このコマンドでは、平文で通信を行うため、現在の実運用ではSSHなどの暗号化された手段が広く利用されている点に注意が必要である。

## コマンドの使用方法

```
rlogin [オプション] リモートホスト名
```

基本的に、指定したリモートホストへ接続を試行し、ログインを行う。

---

## オプション一覧と詳細

- `-l ユーザ名`  
  ログインする際に使用するユーザ名を指定する。  
- `-e 文字`  
  エスケープ文字を変更する。デフォルトは`~`である。  
- `-8`  
  8ビットクリーンモードで動作し、文字コードに関わる制限を緩和する。  
- `-L`  
  出力を端末に表示するのではなく行バッファ（行ごと）ずつ処理する。  
- `-a`  
  自動ログインを行う。.rhostsファイルを利用してパスワードなしログインを試みる。  

---

## ハンズオンと出力例

### 1. 通常のリモートログイン
```bash
rlogin remotehost
```
出力例:
```
Password:
Last login: Thu Feb 27 12:00:00 from localmachine
[remotehost ~]$
```
リモートマシンにログインできると、通常のシェルプロンプトが表示される。

### 2. ユーザを指定してログイン（-l）
```bash
rlogin -l username remotehost
```
出力例:
```
Password for username:
[remotehost ~]$
```
指定したユーザ名でログインできる。

### 3. エスケープ文字の変更（-e）
```bash
rlogin -e ^ remotehost
```
出力例:
```
Password:
...
```
リモート環境で、`^`（この例）をエスケープ文字として扱うようになる。

### 4. 8ビットクリーンモード（-8）
```bash
rlogin -8 remotehost
```
出力例:
```
Password:
...
```
文字コード制限を緩和し、8ビットの文字をそのまま利用できるようになる。

### 5. 行バッファモード（-L）
```bash
rlogin -L remotehost
```
出力例:
```
Password:
...
```
コマンドの出力を行単位で処理するため、大量の出力がある場合に制御がしやすくなる。

### 6. 自動ログインを使用（-a）
```bash
rlogin -a remotehost
```
出力例:
```
[remotehost ~]$
```
`.rhosts`ファイルなどが正しく設定されている場合、パスワード入力なしでログインできる。

