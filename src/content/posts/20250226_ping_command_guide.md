---
title:  【Unix】pingの解説
published: 2025-02-26
description: ""
tags: [ping]
category: Unix
draft: false
---
最終更新：2025-02-26

# 概要

`ping`コマンドは、ネットワーク上にあるホストとの通信疎通を確認するために用いるコマンドである。ICMP (Internet Control Message Protocol) を利用して、相手ホストからの応答の有無や遅延時間を調べることができる。


## コマンドの使用方法

```
ping [オプション] ホスト名またはIPアドレス
```


## オプション一覧と詳細

- `-c count`  
  指定した数だけパケットを送信して終了する。
- `-i interval`  
  パケット送信の間隔を秒単位で指定する。
- `-w deadline`  
  指定した秒数後に自動的に終了する。
- `-W timeout`  
  応答待ちのタイムアウトを秒単位で指定する。
- `-q`  
  統計結果のみ出力（quietモード）。
- `-t ttl`  
  IPヘッダのTTL(Time To Live)を指定する。TTLとは、IPパケットがルータを通過する回数を制限する値である。
- `-a`  
  応答を受け取るたびに端末でベルを鳴らす。
- `-f`  
  大量のICMPパケットを送信するfloodモード。ネットワーク負荷に注意が必要。
- `-D`  
  各応答行の先頭にUNIXタイムスタンプを表示する。



## ハンズオン（全オプション）

### 1. -cオプション（送信回数を指定）
```bash
ping -c 4 example.com
```
出力例:
```
PING example.com (93.184.216.34) 56(84) bytes of data.
64 bytes from 93.184.216.34: icmp_seq=1 ttl=57 time=15.3 ms
64 bytes from 93.184.216.34: icmp_seq=2 ttl=57 time=15.1 ms
...
--- example.com ping statistics ---
4 packets transmitted, 4 received, 0% packet loss, time 3004ms
rtt min/avg/max/mdev = 15.1/15.2/15.3/0.0 ms
```

### 2. -iオプション（送信間隔を指定）
```bash
ping -i 2 example.com
```
2秒ごとにICMPエコーを送信し、応答を受け取る。

### 3. -wオプション（終了までの時間を指定）
```bash
ping -w 5 example.com
```
5秒経過すると自動的に終了する。

### 4. -Wオプション（タイムアウトを指定）
```bash
ping -c 4 -W 1 example.com
```
パケットごとの応答待ち時間を1秒に制限する。

### 5. -qオプション（統計のみ表示）
```bash
ping -c 4 -q example.com
```
実行中の詳細を出力せず、最後の統計のみ表示する。

### 6. -tオプション（TTLを指定）
```bash
ping -t 64 example.com
```
TTLを64に設定し、相手ホストまでの経路上のルータ数を制限する。

### 7. -aオプション（応答時にベルを鳴らす）
```bash
ping -a example.com
```
応答を受け取るたびに端末にベル（ビープ音）を鳴らす。

### 8. -fオプション（floodモード）
```bash
ping -f example.com
```
極めて大量のパケットを連続送信するモードであり、ネットワーク負荷が高いため注意が必要である。

### 9. -Dオプション（タイムスタンプを表示）
```bash
ping -D example.com
```
各応答行の先頭にUNIXタイムスタンプ（エポック秒）を表示する。

### 参考  
- https://eng-entrance.com/linux-command-ping  
- https://webkaru.net/linux/ping-command/  
- https://www.server-memo.net/tips/command/ping/ping.html  

