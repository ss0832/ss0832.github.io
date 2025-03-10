---
title:  【Unix】ifconfigの解説
published: 2025-02-27
description: ""
tags: [ifconfig]
category: Unix
draft: false
---
最終更新：2025-02-27

# 概要

`ifconfig`コマンドは、ネットワークインタフェースの設定や状態確認を行うためのコマンドである。IPアドレスの設定やインタフェースの有効化・無効化などの操作が可能である。

## コマンドの使用方法

```
ifconfig [インタフェース名] [オプション] [...]
```
特定のインタフェースを指定しない場合は、現在アクティブなすべてのネットワークインタフェース情報を表示する。

---

## オプション一覧と詳細

- (インタフェース名未指定)  
  システム上のインタフェース情報をまとめて表示する。  
- (インタフェース名のみ指定)  
  指定したインタフェースの詳細情報を表示する。  
- `up`  
  ネットワークインタフェースを有効化する。  
- `down`  
  ネットワークインタフェースを無効化する。  
- `inet <IPアドレス>`  
  IPv4アドレスを設定する。  
- `netmask <サブネットマスク>`  
  ネットマスクを設定する。  
- `broadcast <ブロードキャストアドレス>`  
  ブロードキャストアドレスを設定する。  
- `mtu <サイズ>`  
  MTU（最大転送単位）を指定する。

---

## ハンズオンと出力例

### 1. すべてのインタフェース情報を表示する（インタフェース名未指定）
```bash
ifconfig
```
出力例:
```
eth0: flags=4163<UP,BROADCAST,RUNNING,MULTICAST>  mtu 1500
        inet 192.168.0.10  netmask 255.255.255.0  broadcast 192.168.0.255
        ...
lo:    flags=73<UP,LOOPBACK,RUNNING>  mtu 65536
        inet 127.0.0.1  netmask 255.0.0.0
        ...
```
現在アクティブなすべてのインタフェースが一覧表示される。

### 2. 特定のインタフェース情報を表示
```bash
ifconfig eth0
```
出力例:
```
eth0: flags=4163<UP,BROADCAST,RUNNING,MULTICAST>  mtu 1500
        inet 192.168.0.10  netmask 255.255.255.0  broadcast 192.168.0.255
        ...
```
指定したインタフェース `eth0` の詳細設定が表示される。

### 3. インタフェースを有効化（up）・無効化（down）
```bash
ifconfig eth0 up
ifconfig eth0 down
```
出力例（通常出力はなし）:
```
(コマンド実行後、ifconfigで状態を確認するとUP/DOWNが変化している)
```

### 4. IPアドレスなどの設定変更
```bash
ifconfig eth0 inet 192.168.0.20 netmask 255.255.255.0 broadcast 192.168.0.255
```
出力例（通常出力はなし）:
```
(設定変更が完了し、ifconfigで確認すると更新されたことがわかる)
```

### 5. MTUの設定
```bash
ifconfig eth0 mtu 1400
```
出力例（通常出力はなし）:
```
(設定変更後、ifconfig eth0 でMTUが1400と表示される)
```

## 補足

- `ifconfig`はネットワークインタフェースの情報確認や簡易的な設定に用いられるが、最新のLinuxディストリビューションでは`ip`コマンド（iproute2ツール）が推奨されることが多い点に留意するとよい。