---
title:  【Unix】nslookupの解説
published: 2025-02-26
description: ""
tags: [nslookup]
category: Unix
draft: false
---
最終更新：2025-02-26

# 概要
`nslookup`コマンドは、DNSサーバへ問い合わせを行い、指定したホスト名やIPアドレスに関するDNSレコードを取得するためのツールである。ネットワークトラブルシューティングやDNS設定の確認に有効である。

## コマンドの使用方法
```
nslookup [オプション] [ホスト名/IPアドレス]
```
オプションを指定せずにホスト名やIPアドレスを引数に与えるだけでも、簡単にDNSクエリを行うことができる。


## オプション一覧と詳細

- (オプション未指定): 指定したホスト名やIPアドレスのデフォルト情報を表示する。  
- `-type=レコードタイプ`: クエリするDNSレコードの種類を指定する。たとえば、以下のようなレコードを指定できる。  
  - `A`: IPv4アドレス  
  - `AAAA`: IPv6アドレス  
  - `CNAME`: 別名レコード  
  - `MX`: メールエクスチェンジャー  
  - `NS`: ネームサーバ  
- `-server=DNSサーバ`: 使用するDNSサーバを指定する。 (例: `-server=8.8.8.8`)  
- `-timeout=秒数`: 応答を待機するタイムアウト時間を指定する。  
- `-retry=回数`: 再試行回数を指定する。 (環境によっては非対応の場合がある)  
- 対話モード: `nslookup`だけを打ち込むと対話モードに入り、`server [DNSサーバ]`などのサブコマンドを対話的に入力できる。


## ハンズオン（オプションと出力例）

### 1. 基本的な使い方（オプションなし）
```
nslookup example.com
```
出力例:
```
Server:         8.8.8.8
Address:        8.8.8.8#53

Non-authoritative answer:
Name:   example.com
Address: 93.184.216.34
```
DNSサーバのアドレスと、クエリしたホストのIPアドレスを確認できる。

### 2. -type=レコードタイプを指定してクエリする
よく使うレコードタイプの例として、以下のものがある。  
- `A` (IPv4アドレス)  
- `AAAA` (IPv6アドレス)  
- `CNAME` (別名レコード)  
- `MX` (メールサーバ情報)  
- `NS` (ネームサーバ)

それぞれの使用例を示す。

#### (a) Aレコード
```
nslookup -type=A example.com
```
出力例:
```
Non-authoritative answer:
Name:   example.com
Address: 93.184.216.34
```

#### (b) AAAAレコード
```
nslookup -type=AAAA example.com
```
出力例（IPv6アドレスが割り当てられていない場合もある）:
```
Non-authoritative answer:
Name:   example.com
Address: 2606:2800:220:1:248:1893:25c8:1946
```

#### (c) CNAMEレコード
```
nslookup -type=CNAME www.example.com
```
出力例:
```
Non-authoritative answer:
www.example.com     canonical name = example.com.
Name: example.com
Address: 93.184.216.34
```
`example.com`がCNAME（別名）として`www.example.com`を指し示している場合に表示される。

#### (d) MXレコード
```
nslookup -type=MX example.com
```
出力例:
```
Non-authoritative answer:
example.com  mail exchanger = 10 mail.example.com.
```
メール受信に使用するサーバ情報が表示される。

#### (e) NSレコード
```
nslookup -type=NS example.com
```
出力例:
```
Non-authoritative answer:
example.com   nameserver = ns1.example.com.
example.com   nameserver = ns2.example.com.
```
ネームサーバのホスト名が表示される。

### 3. -serverオプション（DNSサーバを指定する）
```
nslookup -server=8.8.8.8 example.com
```
出力例:
```
Server:         8.8.8.8
Address:        8.8.8.8#53

Non-authoritative answer:
Name:   example.com
Address: 93.184.216.34
```
ここではGoogleの公共DNSサーバ`8.8.8.8`を使用している。

### 4. -timeoutオプション（応答待ち時間の設定）
```
nslookup -timeout=2 example.com
```
応答を2秒間だけ待機し、タイムアウトした場合はエラーを表示する。

### 5. -retryオプション（再試行回数の指定）
```
nslookup -retry=3 example.com
```
応答が得られない場合、最大3回再試行する。（環境によってはこのオプションに対応していない場合がある。）

### 6. 対話モード
```
nslookup
> set type=MX
> example.com
> server 8.8.4.4
> set type=NS
> example.com
> exit
```
対話モードに入り、DNSサーバやレコードタイプを変更しながら柔軟にDNS情報を取得できる。

