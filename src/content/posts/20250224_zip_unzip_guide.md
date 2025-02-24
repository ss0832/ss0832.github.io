---
title:  【Unix】zipとunzipの解説
published: 2025-02-24
description: ""
tags: [zip, unzip]
category: Unix
draft: false
---
最終更新：2025-02-24

# 概要

`zip`コマンドはファイルを圧縮し、`unzip`コマンドは圧縮されたファイルを解凍するためのツールである。テキストファイルや複数のディレクトリをスムーズにやり取りする際に有用である。

# 導入方法

Unix系ディストリビューションでインストールがされていない場合は、以下のようにパッケージマネージャを使用してインストールすることが多い。

例（Debian/Ubuntu系の場合）:
```bash
sudo apt-get update
sudo apt-get install zip unzip
```

例（CentOS/Fedora系の場合）:
```bash
sudo yum install zip unzip
```

---

# 圧縮（zipコマンド）

使用方法:  
```
zip [オプション] アーカイブ名.zip ファイル...
```

## オプション一覧とハンズオン

### -r (ディレクトリを再帰的に圧縮)

・機能: 指定したディレクトリ以下をすべて圧縮する。  
・例:  
```bash
zip -r archive.zip mydir
```
出力例:  
```
  adding: mydir/ (stored 0%)
  adding: mydir/file1.txt (deflated 70%)
  adding: mydir/subdir/ (stored 0%)
  ...
```

### -9 (最高圧縮率)

・機能: 圧縮率を最高に設定する。  
・例:  
```bash
zip -9 archive.zip largefile.dat
```
出力例:  
```
  adding: largefile.dat (deflated 85%)
```

### -q (出力抑制)

・機能: 処理中の進行メッセージを表示しない。  
・例:
```bash
zip -rq archive.zip mydir
```
出力例:
（特になし。エラーがない限り何も表示されず終了する）

### -e (パスワード保護)

・機能: パスワードを指定して暗号化する。  
・例:
```bash
zip -e archive.zip secrets.txt
```
出力例:
```
Enter password:
Verify password:
  adding: secrets.txt (stored 0%)
```

---

# 解凍（unzipコマンド）

使用方法:  
```
unzip [オプション] アーカイブ名.zip
```

## オプション一覧とハンズオン

### （オプションなし）

・機能: 指定したzipファイルを現在のディレクトリに解凍する。  
・例:
```bash
unzip archive.zip
```
出力例:
```
Archive:  archive.zip
  inflating: file1.txt
  inflating: file2.txt
```

### -d (解凍先ディレクトリ指定)

・機能: 解凍先のディレクトリを指定する。  
・例:
```bash
unzip archive.zip -d /path/to/unzip/
```
出力例:
```
Archive:  archive.zip
  inflating: /path/to/unzip/file1.txt
  inflating: /path/to/unzip/file2.txt
```

### -o (上書きモード)

・機能: 同名ファイルが存在する場合に上書きを行う。  
・例:
```bash
unzip -o archive.zip
```
出力例:
```
Archive:  archive.zip
replace file1.txt?  [y]es/[n]o/[A]ll
  inflating: file1.txt
  inflating: file2.txt
```
（オプションでアラートを抑制しそのまま上書きすることも可能）

### -P (パスワード指定解凍)

・機能: 暗号化されたzipファイルを解凍する際、パスワードを指定する。  
・例:
```bash
unzip -P mysecret archive.zip
```
出力例:
```
Archive:  archive.zip
  inflating: secrets.txt
```

