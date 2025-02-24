---
title:  【Unix】tarの解説
published: 2025-02-24
description: ""
tags: [tar]
category: Unix
draft: false
---
最終更新：2025-02-24

# 概要

`tar`コマンドは、ファイルやディレクトリをまとめたり（アーカイブ化）、まとめたものを取り出したりするためのツールである。圧縮ツールと組み合わせることで、一括で圧縮・解凍ができる便利なコマンドである。

使用方法は以下の通りである。  
```
tar [オプション] [アーカイブファイル] [ファイル...]
```

# オプション

# 出力関係

## -v オプション

・機能: コマンドの実行過程が表示される。

# 圧縮

## -c オプション（アーカイブ作成）

・機能: ファイルやディレクトリをアーカイブファイルにまとめる。  
・使用例:

```bash
tar -cvf archive.tar file1 file2
```

上記では`file1`および`file2`をまとめ、`archive.tar`として保存する。  
出力例（-vを指定した場合）:

```
file1
file2
```

## -z オプション（gzip圧縮を同時に行う）

・機能: gzip形式で圧縮する。  
・使用例:

```bash
tar -czvf archive.tar.gz file1 file2
```

`archive.tar.gz`が生成され、中身はgzip形式で圧縮されている。  
出力例:

```
file1
file2
```

## -j オプション（bzip2圧縮を同時に行う）

・機能: bzip2形式で圧縮する。  
・使用例:

```bash
tar -cjvf archive.tar.bz2 dir1
```

`archive.tar.bz2`が生成され、`dir1`の内容がbzip2形式でまとめられる。  
出力例:

```
dir1/
dir1/file1
dir1/file2
```

# 解凍（展開・取り出し）

## -x オプション（アーカイブの展開）

・機能: アーカイブファイルを展開する。  
・使用例:

```bash
tar -xvf archive.tar
```

`archive.tar`に含まれるファイル・ディレクトリをカレントディレクトリに展開する。  
出力例:

```
file1
file2
```

## -z オプション（gzip圧縮アーカイブの解凍）

・解説: gzip形式に圧縮されたアーカイブを解凍しながら展開する。  
・使用例:

```bash
tar -xzvf archive.tar.gz
```

gzip圧縮された`archive.tar.gz`を解凍し、ファイルを取り出す。  
出力例:

```
file1
file2
```

## -j オプション（bzip2圧縮アーカイブの解凍）

・解説: bzip2形式に圧縮されたアーカイブを解凍しながら展開する。  
・使用例:

```bash
tar -xjvf archive.tar.bz2
```

bzip2形式で圧縮された`archive.tar.bz2`を解凍し、ファイルを取り出す。  
出力例:

```
dir1/
dir1/file1
dir1/file2
```

# その他のオプション

## -t（アーカイブ内容の一覧表示）

・機能: アーカイブファイルの内容を確認する。  
・使用例:

```bash
tar -tvf archive.tar
```

出力例:

```
-rw-r--r-- user/group    1024 date file1
-rw-r--r-- user/group    2048 date file2
```

## -f（アーカイブファイル名を指定）

・機能: 作成・展開時にアーカイブファイルを明示的に指定する。  
・解説: `-f`の次にアーカイブファイル名を指定する。省略すると標準入出力を利用する。

