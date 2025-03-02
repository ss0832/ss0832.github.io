---
title:  【Unix】rmdirの解説
published: 2025-03-02
description: ""
tags: [rmdir]
category: Unix
draft: false
---
最終更新：2025-03-02

# 概要

`rmdir`コマンドは、空のディレクトリを削除するためのUnixコマンドである。ディレクトリ内にファイルや他のディレクトリが存在する場合は削除に失敗する。この性質によって、誤って重要なファイルが含まれたディレクトリを削除してしまうことを防ぐことができる。

## コマンドの使用方法

```
rmdir [オプション] ディレクトリ名...
```

複数のディレクトリ名を指定することで、一度に複数の空ディレクトリを削除できる。


## オプション一覧と詳細

- `-p`, `--parents`  
  指定したパスの親ディレクトリも、それが空になった場合に連続して削除する。  
- `-v`, `--verbose`  
  削除の処理内容を詳細に表示する。  
- `--ignore-fail-on-non-empty`  
  削除対象のディレクトリが空でない場合でもエラーメッセージを表示せず処理を続行する。  
- `--help`  
  ヘルプメッセージを表示して終了する。  
- `--version`  
  バージョン情報を表示して終了する。  


## ハンズオン（オプションと出力例）

### 1. 基本的な使用方法

空のディレクトリを削除する一番シンプルな形式である。

```bash
mkdir test_dir
rmdir test_dir
```
出力例:
```
（何も表示されません）
```
通常、Unixコマンドは成功時に何も出力しないことが多い。

### 2. -p オプション（親ディレクトリも削除）

ネストしたディレクトリ構造を一度に削除する場合に便利である。

```bash
mkdir -p parent/child/grandchild
rmdir -p parent/child/grandchild
```
出力例:
```
（何も表示されません）
```
このコマンドは、`grandchild`ディレクトリを削除し、その後`child`、`parent`ディレクトリが空になれば順次削除する。

### 3. -v オプション（詳細表示モード）

削除した内容を確認したい場合に役立つ。

```bash
mkdir verbose_test
rmdir -v verbose_test
```
出力例:
```
rmdir: removing directory, 'verbose_test'
```

### 4. --ignore-fail-on-non-empty オプション

非空ディレクトリを削除しようとした場合のエラーメッセージを抑制する。

```bash
mkdir non_empty_dir
touch non_empty_dir/file.txt
rmdir --ignore-fail-on-non-empty non_empty_dir
```
出力例:
```
（何も表示されません）
```
ディレクトリ内にファイルがあるため削除は失敗するが、エラーメッセージは表示されない。

### 5. -p と -v の組み合わせ

親ディレクトリも削除し、詳細も表示する。

```bash
mkdir -p combo/test/dir
rmdir -pv combo/test/dir
```
出力例:
```
rmdir: removing directory, 'combo/test/dir'
rmdir: removing directory, 'combo/test'
rmdir: removing directory, 'combo'
```
ディレクトリが削除される過程が順に表示される。

### 6. 複数のディレクトリを一度に削除

```bash
mkdir dir1 dir2 dir3
rmdir dir1 dir2 dir3
```
出力例:
```
（何も表示されません）
```

### 7. 空でないディレクトリを削除しようとした場合

```bash
mkdir non_empty
touch non_empty/sample.txt
rmdir non_empty
```
出力例:
```
rmdir: failed to remove 'non_empty': Directory not empty
```
このエラーメッセージは、ディレクトリが空でないため削除できないことを示している。


`rmdir`コマンドは、特にシェルスクリプトで一時ディレクトリを安全に削除したい場合や、ディレクトリ構造を整理する際に役立つ。より強力な削除には`rm -r`コマンドがあるが、`rmdir`は空ディレクトリだけを対象とする安全なコマンドである点が特徴となっている。