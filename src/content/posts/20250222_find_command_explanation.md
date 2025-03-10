---
title:  【Unix】findコマンドの解説
published: 2025-02-22
description: ""
tags: [find]
category: Unix
draft: false
---
最終更新：2025-02-22

`find`コマンドは、指定したディレクトリ以下のファイルやディレクトリを再帰的に検索するためのユーティリティである。
オプションを組み合わせることで、名前や更新日時、パーミッションなどさまざまな条件を指定してファイルを絞り込むことができる。


## 使用方法

#### 1. `-name <パターン>`
ファイル名が指定したパターンに合致するかを検索する式である。ワイルドカード（`*`など）を使って拡張子やファイル名部分を指定できる。

例:  
```bash
find /path/to/search -name "*.txt"
```
「`.txt`」で終わるファイルを再帰的に検索する。

#### 2. `-type <種類>`
ファイルの種類を指定する。最もよく使われる例は`-type f`（通常ファイル）と`-type d`（ディレクトリ）である。

例:  
```bash
find /path/to/search -type d
```
ディレクトリのみを検索対象とし、一覧表示する。

#### 3. `-mtime <日数>`
ファイルの最終更新日時が指定した日数に基づいて一致するファイルを検索する。`-mtime +N`は「N日より前に更新された」、`-mtime -N`は「N日以内に更新された」を意味する。

例:  
```bash
# 3日以内に更新されたファイルを検索
find /path/to/search -mtime -3
```

#### 4. `-maxdepth <階層>` / `-mindepth <階層>`
検索の深さを制限するオプションである。`-maxdepth 1`を指定すると、直下の階層のみを検索する。`-mindepth 1`であれば、カレントディレクトリそのものは検索しないで、それ以下の階層を対象とする。

例:  
```bash
# 最深1階層だけ検索する
find /path/to/search -maxdepth 1 -name "*.log"

# 2階層目以降だけを検索する
find /path/to/search -mindepth 2 -type f
```

#### 5. `-exec <コマンド> {} \;`
検索条件に合致したファイルに対して、指定したコマンドを実行する。`{}`がファイル名のプレースホルダとして使用される。「`\;`」は、`-exec`オプションの終端を示す。

例:  
```bash
# 拡張子が tmp のファイルを削除する
find /path/to/search -name "*.tmp" -exec rm {} \;
```
ファイルを削除する場合など、必要性を慎重に確認してから実行することが望ましい。


## ハンズオン

以下では、実際に`find`を使った簡単な操作を紹介する。

1. **指定ディレクトリにある「.txt」ファイルを検索する**  
   ```bash
   find /path/to/search -name "*.txt" 
   ```


2. **3日以内に更新された「.log」ファイルを削除する**  
   ```bash
   find /path/to/search -name "*.log" -mtime -3 -exec rm {} \;
   ```
   更新が新しいログファイルだけを削除するので、古いログファイルには影響しない。ただし、誤削除を防ぐためにも実行前に必ず確認することが望ましい。

3. **ディレクトリのみを検索してパスを表示する**  
   ```bash
   find /path/to/search -type d
   ```
   ディレクトリのみを対象とするため、ファイルは結果に表示されない。

4. **サブディレクトリを一切たどらずに検索する (`-maxdepth`)**  
   ```bash
   find /path/to/search -maxdepth 1 -name "*.sh"
   ```
   `-maxdepth 1`を指定しているため、指定ディレクトリの直下だけが検索対象となる。

5. **シェルスクリプトをバックアップするコマンドをまとめて実行する (`-exec`)**  
   ```bash
   find /path/to/search -name "*.sh" -exec cp {} {}.bak \;
   ```
   検索でヒットしたすべてのシェルスクリプトに対して「`.bak`」の拡張子を付けたバックアップファイルを作成する。




