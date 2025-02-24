---
title:  【Unix】sedの解説
published: 2025-02-24
description: ""
tags: [sed]
category: Unix
draft: false
---
最終更新：2025-02-24

# 概要

`sed`コマンドは、テキストストリームを簡易的に編集・置換するストリームエディタである。ファイルやパイプからの入力に対して、指定したパターンに基づく変換を行い、結果を出力する。主な使用例は以下の通り。

```
sed [オプション] [スクリプト] [ファイル...]
```

---

## オプション一覧とハンズオン

### 1. -n（出力抑制）
- 機能: 自動出力を抑制し、表示指示があった行のみを出力する。  
- 使用例:  
  ```
  sed -n '/banana/p' sample.txt
  ```
- 出力例:
  ```
  banana
  ```

### 2. -e（複数スクリプト指定）
- 機能: 複数のスクリプト（置換や操作）を連続で実行する。  
- 使用例:  
  ```
  sed -e 's/apple/orange/g' -e 's/banana/grape/g' sample.txt
  ```
- 出力例:
  ```
  orange
  grape
  orange
  grape
  ```

### 3. -f（スクリプトファイル）
- 機能: 外部ファイルに記述したsedコマンド群を実行する。  
- 使用例:
  ```
  sed -f script.sed sample.txt
  ```
- 出力例:
  ```
  orange
  grape
  orange
  grape
  ```

### 4. -i（ファイル直接編集）
- 機能: 対象ファイルを上書き編集する。  
- 使用例:
  ```
  sed -i 's/apple/orange/g' sample.txt
  ```
- 出力例（実際のファイル内容が変化）:
  ```
  orange
  banana
  orange
  grape
  ```
- バックアップ:
  ```
  sed -i.bak 's/apple/orange/g' sample.txt
  ```

### 5. -r / -E（拡張正規表現）
- 機能: 括弧やパイプなどの拡張正規表現を簡易に使用できる。  
- 使用例:
  ```
  sed -E 's/(apple|banana)/fruit/g' sample.txt
  ```
- 出力例:
  ```
  fruit
  fruit
  orange
  grape
  ```

### 6. -l 数値（行幅指定）
- 機能: 出力を指定文字数で折り返して表示する。  
- 使用例:
  ```
  sed -l 20 's/a/A/g' longtext.txt
  ```
- 出力例: 長文を20文字単位で改行するため、実際の表示は文章次第。

### 7. -z（NUL区切り入力）
- 機能: NUL文字で区切られたテキストを行単位で扱う（環境依存）。  
- 使用例:
  ```
  sed -z 's/apple/orange/g' binarydata.txt
  ```
- 出力例: テキスト表現は変化が見えにくいが、NUL区切りで置換が行われる。

---

## パイプ処理でのハンズオン

他コマンドからの出力を`sed`で変換する例を示す。

```
ls | grep '\.txt$' | sed 's/\.txt$/.md/'
```

- `ls` の出力から `.txt` ファイルのみ抽出し、それらの拡張子を `.md` に変換して表示する。  
- 出力例:
  ```
  file1.md
  file2.md
  memo.md
  ```

