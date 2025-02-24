---
title:  【Unix】awkの解説
published: 2025-02-24
description: ""
tags: [awk]
category: Unix
draft: false
---
最終更新：2025-02-24

# 概要

`awk`コマンドは、テキスト処理やパターンマッチングに特化したスクリプト言語である。行単位でデータを読み込み、フィールド分割、計算、条件分岐などを柔軟に行うことができる。

使用方法は以下の通りである。  
```
awk [オプション] '条件 {アクション}' ファイル...
```

---

## オプション一覧とハンズオン

### 1. -F (フィールド区切り指定)
- 機能: 入力データを区切り文字で分割する際に使用する。  
- 使用例:
  ```bash
  awk -F, '{ print $3 }' data.csv
  ```
  上記例では、CSVファイルをカンマ区切りで分割し、3列目を表示する。

- 出力例:
  ```
  80
  90
  ```

### 2. -v (変数の事前設定)
- 機能: awkスクリプト内で使用する変数を事前に指定する。  
- 使用例:
  ```bash
  awk -v threshold=10 '{ if ($1 > threshold) print $0 }' numbers.txt
  ```
  上記例では、1列目の数値が`threshold`（ここでは10）を超える行を表示する。

- 出力例:
  ```
  15 25 35
  ```

### 3. -f (スクリプトファイルの指定)
- 機能: 外部ファイルに記述したawkスクリプトを適用する。  
- 使用例:
  ```bash
  awk -f script.awk sample.txt
  ```
  上記例では、`script.awk`に記載されたawkコマンドを`sample.txt`に適用する。

- 出力例:
  ```
  orange
  grape
  orange
  grape
  ```

### 4. -W (警告やモード設定)
- 機能: awkの動作モードや警告表示を制御する。  
- 使用例:
  ```bash
  awk -W compat -F, '{ print $2 }' data.csv
  ```
  上記例では、互換モードで動作し、カンマ区切りの2列目を表示する。

- 出力例:
  ```
  25
  22
  ```

---

## ハンズオン例

### 演習1: フィールド分割と集計
1. `numbers.txt`が以下の内容だとする。スペース区切りで3列ある。
   ```
   10 20 30
   15 25 35
   5  10 15
   ```
2. 各行の合計値を表示する例。
   ```bash
   awk '{ sum = $1 + $2 + $3; print sum }' numbers.txt
   ```
3. 出力例:
   ```
   60
   75
   30
   ```

### 演習2: CSVファイルの特定列抽出
1. `data.csv`が以下の内容でカンマ区切りだとする。
   ```
   name,age,score
   alice,25,80
   bob,22,90
   ```
2. `awk`で実行する。
   ```bash
   awk -F, '{ print $1 " is " $2 " years old with score " $3 }' data.csv
   ```
3. 出力例:
   ```
   name is age years old with score score
   alice is 25 years old with score 80
   bob is 22 years old with score 90
   ```

### 演習3: 変数を用いた条件処理
1. `numbers.txt` 内で閾値を変数として設定する例。
   ```bash
   awk -v limit=50 '{ if(($1 + $2 + $3) > limit) print $0 }' numbers.txt
   ```
2. 出力例（`numbers.txt`の合計が50を超える行のみ表示）:
   ```
   15 25 35
   ```

---

## パイプ処理でのハンズオン

awkの前段でsedを利用し、テキストを置換してからawkで処理する例。

```bash
cat sample.txt | sed 's/apple/orange/g' | awk '{print $1, $2}'
```

- 上記例では、`sed`で`apple`を`orange`に置換した後、`awk`で1列目と2列目を表示する。  
- 出力例:
  ```
  orange banana
  orange grape
  ```

---

## アクション部分に関する解説

awkのアクション部分では、以下のような入力を使用することができる。入力に使える一例を以下に示す。
- `$1, $2, ...`: 各フィールドを示す。`$0`は全体行を示す。
- `print`: 標準出力に出力する。
- `if, else`: 条件分岐を行う。
- `sum = $1 + $2`: 計算を行い、変数に格納する。

これらを組み合わせることで、複雑なデータ処理や条件抽出を柔軟に行うことができる。

