---
title:  【Unix】wgetの解説
published: 2025-02-24
description: ""
tags: [wget]
category: Unix
draft: false
---
最終更新：2025-02-24

# 概要

`wget`コマンドは、Webからファイルをダウンロードするための非対話型ネットワークダウンロードツールである。HTTP、HTTPS、FTPなどのプロトコルをサポートしており、再帰的なダウンロードやリジューム機能も備えている。

使用方法は以下の通りである。
```
wget [オプション] [URL]
```

---

## オプション一覧とハンズオン

### 1. -O (出力ファイル名を指定)
- 機能: ダウンロードしたデータを指定したファイル名で保存する。  
- 使用例:
  ```bash
  wget -O output.html http://example.com
  ```
  上記例では、`http://example.com`の内容を`output.html`として保存する。

- 出力例:
  ファイル`output.html`に取得したHTMLデータが保存される。

### 2. -c (ダウンロードの継続)
- 機能: 中断されたダウンロードを再開する。  
- 使用例:
  ```bash
  wget -c http://example.com/largefile.zip
  ```
  上記例では、`largefile.zip`のダウンロードを中断後に再開する。

- 出力例:
  ダウンロードが中断された位置から再開される。

### 3. -r (再帰的ダウンロード)
- 機能: 指定したURL以下のリンクを再帰的に辿ってダウンロードする。  
- 使用例:
  ```bash
  wget -r http://example.com/dir/
  ```
  上記例では、`http://example.com/dir/`以下のリンクを再帰的にダウンロードする。

- 出力例:
  ディレクトリ構造を保ったまま、リンク先のファイルがすべてダウンロードされる。

### 4. -P (保存ディレクトリを指定)
- 機能: ダウンロードしたファイルを指定したディレクトリに保存する。  
- 使用例:
  ```bash
  wget -P /path/to/directory http://example.com/file.txt
  ```
  上記例では、`file.txt`を指定したディレクトリに保存する。

- 出力例:
  ファイルが指定されたディレクトリに保存される。

### 5. -q (静かな出力)
- 機能: ダウンロード中の出力を抑制する。  
- 使用例:
  ```bash
  wget -q http://example.com/file.txt
  ```
  上記例では、ダウンロード中の出力を抑制して`file.txt`をダウンロードする。

- 出力例:
  進捗状況などのメッセージが表示されずにダウンロードが完了する。

### 6. --limit-rate (ダウンロード速度の制限)
- 機能: ダウンロード速度を指定したレートに制限する。  
- 使用例:
  ```bash
  wget --limit-rate=100k http://example.com/largefile.zip
  ```
  上記例では、ダウンロード速度を毎秒100キロバイトに制限してダウンロードを行う。

- 出力例:
  指定した速度でダウンロードが行われる。

### 7. --user, --password (ユーザ認証)
- 機能: HTTPまたはFTPのユーザ認証情報を指定する。  
- 使用例:
  ```bash
  wget --user=username --password=secret http://example.com/protectedfile
  ```
  上記例では、指定したユーザ名とパスワードで認証を行い、ファイルをダウンロードする。

- 出力例:
  認証が成功すれば、ファイルがダウンロードされる。

---

## パイプ処理でのハンズオン

`wget`コマンドを使ったパイプ処理の例を示す。`wget`で取得したデータを`grep`でフィルタリングする。

```bash
wget -qO- http://example.com | grep "Example"
```

- 上記例では、`wget`で取得したHTMLデータから`Example`という単語を含む行を表示する。  
- 出力例:
  ```
  <title>Example Domain</title>
  <h1>Example Domain</h1>
  ```

