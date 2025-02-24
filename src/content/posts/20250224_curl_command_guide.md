---
title:  【Unix】curlの解説
published: 2025-02-24
description: ""
tags: [curl]
category: Unix
draft: false
---
最終更新：2025-02-24

# 概要

`curl`コマンドは、インターネット上のデータを取得したり送信したりするためのツールである。HTTP、HTTPS、FTPなど、様々なプロトコルに対応しており、APIの利用やファイルダウンロードに便利である。

使用方法は以下の通りである。
```
curl [オプション] [URL]
```

---

## オプション一覧とハンズオン

### 1. -o (出力ファイルを指定)
- 機能: ダウンロードしたデータを指定したファイルに保存する。  
- 使用例:
  ```bash
  curl -o output.html http://example.com
  ```
  上記例では、`http://example.com`の内容を`output.html`に保存する。

- 出力例:
  ファイル`output.html`に取得したHTMLデータが保存される。

### 2. -O (URLからファイル名を自動取得)
- 機能: ダウンロードしたデータをURLから自動的に取得したファイル名で保存する。  
- 使用例:
  ```bash
  curl -O http://example.com/file.txt
  ```
  上記例では、`http://example.com/file.txt`の内容をローカルに`file.txt`という名前で保存する。

- 出力例:
  ファイル`file.txt`に取得したデータが保存される。

### 3. -L (リダイレクトを追跡)
- 機能: HTTPリダイレクトを自動的に追跡する。  
- 使用例:
  ```bash
  curl -L http://example.com
  ```
  上記例では、リダイレクトが発生しても最終的なURLの内容を取得する。

- 出力例:
  リダイレクト先の最終ページのデータが表示される。

### 4. -I (ヘッダ情報のみ取得)
- 機能: HTTPヘッダ情報のみを取得する。  
- 使用例:
  ```bash
  curl -I http://example.com
  ```
  上記例では、`http://example.com`のヘッダ情報のみを表示する。

- 出力例:
  ```
  HTTP/1.1 200 OK
  Date: Mon, 24 Feb 2025 03:05:10 GMT
  Content-Type: text/html; charset=UTF-8
  ...
  ```

### 5. -d (データを送信)
- 機能: 指定したデータをPOSTリクエストとして送信する。  
- 使用例:
  ```bash
  curl -d "name=John&age=30" http://example.com/form
  ```
  上記例では、`name`と`age`のデータを`POST`として送信する。

- 出力例:
  サーバからのレスポンスが表示される。

### 6. -H (ヘッダを追加)
- 機能: リクエストにカスタムヘッダを追加する。  
- 使用例:
  ```bash
  curl -H "Authorization: Bearer token" http://example.com
  ```
  上記例では、`Authorization`ヘッダを追加してリクエストを送信する。

- 出力例:
  サーバからのレスポンスが表示される。

### 7. -u (ユーザ認証)
- 機能: HTTPの基本認証を行う。  
- 使用例:
  ```bash
  curl -u username:password http://example.com
  ```
  上記例では、`username`と`password`を用いて認証を行う。

- 出力例:
  認証が成功すれば、サーバからのレスポンスが表示される。

---

## パイプ処理でのハンズオン

curlコマンドを使ったパイプ処理の例を示す。curlで取得したデータをgrepでフィルタリングする。

```bash
curl -s http://example.com | grep "Example"
```

- 上記例では、`curl`で取得したHTMLデータから`Example`という単語を含む行を表示する。  
- 出力例:
  ```
  <title>Example Domain</title>
  <h1>Example Domain</h1>
  ```

