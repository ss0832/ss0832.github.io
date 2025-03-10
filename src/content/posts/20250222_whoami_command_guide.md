---
title:  【Unix】whoamiコマンドの解説
published: 2025-02-22
description: ""
tags: [whoami]
category: Unix
draft: false
---
最終更新：2025-02-22


`whoami`コマンドは、現在のユーザー名を表示するためのコマンドである。自分がどのアカウントとしてログインしているかを確認するために利用される。コマンドを実行すると、そのセッションにおける自分のユーザー名が標準出力に表示される。

## 使い方

```
whoami
```

実行結果は、例えば以下のようになる:

```
you
```

上記の結果は、現在のセッションでログインしているユーザーが`you`であることを示す。

## ハンズオン

以下に、`whoami`コマンドの使用例を示す。

### 1. 基本的な使用例

現在のユーザー名を確認するため、ターミナル上で以下のコマンドを実行せよ。

```bash
whoami
```

このコマンドは、シンプルに現在のユーザー名を返すため、スクリプトの中でユーザー認証やログの記録に活用される場合がある。

### 2. 環境変数との組み合わせ

シェルスクリプト内において、`whoami`を用いて変数に現在のユーザー名を格納し、後続の処理で利用する例である。

```bash
#!/bin/bash
# 現在のユーザー名を取得する
current_user=$(whoami)

# ユーザー名を表示する
echo "現在のユーザーは: ${current_user}"
```

上記のスクリプトを実行すると、「現在のユーザーは: you」のように表示される。

