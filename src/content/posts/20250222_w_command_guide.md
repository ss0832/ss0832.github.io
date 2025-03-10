---
title:  【Unix】wコマンドの解説
published: 2025-02-22
description: ""
tags: [w]
category: Unix
draft: false
---
最終更新：2025-02-22


`w`コマンドは、Unix系システムにおいて現在ログインしているユーザーの情報や、各ユーザーが実行中のプロセスの概要、システムの稼働状況などを表示するためのコマンドである。管理者はこのコマンドを用いて、システムの利用状況、ユーザーの活動状況、システム負荷などを把握することができる。


## 出力結果の例

以下は、`w`コマンドを実行した際の出力例である。

```bash
 04:49:25 up 10 days,  3:12, 4 users, load average: 0.15, 0.20, 0.25
USER     TTY      FROM             LOGIN@   IDLE   JCPU   PCPU WHAT
alice    pts/0    192.168.0.1      04:45    2:00   0.05s  0.05s -bash
bob      pts/1    192.168.0.2      04:44    1:30   0.04s  0.04s -bash
charlie  pts/2    192.168.0.3      04:40    1:10   0.03s  0.03s -bash
dave     pts/3    192.168.0.4      04:41    15:00  0.10s  0.10s -bash
```

## 出力結果の解釈

上記の出力結果は、以下の各部分に分けて解釈することができる。

- **最上行**  
  ```
   04:49:25 up 10 days,  3:12, 4 users, load average: 0.15, 0.20, 0.25
  ```
  - `04:49:25`: 現在のシステム時刻である。  
  - `up 10 days,  3:12`: システムが直近の再起動から10日と3時間12分経過していることを示す。  
  - `4 users`: 現在4人のユーザーがログインしている。  
  - `load average: 0.15, 0.20, 0.25`: 過去1分、5分、15分のロードアベレージがそれぞれ0.15, 0.20, 0.25である。これにより、システムの負荷の傾向を把握できる。

- **ヘッダー行**  
  ```
  USER     TTY      FROM             LOGIN@   IDLE   JCPU   PCPU WHAT
  ```
  この行は、以降に続く各ユーザーの情報の列見出しである。  
  - `USER`: ログインしているユーザー名。  
  - `TTY`: 各ユーザーが使用している端末（TTY）の名称。  
  - `FROM`: リモートホストからログインしている場合、そのホストのIPアドレスまたはホスト名。  
  - `LOGIN@`: ログインを開始した時刻。  
  - `IDLE`: ユーザーが最後に操作したからの経過時間。  
  - `JCPU`: 該当TTYで実行中の全プロセスの合計CPU使用時間。  
  - `PCPU`: 現在実行中のプロセスのCPU使用時間。  
  - `WHAT`: 現在ユーザーが実行中のコマンドやシェルの名称。

- **各ユーザーの行**  
  例として、以下の1行を見てみる。
  ```
  alice    pts/0    192.168.0.1      04:45    2:00   0.05s  0.05s -bash
  ```
  - `alice`: ユーザー名。  
  - `pts/0`: 使用中の端末。  
  - `192.168.0.1`: リモートホストのIPアドレス、またはローカルの場合は省略となる。  
  - `04:45`: ユーザーがログインした時刻。  
  - `2:00`: 最後の操作から2分間アイドル状態である。  
  - `0.05s`: `JCPU`値、端末内で使用された全コマンドのCPU時間の合計。  
  - `0.05s`: `PCPU`値、現在実行中のプロセスのCPU使用時間。  
  - `-bash`: 現在実行中のコマンドで、ここではbashシェルが実行されていることを示す。

## ハンズオン

以下に、`w`コマンドを使用してシステムの稼働状況を把握する例を示す。

### 1. `w`コマンドの基本使用例

ターミナルで以下のコマンドを実行せよ。

```bash
w
```

このコマンドは、上記のような出力結果を表示する。出力結果に含まれる各情報をもとに、システムの稼働状況やユーザーの活動状況を確認することができる。

### 2. 出力結果を解析し、システムの状態をチェックする

出力結果の最上行をよく確認し、以下の点について検証せよ。

- システム時刻とシステムの稼働時間  
- 現在ログインしているユーザー数  
- ロードアベレージを見て、瞬間的なシステム負荷を把握する

また、各ユーザー行についても、どのユーザーがどの端末で作業しているか、アイドル時間が長いユーザーがいないかなどをチェックすることが重要である。

---

このように、`w`コマンドの出力はシステムの全体的な利用状況および各ユーザーの活動状況を即座に把握するために有用である。定期的に実行し、システム管理やトラブルシューティングの参考とすることが望ましい。