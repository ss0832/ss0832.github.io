---
title:  【Unix】killコマンドの解説
published: 2025-02-22
description: ""
tags: [kill]
category: Unix
draft: false
---
最終更新：2025-02-22


`kill`コマンドは、Unix系システムにおいてプロセスにシグナルを送信し、プロセスを終了させる、または特定の動作を実行させるための基本的なツールである。プロセス管理の一環として、必要なプロセスを終了させたり、デバッグやシステム制御のためのシグナル送信に用いられる。以下に、`kill`コマンドの主要なオプションとその具体的な意味、ならびにハンズオンの例を詳述する。

---

# 使用方法

単純に `kill PID` と実行すると、デフォルトのシグナルである `SIGTERM` (15) がプロセスに送信される。例えば、

```bash
kill 1234
```

は、プロセスID `1234` のプロセスに対して `SIGTERM` を送信する。

---

# 主なオプションの詳細

- **-s, --signal \<シグナル名または番号\>**  
  送信するシグナルを明示的に指定するためのオプションである。  
  例:
  ```bash
  kill -s SIGKILL 1234
  ```
  この例では、`SIGKILL` (9) を使用してプロセスを強制終了する。シグナルは数字でも指定可能である:
  ```bash
  kill -s 9 1234
  ```

- **-l, --list**  
  利用可能なシグナルのリストを表示する。また、シグナル番号から対応するシグナル名を参照することもできる。  
  例:
  ```bash
  kill -l
  ```
  このコマンドは、すべてのシグナル名が一覧表示される。特定の番号に対してシグナル名を調べる場合は、番号をパラメータとして指定する。
  ```bash
  kill -l 9
  ```
  これにより、シグナル番号 9 に対応するシグナル名 `KILL` が表示される。

- **-q, --quiet** *(一部の実装のみ)*  
  シグナル送信時のエラーメッセージを抑制する。エラー出力を減らしたい場合に利用される。

- **-w, --wait** *(一部の実装のみ)*  
  シグナル送信後、プロセス終了を待機する。プロセスが終了するまで待機し、その結果を返す。特定の状況下で、送信したシグナルの効果を確認したい場合に有用である。

- **デフォルトシグナル (`SIGTERM`)**  
  オプションを指定しない場合、`kill`はデフォルトでシグナル `SIGTERM` (15) を送信する。`SIGTERM` は、優雅にプロセスを終了させるためのシグナルであるため、まずはこのシグナルを使用するのが一般的である。

- **強制終了シグナル (`SIGKILL`)**  
  プロセスが通常の終了シグナルに反応しない場合、`kill -9 PID` や `kill -s SIGKILL PID` を使用して、プロセスを強制終了することができる。このシグナルは無条件にプロセスを即時終了させるため、データの損失やリソースの適切な解放がなされない可能性がある点に注意が必要である。

---

# ハンズオン

以下に、`kill`コマンドの具体的な使用例と動作確認手順を示す。

## 1. プロセスの終了

1. ターミナルを開く。
2. まず、プロセス一覧を表示して終了させたいプロセスのPIDを確認する。例えば、`ps` コマンドを使用する:
   
   ```bash
   ps aux | grep myprocess
   ```

3. 終了させるプロセスのPIDが `1234` だとして、通常終了シグナル（`SIGTERM`）を送信する:

   ```bash
   kill 1234
   ```

   これにより、プロセスが正常に終了するよう試みる。

## 2. 強制終了

プロセスが `SIGTERM` に反応しない場合、強制終了シグナル `SIGKILL` を送信する:

```bash
kill -s SIGKILL 1234
```

または、短縮形として:

```bash
kill -9 1234
```

これにより、プロセスは無条件に終了される。

## 3. シグナル一覧の確認

利用可能なシグナルの一覧を確認し、シグナル名と番号の対応を知るには、以下のコマンドを実行する:

```bash
kill -l
```

さらに、特定の番号に対するシグナル名が知りたい場合、例えばシグナル番号 `15` なら:

```bash
kill -l 15
```

と実行し、`TERM` などのシグナル名が出力されることを確認する。

## 4. 複数プロセスへのシグナル送信

複数のプロセスに同じシグナルを同時に送信する場合は、PIDをスペースで区切って指定する:

```bash
kill 1234 5678 9012
```

これにより、指定されたすべてのプロセスに `SIGTERM` が送信される。

---
