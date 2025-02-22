---
title: 【Unix】aliasコマンドの解説
published: 2025-02-22
description: ""
tags: [alias]
category: Unix
draft: false
---
最終更新：2025-02-22

`alias`は、任意の文字列を特定のコマンド（群）と結びつける（エイリアスにする）ことで、コマンド実行時の入力を省略できるようにするシェルの機能である。

windows OSのショートカット機能と同様の機能である。

使うシェルによって設定ファイルや書き方が若干異なるため、以下に代表的な4つのシェルにおけるエイリアス設定例を示す。

alias = 偽名や仮名という意味。

## シェルごとの設定例

### bash / zsh

- 設定ファイル: 主に `~/.bashrc` または `~/.bash_profile`, `~/.zshrc`
- エイリアスの書式:  
  ```bash
  alias コマンド名='実際に実行するコマンド'
  ```
- 例:  
  ```bash
  alias ll='ls -alF'
  alias gs='git status'
  ```
- 設定ファイルを読み込むには以下を実行する。  
  ```bash
  source ~/.bashrc
  ```
  または、ログインし直すなどで反映させることができる。

  ```

### tcsh / csh

- 設定ファイル: `~/.tcshrc` や `~/.cshrc`  
- エイリアスの書式は csh と同様:  
  ```tcsh
  alias コマンド名 "実際に実行するコマンド"
  ```
- 例:  
  ```tcsh
  alias ll "ls -alF"
  alias gs "git status"
  ```
- 設定を反映させるには、以下を実行するか、シェルに再起動を行う。  
  ```tcsh
  source ~/.tcshrc
  ```
  
## ハンズオン

1. **一時的なエイリアス設定**  
   現在のシェルセッションのみ有効にする場合、次のようなコマンドを実行する。  
   ```bash
   alias ll='ls -lF'
   ```
   このエイリアスを設定することで、`ll`という短いコマンドで`ls -lF`が実行できるようになる。セッションを終了するとエイリアスは消える。

2. **設定ファイルに書き込む（永続化）**  
   bashを例にすると、`~/.bashrc`などに以下の行を追加し、反映させる。  
   ```bash
   echo "alias ll='ls -alF'" >> ~/.bashrc
   source ~/.bashrc
   ```
   これで、ターミナルを開き直してもエイリアスが自動的に読み込まれるようになる。Zshやcsh系シェルでも同様に設定ファイルを編集する。

3. **複数のコマンドをまとめる**  
   以下のようにすると、一つのエイリアスが複数のコマンドを実行する。  
   ```bash
   alias update-all='sudo apt update && sudo apt upgrade -y'
   ```
   即ち、`update-all`を実行することで、`sudo apt update`、`sudo apt upgrade -y`の順にコマンドが実行される。
   bash、zsh、csh、tcsh いずれの場合も設定ファイルへの記述方法は同様だが、cshとtcshではダブルクォートを用いる点に注意する。  

## 補足
- 筆者は、コマンド名が大文字であり、入力の手間がかかるときに`alias`を用いて小文字にすることがある。