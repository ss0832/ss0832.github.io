---
title: 【PowerShell】Get-Processコマンドの基本的な使い方
published: 2025-03-26
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-26

## 概要

Get-Processコマンドは、Windows PowerShellで稼働中のプロセス情報を取得するための基本的なコマンドレットである。プロセス名やID、メモリ使用量などを簡単に確認できるため、システムの状態把握や不要なプロセスの確認などに役立つ。PowerShell 2.0以降、またCore版を含むPowerShell 7.xでも利用できる。

## 基本的な使い方と出力例

### 1. すべてのプロセスを表示する
```powershell
Get-Process
```
出力例:
```
Handles  NPM(K)    PM(K)      WS(K)     CPU(s)     Id  SI ProcessName
-------  ------    -----      -----     ------     --  -- -----------
    422      24    13496      25380       0.75   7840   1 ApplicationFrameHost
    169      12     2872       7668       0.06  12296   1 ctfmon
    ...
```
- ProcessName: プロセスの名前  
- Id: プロセスID  
- CPU(s): CPU時間（秒）  
- WS(K): ワーキングセットのサイズ（KB）  

### 2. 特定のプロセスを名前で検索する
```powershell
Get-Process -Name chrome
```
出力例:
```
Handles  NPM(K)    PM(K)      WS(K)     CPU(s)     Id  SI ProcessName
-------  ------    -----      -----     ------     --  -- -----------
    510      26    39844      60520      14.70   9568   1 chrome
    ...
```
- -Name: 指定したプロセス名のプロセスのみを取得するオプションである

### 3. ワイルドカードでプロセスを検索する
```powershell
Get-Process -Name "*host*"
```
出力例:
```
Handles  NPM(K)    PM(K)      WS(K)     CPU(s)     Id  SI ProcessName
-------  ------    -----      -----     ------     --  -- -----------
    422      24    13496      25380       0.75   7840   1 ApplicationFrameHost
    ...
```
- "*host*"のようにアスタリスクを用いると、プロセス名にhostを含むものを取得できる

### 4. プロセスIDを指定して検索する
```powershell
Get-Process -Id 7840
```
出力例:
```
Handles  NPM(K)    PM(K)      WS(K)     CPU(s)     Id  SI ProcessName
-------  ------    -----      -----     ------     --  -- -----------
    422      24    13496      25380       0.75   7840   1 ApplicationFrameHost
```

### 5. 特定のプロパティのみ選択して表示する
```powershell
Get-Process -Name chrome | Select-Object Name, Id, CPU, WorkingSet
```
出力例:
```
Name    Id     CPU WorkingSet
----    --     --- ----------
chrome  9568  14.7    61972480
...
```
- Select-Object: 必要なプロパティだけを表示するときに便利である

## すべてのオプションに対するハンズオン

以下の手順を順番に実行することで、Get-Processコマンドのさまざまなオプションを試すことができる。

```powershell
# 1. PowerShellのバージョンを確認
$PSVersionTable.PSVersion

# 2. すべてのプロセスを表示
Get-Process

# 3. 特定のプロセスを名前指定で表示 (例: chrome)
Get-Process -Name chrome

# 4. 複数のプロセスを同時に表示 (例: chrome, explorer)
Get-Process -Name chrome, explorer

# 5. プロセスIDを指定して表示 (例: 1234)
Get-Process -Id 1234

# 6. ワイルドカードを使ってプロセスを検索 (例: *host*)
Get-Process -Name "*host*"

# 7. 特定のプロパティを表示
Get-Process -Name chrome | Select-Object Name, Id, CPU, WorkingSet

# 8. CPU使用率の高いプロセス順に並べて上位5つを表示
Get-Process | Sort-Object -Property CPU -Descending | Select-Object -First 5

# 9. メモリ使用量で並べて上位5つを表示
Get-Process | Sort-Object -Property WorkingSet -Descending | Select-Object -First 5

# 10. 詳細情報をリスト形式で表示 (例: explorer)
Get-Process -Name explorer | Format-List *

# 11. 応答なしプロセスの検索
Get-Process | Where-Object { $_.Responding -eq $false }
```

## 出力例の見方

- Handles: OSリソースのハンドル数  
- NPM(K): カーネルメモリ(Non-paged)使用量(KB)  
- PM(K): メインメモリ(Paged)使用量(KB)  
- WS(K): 実際に使用中のメモリ量(KB)  
- CPU(s): 累積CPU時間(秒)  
- Id: プロセスID  
- SI: セッションID  
- ProcessName: プロセス名  

## 対応PowerShellバージョン

- Windows PowerShell 2.0～5.1  
- PowerShell 6 (Core)・7.x (クロスプラットフォーム対応)

これらのバージョンでGet-Processコマンドは同様の利用が可能である。システム運用の基本として、稼働しているプロセスの状況を素早く把握するために活用するとよい。