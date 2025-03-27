---
title:  【PowerShell】Stop-Processコマンドの基本と応用
published: 2025-03-27
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-27

## 概要

PowerShellの`Stop-Process`コマンドレットは、実行中のプロセスを強制的に停止（終了）させるためのコマンドである。タスクマネージャーで「タスクの終了」を行うのと同様の機能をコマンドラインから実行できる。プロセスの一括終了や条件に基づく終了など、管理タスクやトラブルシューティングで非常に役立つツールである。

## 基本的な使い方

### プロセスIDによる終了

```powershell
# プロセスID 1234 のプロセスを終了する
Stop-Process -Id 1234
```

### プロセス名による終了

```powershell
# notepad という名前のプロセスをすべて終了する
Stop-Process -Name notepad
```

出力例:
```
# 通常は成功時に出力はない
```

### 確認を求める

```powershell
# 確認メッセージを表示してから終了する
Stop-Process -Name chrome -Confirm
```

出力例:
```
本当に Stop-Process コマンドレットを実行しますか?
chrome (1234) を停止します
[Y] はい(Y)  [A] すべて続行(A)  [N] いいえ(N)  [L] すべて無視(L)  [S] 中断(S)  [?] ヘルプ (既定値は "Y"): y
```

## 高度な使用方法と各オプションのハンズオン

### -Id

プロセスIDを指定してプロセスを終了する。

```powershell
# 複数のプロセスIDを指定する場合
Stop-Process -Id 1234, 5678
```

出力例:
```
# 成功時は出力なし
```

### -Name

プロセス名でプロセスを終了する。

```powershell
# 複数のプロセス名を指定する場合
Stop-Process -Name notepad, chrome
```

出力例:
```
# 成功時は出力なし
```

### -InputObject

パイプラインからプロセスオブジェクトを受け取り終了する。

```powershell
# メモリ使用量が100MBを超えるプロセスを終了する
Get-Process | Where-Object { $_.WorkingSet -gt 100MB } | Stop-Process -WhatIf
```

出力例:
```
WhatIf: chrome (1234) を停止します
WhatIf: firefox (5678) を停止します
```

### -PassThru

終了したプロセスオブジェクトを出力に返す。

```powershell
Stop-Process -Name notepad -PassThru
```

出力例:
```
Handles  NPM(K)    PM(K)      WS(K)     CPU(s)     Id  SI ProcessName
-------  ------    -----      -----     ------     --  -- -----------
    189      11    18308      31716       0.45   1234   1 notepad
```

### -Force

終了に抵抗するプロセスを強制的に終了する。

```powershell
Stop-Process -Name "stuck_app" -Force
```

出力例:
```
# 強制終了時も通常は出力なし
```

### -WhatIf

実際には終了せず、何が起こるかをシミュレーションする。

```powershell
Stop-Process -Name chrome -WhatIf
```

出力例:
```
WhatIf: chrome (1234) を停止します
WhatIf: chrome (5678) を停止します
```

### -Confirm

実行前に確認を求める。

```powershell
Stop-Process -Name firefox -Confirm
```

出力例:
```
本当に Stop-Process コマンドレットを実行しますか?
firefox (1234) を停止します
[Y] はい(Y)  [A] すべて続行(A)  [N] いいえ(N)  [L] すべて無視(L)  [S] 中断(S)  [?] ヘルプ (既定値は "Y"): 
```

### プロセスフィルタリングと組み合わせ

```powershell
# 特定のユーザーが実行しているプロセスを終了
$processesToKill = Get-Process | Where-Object { $_.Username -like "*user1*" }
$processesToKill | Stop-Process -Confirm
```

出力例:
```
本当に Stop-Process コマンドレットを実行しますか?
notepad (1234) を停止します
[Y] はい(Y)  [A] すべて続行(A)  [N] いいえ(N)  [L] すべて無視(L)  [S] 中断(S)  [?] ヘルプ (既定値は "Y"):
```

## 出力例の見方

1. **出力なし**: 通常、`Stop-Process`は成功時に何も出力しない。これは正常な動作である。

2. **-PassThru オプション使用時**:
   - `Handles`: プロセスが使用しているハンドル数
   - `NPM(K)`: ページングされない（Non-Paged）メモリ（KB単位）
   - `PM(K)`: ページングされるメモリ（KB単位）
   - `WS(K)`: ワーキングセット（実メモリ使用量、KB単位）
   - `CPU(s)`: プロセスの累積CPU時間（秒）
   - `Id`: プロセスID
   - `SI`: セッションID
   - `ProcessName`: プロセスの名前

3. **-WhatIf または -Confirm オプション使用時**:
   - 表示される情報はプロセス名とプロセスIDであり、何をしようとしているかの確認メッセージである。

4. **エラーメッセージ**:
   - 「指定されたプロセスが見つかりません」: 指定したプロセス名またはIDのプロセスが存在しない
   - 「アクセスが拒否されました」: 管理者権限が必要なプロセスを終了しようとした場合

## 使用可能なPowerShellバージョン

`Stop-Process`コマンドレットは以下のPowerShellバージョンで使用可能である：

- PowerShell 1.0以降（すべてのバージョン）
- Windows PowerShell 5.1
- PowerShell Core 6.0
- PowerShell 7.0以降

PowerShell 5.1以降では、すべてのオプションが利用可能である。PowerShell 7.0では、特に`-PassThru`オプションの出力形式が改善されている。

各バージョンでのコマンドの確認方法：

```powershell
# PowerShellバージョンの確認
$PSVersionTable.PSVersion

# コマンドのヘルプを表示
Get-Help Stop-Process -Full
```

出力例:
```
Major  Minor  Build  Revision
-----  -----  -----  --------
5      1      19041  1682

名前
    Stop-Process

概要
    実行中のプロセスを停止します。
...
```