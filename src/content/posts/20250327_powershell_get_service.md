---
title:  【PowerShell】Get-Serviceコマンドの基本と応用
published: 2025-03-27
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-27

## 概要

PowerShellの`Get-Service`コマンドレットは、ローカルおよびリモートコンピュータ上のWindowsサービスに関する情報を取得するためのコマンドである。サービスの状態、スタートアップの種類、実行中プロセスなどの情報を確認できる。システム管理者やトラブルシューティングを行うユーザーにとって非常に役立つツールである。

## 基本的な使い方

### すべてのサービスを取得する

```powershell
# コンピューターで実行されているすべてのサービスを取得
Get-Service
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Stopped  AarSvc_828cab      Agent Activation Runtime_828cab
Running  AdobeARMservice    Adobe Acrobat Update Service
Stopped  AJRouter           AllJoyn Router Service
Running  AppIDSvc           Application Identity
Running  Appinfo            Application Information
...
```

### 特定のサービスを名前で取得する

```powershell
# 名前でサービスを取得
Get-Service -Name "wuauserv"
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  wuauserv           Windows Update
```

### 表示名でサービスを取得する

```powershell
# 表示名でサービスを取得
Get-Service -DisplayName "*update*"
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  wuauserv           Windows Update
Running  AdobeARMservice    Adobe Acrobat Update Service
Stopped  edgeupdate         Microsoft Edge Update Service
```

## 全オプションに対するハンズオン

### -Name

サービスの名前を指定して情報を取得する。ワイルドカード文字も使用可能。

```powershell
# 複数のサービス名を指定
Get-Service -Name "wuauserv", "bits"

# ワイルドカードを使用
Get-Service -Name "w*"
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  wuauserv           Windows Update
Running  bits               Background Intelligent Transfer Service

Status   Name               DisplayName
------   ----               -----------
Running  W32Time            Windows Time
Running  WaaSMedicSvc       Windows Update Medic Service
Running  WalletService      WalletService
...
```

### -DisplayName

サービスの表示名を指定して情報を取得する。ワイルドカード文字も使用可能。

```powershell
# 表示名の一部を指定
Get-Service -DisplayName "*network*"
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  lmhosts            TCP/IP NetBIOS Helper
Running  nsi                Network Store Interface Service
Running  netprofm           Network List Service
...
```

### -DependentServices

指定したサービスに依存する他のサービスを取得する。

```powershell
# RpcSsサービスに依存するサービスを取得
Get-Service -Name "RpcSs" -DependentServices
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  AudioSrv           Windows Audio
Running  BITS               Background Intelligent Transfer Service
Running  CryptSvc           Cryptographic Services
...
```

### -RequiredServices

指定したサービスが依存する他のサービスを取得する。

```powershell
# Windefendサービスが依存するサービスを取得
Get-Service -Name "Windefend" -RequiredServices
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  RpcSs              Remote Procedure Call (RPC)
Running  BFE                Base Filtering Engine
...
```

### -Include

特定のサービスを含める場合に使用する。基本的なフィルタリングの後に適用される。

```powershell
# すべてのサービスを取得し、その中から"win"で始まるサービスのみを含める
Get-Service -Include "win*"
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  WinDefend          Windows Defender Service
Running  WinHttpAutoProxySvc WinHTTP Web Proxy Auto-Discovery Service
...
```

### -Exclude

特定のサービスを除外する場合に使用する。基本的なフィルタリングの後に適用される。

```powershell
# すべてのサービスを取得し、その中から"w"で始まるサービスを除外
Get-Service -Exclude "w*"
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  AdobeARMservice    Adobe Acrobat Update Service
Running  AppIDSvc           Application Identity
...
```

### -ComputerName

リモートコンピュータ上のサービスを取得する。

```powershell
# リモートコンピュータのサービスを取得
Get-Service -ComputerName "Server01" -Name "spooler"
```

出力例:
```
Status   Name               DisplayName                            PSComputerName
------   ----               -----------                            --------------
Running  spooler            Print Spooler                          Server01
```

### -InputObject

パイプラインから渡されたサービスオブジェクトを処理する。

```powershell
# 実行中のサービスだけを取得し、その中からWで始まるサービスをフィルタリング
Get-Service | Where-Object {$_.Status -eq "Running"} | Where-Object {$_.Name -like "W*"}
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  W32Time            Windows Time
Running  WaaSMedicSvc       Windows Update Medic Service
Running  WinDefend          Windows Defender Service
...
```

## フィルタリングと並べ替え

### 実行中のサービスのみ表示

```powershell
# Statusプロパティでフィルタリング
Get-Service | Where-Object {$_.Status -eq "Running"}
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  AdobeARMservice    Adobe Acrobat Update Service
Running  AppIDSvc           Application Identity
Running  Appinfo            Application Information
...
```

### サービスの状態で並べ替え

```powershell
# 複数の条件で並べ替え
Get-Service | Sort-Object -Property Status, Name
```

出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  AdobeARMservice    Adobe Acrobat Update Service
Running  AppIDSvc           Application Identity
...
Stopped  AarSvc_828cab      Agent Activation Runtime_828cab
Stopped  AJRouter           AllJoyn Router Service
...
```

### サービスの詳細情報を表示

```powershell
# Format-Listで詳細表示
Get-Service -Name "wuauserv" | Format-List *
```

出力例:
```
Name                : wuauserv
RequiredServices    : {rpcss}
CanPauseAndContinue : False
CanShutdown         : False
CanStop             : True
DisplayName         : Windows Update
DependentServices   : {}
MachineName         : .
ServiceName         : wuauserv
ServicesDependedOn  : {rpcss}
ServiceHandle       : SafeServiceHandle
Status              : Running
ServiceType         : Win32OwnProcess
StartType           : Automatic
Site                :
Container           :
```

## 出力例の見方

1. **基本的な出力フォーマット**:
   - `Status`: サービスの状態（Running/Stopped/Paused）
   - `Name`: サービスの名前（プログラム的に参照される内部名）
   - `DisplayName`: サービスの表示名（ユーザーフレンドリーな名前）

2. **Format-Listを使用した詳細出力**:
   - `Name`: サービスの内部名称（スクリプトやコマンドラインから参照するための短い名前）
   - `ServiceName`: サービスの正式な内部名称（`Name`と同じ値を持つ）
   - `DisplayName`: サービスの表示名（ユーザーインターフェース上で表示される説明的な名前）
   - `Status`: サービスの現在の状態（Running/Stopped/Paused/StartPending/StopPending など）
   - `StartType`: サービスの起動種類
     - `Automatic`: システム起動時に自動的に開始
     - `AutomaticDelayedStart`: システム起動後、遅延して自動的に開始
     - `Manual`: 手動または他のサービス・アプリケーションによって開始される必要がある
     - `Disabled`: サービスが無効化されており、開始できない
   - `ServiceType`: サービスの種類
     - `Win32OwnProcess`: 独自のプロセスで実行される Win32 サービス
     - `Win32ShareProcess`: 他のサービスと同じプロセスで実行される Win32 サービス
     - `KernelDriver`: カーネルデバイスドライバー
     - `FileSystemDriver`: ファイルシステムドライバー
     - `InteractiveProcess`: 対話型デスクトップとやり取りできるサービス
   - `RequiredServices`/`ServicesDependedOn`: このサービスが依存している他のサービスのリスト（これらが動作していないと、このサービスは起動できない）
   - `DependentServices`: このサービスに依存している他のサービスのリスト（このサービスが停止すると、これらも停止する可能性がある）
   - `CanPauseAndContinue`: サービスを一時停止して後で再開できるかどうか（True/False）
   - `CanStop`: サービスを停止できるかどうか（True/False）
   - `CanShutdown`: システムシャットダウン時に通知を受け取るかどうか（True/False）
   - `MachineName`: サービスが実行されているコンピュータ名（`.`はローカルコンピュータを示す）
   - `ServiceHandle`: サービスへのハンドル（システム内部で使用）
   - `Site`: Webサイト情報（主にIISサービスで使用）
   - `Container`: サービスが属するコンテナ情報

3. **-ComputerNameオプション使用時**:
   - `PSComputerName`: サービス情報を取得したコンピュータ名

## 使用可能なPowerShellバージョン

`Get-Service`コマンドレットは以下のPowerShellバージョンで使用可能である：

- PowerShell 1.0以降（基本機能）
- PowerShell 2.0以降（-ComputerNameを含むすべてのリモート操作機能）
- PowerShell 3.0以降（-Includeおよび-Excludeパラメータを含むフィルタリング機能強化）
- PowerShell 5.1および7.x（Windows PowerShellとPowerShell Core）

PowerShell 7.x（Core）では、基本的な機能は維持されているが、一部のサービス管理機能はWindowsプラットフォームに限定される点に注意。

バージョン確認方法：

```powershell
# PowerShellバージョンの確認
$PSVersionTable.PSVersion

# コマンドのヘルプ表示
Get-Help Get-Service -Full
```

出力例:
```
Major  Minor  Build  Revision
-----  -----  -----  --------
5      1      19041  1682

名前
    Get-Service

概要
    コンピューター上のサービスを取得します。
...
```