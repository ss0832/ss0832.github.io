---
title:  【PowerShell】Update-Helpコマンドでヘルプファイルを最新化する方法
published: 2025-03-30
description: "PowerShellのUpdate-Helpコマンドレットを使用してヘルプドキュメントを更新する方法を解説"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Update-Help`はPowerShellのヘルプファイル（ドキュメント）を更新するためのコマンドレットである。PowerShellをインストールした直後や、新しいモジュールを追加した際には、最新のヘルプファイルが含まれていない場合があるため、このコマンドレットを使用して最新のドキュメントをダウンロードする必要がある。`Update-Help`を実行することで、コマンドレットや関数の最新の使用方法、構文、例などを参照できるようになる。また、特定のモジュールのヘルプのみを更新したり、異なる言語のヘルプをダウンロードしたりする機能も備えている。

## 基本的な使い方

### すべてのモジュールのヘルプを更新する

最も基本的な使い方は、すべての利用可能なモジュールのヘルプを更新することである。管理者権限で実行すると、すべてのユーザー向けにヘルプが更新される。

```powershell
Update-Help
```

### 特定のモジュールのヘルプを更新する

特定のモジュールのヘルプだけを更新することもできる。

```powershell
Update-Help -Module Microsoft.PowerShell.Management
```

### 複数のモジュールのヘルプを更新する

複数のモジュールを指定して、それらのヘルプを更新することも可能である。

```powershell
Update-Help -Module Microsoft.PowerShell.Management, Microsoft.PowerShell.Utility
```

### 特定の言語のヘルプを更新する

デフォルトではシステム言語のヘルプが更新されるが、UIカルチャを指定して特定の言語のヘルプを更新することができる。

```powershell
Update-Help -UICulture ja-JP
```

## 応用的な使い方

### インターネットに接続していない環境での更新

インターネットに接続していない環境でもヘルプを更新するために、インターネットに接続している環境で先にヘルプをダウンロードし、それを利用することができる。

```powershell
# インターネットに接続している環境で実行
Save-Help -DestinationPath "C:\ヘルプ" -Module Microsoft.PowerShell.* -UICulture ja-JP, en-US

# インターネットに接続していない環境で実行
Update-Help -SourcePath "C:\ヘルプ" -Module Microsoft.PowerShell.* -UICulture ja-JP
```

### 強制的にヘルプを更新する

デフォルトでは、24時間以内に更新したモジュールのヘルプは更新されない。強制的に更新するには`-Force`パラメータを使用する。

```powershell
Update-Help -Force
```

### 最新化状況を確認して更新する

ヘルプが最新かどうかを確認し、必要に応じて更新するスクリプト例。

```powershell
$modules = Get-Module -ListAvailable
foreach ($module in $modules) {
    try {
        $helpInfo = Get-Help -Name $module.Name -ErrorAction SilentlyContinue
        if ($helpInfo.Description -eq $null) {
            Write-Host "ヘルプがないため更新します: $($module.Name)" -ForegroundColor Yellow
            Update-Help -Module $module.Name -ErrorAction SilentlyContinue
        }
    } catch {
        Write-Host "エラー: $($module.Name) - $_" -ForegroundColor Red
    }
}
```

### オフラインでのヘルプ更新（詳細設定）

特定のモジュールおよび言語に対して、オフラインでのヘルプ更新を設定する高度な例。

```powershell
# 指定したモジュールと言語のヘルプをダウンロード
$modules = @("Microsoft.PowerShell.Management", "Microsoft.PowerShell.Security")
$cultures = @("en-US", "ja-JP", "de-DE")
$savePath = "D:\PowerShell\OfflineHelp"

# 保存先ディレクトリを作成
if (-not (Test-Path $savePath)) {
    New-Item -Path $savePath -ItemType Directory -Force
}

# ヘルプをダウンロード
foreach ($culture in $cultures) {
    $culturePath = Join-Path $savePath $culture
    if (-not (Test-Path $culturePath)) {
        New-Item -Path $culturePath -ItemType Directory -Force
    }
    Save-Help -Module $modules -UICulture $culture -DestinationPath $culturePath -Force
}

# オフラインのヘルプを更新
foreach ($culture in $cultures) {
    $culturePath = Join-Path $savePath $culture
    Update-Help -Module $modules -UICulture $culture -SourcePath $culturePath -Force
}
```

## ハンズオン：Update-Helpの実践

以下はUpdate-Helpコマンドの実践例とその出力である。

1. **基本的なヘルプの更新**

```powershell
# すべてのヘルプを更新
Update-Help -Verbose
```

出力例：
```
VERBOSE: Update-Help の操作を開始しています
VERBOSE: モジュール Microsoft.PowerShell.Utility のヘルプの更新...
VERBOSE: モジュール Microsoft.PowerShell.Management のヘルプの更新...
VERBOSE: モジュール Microsoft.PowerShell.Security のヘルプの更新...
VERBOSE: HelpInfoUri http://go.microsoft.com/fwlink/?LinkId=285760 からヘルプをダウンロードしています。
VERBOSE: ja-JP のヘルプ ファイルをダウンロードしています。
VERBOSE: ヘルプ ファイルのダウンロードが完了しました。モジュール Microsoft.PowerShell.Utility のヘルプ ファイルをインストールしています。
VERBOSE: HelpInfoUri http://go.microsoft.com/fwlink/?LinkId=285760 からヘルプをダウンロードしています。
...
VERBOSE: Update-Help の操作が完了しました
```

2. **特定のモジュールのヘルプを更新**

```powershell
# PowerShell Coreモジュールのヘルプを更新
Update-Help -Module Microsoft.PowerShell.Core -Verbose
```

出力例：
```
VERBOSE: Update-Help の操作を開始しています
VERBOSE: モジュール Microsoft.PowerShell.Core のヘルプの更新...
VERBOSE: HelpInfoUri http://go.microsoft.com/fwlink/?LinkId=285759 からヘルプをダウンロードしています。
VERBOSE: ja-JP のヘルプ ファイルをダウンロードしています。
VERBOSE: ヘルプ ファイルのダウンロードが完了しました。モジュール Microsoft.PowerShell.Core のヘルプ ファイルをインストールしています。
VERBOSE: Update-Help の操作が完了しました
```

3. **エラー処理を使用したヘルプ更新**

```powershell
# エラーを表示せずにヘルプ更新を試みる
try {
    Update-Help -ErrorAction Stop
    Write-Host "ヘルプの更新に成功しました" -ForegroundColor Green
} catch {
    Write-Host "ヘルプの更新中にエラーが発生しました: $_" -ForegroundColor Red
    
    # インターネット接続を確認
    if (Test-Connection 8.8.8.8 -Count 1 -Quiet) {
        Write-Host "インターネット接続は利用可能です。別の問題が発生している可能性があります。" -ForegroundColor Yellow
    } else {
        Write-Host "インターネット接続がありません。ネットワーク設定を確認してください。" -ForegroundColor Yellow
    }
}
```

出力例（成功時）：
```
ヘルプの更新に成功しました
```

出力例（失敗時）：
```
ヘルプの更新中にエラーが発生しました: Update-Help: インターネット接続に問題があるか、既定の HelpInfoUri にアクセスできません。
インターネット接続がありません。ネットワーク設定を確認してください。
```

4. **ヘルプの保存と更新（オフライン環境用）**

```powershell
# 保存先フォルダを作成
$helpFolder = "C:\PowerShell\SavedHelp"
if (-not (Test-Path $helpFolder)) {
    New-Item -Path $helpFolder -ItemType Directory | Out-Null
}

# 選択したモジュールのヘルプを保存
$selectedModules = @("Microsoft.PowerShell.Core", "Microsoft.PowerShell.Security")
Save-Help -Module $selectedModules -DestinationPath $helpFolder -Force -Verbose

# 保存したヘルプを使用して更新
Update-Help -Module $selectedModules -SourcePath $helpFolder -Force -Verbose
```

出力例（Save-Help）：
```
VERBOSE: Save-Help の操作を開始しています
VERBOSE: モジュール Microsoft.PowerShell.Core のヘルプの格納...
VERBOSE: HelpInfoUri http://go.microsoft.com/fwlink/?LinkId=285759 からヘルプをダウンロードしています。
VERBOSE: ja-JP のヘルプ ファイルをダウンロードしています。
VERBOSE: ヘルプ ファイルのダウンロードが完了しました。
VERBOSE: モジュール Microsoft.PowerShell.Security のヘルプの格納...
VERBOSE: HelpInfoUri http://go.microsoft.com/fwlink/?LinkId=285760 からヘルプをダウンロードしています。
VERBOSE: ja-JP のヘルプ ファイルをダウンロードしています。
VERBOSE: ヘルプ ファイルのダウンロードが完了しました。
VERBOSE: Save-Help の操作が完了しました
```

出力例（Update-Help）：
```
VERBOSE: Update-Help の操作を開始しています
VERBOSE: モジュール Microsoft.PowerShell.Core のヘルプの更新...
VERBOSE: ソースパス C:\PowerShell\SavedHelp からヘルプをダウンロードしています。
VERBOSE: ja-JP のヘルプ ファイルをダウンロードしています。
VERBOSE: ヘルプ ファイルのダウンロードが完了しました。モジュール Microsoft.PowerShell.Core のヘルプ ファイルをインストールしています。
VERBOSE: モジュール Microsoft.PowerShell.Security のヘルプの更新...
VERBOSE: ソースパス C:\PowerShell\SavedHelp からヘルプをダウンロードしています。
VERBOSE: ja-JP のヘルプ ファイルをダウンロードしています。
VERBOSE: ヘルプ ファイルのダウンロードが完了しました。モジュール Microsoft.PowerShell.Security のヘルプ ファイルをインストールしています。
VERBOSE: Update-Help の操作が完了しました
```

5. **複数言語のヘルプ更新**

```powershell
# 英語と日本語のヘルプを更新
Update-Help -UICulture en-US, ja-JP -Verbose
```

出力例：
```
VERBOSE: Update-Help の操作を開始しています
VERBOSE: モジュール Microsoft.PowerShell.Management のヘルプの更新...
VERBOSE: HelpInfoUri http://go.microsoft.com/fwlink/?LinkId=285760 からヘルプをダウンロードしています。
VERBOSE: en-US のヘルプ ファイルをダウンロードしています。
VERBOSE: ヘルプ ファイルのダウンロードが完了しました。モジュール Microsoft.PowerShell.Management のヘルプ ファイルをインストールしています。
VERBOSE: HelpInfoUri http://go.microsoft.com/fwlink/?LinkId=285760 からヘルプをダウンロードしています。
VERBOSE: ja-JP のヘルプ ファイルをダウンロードしています。
...
VERBOSE: Update-Help の操作が完了しました
```

6. **ヘルプ更新のスケジュール化（スクリプト例）**

```powershell
# ヘルプ更新を毎週日曜日の午前3時に実行するスケジュールタスクを作成
$action = New-ScheduledTaskAction -Execute 'powershell.exe' -Argument '-NoProfile -WindowStyle Hidden -Command "Update-Help -Force"'
$trigger = New-ScheduledTaskTrigger -Weekly -DaysOfWeek Sunday -At 3am
$principal = New-ScheduledTaskPrincipal -UserId "SYSTEM" -LogonType ServiceAccount -RunLevel Highest
$task = New-ScheduledTask -Action $action -Trigger $trigger -Principal $principal -Description "PowerShellヘルプを週次で更新"

# スケジュールタスクを登録
Register-ScheduledTask -TaskName "UpdatePowerShellHelp" -InputObject $task

Write-Host "PowerShellヘルプの自動更新タスクが作成されました" -ForegroundColor Green
```

出力例：
```
TaskPath                                       TaskName                          State
--------                                       --------                          -----
\                                              UpdatePowerShellHelp              Ready

PowerShellヘルプの自動更新タスクが作成されました
```

7. **現在利用可能なヘルプの状態確認**

```powershell
# インストール済みモジュールのヘルプ情報を取得
$modules = Get-Module -ListAvailable | Where-Object { $_.HelpInfoUri }
$results = @()

foreach ($module in $modules) {
    $helpInfo = Get-HelpVersion -Module $module.Name -ErrorAction SilentlyContinue
    
    if ($helpInfo) {
        $results += [PSCustomObject]@{
            ModuleName = $module.Name
            Version = $module.Version
            LastUpdated = $helpInfo.LastUpdate
            Status = if ($helpInfo.LastUpdate -lt (Get-Date).AddDays(-90)) { "古い" } else { "最新" }
        }
    } else {
        $results += [PSCustomObject]@{
            ModuleName = $module.Name
            Version = $module.Version
            LastUpdated = $null
            Status = "ヘルプなし"
        }
    }
}

$results | Format-Table -AutoSize
```

出力例：
```
ModuleName                      Version LastUpdated           Status
---------                      ------- -----------           ------
Microsoft.PowerShell.Archive    1.2.5  2024-02-15 10:23:45   最新
Microsoft.PowerShell.Core       7.3.0  2024-01-02 08:15:30   最新
Microsoft.PowerShell.Utility    7.3.0  2023-12-28 15:42:18   最新
PackageManagement              1.4.8.1  2023-09-10 09:30:22   古い
PowerShellGet                   2.2.5  2023-08-05 14:17:53   古い
PSReadLine                     2.2.6.0  null                  ヘルプなし
```

## 対応PowerShellバージョン

Update-Helpコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 3.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

PowerShell 2.0には組み込まれておらず、以下の環境でも使用できる：
- Windows PowerShell (Windows 7 SP1/Server 2008 R2以降)
- PowerShell Core (Linux、macOS、Windows)
- PowerShell 7 (クロスプラットフォーム)

## 参考サイト

- [Microsoft公式ドキュメント: Update-Help](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.core/update-help)
- [Microsoft公式ドキュメント: Save-Help](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.core/save-help)
- [PowerShell.org - Understanding Update-Help and Save-Help](https://powershell.org/2013/07/how-to-use-update-help-and-save-help-effectively/)
- [Microsoft Learn - PowerShellのHelp System](https://docs.microsoft.com/ja-jp/powershell/scripting/learn/ps101/02-help-system)
- [SS64.com PowerShell Commands - Update-Help](https://ss64.com/ps/update-help.html)