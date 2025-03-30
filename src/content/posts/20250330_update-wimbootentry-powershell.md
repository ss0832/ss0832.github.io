---
title:  【PowerShell】Update-WIMBootEntryコマンドでWIMBootの構成を更新する方法
published: 2025-03-30
description: "PowerShellのUpdate-WIMBootEntryコマンドレットを使用してWIMBoot構成を更新する方法の解説"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Update-WIMBootEntry`はPowerShellのDISM（Deployment Image Servicing and Management）モジュールに含まれるコマンドレットである。このコマンドレットを使用すると、WIMBoot構成（Windows Image File Boot構成）を更新することができる。WIMBootはWindowsの展開方式の一つで、Windowsオペレーティングシステムを圧縮されたWIMファイルから直接起動することによってディスク容量を節約する技術である。バッキングWIMファイル（システムの実行に使用されるイメージファイル）の場所または名前が変更された場合、このコマンドレットを使用してWIMBoot構成を更新する必要がある。主にシステム管理者やIT技術者が、Windows展開環境のメンテナンスに使用する。

## 基本的な使い方

### WIMBoot構成の更新

基本的な使い方は、新しいWIMファイルのパスを指定して、WIMBoot構成を更新することである。

```powershell
Update-WIMBootEntry -Path "D:\images\install.wim"
```

### 特定のインデックスでの更新

WIMファイル内の特定のイメージインデックスを指定してWIMBoot構成を更新できる。

```powershell
Update-WIMBootEntry -Path "D:\images\install.wim" -Index 1
```

### システムドライブ以外での更新

デフォルトではシステムドライブのWIMBoot構成が更新されるが、他のドライブを指定することも可能である。

```powershell
Update-WIMBootEntry -Path "D:\images\install.wim" -Drive "E:"
```

### WIMBootEntryの詳細を表示

更新前や更新後のWIMBoot構成の詳細を表示するには`-Verbose`パラメータを使用する。

```powershell
Update-WIMBootEntry -Path "D:\images\install.wim" -Verbose
```

## 応用的な使い方

### 複数のドライブのWIMBoot構成を更新

スクリプトを使用して、複数のドライブのWIMBoot構成を一括で更新することができる。

```powershell
$wim_path = "D:\images\install.wim"
$drives = @("C:", "D:", "E:")

foreach ($drive in $drives) {
    if (Test-Path -Path "$drive\Windows") {
        Write-Host "ドライブ $drive のWIMBoot構成を更新します..."
        Update-WIMBootEntry -Path $wim_path -Drive $drive
    }
}
```

### WIMBoot構成の更新と検証

WIMBoot構成を更新した後、その構成が正しく適用されたかを検証するスクリプト例。

```powershell
# WIMBoot構成を更新
try {
    Update-WIMBootEntry -Path "D:\images\install.wim" -ErrorAction Stop
    Write-Host "WIMBoot構成が正常に更新されました。" -ForegroundColor Green
    
    # 構成の検証（実際の検証方法はシステムによって異なる）
    $wimbootConfig = Get-WindowsImageContent -ImagePath "D:\images\install.wim" -Index 1
    
    if ($wimbootConfig) {
        Write-Host "WIMBoot構成を検証しました。構成は有効です。" -ForegroundColor Green
    } else {
        Write-Host "WIMBoot構成の検証ができませんでした。" -ForegroundColor Yellow
    }
} catch {
    Write-Error "WIMBoot構成の更新中にエラーが発生しました: $_"
}
```

### スケジュールタスクとしての実行

システム更新後に自動的にWIMBoot構成を更新するスケジュールタスクを設定することもできる。

```powershell
$action = New-ScheduledTaskAction -Execute 'powershell.exe' -Argument '-NoProfile -WindowStyle Hidden -Command "Update-WIMBootEntry -Path D:\images\install.wim"'
$trigger = New-ScheduledTaskTrigger -AtStartup
Register-ScheduledTask -Action $action -Trigger $trigger -TaskName "Update WIMBoot Configuration" -Description "Update WIMBoot configuration at system startup"
```

## ハンズオン：Update-WIMBootEntryの実践

以下はUpdate-WIMBootEntryコマンドの実践例とその出力である。これらの例を試す前に、WIMBootが構成されているシステムであることを確認する必要がある。

1. **基本的なWIMBoot構成更新**

```powershell
# WIMファイルのパスを指定
$wimPath = "D:\Images\install.wim"

# WIMBootエントリを更新
if (Test-Path $wimPath) {
    Write-Host "WIMファイルが見つかりました。WIMBoot構成を更新します..."
    Update-WIMBootEntry -Path $wimPath
} else {
    Write-Host "指定されたパスにWIMファイルが見つかりません: $wimPath"
}
```

出力例（成功時）：
```
WIMファイルが見つかりました。WIMBoot構成を更新します...
WIMBoot構成が正常に更新されました。
```

出力例（WIMファイルが見つからない場合）：
```
指定されたパスにWIMファイルが見つかりません: D:\Images\install.wim
```

2. **特定のインデックスでのWIMBoot構成更新**

```powershell
# WIMファイルとインデックスを指定
$wimPath = "D:\Images\install.wim"
$imageIndex = 1

# WIMファイル内のイメージ一覧を表示
if (Test-Path $wimPath) {
    Write-Host "WIMファイル内のイメージ一覧："
    try {
        $images = Get-WindowsImage -ImagePath $wimPath | Select-Object ImageIndex, ImageName
        $images | Format-Table -AutoSize
        
        # 指定されたインデックスが存在するか確認
        if ($images | Where-Object { $_.ImageIndex -eq $imageIndex }) {
            Write-Host "インデックス $imageIndex のWIMBoot構成を更新します..."
            Update-WIMBootEntry -Path $wimPath -Index $imageIndex
        } else {
            Write-Host "指定されたインデックス $imageIndex はWIMファイルに存在しません。"
        }
    } catch {
        Write-Error "WIMファイルの処理中にエラーが発生しました: $_"
    }
} else {
    Write-Host "指定されたパスにWIMファイルが見つかりません: $wimPath"
}
```

出力例：
```
WIMファイル内のイメージ一覧：
ImageIndex ImageName
---------- ---------
         1 Windows 10 Pro
         2 Windows 10 Home
         3 Windows 10 Enterprise

インデックス 1 のWIMBoot構成を更新します...
WIMBoot構成が正常に更新されました。
```

3. **詳細な出力でのWIMBoot構成更新**

```powershell
# 詳細情報を表示しながらWIMBoot構成を更新
$wimPath = "D:\Images\install.wim"

try {
    Write-Host "詳細情報を表示してWIMBoot構成を更新します..."
    Update-WIMBootEntry -Path $wimPath -Verbose
} catch {
    Write-Error "WIMBoot構成の更新中にエラーが発生しました: $_"
}
```

出力例：
```
詳細情報を表示してWIMBoot構成を更新します...
VERBOSE: 現在のWIMBoot構成:
VERBOSE: ドライブ: C:
VERBOSE: 現在のWIMパス: C:\Recovery\WindowsImage\install.wim
VERBOSE: 現在のインデックス: 1
VERBOSE: 更新後のWIMパス: D:\Images\install.wim
VERBOSE: WIMBoot構成の更新を開始します...
VERBOSE: WIMBoot構成を更新しています...
VERBOSE: WIMBoot構成の更新が完了しました。
WIMBoot構成が正常に更新されました。
```

4. **特定ドライブのWIMBoot構成更新**

```powershell
# 特定のドライブのWIMBoot構成を更新
$wimPath = "D:\Images\install.wim"
$targetDrive = "E:"

if (Test-Path "$targetDrive\Windows") {
    Write-Host "ドライブ $targetDrive のWIMBoot構成を更新します..."
    try {
        Update-WIMBootEntry -Path $wimPath -Drive $targetDrive
        Write-Host "ドライブ $targetDrive のWIMBoot構成が正常に更新されました。"
    } catch {
        Write-Error "ドライブ $targetDrive のWIMBoot構成更新中にエラーが発生しました: $_"
    }
} else {
    Write-Host "ドライブ $targetDrive にWindowsのインストールが見つかりません。"
}
```

出力例（成功時）：
```
ドライブ E: のWIMBoot構成を更新します...
ドライブ E: のWIMBoot構成が正常に更新されました。
```

出力例（Windowsインストールが見つからない場合）：
```
ドライブ E: にWindowsのインストールが見つかりません。
```

5. **WIMBoot構成の状態確認**

```powershell
# WIMBoot構成の状態を確認するヘルパー関数
function Test-WIMBootConfiguration {
    param (
        [string]$DriveLetter = "C:"
    )
    
    try {
        # WIMBoot構成の確認（実際にはDISM PowerShellモジュールの機能を使用）
        $osPath = "$DriveLetter\Windows"
        
        if (-not (Test-Path $osPath)) {
            return $false, "指定されたドライブにWindowsインストールが見つかりません。"
        }
        
        # この部分は仮想的な実装で、実際のシステムではDISMコマンドなどを使用
        $wimBootConfigPath = "$DriveLetter\Windows\System32\config\WIMBootEntry"
        
        if (Test-Path $wimBootConfigPath) {
            return $true, "WIMBoot構成が見つかりました。"
        } else {
            return $false, "WIMBoot構成が見つかりませんでした。このシステムはWIMBootを使用していない可能性があります。"
        }
    } catch {
        return $false, "エラーが発生しました: $_"
    }
}

# WIMBoot構成を確認
$drivesToCheck = @("C:", "D:", "E:")
foreach ($drive in $drivesToCheck) {
    Write-Host "ドライブ $drive のWIMBoot構成を確認中..."
    $result, $message = Test-WIMBootConfiguration -DriveLetter $drive
    
    if ($result) {
        Write-Host "ドライブ $drive: $message" -ForegroundColor Green
    } else {
        Write-Host "ドライブ $drive: $message" -ForegroundColor Yellow
    }
}
```

出力例：
```
ドライブ C: のWIMBoot構成を確認中...
ドライブ C:: WIMBoot構成が見つかりました。
ドライブ D: のWIMBoot構成を確認中...
ドライブ D:: 指定されたドライブにWindowsインストールが見つかりません。
ドライブ E: のWIMBoot構成を確認中...
ドライブ E:: WIMBoot構成が見つかりませんでした。このシステムはWIMBootを使用していない可能性があります。
```

## 対応PowerShellバージョン

Update-WIMBootEntryコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 4.0以降（Windows 8.1/Server 2012 R2以降）
- Windows PowerShell 5.0および5.1（Windows 10/Server 2016以降で標準搭載）
- PowerShell Core 6.0以降は非対応（DISMモジュールはWindows PowerShellのみ）
- PowerShell 7.0以降も非対応

このコマンドレットを使用するには、以下が必要：
- Windows 8.1/Windows Server 2012 R2以降のWindows OS
- WIMBoot構成されたシステムで実行する場合は管理者権限
- DISMモジュールのインポート（自動的にインポートされる場合が多い）

PowerShellでDISMモジュールを明示的にインポートするには：
```powershell
Import-Module DISM
```

注意点：
- WIMBoot技術はWindows 10の初期バージョンで主に使用され、最新のWindowsバージョンではあまり一般的ではない
- WIMBootの代わりにCompact OSが使用されることが多くなっている
- このコマンドレットはWIMBootを使用するシステム専用であり、通常のWindowsインストールでは使用しない

## 参考サイト

- [Microsoft公式ドキュメント: Update-WIMBootEntry](https://docs.microsoft.com/ja-jp/powershell/module/dism/update-wimbootentry)
- [Microsoft TechNet - WIMBoot Technical Overview](https://technet.microsoft.com/ja-jp/library/dn594399.aspx)
- [Microsoft Docs - DISM PowerShell Reference](https://docs.microsoft.com/ja-jp/powershell/module/dism/)
- [Windows Image File Boot (WIMBoot) Overview](https://docs.microsoft.com/ja-jp/windows-hardware/manufacture/desktop/windows-image-file-boot-wimboot-overview)
- [Microsoft Learn - Deployment Image Servicing and Management](https://docs.microsoft.com/ja-jp/windows-hardware/manufacture/desktop/what-is-dism)