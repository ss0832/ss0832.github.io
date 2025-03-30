---
title:  【PowerShell】Use-WindowsUnattendコマンドで無人応答ファイルを適用する方法
published: 2025-03-30
description: "PowerShellのUse-WindowsUnattendコマンドレットを使ったWindows無人セットアップの基本と活用法"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Use-WindowsUnattend`はWindows PowerShellのDISM（Deployment Image Servicing and Management）モジュールに含まれるコマンドレットである。このコマンドレットを使用すると、オフライン/オンラインのWindows環境に無人応答ファイル（unattend.xml）を適用することができる。無人応答ファイルはWindows展開時の設定を自動化するためのXMLフォーマットの設定ファイルであり、`Use-WindowsUnattend`を使うことで、オペレーティングシステムの展開段階だけでなく、既にインストール済みのシステムへも設定を適用できる。主にシステム管理者がWindows環境の構成を一貫して行うために使用される。

## 基本的な使い方

### オンラインWindows環境への無人応答ファイルの適用

現在実行中のWindows環境に無人応答ファイルを適用する最も基本的な使い方である。

```powershell
Use-WindowsUnattend -Path C:\unattend\unattend.xml -Online
```

### オフラインWindows環境への無人応答ファイルの適用

マウントされたオフラインのWindowsイメージに無人応答ファイルを適用できる。

```powershell
Use-WindowsUnattend -Path C:\unattend\unattend.xml -WindowsPath D:\mounted_image
```

### 特定のパスにある無人応答ファイルの適用

無人応答ファイルを適用する際に、特定のパスを指定することができる。

```powershell
Use-WindowsUnattend -Path C:\unattend\specialized.xml -WindowsPath E:\windows_image -ScratchDirectory C:\temp
```

### 進行状況の表示

`-NoRestart`パラメータを使用して、コマンド完了後に自動的に再起動しないようにする。

```powershell
Use-WindowsUnattend -Path C:\unattend\unattend.xml -Online -NoRestart
```

## 応用的な使い方

### 複数の無人応答ファイルの適用

複数の無人応答ファイルを順番に適用して、異なる構成を組み合わせることができる。

```powershell
Use-WindowsUnattend -Path C:\unattend\base.xml -Online
Use-WindowsUnattend -Path C:\unattend\apps.xml -Online
Use-WindowsUnattend -Path C:\unattend\network.xml -Online -NoRestart
```

### Windows PEでのオフラインイメージへの適用

Windows PE環境からオフラインのWindowsイメージに無人応答ファイルを適用する例である。

```powershell
# まずWindowsイメージをマウント
Mount-WindowsImage -ImagePath D:\sources\install.wim -Index 1 -Path E:\mounted_windows

# 無人応答ファイルを適用
Use-WindowsUnattend -Path C:\unattend\sysprep.xml -WindowsPath E:\mounted_windows

# 変更を保存
Dismount-WindowsImage -Path E:\mounted_windows -Save
```

### LogPathの指定

ログファイルの出力先を指定することができる。

```powershell
Use-WindowsUnattend -Path C:\unattend\unattend.xml -Online -LogPath C:\logs\unattend_application.log
```

### システム言語の設定例

言語やキーボードレイアウトなどを設定した無人応答ファイルを適用する例。

```powershell
Use-WindowsUnattend -Path C:\unattend\language_settings.xml -Online -NoRestart
```

## ハンズオン：Use-WindowsUnattendの実践

以下はUse-WindowsUnattendコマンドの実践例とその出力である。

1. **基本的な無人応答ファイルの適用（オンライン）**

まずは基本的な無人応答ファイルを作成する。

```powershell
# サンプルの無人応答ファイルを作成
$unattendXml = @"
<?xml version="1.0" encoding="utf-8"?>
<unattend xmlns="urn:schemas-microsoft-com:unattend">
    <settings pass="oobeSystem">
        <component name="Microsoft-Windows-Shell-Setup" processorArchitecture="amd64" publicKeyToken="31bf3856ad364e35" language="neutral" versionScope="nonSxS" xmlns:wcm="http://schemas.microsoft.com/WMIConfig/2002/State">
            <OOBE>
                <HideEULAPage>true</HideEULAPage>
                <HideWirelessSetupInOOBE>true</HideWirelessSetupInOOBE>
                <ProtectYourPC>3</ProtectYourPC>
            </OOBE>
            <TimeZone>Tokyo Standard Time</TimeZone>
        </component>
    </settings>
</unattend>
"@

# ファイルとして保存
$unattendPath = "C:\temp\basic_unattend.xml"
New-Item -Path $unattendPath -ItemType File -Force
$unattendXml | Out-File -FilePath $unattendPath -Encoding utf8

# 無人応答ファイルを適用
Use-WindowsUnattend -Path $unattendPath -Online -NoRestart
```

出力例：
```
Path          : C:\temp\basic_unattend.xml
Online        : True
NoRestart     : True
LogPath       : C:\Windows\Logs\Dism\dism.log
ScratchDir    : 
LogLevel      : 
RestartNeeded : False
```

2. **特定のコンポーネントを対象とした無人応答ファイルの適用**

```powershell
# 特定コンポーネントの無人応答ファイル
$componentUnattendXml = @"
<?xml version="1.0" encoding="utf-8"?>
<unattend xmlns="urn:schemas-microsoft-com:unattend">
    <settings pass="specialize">
        <component name="Microsoft-Windows-IE-ESC" processorArchitecture="amd64" publicKeyToken="31bf3856ad364e35" language="neutral" versionScope="nonSxS">
            <IEHardenAdmin>false</IEHardenAdmin>
            <IEHardenUser>false</IEHardenUser>
        </component>
    </settings>
</unattend>
"@

# ファイルとして保存
$componentUnattendPath = "C:\temp\ie_esc_unattend.xml"
New-Item -Path $componentUnattendPath -ItemType File -Force
$componentUnattendXml | Out-File -FilePath $componentUnattendPath -Encoding utf8

# 無人応答ファイルを適用
$logPath = "C:\temp\unattend_application.log"
Use-WindowsUnattend -Path $componentUnattendPath -Online -NoRestart -LogPath $logPath
```

出力例：
```
Path          : C:\temp\ie_esc_unattend.xml
Online        : True
NoRestart     : True
LogPath       : C:\temp\unattend_application.log
ScratchDir    : 
LogLevel      : 
RestartNeeded : False
```

3. **ログレベルを指定した無人応答ファイルの適用**

```powershell
# ログレベルを指定して無人応答ファイルを適用
Use-WindowsUnattend -Path $unattendPath -Online -NoRestart -LogLevel 4
```

出力例：
```
Path          : C:\temp\basic_unattend.xml
Online        : True
NoRestart     : True
LogPath       : C:\Windows\Logs\Dism\dism.log
ScratchDir    : 
LogLevel      : 4
RestartNeeded : False
```

4. **スクラッチディレクトリを指定した無人応答ファイルの適用**

```powershell
# テンポラリディレクトリを作成
$scratchDir = "C:\temp\scratch"
New-Item -Path $scratchDir -ItemType Directory -Force

# スクラッチディレクトリを指定して無人応答ファイルを適用
Use-WindowsUnattend -Path $unattendPath -Online -NoRestart -ScratchDirectory $scratchDir
```

出力例：
```
Path          : C:\temp\basic_unattend.xml
Online        : True
NoRestart     : True
LogPath       : C:\Windows\Logs\Dism\dism.log
ScratchDir    : C:\temp\scratch
LogLevel      : 
RestartNeeded : False
```

5. **オフラインWindowsイメージへの適用（管理者権限必要）**

このサンプルは管理者権限が必要で、Windows ISOイメージまたはWIMファイルが必要。

```powershell
# 一時マウントポイントを作成
$mountPath = "C:\temp\mount"
New-Item -Path $mountPath -ItemType Directory -Force -ErrorAction SilentlyContinue

# Windowsイメージをマウント（この例ではWIMファイルを使用）
try {
    # WIMファイルが存在する場合のみ実行
    if (Test-Path "D:\sources\install.wim") {
        Mount-WindowsImage -ImagePath "D:\sources\install.wim" -Index 1 -Path $mountPath

        # 無人応答ファイルを適用
        Use-WindowsUnattend -Path $unattendPath -WindowsPath $mountPath

        # 変更を保存してマウント解除
        Dismount-WindowsImage -Path $mountPath -Save
        
        Write-Host "オフラインイメージに無人応答ファイルを適用しました。" -ForegroundColor Green
    } else {
        Write-Host "WIMファイルが見つかりません。このサンプルはスキップします。" -ForegroundColor Yellow
    }
} catch {
    Write-Host "エラーが発生しました: $_" -ForegroundColor Red
    
    # エラー発生時はマウント解除を試みる
    if (Test-Path $mountPath) {
        try {
            Dismount-WindowsImage -Path $mountPath -Discard
        } catch {
            Write-Host "マウント解除中にエラーが発生しました。" -ForegroundColor Red
        }
    }
}
```

出力例（WIMファイルがある場合）：
```
マウント操作は正常に完了しました。
Path          : C:\temp\basic_unattend.xml
WindowsPath   : C:\temp\mount
LogPath       : C:\Windows\Logs\Dism\dism.log
ScratchDir    : 
LogLevel      : 
RestartNeeded : False
マウント解除操作は正常に完了しました。
オフラインイメージに無人応答ファイルを適用しました。
```

6. **無人応答ファイルの言語設定を適用**

```powershell
# 言語設定用の無人応答ファイル
$languageUnattendXml = @"
<?xml version="1.0" encoding="utf-8"?>
<unattend xmlns="urn:schemas-microsoft-com:unattend">
    <settings pass="oobeSystem">
        <component name="Microsoft-Windows-International-Core" processorArchitecture="amd64" publicKeyToken="31bf3856ad364e35" language="neutral" versionScope="nonSxS">
            <InputLocale>ja-JP</InputLocale>
            <SystemLocale>ja-JP</SystemLocale>
            <UILanguage>ja-JP</UILanguage>
            <UserLocale>ja-JP</UserLocale>
        </component>
    </settings>
</unattend>
"@

# ファイルとして保存
$languageUnattendPath = "C:\temp\language_unattend.xml"
New-Item -Path $languageUnattendPath -ItemType File -Force
$languageUnattendXml | Out-File -FilePath $languageUnattendPath -Encoding utf8

# 無人応答ファイルを適用
Use-WindowsUnattend -Path $languageUnattendPath -Online -NoRestart
```

出力例：
```
Path          : C:\temp\language_unattend.xml
Online        : True
NoRestart     : True
LogPath       : C:\Windows\Logs\Dism\dism.log
ScratchDir    : 
LogLevel      : 
RestartNeeded : False
```

7. **無人応答ファイルの適用状態の確認**

```powershell
# ログファイルの確認
if (Test-Path "C:\Windows\Logs\Dism\dism.log") {
    $logContent = Get-Content -Path "C:\Windows\Logs\Dism\dism.log" -Tail 20
    Write-Host "DISMログの最終20行:"
    $logContent
} else {
    Write-Host "DISMログファイルが見つかりません。"
}

# 適用された設定を確認（例: タイムゾーン）
$currentTimeZone = Get-TimeZone
Write-Host "現在のタイムゾーン設定: $($currentTimeZone.Id)"
```

出力例：
```
DISMログの最終20行:
2025-03-30 02:10:15, Info                  DISM   DISM Provider Store: PID=1234 Loading provider from location C:\Windows\System32\Dism\UnattendProvider.dll
2025-03-30 02:10:15, Info                  DISM   DISM Provider Store: PID=1234 Successfully loaded the provider
2025-03-30 02:10:15, Info                  DISM   DISM Unattend Provider: PID=1234 Processing XML content
...
2025-03-30 02:10:16, Info                  DISM   DISM Provider Store: PID=1234 Finalizing session
2025-03-30 02:10:16, Info                  DISM   DISM.EXE: Image session has been closed. Reboot required=no.
2025-03-30 02:10:16, Info                  DISM   DISM.EXE: Process complete with exit code 0

現在のタイムゾーン設定: Tokyo Standard Time
```

## 対応PowerShellバージョン

Use-WindowsUnattendコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 3.0以降（Windows 8以降、Windows Server 2012以降）
- Windows 10またはWindows Server 2016以降では標準搭載
- PowerShell Coreでは利用不可（Windows PowerShellの機能）

このコマンドレットを使用するには、以下が必要：
- 管理者権限でのPowerShell実行
- DISMモジュールのインポート（Windows標準インストールに含まれる）

PowerShellでDISMモジュールを明示的にインポートするには：
```powershell
Import-Module DISM
```

## 参考サイト

- [Microsoft公式ドキュメント: Use-WindowsUnattend](https://docs.microsoft.com/ja-jp/powershell/module/dism/use-windowsunattend)
- [Microsoft Learn - Windows の無人インストール](https://docs.microsoft.com/ja-jp/windows-hardware/manufacture/desktop/windows-setup-automation-overview)
- [Microsoft Docs - Unattended Windows Setup Reference](https://docs.microsoft.com/ja-jp/windows-hardware/customize/desktop/unattend/)
- [Microsoft TechNet - Deployment Image Servicing and Management (DISM) Technical Reference](https://technet.microsoft.com/ja-jp/library/dd744566(v=ws.10).aspx)