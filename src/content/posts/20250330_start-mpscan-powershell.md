---
title:  【PowerShell】Start-MpScanコマンドでWindows Defenderスキャンを実行する方法
published: 2025-03-30
description: "PowerShellのStart-MpScanコマンドレットを使ったWindows Defenderウイルススキャンの実行方法"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Start-MpScan`はPowerShellのWindowsDefender モジュールに含まれるコマンドレットであり、Windows Defenderアンチマルウェアスキャンをコマンドラインから実行するために使用される。このコマンドレットを使用することで、クイックスキャン、フルスキャン、カスタムスキャンなど様々な種類のスキャンを実行でき、スクリプトや自動化タスクにセキュリティスキャンを組み込むことができる。また、特定のドライブやフォルダに対象を絞ったスキャンや、プロセスメモリのスキャンなどの高度なオプションも提供されている。システム管理者やセキュリティ担当者にとって、定期的なセキュリティチェックやインシデント対応時の調査に役立つツールである。

## 基本的な使い方

### クイックスキャンの実行

最も基本的な使い方は、システムのクイックスキャンを実行することである。

```powershell
Start-MpScan -ScanType QuickScan
```

### フルスキャンの実行

システム全体の完全なスキャンを実行するには、フルスキャンを指定する。

```powershell
Start-MpScan -ScanType FullScan
```

### カスタムスキャンの実行

特定のフォルダやファイルに対象を絞ったスキャンを実行できる。

```powershell
Start-MpScan -ScanType CustomScan -ScanPath "C:\Users\Username\Downloads"
```

### ブートセクタスキャンの実行

ブートセクターのマルウェアをスキャンするには、以下のコマンドを使用する。

```powershell
Start-MpScan -ScanType BootSectorScan
```

## 応用的な使い方

### 複数のパスを指定したカスタムスキャン

複数のパスを指定してカスタムスキャンを実行することが可能である。

```powershell
Start-MpScan -ScanType CustomScan -ScanPath "C:\Users\Username\Downloads", "D:\Data"
```

### 低優先度でのスキャン実行

バックグラウンドでスキャンを実行し、システムパフォーマンスへの影響を最小限に抑える。

```powershell
Start-MpScan -ScanType FullScan -AsJob -ScanParameters LowCPU
```

### ネットワークドライブのスキャン

ネットワークドライブをスキャン対象に含めることも可能である。

```powershell
Start-MpScan -ScanType CustomScan -ScanPath "\\Server\Share" -ScanParameters DisableRemediation
```

### スケジュールタスクとしてのスキャン設定

特定の時間にスキャンを実行するためのスケジュールタスクを設定することができる。

```powershell
$action = New-ScheduledTaskAction -Execute 'powershell.exe' -Argument '-NoProfile -WindowStyle Hidden -Command "Start-MpScan -ScanType QuickScan"'
$trigger = New-ScheduledTaskTrigger -Daily -At 2am
Register-ScheduledTask -Action $action -Trigger $trigger -TaskName "Daily Defender Scan" -Description "Daily Windows Defender Quick Scan"
```

## ハンズオン：Start-MpScanの実践

以下はStart-MpScanコマンドの実践例とその出力である。

1. **クイックスキャンの実行**

```powershell
# クイックスキャンの実行
Write-Host "Windows Defender クイックスキャンを開始します..."
Start-MpScan -ScanType QuickScan
```

出力例：
```
Windows Defender クイックスキャンを開始します...

Scanning: C:\Windows\System32\drivers\etc\hosts
Scanning: C:\ProgramData\Microsoft\Windows Defender\Scans\mpenginedb.db
...
Scan ID: {0A1B2C3D-4E5F-6789-0A1B-2C3D4E5F6789}
Scan Type: Quick Scan
Scan Parameters: Default
Scan Resources:
        Standard
Scan Status: Complete
Number of Scanned Objects: 12345
Time Elapsed: 00:01:23

Scan finished successfully.
```

2. **特定フォルダのカスタムスキャン**

```powershell
# ダウンロードフォルダのカスタムスキャン
$downloadFolder = "$env:USERPROFILE\Downloads"
Write-Host "ダウンロードフォルダ: $downloadFolder のスキャンを開始します..."
Start-MpScan -ScanType CustomScan -ScanPath $downloadFolder
```

出力例：
```
ダウンロードフォルダ: C:\Users\Username\Downloads のスキャンを開始します...

Scanning: C:\Users\Username\Downloads\report.pdf
Scanning: C:\Users\Username\Downloads\setup.exe
...
Scan ID: {9A8B7C6D-5E4F-3210-9A8B-7C6D5E4F3210}
Scan Type: Custom Scan
Scan Parameters: Default
Scan Resources:
        C:\Users\Username\Downloads
Scan Status: Complete
Number of Scanned Objects: 42
Time Elapsed: 00:00:17

Scan finished successfully.
```

3. **バックグラウンドでのフルスキャン**

```powershell
# バックグラウンドでフルスキャンを実行
$job = Start-MpScan -ScanType FullScan -AsJob
Write-Host "フルスキャンをバックグラウンドで開始しました。ジョブID: $($job.Id)"

# ジョブの状態を確認
Start-Sleep -Seconds 5
$jobStatus = Get-Job -Id $job.Id
Write-Host "現在のスキャン状態: $($jobStatus.State)"
```

出力例：
```
フルスキャンをバックグラウンドで開始しました。ジョブID: 3

現在のスキャン状態: Running
```

4. **スキャン結果の確認とジョブの監視**

```powershell
# 前のステップで開始したジョブの継続監視（実際のスクリプトではもう少し洗練された方法を使用する）
for ($i = 0; $i -lt 5; $i++) {
    $jobStatus = Get-Job -Id $job.Id
    Write-Host "ジョブ状態 ($i): $($jobStatus.State)"
    if ($jobStatus.State -eq "Completed") {
        break
    }
    Start-Sleep -Seconds 10
}

# ジョブが終了したかを確認
if ((Get-Job -Id $job.Id).State -eq "Completed") {
    Write-Host "スキャンが完了しました。結果を取得します..."
    Receive-Job -Id $job.Id
    Remove-Job -Id $job.Id
} else {
    Write-Host "スキャンはまだ実行中です。バックグラウンドで継続します。"
    Write-Host "後で 'Receive-Job -Id $($job.Id)' を実行して結果を確認してください。"
}
```

出力例：
```
ジョブ状態 (0): Running
ジョブ状態 (1): Running
ジョブ状態 (2): Running
ジョブ状態 (3): Running
ジョブ状態 (4): Running
スキャンはまだ実行中です。バックグラウンドで継続します。
後で 'Receive-Job -Id 3' を実行して結果を確認してください。
```

5. **スキャンの詳細設定とリアルタイム保護の一時停止/再開**

```powershell
# リアルタイム保護の状態を確認
$rtProtection = Get-MpPreference | Select-Object DisableRealtimeMonitoring
Write-Host "リアルタイム保護の無効化状態: $($rtProtection.DisableRealtimeMonitoring)"

# スキャン前にリアルタイム保護を一時的に無効化（注意: セキュリティリスクがあります）
Write-Host "リアルタイム保護を一時的に無効化します..."
Set-MpPreference -DisableRealtimeMonitoring $true

# 高度な設定でスキャンを実行
Write-Host "アーカイブスキャンを有効にして、カスタムスキャンを開始します..."
Start-MpScan -ScanType CustomScan -ScanPath "C:\Test" -ScanParameters DisableRemovableDriveScanning,IgnoreSystemBackupFiles

# スキャン後にリアルタイム保護を再度有効化
Write-Host "リアルタイム保護を再度有効化します..."
Set-MpPreference -DisableRealtimeMonitoring $false

# 最終状態を確認
$rtProtection = Get-MpPreference | Select-Object DisableRealtimeMonitoring
Write-Host "リアルタイム保護の無効化状態: $($rtProtection.DisableRealtimeMonitoring)"
```

出力例：
```
リアルタイム保護の無効化状態: False
リアルタイム保護を一時的に無効化します...
アーカイブスキャンを有効にして、カスタムスキャンを開始します...

Scanning: C:\Test\document1.docx
Scanning: C:\Test\archive.zip
...
Scan ID: {5E4F3210-9A8B-7C6D-5E4F-3210-9A8B}
Scan Type: Custom Scan
Scan Parameters: DisableRemovableDriveScanning,IgnoreSystemBackupFiles
Scan Resources:
        C:\Test
Scan Status: Complete
Number of Scanned Objects: 15
Time Elapsed: 00:00:08

Scan finished successfully.

リアルタイム保護を再度有効化します...
リアルタイム保護の無効化状態: False
```

6. **スキャン履歴の確認**

```powershell
# 最近のスキャン履歴を確認
Write-Host "Windows Defenderスキャン履歴の取得を試みます..."
try {
    $history = Get-MpThreatDetection | Select-Object -First 5
    if ($history) {
        Write-Host "直近の5件のスキャン/検出履歴:"
        $history | Format-Table ThreatID, ThreatName, Resources, InitialDetectionTime -AutoSize
    } else {
        Write-Host "直近の検出履歴はありません。"
    }
    
    # 代わりにWindows Defenderログも確認
    Write-Host "Windows Defenderイベントログを確認中..."
    $defenderLogs = Get-WinEvent -LogName "Microsoft-Windows-Windows Defender/Operational" -MaxEvents 5 -ErrorAction SilentlyContinue
    if ($defenderLogs) {
        $defenderLogs | Format-Table TimeCreated, Id, Message -AutoSize -Wrap
    } else {
        Write-Host "Windows Defenderのイベントログにアクセスできませんでした。"
    }
} catch {
    Write-Host "履歴の取得中にエラーが発生しました: $_"
}
```

出力例：
```
Windows Defenderスキャン履歴の取得を試みます...
直近の検出履歴はありません。
Windows Defenderイベントログを確認中...

TimeCreated          Id Message                                                
-----------          -- -------                                                
3/30/2025 2:45:12 AM 1001 Windows Defenderスキャンが完了しました。スキャンタイプ: クイックスキャン。脅威は検出されませんでした。
3/30/2025 2:40:25 AM 1000 Windows Defenderスキャンが開始されました。スキャンタイプ: クイックスキャン
3/29/2025 10:15:33 PM 5001 リアルタイム保護が無効に設定されています。
3/29/2025 10:15:05 PM 5000 Windows Defender機能が変更されました。詳細: リアルタイム保護 - 有効
3/29/2025 9:00:00 PM 1000 Windows Defenderスキャンが開始されました。スキャンタイプ: スケジュールされたスキャン
```

7. **脅威に対するアクションの実行**

```powershell
# 検出された脅威の一覧を表示
$threats = Get-MpThreatDetection
if ($threats) {
    Write-Host "検出された脅威の一覧:"
    $threats | Format-Table ThreatID, ThreatName, Resources -AutoSize
    
    # 特定の脅威に対してアクションを実行
    $threatId = $threats[0].ThreatID
    Write-Host "脅威ID: $threatId に対して削除アクションを実行します..."
    Remove-MpThreat -ThreatID $threatId
} else {
    Write-Host "検出された脅威はありません。テスト用のICSAテストファイルを使用してシミュレーションできます。"
    Write-Host "注意: テスト用ファイルは本当のマルウェアではありませんが、アンチウイルスソフトによって検出されます。"
    Write-Host "テストファイルの詳細はMicrosoft公式サイトで確認できます。"
}
```

出力例（脅威が検出されていない場合）：
```
検出された脅威はありません。テスト用のICSAテストファイルを使用してシミュレーションできます。
注意: テスト用ファイルは本当のマルウェアではありませんが、アンチウイルスソフトによって検出されます。
テストファイルの詳細はMicrosoft公式サイトで確認できます。
```

出力例（脅威が検出された場合）：
```
検出された脅威の一覧:

ThreatID                               ThreatName                 Resources
--------                               ----------                 ---------
2147519003                             Trojan:Win32/Test.EICAR    C:\Test\eicar.com
2147519004                             PUA:Win32/OperativeExample D:\Downloads\example.exe

脅威ID: 2147519003 に対して削除アクションを実行します...
```

## 対応PowerShellバージョン

Start-MpScanコマンドレットは以下のPowerShellバージョンと環境で利用可能である：
- Windows PowerShell 4.0以降（Windows 8.1/Server 2012 R2以降）
- Windows PowerShell 5.0および5.1（Windows 10/Server 2016以降で標準搭載）
- PowerShell Core 6.0以降は非対応（WindowsDefenderモジュールはWindows PowerShellのみ）
- PowerShell 7.0以降は、Windows OSで実行している場合のみWindowsCompatibilityモジュールを通じて利用可能

このコマンドレットを使用するには、以下が必要：
- Windows 8.1/Windows Server 2012 R2以降のWindows OS
- Windows Defenderがインストールされ有効になっていること
- 管理者権限での実行（一部の機能に必要）
- WindowsDefenderモジュールのインポート（自動的にインポートされる場合が多い）

PowerShellでWindowsDefenderモジュールを明示的にインポートするには：
```powershell
Import-Module Defender
```

## 参考サイト

- [Microsoft公式ドキュメント: Start-MpScan](https://docs.microsoft.com/ja-jp/powershell/module/defender/start-mpscan)
- [Microsoft Learn - Windows Defenderを使用したセキュリティ管理](https://docs.microsoft.com/ja-jp/learn/modules/implement-windows-defender-advanced-threat-protection/)
- [Microsoft Security Blog - Using PowerShell to manage Windows Defender](https://www.microsoft.com/security/blog/2018/05/14/powershell-for-windows-defender-atp/)
- [Windows Defenderセキュリティインテリジェンス更新プログラム](https://www.microsoft.com/en-us/wdsi/)
- [Microsoft Docs - Windows Defenderコマンドライン管理](https://docs.microsoft.com/ja-jp/windows/security/threat-protection/microsoft-defender-antivirus/command-line-arguments-microsoft-defender-antivirus)