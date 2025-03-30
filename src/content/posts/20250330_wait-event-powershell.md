---
title:  【PowerShell】Wait-Eventコマンドでイベント発生を待機する方法
published: 2025-03-30
description: "PowerShellのWait-Eventコマンドレットを使用したイベント処理の基本と応用"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Wait-Event`はPowerShellにおいて、特定のイベントが発生するまでスクリプトの実行を一時停止するためのコマンドレットである。非同期処理やイベント駆動型のスクリプトを作成する際に役立ち、PowerShellのイベントシステムと組み合わせることで、ファイル変更の検出、プロセスの終了、タイマーイベントなど様々な状況に対応できる。イベントが発生するまでブロックする同期的な待機メカニズムを提供することで、イベントベースのプログラミングを可能にしている。

## 基本的な使い方

### イベントの登録と待機

基本的な使い方は、まずイベントを登録し、そのイベントが発生するのを待機することである。

```powershell
# タイマーイベントを登録
Register-ObjectEvent -InputObject ([System.Timers.Timer]::new(3000)) -EventName Elapsed -SourceIdentifier "Timer"
# タイマーを開始
$Timer = Get-EventSubscriber -SourceIdentifier "Timer" | ForEach-Object {$_.Action = {Write-Host "タイマー発火"}; $_.MessageData -as [System.Timers.Timer];}
$Timer.Start()
# イベントを待機
Wait-Event -SourceIdentifier "Timer"
```

### 特定のイベント名での待機

イベント名を指定して待機することもできる。

```powershell
Register-EngineEvent -SourceIdentifier "CustomEvent" -Action { Write-Host "カスタムイベントが発生しました" }
# 別のセッションまたは並行処理で以下を実行：
# New-Event -SourceIdentifier "CustomEvent"
Wait-Event -SourceIdentifier "CustomEvent"
```

### タイムアウトの設定

指定した時間だけ待機し、その時間内にイベントが発生しない場合は処理を続行する。

```powershell
# 最大5秒間イベントを待機
Wait-Event -SourceIdentifier "TimedEvent" -Timeout 5
```

## 応用的な使い方

### ファイル変更イベントの待機

ファイルシステムの変更を監視し、変更があったら処理を行う例である。

```powershell
$watcher = New-Object System.IO.FileSystemWatcher
$watcher.Path = "C:\Temp"
$watcher.Filter = "*.txt"
$watcher.IncludeSubdirectories = $false
$watcher.EnableRaisingEvents = $true

Register-ObjectEvent -InputObject $watcher -EventName Changed -SourceIdentifier "FileChanged" -Action {
    Write-Host "ファイルが変更されました: $($Event.SourceEventArgs.Name)"
}

# ファイル変更イベントを待機
Wait-Event -SourceIdentifier "FileChanged"
```

### WMIイベントの待機

WMIイベントを監視し、特定のシステムイベントを待機する例である。

```powershell
Register-WmiEvent -Query "SELECT * FROM Win32_ProcessStartTrace" -SourceIdentifier "ProcessStarted" -Action {
    Write-Host "新しいプロセスが起動されました: $($Event.SourceEventArgs.NewEvent.ProcessName)"
}

# プロセス起動イベントを待機
Wait-Event -SourceIdentifier "ProcessStarted"
```

### 複数イベントの処理

複数のイベントを登録し、それらのいずれかを待機する例である。

```powershell
# 複数のイベントを登録
Register-ObjectEvent -InputObject ([System.Timers.Timer]::new(5000)) -EventName Elapsed -SourceIdentifier "Timer1"
Register-ObjectEvent -InputObject ([System.Timers.Timer]::new(8000)) -EventName Elapsed -SourceIdentifier "Timer2"

# タイマーを開始
$Timer1 = Get-EventSubscriber -SourceIdentifier "Timer1" | ForEach-Object {$_.Action = {Write-Host "タイマー1発火"}; $_.MessageData -as [System.Timers.Timer];}
$Timer2 = Get-EventSubscriber -SourceIdentifier "Timer2" | ForEach-Object {$_.Action = {Write-Host "タイマー2発火"}; $_.MessageData -as [System.Timers.Timer];}
$Timer1.Start()
$Timer2.Start()

# 最初に発生したイベントを待機
$event = Wait-Event -SourceIdentifier "Timer1", "Timer2"
Write-Host "イベントが発生しました: $($event.SourceIdentifier)"
```

## ハンズオン：Wait-Eventの実践

以下はWait-Eventコマンドの実践例とその出力である。

1. **タイマーイベントの待機**

```powershell
# タイマーオブジェクトを作成
$timer = New-Object System.Timers.Timer
$timer.Interval = 3000  # 3秒
$timer.AutoReset = $false  # 一度だけ発火

# イベントを登録
Register-ObjectEvent -InputObject $timer -EventName Elapsed -SourceIdentifier "SimpleTimer"

Write-Host "タイマーを開始します。3秒後にイベントが発生します - $(Get-Date)"
$timer.Start()

# イベントを待機
Write-Host "イベントの発生を待機中..."
$event = Wait-Event -SourceIdentifier "SimpleTimer"

Write-Host "イベントが発生しました! - $(Get-Date)"
Write-Host "イベント情報: $($event.SourceIdentifier), 時間: $($event.TimeGenerated)"

# クリーンアップ
Remove-Event -SourceIdentifier "SimpleTimer"
Unregister-Event -SourceIdentifier "SimpleTimer"
```

出力例：
```
タイマーを開始します。3秒後にイベントが発生します - 2025/03/30 02:10:00
イベントの発生を待機中...
イベントが発生しました! - 2025/03/30 02:10:03
イベント情報: SimpleTimer, 時間: 2025/03/30 02:10:03
```

2. **タイムアウト付きイベント待機**

```powershell
# ソースIDだけでイベントを登録（実際には発生しない）
Register-EngineEvent -SourceIdentifier "ManualEvent" -Action { Write-Host "このイベントは発生しません" }

Write-Host "イベントを5秒間だけ待機します - $(Get-Date)"
$event = Wait-Event -SourceIdentifier "ManualEvent" -Timeout 5

if ($event -eq $null) {
    Write-Host "タイムアウトしました。イベントは発生しませんでした - $(Get-Date)"
} else {
    Write-Host "イベントが発生しました: $($event.SourceIdentifier)"
}

# クリーンアップ
Unregister-Event -SourceIdentifier "ManualEvent"
```

出力例：
```
イベントを5秒間だけ待機します - 2025/03/30 02:15:00
タイムアウトしました。イベントは発生しませんでした - 2025/03/30 02:15:05
```

3. **カスタムイベントの生成と待機**

```powershell
# バックグラウンドジョブでイベントを発生させる
$job = Start-Job -ScriptBlock {
    Start-Sleep -Seconds 2
    New-Event -SourceIdentifier "CustomNotification" -MessageData "カスタムメッセージ"
}

Write-Host "カスタムイベントを待機します - $(Get-Date)"
$event = Wait-Event -SourceIdentifier "CustomNotification"

Write-Host "イベントが発生しました! - $(Get-Date)"
Write-Host "イベントのメッセージ: $($event.MessageData)"

# クリーンアップ
Remove-Event -SourceIdentifier "CustomNotification"
Remove-Job -Job $job -Force
```

出力例：
```
カスタムイベントを待機します - 2025/03/30 02:20:00
イベントが発生しました! - 2025/03/30 02:20:02
イベントのメッセージ: カスタムメッセージ
```

4. **ファイルシステム監視イベント**

```powershell
# 一時ディレクトリパスを作成
$tempDir = Join-Path $env:TEMP "EventTest"
if (-not (Test-Path $tempDir)) {
    New-Item -Path $tempDir -ItemType Directory | Out-Null
}
$filePath = Join-Path $tempDir "test.txt"

# ファイル監視を設定
$watcher = New-Object System.IO.FileSystemWatcher
$watcher.Path = $tempDir
$watcher.Filter = "test.txt"
$watcher.EnableRaisingEvents = $true

# 変更イベントを登録
Register-ObjectEvent -InputObject $watcher -EventName Changed -SourceIdentifier "FileChangeEvent"

# 別のプロセスでファイルを作成・変更するスクリプトを開始
$job = Start-Job -ScriptBlock {
    param($path)
    Start-Sleep -Seconds 3
    Set-Content -Path $path -Value "テストファイルの内容" -Force
} -ArgumentList $filePath

Write-Host "ファイル変更イベントを待機しています..."
$fileEvent = Wait-Event -SourceIdentifier "FileChangeEvent"

Write-Host "ファイルが変更されました！"
Write-Host "変更されたファイル: $($fileEvent.SourceEventArgs.FullPath)"
Write-Host "変更タイプ: $($fileEvent.SourceEventArgs.ChangeType)"

# クリーンアップ
Remove-Event -SourceIdentifier "FileChangeEvent"
Unregister-Event -SourceIdentifier "FileChangeEvent"
Remove-Job -Job $job -Force
Remove-Item -Path $tempDir -Recurse -Force
```

出力例：
```
ファイル変更イベントを待機しています...
ファイルが変更されました！
変更されたファイル: C:\Users\User\AppData\Local\Temp\EventTest\test.txt
変更タイプ: Changed
```

5. **複数のイベントソースの監視**

```powershell
# 2つのタイマーを設定
$timer1 = New-Object System.Timers.Timer
$timer1.Interval = 4000  # 4秒
$timer1.AutoReset = $false

$timer2 = New-Object System.Timers.Timer
$timer2.Interval = 2000  # 2秒
$timer2.AutoReset = $false

# イベントを登録
Register-ObjectEvent -InputObject $timer1 -EventName Elapsed -SourceIdentifier "Timer1"
Register-ObjectEvent -InputObject $timer2 -EventName Elapsed -SourceIdentifier "Timer2"

# タイマーを開始
Write-Host "両方のタイマーを開始します - $(Get-Date)"
$timer1.Start()
$timer2.Start()

# いずれかのイベントを待機
Write-Host "最初に発生するイベントを待機します..."
$firstEvent = Wait-Event -SourceIdentifier "Timer1", "Timer2"

Write-Host "最初に発生したイベント: $($firstEvent.SourceIdentifier) - $(Get-Date)"

# 残りのイベントを待機
Write-Host "残りのイベントを待機します..."
$remainingEvent = Wait-Event -SourceIdentifier @("Timer1", "Timer2") | Where-Object { $_.SourceIdentifier -ne $firstEvent.SourceIdentifier }

Write-Host "2番目のイベント: $($remainingEvent.SourceIdentifier) - $(Get-Date)"

# クリーンアップ
Get-EventSubscriber | Unregister-Event
Get-Event | Remove-Event
```

出力例：
```
両方のタイマーを開始します - 2025/03/30 02:30:00
最初に発生するイベントを待機します...
最初に発生したイベント: Timer2 - 2025/03/30 02:30:02
残りのイベントを待機します...
2番目のイベント: Timer1 - 2025/03/30 02:30:04
```

## 対応PowerShellバージョン

Wait-Eventコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 2.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

## 参考サイト

- [Microsoft公式ドキュメント: Wait-Event](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.utility/wait-event)
- [Microsoft Learn - PowerShellのイベントシステム](https://docs.microsoft.com/ja-jp/powershell/scripting/learn/deep-dives/event-driven-programming)
- [SS64.com PowerShell Commands - Wait-Event](https://ss64.com/ps/wait-event.html)
- [PowerShell.org - イベント処理のベストプラクティス](https://powershell.org/2018/10/powershell-event-processing/)