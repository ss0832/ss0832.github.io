---
title:  【PowerShell】Write-EventLogコマンドでWindowsイベントログを操作する方法
published: 2025-03-29
description: "PowerShellのWrite-EventLogコマンドを使ったWindowsイベントログへの書き込み方法を解説"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-29

## 概要

`Write-EventLog`はPowerShellでWindowsイベントログにエントリを書き込むためのコマンドレットである。アプリケーションの動作記録、エラー追跡、セキュリティ監視など、システム管理やトラブルシューティングに役立つログを残すことができる。このコマンドレットを使用することで、スクリプトやアプリケーションの実行状態をWindowsの標準的なログ機能を通じて管理することが可能になる。

## 基本的な使い方

### イベントログへの基本的な書き込み

既存のログソースを使用してイベントログに情報を書き込む最も基本的な例である。

```powershell
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 100 -Message "テストメッセージ"
```

### イベントの種類を指定

イベントの種類（情報、警告、エラーなど）を指定することができる。

```powershell
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 101 -EntryType Warning -Message "警告：ディスク容量が少なくなっています"
```

### カテゴリの指定

カテゴリ番号を指定してイベントを分類できる。

```powershell
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 102 -Category 0 -Message "カテゴリ付きメッセージ"
```

### 新しいイベントソースの登録

新しいイベントソースを使用する場合は、事前に登録が必要である（管理者権限が必要）。

```powershell
if (-not [System.Diagnostics.EventLog]::SourceExists("MyScript")) {
    New-EventLog -LogName "Application" -Source "MyScript"
}
Write-EventLog -LogName "Application" -Source "MyScript" -EventId 1000 -Message "新しいソースからのメッセージ"
```

## 応用的な使い方

### Try-Catchブロックでのエラーログ

エラーハンドリングと組み合わせて使用する例である。

```powershell
try {
    # 何らかの処理
    $result = 10 / 0  # ゼロ除算エラー
} catch {
    $errorMessage = $_.Exception.Message
    Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 900 -EntryType Error -Message "エラーが発生しました: $errorMessage"
}
```

### リモートコンピュータのイベントログへの書き込み

リモートコンピュータのイベントログに書き込む場合は`-ComputerName`パラメータを使用する。

```powershell
Write-EventLog -ComputerName "Server01" -LogName "Application" -Source "PowerShell" -EventId 200 -Message "リモートサーバーへのログメッセージ"
```

### バイナリデータの追加

バイナリデータをイベントログに含めることもできる。

```powershell
$binaryData = [System.BitConverter]::GetBytes(1234)
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 300 -Message "バイナリデータ付きメッセージ" -RawData $binaryData
```

### スクリプト実行のログ記録

スクリプトの開始と終了をログに記録する例である。

```powershell
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 5000 -EntryType Information -Message "スクリプト実行開始"
# スクリプトの処理
Start-Sleep -Seconds 2  # 処理の代わりに2秒待機
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 5001 -EntryType Information -Message "スクリプト実行完了"
```

## ハンズオン：Write-EventLogの実践

以下は一連のWrite-EventLogコマンドとその確認方法である。

1. **基本的なイベントログの書き込みとその確認**

```powershell
# イベントログの書き込み
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 1 -Message "基本的なログメッセージ"

# 書き込んだログの確認
Get-EventLog -LogName "Application" -Newest 1 -Source "PowerShell"
```

出力例：
```
Index Time          EntryType   Source   InstanceID Message
----- ----          ---------   ------   ---------- -------
12345 Mar 29 01:55  Information PowerShell        1 基本的なログメッセージ
```

2. **異なる種類のイベントログの書き込み**

```powershell
# エラーイベントの書き込み
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 2 -EntryType Error -Message "これはエラーメッセージです"

# 警告イベントの書き込み
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 3 -EntryType Warning -Message "これは警告メッセージです"

# 直近2件のログを表示して確認
Get-EventLog -LogName "Application" -Newest 2 -Source "PowerShell"
```

出力例：
```
Index Time          EntryType   Source   InstanceID Message
----- ----          ---------   ------   ---------- -------
12347 Mar 29 01:56  Warning     PowerShell        3 これは警告メッセージです
12346 Mar 29 01:56  Error       PowerShell        2 これはエラーメッセージです
```

3. **新しいイベントソースの作成と使用**

```powershell
# 管理者権限で実行する必要があります
if (-not [System.Diagnostics.EventLog]::SourceExists("MyCustomApp")) {
    New-EventLog -LogName "Application" -Source "MyCustomApp"
    Write-Output "新しいイベントソース 'MyCustomApp' を作成しました"
}

# 新しいソースを使ってログを書き込む
Write-EventLog -LogName "Application" -Source "MyCustomApp" -EventId 100 -Message "カスタムアプリからのメッセージ"

# 新しいソースからのログを確認
Get-EventLog -LogName "Application" -Newest 1 -Source "MyCustomApp"
```

出力例：
```
新しいイベントソース 'MyCustomApp' を作成しました

Index Time          EntryType   Source       InstanceID Message
----- ----          ---------   ------       ---------- -------
12348 Mar 29 01:57  Information MyCustomApp        100 カスタムアプリからのメッセージ
```

4. **Try-Catchブロックでのエラーログ**

```powershell
try {
    # 意図的なエラーを発生させる
    $null.Property  # NullReferenceException
} catch {
    $errorMessage = $_.Exception.Message
    Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 500 -EntryType Error -Message "エラーが発生しました: $errorMessage"
    Write-Output "エラーをイベントログに記録しました"
}

# エラーログを確認
Get-EventLog -LogName "Application" -Newest 1 -Source "PowerShell"
```

出力例：
```
エラーをイベントログに記録しました

Index Time          EntryType   Source     InstanceID Message
----- ----          ---------   ------     ---------- -------
12349 Mar 29 01:58  Error       PowerShell        500 エラーが発生しました: null 参照の例外がスローされました。
```

5. **カテゴリ付きイベントの書き込み**

```powershell
# カテゴリ1のイベント
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 600 -Category 1 -Message "ネットワーク関連イベント"

# カテゴリ2のイベント
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 601 -Category 2 -Message "セキュリティ関連イベント"

# カテゴリで絞り込んでログを確認
Get-EventLog -LogName "Application" -Newest 10 -Source "PowerShell" | Where-Object { $_.CategoryNumber -eq 1 }
```

出力例：
```
Index Time          EntryType   Source     InstanceID Message
----- ----          ---------   ------     ---------- -------
12350 Mar 29 01:59  Information PowerShell        600 ネットワーク関連イベント
```

6. **バイナリデータ付きのイベント**

```powershell
# バイナリデータを準備
$binaryData = [System.Text.Encoding]::ASCII.GetBytes("Sample data")

# バイナリデータ付きでイベントを書き込み
Write-EventLog -LogName "Application" -Source "PowerShell" -EventId 700 -Message "バイナリデータ付きイベント" -RawData $binaryData

# 書き込んだログを確認
$event = Get-EventLog -LogName "Application" -Newest 1 -Source "PowerShell"
Write-Output "イベントID: $($event.EventID)"
Write-Output "メッセージ: $($event.Message)"
```

出力例：
```
イベントID: 700
メッセージ: バイナリデータ付きイベント
```

## 対応PowerShellバージョン

Write-EventLogコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 1.0以降
- Windows PowerShellに限定されたコマンドレット（PowerShell Core 6.0以降ではサポートされていない）

**注意**: PowerShell Core 6.0およびPowerShell 7.0以降では、このコマンドレットはWindows互換性モジュールの一部として利用可能だが、全ての環境で動作するわけではない。

## 参考サイト

- [Microsoft公式ドキュメント: Write-EventLog](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.management/write-eventlog)
- [Windows PowerShell Cookbook - Event Logging](https://docs.microsoft.com/en-us/previous-versions/windows/it-pro/windows-powershell-1.0/ff730949(v=technet.10))
- [Learn Windows PowerShell in a Month of Lunches - Working with EventLogs](https://www.manning.com/books/learn-windows-powershell-in-a-month-of-lunches-third-edition)
- [SS64.com PowerShell Commands - Write-EventLog](https://ss64.com/ps/write-eventlog.html)