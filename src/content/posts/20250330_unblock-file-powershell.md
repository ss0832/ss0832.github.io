---
title:  【PowerShell】Unblock-Fileコマンドでセキュリティブロックを解除する方法
published: 2025-03-30
description: "PowerShellのUnblock-Fileコマンドレットを使ったダウンロードファイルのセキュリティブロック解除方法"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Unblock-File`はPowerShellにおいて、インターネットからダウンロードしたファイルに適用されるセキュリティブロックを解除するためのコマンドレットである。Windowsはセキュリティ上の理由から、インターネットからダウンロードしたファイルに対し「Zone.Identifier」と呼ばれる代替データストリーム（ADS）を追加し、そのファイルが信頼できない場所から取得されたことを示す。このブロックが適用されていると、ファイルを実行しようとした際に「このファイルは信頼できない可能性があります」といった警告が表示される。`Unblock-File`コマンドレットを使用すると、この代替データストリームを削除し、警告なしでファイルを使用できるようになる。

## 基本的な使い方

### 単一ファイルのブロック解除

特定のファイルのセキュリティブロックを解除するには、以下のコマンドを使用する。

```powershell
Unblock-File -Path "C:\Downloads\example.ps1"
```

### パスをパイプラインから渡す

ファイルパスをパイプライン経由で渡すこともできる。

```powershell
"C:\Downloads\example.ps1" | Unblock-File
```

### 複数ファイルの一括ブロック解除

複数のファイルを同時に処理するには、配列やワイルドカードを使用する。

```powershell
Unblock-File -Path "C:\Downloads\*.ps1"
```

### Get-ChildItemと組み合わせたブロック解除

Get-ChildItemコマンドレットと組み合わせて、特定の条件に合うファイルのブロックを解除できる。

```powershell
Get-ChildItem -Path "C:\Downloads" -Filter "*.zip" | Unblock-File
```

## 応用的な使い方

### フォルダ内のすべてのファイルを再帰的にブロック解除

フォルダとそのサブフォルダ内のすべてのファイルのブロックを解除するには、再帰的な検索を使用する。

```powershell
Get-ChildItem -Path "C:\Downloads" -Recurse | Unblock-File
```

### WhatIfパラメータを使った動作確認

実際にブロックを解除する前に、どのファイルが対象になるかを確認できる。

```powershell
Get-ChildItem -Path "C:\Downloads\*.ps1" | Unblock-File -WhatIf
```

### ブロックされたファイルのみを対象にする

Zone.Identifierを持つファイルのみを検出してブロック解除する。

```powershell
Get-ChildItem -Path "C:\Downloads" -Recurse | Where-Object { 
    (Get-Item $_.FullName -Stream Zone.Identifier -ErrorAction SilentlyContinue) 
} | Unblock-File
```

### 特定の拡張子を持つファイルだけをブロック解除

特定の拡張子（例：.ps1, .exe, .dll）を持つファイルのみをブロック解除する。

```powershell
$extensions = @(".ps1", ".exe", ".dll")
Get-ChildItem -Path "C:\Downloads" -Recurse | Where-Object { 
    $extensions -contains $_.Extension 
} | Unblock-File
```

## ハンズオン：Unblock-Fileの実践

以下はUnblock-Fileコマンドの実践例とその出力である。

1. **ファイルのブロック状態を確認する**

まず、ファイルがブロックされているかどうかを確認する。

```powershell
# テストファイルを作成（実際にはダウンロードされたファイルを使用）
$testFile = "C:\Temp\blockedScript.ps1"
Set-Content -Path $testFile -Value 'Write-Host "This is a test script."'

# 手動でZone.Identifierを追加してブロック状態をシミュレート
Set-Content -Path $testFile -Stream Zone.Identifier -Value '[ZoneTransfer]
ZoneId=3
ReferrerUrl=http://example.com
HostUrl=http://example.com/script.ps1'

# ブロック状態の確認
Get-Item -Path $testFile -Stream Zone.Identifier -ErrorAction SilentlyContinue
```

出力例：
```
FileName: C:\Temp\blockedScript.ps1

Stream          Length
------          ------
Zone.Identifier    82
```

2. **単一ファイルのブロック解除**

```powershell
# ファイルのブロックを解除
Unblock-File -Path $testFile

# ブロック解除後の状態を確認
Get-Item -Path $testFile -Stream Zone.Identifier -ErrorAction SilentlyContinue
```

出力例（ブロック解除後）：
```
# 出力なし（Zone.Identifierストリームが削除されたため）
```

3. **複数のファイルを一括でブロック解除**

```powershell
# テスト用に複数のブロックされたファイルを作成
$testDir = "C:\Temp\BlockedFiles"
New-Item -Path $testDir -ItemType Directory -Force | Out-Null

# 複数のテストファイルを作成
1..5 | ForEach-Object {
    $file = Join-Path $testDir "blocked$_.ps1"
    Set-Content -Path $file -Value "Write-Host 'Test script $_'"
    Set-Content -Path $file -Stream Zone.Identifier -Value '[ZoneTransfer]
ZoneId=3'
}

# ディレクトリ内のファイル数を確認
$blockedFiles = Get-ChildItem -Path $testDir
Write-Host "ブロックされたファイル数: $($blockedFiles.Count)"

# すべてのファイルのブロックを解除
Get-ChildItem -Path $testDir | Unblock-File

# ブロック解除後に再確認
$stillBlocked = Get-ChildItem -Path $testDir | Where-Object {
    Get-Item $_.FullName -Stream Zone.Identifier -ErrorAction SilentlyContinue
}

Write-Host "ブロック解除後、まだブロックされているファイル数: $($stillBlocked.Count)"
```

出力例：
```
ブロックされたファイル数: 5
ブロック解除後、まだブロックされているファイル数: 0
```

4. **WhatIfパラメータの使用**

```powershell
# 新しいテストファイルを作成
$newTestDir = "C:\Temp\MoreBlockedFiles"
New-Item -Path $newTestDir -ItemType Directory -Force | Out-Null

# 複数の種類のファイルを作成
$extensions = @(".ps1", ".txt", ".exe", ".bat")
foreach ($ext in $extensions) {
    $file = Join-Path $newTestDir "test$ext"
    Set-Content -Path $file -Value "Test content"
    Set-Content -Path $file -Stream Zone.Identifier -Value '[ZoneTransfer]
ZoneId=3'
}

# WhatIf モードでブロック解除をシミュレート
Write-Host "WhatIfモードでのブロック解除シミュレーション："
Get-ChildItem -Path $newTestDir -Filter "*.ps1" | Unblock-File -WhatIf
```

出力例：
```
WhatIfモードでのブロック解除シミュレーション：
What if: C:\Temp\MoreBlockedFiles\test.ps1 のブロックを解除します
```

5. **特定の条件に合うファイルのみをブロック解除**

```powershell
# 条件に基づいてファイルをブロック解除
Write-Host "実行ファイルのみをブロック解除："
$executableExtensions = @(".ps1", ".exe", ".bat")
Get-ChildItem -Path $newTestDir | Where-Object { 
    $executableExtensions -contains $_.Extension 
} | ForEach-Object {
    Write-Host "ブロック解除: $($_.Name)"
    Unblock-File -Path $_.FullName
}

# 結果を確認
$remainingBlocked = Get-ChildItem -Path $newTestDir | Where-Object {
    Get-Item $_.FullName -Stream Zone.Identifier -ErrorAction SilentlyContinue
}

Write-Host "`n残りのブロックされているファイル:"
$remainingBlocked | Select-Object Name
```

出力例：
```
実行ファイルのみをブロック解除：
ブロック解除: test.ps1
ブロック解除: test.exe
ブロック解除: test.bat

残りのブロックされているファイル:
Name
----
test.txt
```

6. **クリーンアップ**

```powershell
# テスト用ディレクトリとファイルを削除
Remove-Item -Path $testFile -Force -ErrorAction SilentlyContinue
Remove-Item -Path $testDir -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item -Path $newTestDir -Recurse -Force -ErrorAction SilentlyContinue
Write-Host "テストファイルとディレクトリをクリーンアップしました。"
```

出力例：
```
テストファイルとディレクトリをクリーンアップしました。
```

## 対応PowerShellバージョン

Unblock-Fileコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 3.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

ただし、Windows以外のオペレーティングシステム（Linux、macOS）ではZone.Identifierの概念がないため、これらのシステムではUnblock-Fileは実質的に何も行わない。

## 参考サイト

- [Microsoft公式ドキュメント: Unblock-File](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.utility/unblock-file)
- [Microsoft Learn - PowerShellセキュリティ機能](https://docs.microsoft.com/ja-jp/powershell/scripting/learn/ps101/09-functions)
- [Alternate Data Streams in NTFS](https://blogs.technet.microsoft.com/askcore/2013/03/24/alternate-data-streams-in-ntfs/)
- [PowerShellマガジン - セキュリティブロックの仕組みとUnblock-File](https://powershellmagazine.com/2013/07/19/pstip-how-to-check-if-a-file-has-been-blocked/)
- [SS64.com PowerShell Commands - Unblock-File](https://ss64.com/ps/unblock-file.html)