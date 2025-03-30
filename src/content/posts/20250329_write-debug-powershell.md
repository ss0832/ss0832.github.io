---
title:  【PowerShell】Write-Debugコマンドでデバッグメッセージを出力する方法
published: 2025-03-29
description: "PowerShellのWrite-Debugコマンドレットの基本と活用テクニック"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-29

## 概要

`Write-Debug`はPowerShellでデバッグメッセージを出力するためのコマンドレットである。このコマンドレットは主にスクリプト開発時の問題解決やトラブルシューティングを支援するために使用される。デバッグメッセージは通常、スクリプトの実行中に変数の値や処理の進行状況を確認するために使用され、実行環境に応じて表示・非表示を切り替えることが可能である。

## 基本的な使い方

### 単純なデバッグメッセージの出力

最も基本的な使い方は、単にデバッグメッセージを出力することである。

```powershell
Write-Debug "これはデバッグメッセージです"
```

ただし、デフォルトでは`$DebugPreference`が`SilentlyContinue`に設定されているため、このメッセージは表示されない。

### デバッグメッセージを表示するための設定

デバッグメッセージを表示するには、`$DebugPreference`を`Continue`に変更するか、`-Debug`パラメータを使用する必要がある。

```powershell
$DebugPreference = "Continue"
Write-Debug "このデバッグメッセージは表示されます"
```

または

```powershell
Write-Debug "このデバッグメッセージは-Debugパラメータで表示されます" -Debug
```

### デバッグアクションの設定

デバッグメッセージの表示方法は、`-DebugAction`パラメータで制御することもできる。

```powershell
Write-Debug "カスタムアクションでのデバッグメッセージ" -Debug:$true -DebugAction Inquire
```

## 応用的な使い方

### 条件付きデバッグ出力

条件に基づいてデバッグメッセージを出力することが可能である。

```powershell
$value = 10
if ($value -gt 5) {
    Write-Debug "値が5より大きいです: $value"
}
```

### 関数内でのデバッグ情報

関数内で処理の流れや変数の状態を確認するためのデバッグ情報を出力できる。

```powershell
function Test-DebugFunction {
    [CmdletBinding()]
    param (
        [string]$Name
    )
    
    Write-Debug "関数の開始"
    Write-Debug "受け取った名前: $Name"
    
    $result = "Hello, $Name!"
    Write-Debug "結果: $result"
    
    Write-Debug "関数の終了"
    return $result
}

Test-DebugFunction -Name "PowerShell" -Debug
```

### スコープごとのデバッグ設定

スクリプトの特定の部分だけデバッグメッセージを表示することもできる。

```powershell
# グローバルスコープでは非表示
$DebugPreference = "SilentlyContinue"
Write-Debug "このメッセージは表示されません"

# 一時的にデバッグメッセージを有効にする
$originalDebugPreference = $DebugPreference
$DebugPreference = "Continue"
Write-Debug "このメッセージは表示されます"

# 元の設定に戻す
$DebugPreference = $originalDebugPreference
Write-Debug "このメッセージは再び表示されません"
```

### 高度なデバッグ情報の表示

複雑なオブジェクトの状態をデバッグ出力に含めることができる。

```powershell
$process = Get-Process -Id $PID
Write-Debug "現在のプロセス: $($process.Name), ID: $($process.Id), メモリ: $($process.WorkingSet64 / 1MB) MB"
```

## ハンズオン：Write-Debugの実践

以下はWrite-Debugコマンドの実践例とその出力である。

1. **基本的なデバッグメッセージ（表示されない）**

```powershell
# デフォルトの設定では表示されない
$DebugPreference
Write-Debug "このデバッグメッセージは表示されません"
```

出力例：
```
SilentlyContinue
```

2. **デバッグメッセージを表示する**

```powershell
# $DebugPreference を変更する
$DebugPreference = "Continue"
Write-Debug "このデバッグメッセージは表示されます"
```

出力例：
```
DEBUG: このデバッグメッセージは表示されます
```

3. **-Debugパラメータの使用**

```powershell
# $DebugPreference をデフォルトに戻す
$DebugPreference = "SilentlyContinue"
Write-Debug "このメッセージは-Debugフラグで表示されます" -Debug
```

出力例：
```
DEBUG: このメッセージは-Debugフラグで表示されます
```

4. **変数の値をデバッグ出力**

```powershell
$userName = "PowerShell User"
$loginCount = 42
Write-Debug "ユーザー: $userName, ログイン回数: $loginCount" -Debug
```

出力例：
```
DEBUG: ユーザー: PowerShell User, ログイン回数: 42
```

5. **Inquireアクションの使用**

```powershell
Write-Debug "処理を続行しますか？" -Debug:$true -DebugAction Inquire
```

出力例（インタラクティブ）：
```
DEBUG: 処理を続行しますか？
[Y] はい(Y)  [A] すべて続行(A)  [H] 保留(H)  [S] 中断(S)  [?] ヘルプ (既定値は "Y"): 
```

6. **関数内でのデバッグ情報**

```powershell
function Process-Item {
    [CmdletBinding()]
    param (
        [string]$Item
    )
    
    Write-Debug "アイテム処理開始: $Item"
    
    # 処理のシミュレーション
    $result = $Item.ToUpper()
    Write-Debug "処理結果: $result"
    
    return $result
}

$output = Process-Item -Item "test-data" -Debug
Write-Output "関数の戻り値: $output"
```

出力例：
```
DEBUG: アイテム処理開始: test-data
DEBUG: 処理結果: TEST-DATA
関数の戻り値: TEST-DATA
```

7. **配列処理でのデバッグ**

```powershell
$items = @("項目1", "項目2", "項目3")
foreach ($item in $items) {
    Write-Debug "処理中の項目: $item" -Debug
    # 何らかの処理
    Start-Sleep -Milliseconds 500
}
```

出力例：
```
DEBUG: 処理中の項目: 項目1
DEBUG: 処理中の項目: 項目2
DEBUG: 処理中の項目: 項目3
```

8. **オブジェクトの内容をデバッグ表示**

```powershell
$computerInfo = [PSCustomObject]@{
    Name = $env:COMPUTERNAME
    OS = (Get-CimInstance Win32_OperatingSystem).Caption
    Memory = [math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB, 2)
}

Write-Debug "コンピュータ情報: $($computerInfo | ConvertTo-Json -Compress)" -Debug
```

出力例：
```
DEBUG: コンピュータ情報: {"Name":"DESKTOP-ABC123","OS":"Microsoft Windows 10 Pro","Memory":16.0}
```

## 対応PowerShellバージョン

Write-Debugコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 2.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

## 参考サイト

- [Microsoft公式ドキュメント: Write-Debug](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.utility/write-debug)
- [PowerShell.org - Advanced Functions and Debugging](https://powershell.org/2018/08/advanced-functions-and-debugging/)
- [SS64.com PowerShell Commands - Write-Debug](https://ss64.com/ps/write-debug.html)
- [4sysops - PowerShell Write-Debug and debugging basics](https://4sysops.com/archives/powershell-write-debug-and-debugging-basics/)