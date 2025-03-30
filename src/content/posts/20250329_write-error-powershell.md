---
title:  【PowerShell】Write-Errorコマンドによるエラー処理の実装
published: 2025-03-29
description: "PowerShellのWrite-Errorコマンドの基本的な使い方と応用例を解説"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-29

## 概要

`Write-Error`はPowerShellでエラーメッセージを生成するためのコマンドレットである。スクリプトやコマンドの実行中に問題が発生した場合に、エラーストリームにメッセージを出力して利用者に通知することができる。通常のコンソール出力と異なり、エラーとして明示的に扱われるため、エラーハンドリングやログ記録などの自動化処理に適している。

## 基本的な使い方

### 単純なエラーメッセージの出力

最も基本的な使い方は、エラーメッセージを指定するだけである。

```powershell
Write-Error "設定ファイルが見つかりません。"
```

### エラーカテゴリの指定

エラーの種類をカテゴリとして指定することができる。

```powershell
Write-Error -Message "ファイルへのアクセスが拒否されました。" -Category PermissionDenied
```

### エラーIDの指定

エラーを識別するためのIDを指定することができる。

```powershell
Write-Error -Message "データベース接続エラー" -ErrorId "DB001"
```

### 例外オブジェクトの指定

特定の例外に関連付けたエラーを出力することができる。

```powershell
$exception = [System.InvalidOperationException]::new("操作が無効です。")
Write-Error -Exception $exception
```

## 応用的な使い方

### ターゲットオブジェクトの指定

エラーの原因となったオブジェクトを指定することができる。

```powershell
$problematicObject = Get-Item "C:\temp"
Write-Error -Message "アクセス権限がありません" -TargetObject $problematicObject
```

### カスタムエラーレコードの作成

詳細なエラー情報を持つエラーレコードを作成できる。

```powershell
$errorRecord = [System.Management.Automation.ErrorRecord]::new(
    [System.Exception]::new("カスタムエラー"),
    "CustomError001",
    [System.Management.Automation.ErrorCategory]::InvalidOperation,
    $null
)
Write-Error -ErrorRecord $errorRecord
```

### Try-Catchブロック内でのエラー処理

try-catchブロック内でエラーを適切に処理する例である。

```powershell
try {
    # 何らかの処理
    $result = 10 / 0  # ゼロ除算エラー
} catch {
    Write-Error -Message "計算中にエラーが発生しました: $_" -Category InvalidOperation
}
```

### カテゴリ別のエラー処理

異なるエラーカテゴリを使い分けることで、適切なエラーハンドリングが可能になる。

```powershell
function Test-ErrorCategories {
    param (
        [string]$FilePath,
        [int]$Value
    )
    
    if (-not (Test-Path $FilePath)) {
        Write-Error -Message "ファイル '$FilePath' が存在しません。" -Category ObjectNotFound
        return
    }
    
    if ($Value -lt 0) {
        Write-Error -Message "値 '$Value' は正の数でなければなりません。" -Category InvalidArgument
        return
    }
    
    # 正常処理
}
```

## ハンズオン：Write-Errorの実践

以下は一連のWrite-Errorコマンドとその出力例である。

1. **基本的なエラーメッセージ**

```powershell
Write-Error "これは基本的なエラーメッセージです。"
```

出力例：
```
Write-Error: これは基本的なエラーメッセージです。
    + CategoryInfo          : NotSpecified: (:) [Write-Error], WriteErrorException
    + FullyQualifiedErrorId : Microsoft.PowerShell.Commands.WriteErrorException
```

2. **エラーカテゴリの指定**

```powershell
Write-Error -Message "ファイルが見つかりません" -Category ObjectNotFound
```

出力例：
```
Write-Error: ファイルが見つかりません
    + CategoryInfo          : ObjectNotFound: (:) [Write-Error], WriteErrorException
    + FullyQualifiedErrorId : Microsoft.PowerShell.Commands.WriteErrorException
```

3. **エラーIDの指定**

```powershell
Write-Error -Message "無効な操作です" -ErrorId "InvalidOp001"
```

出力例：
```
Write-Error: 無効な操作です
    + CategoryInfo          : NotSpecified: (:) [Write-Error], WriteErrorException
    + FullyQualifiedErrorId : InvalidOp001,Microsoft.PowerShell.Commands.WriteErrorException
```

4. **例外オブジェクトの使用**

```powershell
$exception = [System.ArgumentNullException]::new("param1")
Write-Error -Message "null引数が検出されました" -Exception $exception
```

出力例：
```
Write-Error: null引数が検出されました
    + CategoryInfo          : NotSpecified: (:) [Write-Error], ArgumentNullException
    + FullyQualifiedErrorId : System.ArgumentNullException,Microsoft.PowerShell.Commands.WriteErrorException
```

5. **ターゲットオブジェクトの指定**

```powershell
$targetObj = "問題のあるデータ"
Write-Error -Message "データの検証に失敗しました" -TargetObject $targetObj
```

出力例：
```
Write-Error: データの検証に失敗しました
    + CategoryInfo          : NotSpecified: (問題のあるデータ:String) [Write-Error], WriteErrorException
    + FullyQualifiedErrorId : Microsoft.PowerShell.Commands.WriteErrorException
```

6. **カスタムエラーレコードの使用**

```powershell
$errorRecord = [System.Management.Automation.ErrorRecord]::new(
    [System.Exception]::new("ネットワーク接続が切断されました"),
    "NetworkError001",
    [System.Management.Automation.ErrorCategory]::ConnectionError,
    $null
)
Write-Error -ErrorRecord $errorRecord
```

出力例：
```
Write-Error: ネットワーク接続が切断されました
    + CategoryInfo          : ConnectionError: (:) [], NetworkError001
    + FullyQualifiedErrorId : NetworkError001
```

7. **エラー変数への格納**

```powershell
Write-Error "エラーを変数に格納" -ErrorVariable myError
$myError | Format-List * -Force
```

出力例（最初のエラーメッセージの後）：
```
PSMessageDetails      : 
Exception             : Microsoft.PowerShell.Commands.WriteErrorException: エラーを変数に格納
TargetObject          : 
CategoryInfo          : NotSpecified: (:) [Write-Error], WriteErrorException
FullyQualifiedErrorId : Microsoft.PowerShell.Commands.WriteErrorException
ErrorDetails          : 
InvocationInfo        : System.Management.Automation.InvocationInfo
ScriptStackTrace      : at <ScriptBlock>, <No file>: line 1
PipelineIterationInfo : {}
```

8. **Try-Catchブロックでのエラーハンドリング**

```powershell
try {
    # エラーを発生させる
    throw "強制的なエラー"
} catch {
    Write-Error -Message "処理中にエラーが発生しました: $_" -Category OperationStopped -ErrorAction Continue
    Write-Host "エラー後の処理を継続します" -ForegroundColor Yellow
}
```

出力例：
```
Write-Error: 処理中にエラーが発生しました: 強制的なエラー
    + CategoryInfo          : OperationStopped: (:) [Write-Error], WriteErrorException
    + FullyQualifiedErrorId : Microsoft.PowerShell.Commands.WriteErrorException

エラー後の処理を継続します
```

## 対応PowerShellバージョン

Write-Errorコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 1.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

## 参考サイト

- [Microsoft公式ドキュメント: Write-Error](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.utility/write-error)
- [PowerShell.org - Understanding Error Handling](https://powershell.org/2019/07/understanding-error-handling-in-powershell/)
- [SS64.com PowerShell Commands - Write-Error](https://ss64.com/ps/write-error.html)
- [Microsoft Learn - Powershell Error Handling Best Practices](https://docs.microsoft.com/en-us/powershell/scripting/learn/deep-dives/everything-about-exceptions)