---
title:  【PowerShell】Wait-Processコマンドでプロセス終了を待機する方法
published: 2025-03-30
description: "PowerShellのWait-Processコマンドレットの基本的な使い方と活用例"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Wait-Process`はPowerShellにおいて、指定したプロセスが終了するまで実行を一時停止するためのコマンドレットである。バッチ処理やスクリプト内で特定のプロセスの完了を待ってから次の処理に進みたい場合に非常に有用である。プロセスはIDまたは名前で指定できるほか、タイムアウト設定も可能である。このコマンドレットを使うことで、プロセス間の依存関係を管理し、スクリプトの実行順序を制御することができる。

## 基本的な使い方

### プロセスIDによる待機

特定のプロセスIDのプロセスが終了するまで待機する。

```powershell
Wait-Process -Id 1234
```

### プロセス名による待機

指定した名前のプロセスが終了するまで待機する。同じ名前の複数のプロセスが存在する場合はすべてのプロセスが終了するまで待機する。

```powershell
Wait-Process -Name "notepad"
```

### タイムアウト設定

指定した時間だけ待機し、その時間が過ぎてもプロセスが終了しない場合は処理を続行する。

```powershell
Wait-Process -Name "calc" -Timeout 10
```

### パイプラインを使った待機

Get-Processの結果をパイプラインでWait-Processに渡すこともできる。

```powershell
Get-Process -Name "notepad" | Wait-Process
```

## 応用的な使い方

### 複数プロセスの待機

複数のプロセスを同時に待機することができる。

```powershell
Wait-Process -Name "notepad", "calc", "mspaint"
```

### エラー処理の追加

プロセスが見つからない場合や待機中にエラーが発生した場合の処理を追加する。

```powershell
try {
    Wait-Process -Name "non_existent_process" -ErrorAction Stop
} catch {
    Write-Host "プロセスが見つからないか、エラーが発生しました: $_"
}
```

### コンソール出力を伴う待機

プロセスの待機中に情報を表示する。

```powershell
$process = Start-Process -FilePath "notepad" -PassThru
Write-Host "プロセス $($process.Id) の終了を待機しています..."
Wait-Process -Id $process.Id
Write-Host "プロセスが終了しました"
```

### 待機後の処理

プロセスが終了した後に特定の処理を行う。

```powershell
$process = Start-Process -FilePath "calc" -PassThru
Wait-Process -Id $process.Id
if ($process.ExitCode -eq 0) {
    Write-Host "プロセスは正常に終了しました"
} else {
    Write-Host "プロセスはエラーコード $($process.ExitCode) で終了しました"
}
```

## ハンズオン：Wait-Processの実践

以下はWait-Processコマンドの実践例とその出力である。

1. **プロセスを起動して終了を待つ基本例**

```powershell
# メモ帳を起動
$notepad = Start-Process -FilePath "notepad" -PassThru
Write-Host "メモ帳を起動しました (PID: $($notepad.Id))。終了を待機しています..."

# プロセスが終了するまで待機
Wait-Process -Id $notepad.Id

Write-Host "メモ帳プロセスが終了しました"
```

出力例（メモ帳を手動で閉じた後）：
```
メモ帳を起動しました (PID: 1234)。終了を待機しています...
メモ帳プロセスが終了しました
```

2. **タイムアウト付きの待機**

```powershell
# 電卓を起動
$calc = Start-Process -FilePath "calc" -PassThru
Write-Host "電卓を起動しました (PID: $($calc.Id))。5秒間待機します..."

try {
    # 5秒間だけ待機
    Wait-Process -Id $calc.Id -Timeout 5 -ErrorAction Stop
    Write-Host "電卓が5秒以内に終了しました"
} catch {
    if ($_.Exception.Message -like "*timed out*") {
        Write-Host "タイムアウトしました。電卓はまだ実行中です"
        Stop-Process -Id $calc.Id
        Write-Host "電卓を強制終了しました"
    } else {
        Write-Host "エラーが発生しました: $_"
    }
}
```

出力例（タイムアウトした場合）：
```
電卓を起動しました (PID: 5678)。5秒間待機します...
タイムアウトしました。電卓はまだ実行中です
電卓を強制終了しました
```

3. **複数プロセスの待機**

```powershell
# 複数のプロセスを起動
$notepad = Start-Process -FilePath "notepad" -PassThru
$mspaint = Start-Process -FilePath "mspaint" -PassThru
Write-Host "メモ帳 (PID: $($notepad.Id)) とペイント (PID: $($mspaint.Id)) を起動しました"

# どちらか一方のプロセスが終了するのを待つ
$processToWait = $notepad.Id
Write-Host "プロセス $processToWait の終了を待機しています..."
Wait-Process -Id $processToWait

# 残りのプロセスを終了
if (Get-Process -Id $mspaint.Id -ErrorAction SilentlyContinue) {
    Stop-Process -Id $mspaint.Id
    Write-Host "ペイントを終了しました"
}

Write-Host "すべてのプロセスが終了しました"
```

出力例（メモ帳を手動で閉じた場合）：
```
メモ帳 (PID: 1234) とペイント (PID: 5678) を起動しました
プロセス 1234 の終了を待機しています...
ペイントを終了しました
すべてのプロセスが終了しました
```

4. **プロセス名による待機**

```powershell
# 複数のメモ帳インスタンスを起動
Start-Process -FilePath "notepad" -WindowStyle Minimized
Start-Process -FilePath "notepad" -WindowStyle Minimized
$count = (Get-Process -Name "notepad").Count
Write-Host "$count 個のメモ帳プロセスが実行中です。全てのメモ帳が閉じられるまで待機します..."

# すべてのメモ帳が閉じられるまで待機
Wait-Process -Name "notepad" -ErrorAction SilentlyContinue

Write-Host "すべてのメモ帳プロセスが終了しました"
```

出力例（すべてのメモ帳を手動で閉じた後）：
```
2 個のメモ帳プロセスが実行中です。全てのメモ帳が閉じられるまで待機します...
すべてのメモ帳プロセスが終了しました
```

5. **パイプラインを使った待機**

```powershell
# wordpadを起動
Start-Process -FilePath "wordpad" -PassThru
Write-Host "Wordpadを起動しました。プロセスを取得して待機します..."

# プロセスをパイプラインで渡して待機
Get-Process -Name "wordpad" | ForEach-Object {
    Write-Host "Wordpad (PID: $($_.Id)) の終了を待機しています..."
} | Wait-Process

Write-Host "Wordpadプロセスが終了しました"
```

出力例（Wordpadを手動で閉じた後）：
```
Wordpadを起動しました。プロセスを取得して待機します...
Wordpad (PID: 9876) の終了を待機しています...
Wordpadプロセスが終了しました
```

6. **存在しないプロセスのエラーハンドリング**

```powershell
try {
    Write-Host "存在しないプロセス名を指定しています..."
    Wait-Process -Name "non_existent_process_name" -ErrorAction Stop
} catch {
    Write-Host "エラーが発生しました: $_"
}
```

出力例：
```
存在しないプロセス名を指定しています...
エラーが発生しました: プロセス non_existent_process_name が見つかりません。
```

7. **終了コードの確認**

```powershell
# 終了コードを返すプログラムの例
$tempScript = New-TemporaryFile | Rename-Item -NewName { $_.Name + ".ps1" } -PassThru
Add-Content -Path $tempScript.FullName -Value "exit 123"

$process = Start-Process -FilePath "powershell" -ArgumentList "-File `"$($tempScript.FullName)`"" -PassThru
Write-Host "プロセス $($process.Id) の終了を待機しています..."
Wait-Process -Id $process.Id
Write-Host "プロセスの終了コード: $($process.ExitCode)"

# 一時ファイルを削除
Remove-Item -Path $tempScript.FullName
```

出力例：
```
プロセス 3456 の終了を待機しています...
プロセスの終了コード: 123
```

## 対応PowerShellバージョン

Wait-Processコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 2.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

## 参考サイト

- [Microsoft公式ドキュメント: Wait-Process](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.management/wait-process)
- [PowerShell.org - Process Management in PowerShell](https://powershell.org/2018/06/process-management-in-powershell/)
- [SS64.com PowerShell Commands - Wait-Process](https://ss64.com/ps/wait-process.html)
- [Microsoft Learn - PowerShellのプロセス管理](https://docs.microsoft.com/ja-jp/powershell/scripting/samples/managing-processes-with-process-cmdlets)