---
title:  【PowerShell】Wait-Jobコマンドでバックグラウンドジョブの完了を待機する方法
published: 2025-03-30
description: "PowerShellのWait-Jobコマンドレットを使用したバックグラウンドジョブ管理の基礎"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Wait-Job`はPowerShellにおいて、バックグラウンドで実行されているジョブの完了を待機するためのコマンドレットである。これを使用することで、ジョブが終了するまでスクリプトの実行を一時停止し、ジョブの完了後に続行することができる。バックグラウンドジョブを使うと複数の処理を並行して実行できるが、後続の処理がジョブの結果に依存する場合などに`Wait-Job`を使用して同期をとることが可能である。

## 基本的な使い方

### ジョブを起動して待機する

基本的な使い方は、Start-Jobでジョブを起動し、それを待機することである。

```powershell
$job = Start-Job -ScriptBlock { Start-Sleep -Seconds 5; "ジョブが完了しました" }
Wait-Job -Job $job
Receive-Job -Job $job
```

### ジョブIDで待機する

ジョブIDを指定して待機することもできる。

```powershell
$job = Start-Job -ScriptBlock { Start-Sleep -Seconds 5; "ジョブが完了しました" }
Wait-Job -Id $job.Id
```

### ジョブ名で待機する

ジョブに名前をつけて、その名前で待機することもできる。

```powershell
Start-Job -Name "長時間処理" -ScriptBlock { Start-Sleep -Seconds 10; "長時間処理完了" }
Wait-Job -Name "長時間処理"
```

### タイムアウトを設定する

指定した時間だけ待機し、その時間内にジョブが終了しない場合は処理を続行する。

```powershell
$job = Start-Job -ScriptBlock { Start-Sleep -Seconds 30; "処理完了" }
Wait-Job -Job $job -Timeout 10
```

## 応用的な使い方

### 複数ジョブのいずれかの完了を待つ

`-Any`パラメータを使用すると、指定した複数のジョブのうち、どれか1つが完了するのを待機できる。

```powershell
$job1 = Start-Job -Name "Job1" -ScriptBlock { Start-Sleep -Seconds 5; "Job1完了" }
$job2 = Start-Job -Name "Job2" -ScriptBlock { Start-Sleep -Seconds 10; "Job2完了" }
Wait-Job -Job $job1, $job2 -Any
```

### フィルターを使用して特定のジョブを待機

`-Filter`パラメータを使用して、特定の条件に一致するジョブを待機できる。

```powershell
Start-Job -Name "重要処理" -ScriptBlock { Start-Sleep -Seconds 5 }
Start-Job -Name "通常処理" -ScriptBlock { Start-Sleep -Seconds 8 }
Wait-Job -Filter @{Name="重要*"}
```

### ジョブの状態を確認しながら待機

ジョブの状態を監視しながら待機することも可能である。

```powershell
$job = Start-Job -ScriptBlock { 
    foreach ($i in 1..10) {
        Write-Progress -Activity "処理中" -Status "ステップ $i/10" -PercentComplete ($i*10)
        Start-Sleep -Seconds 1
    }
    "処理完了"
}

while ($job.State -eq "Running") {
    $jobInfo = Get-Job -Id $job.Id
    Write-Host "ジョブの状態: $($jobInfo.State) - $(Get-Date)"
    Start-Sleep -Seconds 2
    if ($job.State -ne "Running") {
        Wait-Job -Job $job
        break
    }
}
```

### Wait-JobとReceive-Jobの組み合わせ

ジョブの完了を待ってから結果を受け取る典型的なパターンである。

```powershell
$job = Start-Job -ScriptBlock { 
    $result = 0
    for ($i = 1; $i -le 5; $i++) {
        $result += $i
        Start-Sleep -Seconds 1
    }
    return $result
}

Wait-Job -Job $job
$result = Receive-Job -Job $job
Remove-Job -Job $job
Write-Host "計算結果: $result"
```

## ハンズオン：Wait-Jobの実践

以下はWait-Jobコマンドの実践例とその出力である。

1. **基本的なジョブ待機**

```powershell
# シンプルなジョブを起動
$job = Start-Job -ScriptBlock {
    Start-Sleep -Seconds 3
    return "基本的なジョブが完了しました"
}

Write-Host "ジョブを起動しました（ID: $($job.Id)）- $(Get-Date)"
Wait-Job -Job $job | Out-Null
Write-Host "ジョブが完了しました - $(Get-Date)"
$result = Receive-Job -Job $job
Write-Host "ジョブの結果: $result"
Remove-Job -Job $job
```

出力例：
```
ジョブを起動しました（ID: 1） - 2025/03/30 12:00:00
ジョブが完了しました - 2025/03/30 12:00:03
ジョブの結果: 基本的なジョブが完了しました
```

2. **タイムアウトを設定したジョブ待機**

```powershell
$job = Start-Job -ScriptBlock {
    Start-Sleep -Seconds 15
    return "長時間処理が完了しました"
}

Write-Host "ジョブを起動し、5秒間だけ待機します - $(Get-Date)"
$completed = Wait-Job -Job $job -Timeout 5
if ($completed) {
    Write-Host "ジョブはタイムアウト前に完了しました"
    Receive-Job -Job $job
} else {
    Write-Host "ジョブはタイムアウトしました。まだ実行中です - $(Get-Date)"
    Write-Host "ジョブの現在の状態: $($job.State)"
    # 追加で待機する場合
    Write-Host "ジョブの完了を引き続き待機します..."
    Wait-Job -Job $job | Out-Null
    Write-Host "ジョブが完了しました - $(Get-Date)"
    Receive-Job -Job $job
}
Remove-Job -Job $job
```

出力例：
```
ジョブを起動し、5秒間だけ待機します - 2025/03/30 12:05:00
ジョブはタイムアウトしました。まだ実行中です - 2025/03/30 12:05:05
ジョブの現在の状態: Running
ジョブの完了を引き続き待機します...
ジョブが完了しました - 2025/03/30 12:05:15
長時間処理が完了しました
```

3. **複数のジョブを管理**

```powershell
# 複数のジョブを起動
$job1 = Start-Job -Name "Process1" -ScriptBlock {
    Start-Sleep -Seconds 3
    return "プロセス1が完了しました"
}

$job2 = Start-Job -Name "Process2" -ScriptBlock {
    Start-Sleep -Seconds 6
    return "プロセス2が完了しました"
}

$job3 = Start-Job -Name "Process3" -ScriptBlock {
    Start-Sleep -Seconds 9
    return "プロセス3が完了しました"
}

# すべてのジョブを配列に格納
$allJobs = @($job1, $job2, $job3)

Write-Host "3つのジョブが起動されました - $(Get-Date)"
Write-Host "最初のジョブが完了するのを待機します..."

# どれか1つのジョブが完了するまで待機
$completedJob = Wait-Job -Job $allJobs -Any
Write-Host "ジョブ「$($completedJob.Name)」が最初に完了しました - $(Get-Date)"
Write-Host "結果: $(Receive-Job -Job $completedJob)"

# 残りのジョブが完了するのを待機
Write-Host "残りのジョブの完了を待機します..."
Wait-Job -Job $allJobs | Out-Null
Write-Host "すべてのジョブが完了しました - $(Get-Date)"

# すべての結果を表示
foreach ($job in $allJobs) {
    Write-Host "$($job.Name): $(Receive-Job -Job $job -Keep)"
}

# ジョブをクリーンアップ
Remove-Job -Job $allJobs
```

出力例：
```
3つのジョブが起動されました - 2025/03/30 12:10:00
最初のジョブが完了するのを待機します...
ジョブ「Process1」が最初に完了しました - 2025/03/30 12:10:03
結果: プロセス1が完了しました
残りのジョブの完了を待機します...
すべてのジョブが完了しました - 2025/03/30 12:10:09
Process1: プロセス1が完了しました
Process2: プロセス2が完了しました
Process3: プロセス3が完了しました
```

4. **フィルターを使用したジョブ待機**

```powershell
# 複数のジョブを異なる名前で起動
Start-Job -Name "DataProcess_1" -ScriptBlock { Start-Sleep -Seconds 3; return "データ処理1完了" }
Start-Job -Name "DataProcess_2" -ScriptBlock { Start-Sleep -Seconds 5; return "データ処理2完了" }
Start-Job -Name "Maintenance" -ScriptBlock { Start-Sleep -Seconds 2; return "メンテナンス完了" }

Write-Host "3つのジョブを起動しました - $(Get-Date)"
Write-Host "「DataProcess」で始まるジョブの完了を待機します..."

# フィルターでデータ処理ジョブを待機
$dataJobs = Wait-Job -Filter @{Name="DataProcess*"}
Write-Host "データ処理ジョブが完了しました - $(Get-Date)"

# 結果を表示
foreach ($job in $dataJobs) {
    Write-Host "$($job.Name): $(Receive-Job -Job $job)"
}

# すべてのジョブの完了を確認
Wait-Job -State Running
Write-Host "すべてのジョブが完了しました"

# ジョブをクリーンアップ
Get-Job | Remove-Job
```

出力例：
```
3つのジョブを起動しました - 2025/03/30 12:15:00
「DataProcess」で始まるジョブの完了を待機します...
データ処理ジョブが完了しました - 2025/03/30 12:15:05
DataProcess_1: データ処理1完了
DataProcess_2: データ処理2完了
すべてのジョブが完了しました
```

5. **エラー処理を含んだジョブ待機**

```powershell
# 正常なジョブとエラーを含むジョブを起動
$goodJob = Start-Job -Name "GoodJob" -ScriptBlock { 
    Start-Sleep -Seconds 2
    return "正常に完了しました"
}

$errorJob = Start-Job -Name "ErrorJob" -ScriptBlock { 
    Start-Sleep -Seconds 4
    throw "エラーが発生しました"
}

Write-Host "2つのジョブを起動しました - $(Get-Date)"
Write-Host "すべてのジョブの完了を待機します..."

# すべてのジョブを待機
Wait-Job -Job $goodJob, $errorJob | Out-Null
Write-Host "すべてのジョブが完了しました - $(Get-Date)"

# 結果とエラーを処理
foreach ($job in $goodJob, $errorJob) {
    Write-Host "`n$($job.Name) の状態: $($job.State)"
    try {
        $result = Receive-Job -Job $job -ErrorAction Stop
        Write-Host "結果: $result"
    } catch {
        Write-Host "エラー: $_" -ForegroundColor Red
    }
}

# ジョブをクリーンアップ
Remove-Job -Job $goodJob, $errorJob
```

出力例：
```
2つのジョブを起動しました - 2025/03/30 12:20:00
すべてのジョブの完了を待機します...
すべてのジョブが完了しました - 2025/03/30 12:20:04

GoodJob の状態: Completed
結果: 正常に完了しました

ErrorJob の状態: Failed
エラー: エラーが発生しました
```

## 対応PowerShellバージョン

Wait-Jobコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 2.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

## 参考サイト

- [Microsoft公式ドキュメント: Wait-Job](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.core/wait-job)
- [PowerShell.org - Background Jobs in PowerShell](https://powershell.org/2013/04/powershell-background-jobs/)
- [SS64.com PowerShell Commands - Wait-Job](https://ss64.com/ps/wait-job.html)
- [Microsoft Learn - PowerShellのジョブの使用](https://docs.microsoft.com/ja-jp/powershell/scripting/learn/ps101/10-script-modules)