---
title:  【PowerShell】Write-Progressコマンドによる進捗表示の実装方法
published: 2025-03-29
description: "PowerShellのWrite-Progressコマンドの基礎から応用までを解説"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-29

## 概要

`Write-Progress`は、PowerShellで長時間実行されるスクリプトやコマンドの進捗状況をコンソール上に視覚的に表示するためのコマンドレットである。ユーザーに処理の進行状況をリアルタイムで伝えることで、実行中のプロセスが停止しているわけではなく、正常に動作していることを示すことができる。

## 基本的な使い方

### 単純な進捗バーの表示

最も基本的な使い方は、アクティビティの説明と進捗率を指定することである。

```powershell
Write-Progress -Activity "ファイルをコピーしています" -PercentComplete 50
```

### 状態情報の追加

現在の状態を示す情報を追加することもできる。

```powershell
Write-Progress -Activity "ファイルのダウンロード" -Status "25 / 100 ファイル完了" -PercentComplete 25
```

### 現在の処理項目の表示

現在処理中の特定のアイテムを表示することができる。

```powershell
Write-Progress -Activity "ファイルの処理" -Status "ファイル処理中" -CurrentOperation "document.docx" -PercentComplete 30
```

### 進捗バーのループでの使用

進捗バーは通常、繰り返し処理の中で使用される。

```powershell
1..100 | ForEach-Object {
    Write-Progress -Activity "カウントアップ中" -Status "現在の数値: $_" -PercentComplete $_
    Start-Sleep -Milliseconds 100
}
```

## 応用的な使い方

### ネストされた進捗バーの作成

進捗バーはネストして表示することができる。親と子のタスクを視覚的に区別できる。

```powershell
1..10 | ForEach-Object -Process {
    $parentProgress = $_*10
    Write-Progress -Id 0 -Activity "メインタスク" -Status "処理中..." -PercentComplete $parentProgress
    
    1..10 | ForEach-Object -Process {
        $childProgress = $_*10
        Write-Progress -Id 1 -ParentId 0 -Activity "サブタスク" -Status "サブ処理中..." -PercentComplete $childProgress
        Start-Sleep -Milliseconds 100
    }
}
```

### 進捗バーの完了

処理が完了した進捗バーは、`-Completed`パラメータを使って明示的に完了を示すことができる。

```powershell
Write-Progress -Activity "データ処理" -Completed
```

### SecondsRemaining パラメータの使用

残り時間を表示することもできる。

```powershell
Write-Progress -Activity "ファイルのアップロード" -Status "処理中..." -PercentComplete 75 -SecondsRemaining 30
```

## ハンズオン：Write-Progressの実践

以下は一連のWrite-Progressコマンドとその動作例である。

1. **基本的な進捗表示**

```powershell
for ($i = 1; $i -le 100; $i++) {
    Write-Progress -Activity "基本的な進捗バーのデモ" -Status "進行中..." -PercentComplete $i
    Start-Sleep -Milliseconds 50
}
```

出力例：
```
基本的な進捗バーのデモ
進行中...
[===============================================>                      ] 75%
```

2. **状態情報を含む進捗表示**

```powershell
$files = 1..20
foreach ($file in $files) {
    $percentComplete = ($file / $files.Count) * 100
    Write-Progress -Activity "ファイル処理中" -Status "$file / $($files.Count) ファイル完了" -PercentComplete $percentComplete
    Start-Sleep -Milliseconds 200
}
```

出力例：
```
ファイル処理中
10 / 20 ファイル完了
[===========================>                                          ] 50%
```

3. **現在の操作情報を含む進捗表示**

```powershell
$files = @("file1.txt", "image.jpg", "document.docx", "data.xlsx")
foreach ($file in $files) {
    $index = [array]::IndexOf($files, $file)
    $percentComplete = ($index / $files.Count) * 100
    Write-Progress -Activity "ファイル変換" -Status "処理中..." -CurrentOperation "現在のファイル: $file" -PercentComplete $percentComplete
    Start-Sleep -Seconds 1
}
```

出力例：
```
ファイル変換
処理中...
現在のファイル: image.jpg
[=============>                                                        ] 25%
```

4. **ネストされた進捗バーの表示**

```powershell
foreach ($i in 1..3) {
    Write-Progress -Id 0 -Activity "メインプロセス" -Status "ステップ $i/3" -PercentComplete ($i/3*100)
    
    foreach ($j in 1..5) {
        Write-Progress -Id 1 -ParentId 0 -Activity "サブプロセス" -Status "サブステップ $j/5" -PercentComplete ($j/5*100)
        Start-Sleep -Milliseconds 200
    }
}
```

出力例：
```
メインプロセス
ステップ 2/3
[======================>                                               ] 66%

サブプロセス
サブステップ 3/5
[==================>                                                   ] 60%
```

5. **残り時間の表示**

```powershell
$total = 100
for ($i = 1; $i -le $total; $i++) {
    $secondsRemaining = ($total - $i) / 10
    Write-Progress -Activity "データダウンロード" -Status "進行中..." -PercentComplete $i -SecondsRemaining $secondsRemaining
    Start-Sleep -Milliseconds 100
}
```

出力例：
```
データダウンロード
進行中...
[============================================>                         ] 70% 残り3秒
```

6. **進捗バーの完了表示**

```powershell
Write-Progress -Activity "インストール" -Status "プログラムの設定" -PercentComplete 80
Start-Sleep -Seconds 2
Write-Progress -Activity "インストール" -Completed
```

出力例（完了前）：
```
インストール
プログラムの設定
[========================================>                             ] 80%
```

完了後、進捗バーは消えて次のコマンドプロンプトが表示される

## 対応PowerShellバージョン

Write-Progressコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 2.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

## 参考サイト

- [Microsoft公式ドキュメント: Write-Progress](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.utility/write-progress)
- [PowerShell.org - Using Write-Progress](https://powershell.org/2013/08/using-write-progress/)
- [SS64.com PowerShell Commands - Write-Progress](https://ss64.com/ps/write-progress.html)