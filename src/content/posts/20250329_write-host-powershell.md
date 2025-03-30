---
title:  【PowerShell】Write-Hostコマンドの基本と活用法
published: 2025-03-29
description: "PowerShellのWrite-Hostコマンドを使ったコンソール出力の方法とテクニック"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-29

## 概要

`Write-Host`はPowerShellでコンソールに直接テキストや情報を表示するためのコマンドレットである。このコマンドレットの特徴は、パイプラインにデータを送らず、常にホストプログラム（通常はコンソール画面）に直接出力することである。テキストの色や背景色を変更したり、複数のオブジェクトを一行に表示したりすることができ、視覚的なフィードバックを提供する際に役立つ。

## 基本的な使い方

### 単純なテキスト出力

最も基本的な使い方は、単純にテキストを表示することである。

```powershell
Write-Host "Hello, PowerShell!"
```

### カラー表示

`-ForegroundColor`パラメータと`-BackgroundColor`パラメータを使用して、テキストの色と背景色を指定できる。

```powershell
Write-Host "エラーメッセージ" -ForegroundColor Red
Write-Host "成功メッセージ" -ForegroundColor Green
Write-Host "警告メッセージ" -ForegroundColor Yellow -BackgroundColor Black
```

### 複数オブジェクトの表示

複数の項目を一行に表示することができる。

```powershell
Write-Host "ユーザー名:" -NoNewline
Write-Host " Admin" -ForegroundColor Cyan
```

### 区切り文字の指定

`-Separator`パラメータを使用して、複数のオブジェクト間の区切り文字を指定できる。

```powershell
Write-Host "項目1", "項目2", "項目3" -Separator ", "
```

## 応用的な使い方

### 変数と組み合わせて使用

変数の内容をカラー表示することができる。

```powershell
$userName = "PowerShell User"
$status = "オンライン"
Write-Host "現在のユーザー: " -NoNewline
Write-Host $userName -ForegroundColor Cyan -NoNewline
Write-Host " ($status)" -ForegroundColor Green
```

### 条件分岐による色の変更

条件によって表示色を変えることができる。

```powershell
$errorCount = 0
if ($errorCount -gt 0) {
    Write-Host "エラーが発生しています！" -ForegroundColor Red
} else {
    Write-Host "正常に動作しています" -ForegroundColor Green
}
```

### プログレスバーの表現

単純なプログレスバーを表現することも可能である。

```powershell
for ($i = 0; $i -le 100; $i += 10) {
    Write-Host "`rロード中: $i% " -NoNewline
    Write-Host "["
    Write-Host "#" * ($i / 2) -NoNewline -ForegroundColor Green
    Write-Host " " * (50 - $i / 2) -NoNewline
    Write-Host "]" -NoNewline
    Start-Sleep -Milliseconds 200
}
Write-Host "`rロード完了!     " -ForegroundColor Green
```

## ハンズオン：Write-Hostの実践

以下はWrite-Hostコマンドの各種使用例とその出力である。

1. **基本的なテキスト出力**

```powershell
Write-Host "PowerShellへようこそ！"
```

出力例：
```
PowerShellへようこそ！
```

2. **カラーテキストの表示**

```powershell
Write-Host "これは赤色のテキストです" -ForegroundColor Red
Write-Host "これは緑色のテキストです" -ForegroundColor Green
Write-Host "これは白色テキスト＋青色背景です" -ForegroundColor White -BackgroundColor Blue
```

出力例（色付きで表示される）：
```
これは赤色のテキストです
これは緑色のテキストです
これは白色テキスト＋青色背景です
```

3. **改行なしの出力**

```powershell
Write-Host "こんにちは " -NoNewline
Write-Host "PowerShell!" -ForegroundColor Cyan
```

出力例：
```
こんにちは PowerShell!
```

4. **区切り文字を使った複数オブジェクトの表示**

```powershell
Write-Host "りんご", "バナナ", "オレンジ" -Separator " > "
```

出力例：
```
りんご > バナナ > オレンジ
```

5. **変数と組み合わせた使用例**

```powershell
$count = 5
$total = 10
Write-Host "処理状況: " -NoNewline
Write-Host "$count/$total" -ForegroundColor Cyan -NoNewline
Write-Host " 項目完了" -ForegroundColor Green
```

出力例：
```
処理状況: 5/10 項目完了
```

6. **条件分岐と組み合わせた例**

```powershell
$temperature = 35
Write-Host "現在の温度: $temperature°C " -NoNewline
if ($temperature -gt 30) {
    Write-Host "(高温注意)" -ForegroundColor Red
} elseif ($temperature -lt 5) {
    Write-Host "(低温注意)" -ForegroundColor Blue
} else {
    Write-Host "(適温)" -ForegroundColor Green
}
```

出力例：
```
現在の温度: 35°C (高温注意)
```

7. **特殊文字とエスケープシーケンスの使用**

```powershell
Write-Host "タブを使用:\t区切り"
Write-Host "バックティック(``)を使った特殊文字: `n改行"
```

出力例：
```
タブを使用:	区切り
バックティック(`)を使った特殊文字: 
改行
```

8. **オブジェクトの直接表示**

```powershell
$process = Get-Process -Id $PID | Select-Object Name, Id
Write-Host "現在のプロセス: " -NoNewline
Write-Host $process -ForegroundColor Cyan
```

出力例：
```
現在のプロセス: @{Name=powershell; Id=1234}
```

## 対応PowerShellバージョン

Write-Hostコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 1.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

## 参考サイト

- [Microsoft公式ドキュメント: Write-Host](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.utility/write-host)
- [PowerShell.org - Understanding Output in PowerShell](https://powershell.org/2013/09/understanding-output-in-powershell/)
- [SS64.com PowerShell Commands - Write-Host](https://ss64.com/ps/write-host.html)