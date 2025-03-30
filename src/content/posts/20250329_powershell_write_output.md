---
title:  【PowerShell】Write-Outputコマンドの基礎と活用方法
published: 2025-03-29
description: "PowerShellのWrite-Outputコマンドレットの使い方を解説"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-29

## 概要

`Write-Output`は、PowerShellにおいてオブジェクトをパイプラインに出力するための基本的なコマンドレットである。単にコンソールに表示するだけでなく、後続のコマンドレットにデータを渡すためにも使用される。PowerShellでは、`echo`や`out`という別名でも利用可能である。

## 基本的な使い方

### 単純な文字列出力

最も基本的な使い方は、文字列を出力することである。

```powershell
Write-Output "Hello, PowerShell!"
```

### 複数の値を出力

複数の値をカンマ区切りで指定すると、それぞれの値が個別に出力される。

```powershell
Write-Output "値1", "値2", "値3"
```

### 変数の内容を出力

変数に格納された値を出力できる。

```powershell
$name = "PowerShell"
Write-Output "こんにちは、$name!"
```

### 配列の出力

配列を出力すると、各要素が別々に表示される。

```powershell
$array = @("要素1", "要素2", "要素3")
Write-Output $array
```

### オブジェクト出力

オブジェクトを出力すると、PowerShellの既定の形式で表示される。

```powershell
$process = Get-Process -Id $PID
Write-Output $process
```

## 応用的な使い方

### -NoEnumerateパラメータ

`-NoEnumerate`パラメータを使用すると、配列を単一のオブジェクトとして出力する。

```powershell
$array = @("要素1", "要素2", "要素3")
Write-Output -NoEnumerate $array
```

### パイプラインでの活用

`Write-Output`はパイプラインの一部として他のコマンドレットと組み合わせて使用できる。

```powershell
Write-Output "テキスト1", "テキスト2", "テキスト3" | Where-Object { $_ -like "*2*" }
```

### 別名の使用

`Write-Output`は`echo`や`out`という別名でも使用できる。

```powershell
echo "Hello, PowerShell!"
```

### Write-HostとWrite-Outputの違い

`Write-Host`はコンソールに直接出力するのに対し、`Write-Output`はパイプラインにオブジェクトを送る。

```powershell
# これはパイプラインにオブジェクトを送るため、後続の処理に使える
Write-Output "テスト" | Get-Member

# これは単にコンソールに表示するだけで、パイプラインには何も送らない
Write-Host "テスト" | Get-Member  # Get-Memberは何も受け取れない
```

## ハンズオン：Write-Outputの実践

以下は各種Write-Outputコマンドとその出力例である。

1. **基本的な文字列出力**

```powershell
Write-Output "PowerShellを学習中です"
```

出力例：
```
PowerShellを学習中です
```

2. **複数の値の出力**

```powershell
Write-Output "値1", "値2", "値3"
```

出力例：
```
値1
値2
値3
```

3. **変数の出力**

```powershell
$version = $PSVersionTable.PSVersion
Write-Output "PowerShellのバージョン: $version"
```

出力例：
```
PowerShellのバージョン: 7.3.0
```

4. **配列の出力（通常）**

```powershell
$fruits = @("りんご", "バナナ", "オレンジ")
Write-Output $fruits
```

出力例：
```
りんご
バナナ
オレンジ
```

5. **配列の出力（NoEnumerate使用）**

```powershell
$fruits = @("りんご", "バナナ", "オレンジ")
Write-Output -NoEnumerate $fruits
```

出力例：
```
System.Object[] (配列オブジェクトとして表示)
```

6. **パイプラインでのフィルタリング**

```powershell
Write-Output "file1.txt", "file2.doc", "file3.pdf" | Where-Object { $_ -like "*.txt" }
```

出力例：
```
file1.txt
```

7. **別名を使った出力**

```powershell
echo "これはechoの別名を使った出力です"
```

出力例：
```
これはechoの別名を使った出力です
```

8. **オブジェクトのプロパティ選択**

```powershell
$process = Get-Process -Id $PID
Write-Output $process | Select-Object Name, Id, CPU
```

出力例：
```
Name          Id    CPU
----          --    ---
powershell  1234  10.53
```

## 対応PowerShellバージョン

Write-Outputコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 1.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

## 参考サイト

- [Microsoft公式ドキュメント: Write-Output](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.utility/write-output)
- [PowerShell.org - Understanding Output in PowerShell](https://powershell.org/2013/09/understanding-output-in-powershell/)
- [SS64.com PowerShell Commands - Write-Output](https://ss64.com/ps/write-output.html)