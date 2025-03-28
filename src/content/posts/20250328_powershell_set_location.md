---
title:  【PowerShell】Set-Locationコマンドの基本と応用
published: 2025-03-28
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-28

## 概要

PowerShellの`Set-Location`コマンドレット（エイリアス：`cd`, `chdir`, `sl`）は、現在の作業場所（カレントディレクトリ）を移動するためのコマンドである。ファイルシステムだけでなく、レジストリやその他のPowerShellドライブなど、あらゆるPowerShellプロバイダーでの位置の移動に使用できる。このコマンドレットは、Windowsコマンドプロンプトの`cd`コマンドやUNIX系OSの`cd`コマンドと類似の機能を提供するが、PowerShellの環境変数やプロバイダーシステムと統合されているためより柔軟性が高い。

## 基本的な使い方

### 特定のディレクトリへの移動

```powershell
# Cドライブの特定のフォルダへ移動
Set-Location -Path C:\Windows\System32
```

出力例:
```
# 成功時は出力なし（プロンプトが変わるのみ）
PS C:\Windows\System32>
```

### 相対パスでの移動

```powershell
# 現在の場所から相対的に移動
Set-Location -Path ..\Program Files
```

出力例:
```
PS C:\Windows> Set-Location -Path ..\Program Files
PS C:\Program Files>
```

### エイリアスを使った簡潔な記述

```powershell
# cdエイリアスを使用（最も一般的）
cd C:\Users\ss0832\Documents
```

出力例:
```
PS C:\> cd C:\Users\ss0832\Documents
PS C:\Users\ss0832\Documents>
```

## 全オプションに対するハンズオン

### -Path

目的の場所のパスを指定する。

```powershell
# 基本的なディレクトリ移動
Set-Location -Path C:\Temp

# 環境変数を使用したパス指定
Set-Location -Path $env:USERPROFILE

# 相対パス指定
Set-Location -Path ..\..
```

出力例:
```
PS C:\Windows> Set-Location -Path C:\Temp
PS C:\Temp>

PS C:\Temp> Set-Location -Path $env:USERPROFILE
PS C:\Users\ss0832>

PS C:\Users\ss0832> Set-Location -Path ..\..
PS C:\>
```

### -LiteralPath

特殊文字（ワイルドカードなど）を含むパスをそのままのパスとして使用する。

```powershell
# 角括弧を含むフォルダ名（通常はワイルドカードとして解釈される）
Set-Location -LiteralPath "C:\Test[Folder]"
```

出力例:
```
PS C:\> Set-Location -LiteralPath "C:\Test[Folder]"
PS C:\Test[Folder]>
```

### -PassThru

現在の場所を表すオブジェクトをパイプラインに返す。

```powershell
# 移動後の場所を変数に保存
$newLoc = Set-Location -Path C:\Windows -PassThru
$newLoc
```

出力例:
```
PS C:\> $newLoc = Set-Location -Path C:\Windows -PassThru
PS C:\Windows> $newLoc

Path
----
C:\Windows
```

### -StackName

指定された場所スタックを使用する。場所スタックは、以前のディレクトリを記憶するための仕組み。

```powershell
# 新しいスタックを作成し使用
Push-Location C:\Windows -StackName "WindowsNav"
Set-Location C:\Program Files
Set-Location -StackName "WindowsNav" -Path (Pop-Location -StackName "WindowsNav" -PassThru).Path
```

出力例:
```
PS C:\> Push-Location C:\Windows -StackName "WindowsNav"
PS C:\Windows> Set-Location C:\Program Files
PS C:\Program Files> Set-Location -StackName "WindowsNav" -Path (Pop-Location -StackName "WindowsNav" -PassThru).Path
PS C:\Windows>
```

## 特殊な使用方法

### ドライブ間の移動

```powershell
# C:ドライブからD:ドライブに移動
Set-Location D:
```

出力例:
```
PS C:\Users\ss0832> Set-Location D:
PS D:\>
```

### レジストリへの移動

```powershell
# レジストリHKLMドライブに移動
Set-Location HKLM:
```

出力例:
```
PS C:\> Set-Location HKLM:
PS HKLM:\>
```

### 以前の場所に戻る

```powershell
# 一時的に別の場所へ移動し、戻る
Push-Location C:\Windows
Set-Location C:\Program Files
Pop-Location
```

出力例:
```
PS C:\> Push-Location C:\Windows
PS C:\Windows> Set-Location C:\Program Files
PS C:\Program Files> Pop-Location
PS C:\Windows>
```

### 環境変数を使用したパス指定

```powershell
# ユーザープロファイルに移動
Set-Location -Path $env:USERPROFILE
```

出力例:
```
PS C:\Windows> Set-Location -Path $env:USERPROFILE
PS C:\Users\ss0832>
```

### ホームディレクトリへの移動

```powershell
# チルダを使用してホームディレクトリに移動
Set-Location ~
```

出力例:
```
PS C:\Windows> Set-Location ~
PS C:\Users\ss0832>
```

## PowerShellドライブプロバイダーとの連携

PowerShellは、ファイルシステム以外のデータストアにもアクセスできるプロバイダーシステムを持っている。

```powershell
# レジストリHKCUドライブに移動
Set-Location HKCU:
```

出力例:
```
PS C:\> Set-Location HKCU:
PS HKCU:\>
```

```powershell
# PowerShell変数ドライブに移動
Set-Location Variable:
```

出力例:
```
PS C:\> Set-Location Variable:
PS Variable:>
```

## 出力例の見方

`Set-Location`コマンドレットは、基本的に成功した場合は何も出力しない。代わりにPowerShellのプロンプトが変わることで、現在の場所が視覚的に示される。

1. **変更後のプロンプト**:
   - `PS C:\Windows>` のように、「PS」の後に現在のパスが表示される
   - パスの後に「>」が続く

2. **-PassThruオプション使用時の出力**:
   - `Path` プロパティを持つオブジェクトが返される
   - このオブジェクトは現在の場所を示しており、他のコマンドにパイプで渡すことができる

3. **エラーメッセージ**:
   - 指定したパスが存在しない場合: 「指定されたパスが見つかりません。」
   - アクセス権がない場合: 「アクセスが拒否されました。」
   - パスの形式が正しくない場合: 「パスの形式が無効です。」

