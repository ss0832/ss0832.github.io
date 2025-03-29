---
title:  【PowerShell】Write-Warningコマンドの基本と応用
published: 2025-03-29
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-29

## 概要

PowerShellの`Write-Warning`コマンドレットは、警告メッセージをユーザーに表示するための機能である。このコマンドレットを使用すると、黄色のテキスト（デフォルト設定）で警告メッセージが表示され、スクリプトの実行中に潜在的な問題や注意事項をユーザーに通知することができる。警告メッセージは通常のコマンド出力とは異なり、エラーストリームに出力されるが、エラーとは異なりスクリプトの実行は停止しない。また、`$WarningPreference`変数を使用して、警告メッセージの表示方法を制御することも可能である。

自作のpowershell関数を作成した際に警告を表示するときに役に立つ。

## 基本的な使い方

### 単純な警告メッセージの表示

```powershell
# 基本的な警告メッセージを表示
Write-Warning "これは警告メッセージです"
```

出力例:
```
WARNING: これは警告メッセージです
```

### 変数を含む警告メッセージの表示

```powershell
# 変数の値を含めた警告メッセージ
$diskSpace = 15
Write-Warning "ディスク空き容量が $diskSpace GB を下回っています。クリーンアップをお勧めします。"
```

出力例:
```
WARNING: ディスク空き容量が 15 GB を下回っています。クリーンアップをお勧めします。
```

## 全オプションに対するハンズオン

### -Message

表示する警告メッセージを指定する。これはWrite-Warningの最も基本的なパラメータである。

```powershell
# 明示的にMessageパラメータを使用
Write-Warning -Message "設定ファイルが見つかりません。デフォルト設定を使用します。"
```

出力例:
```
WARNING: 設定ファイルが見つかりません。デフォルト設定を使用します。
```

### パイプラインからの入力

テキストデータをパイプラインでWrite-Warningに渡すことができる。

```powershell
# パイプラインを使用した警告メッセージ
"ファイルが更新されていません" | Write-Warning
```

出力例:
```
WARNING: ファイルが更新されていません
```

### $WarningPreference変数の利用

PowerShellの組み込み変数`$WarningPreference`を使用して、警告メッセージの表示動作を制御できる。

```powershell
# 現在の警告表示設定を確認
$WarningPreference

# 警告を非表示にする
$WarningPreference = 'SilentlyContinue'
Write-Warning "この警告メッセージは表示されません"

# 警告をエラーとして扱う
$WarningPreference = 'Stop'
Write-Warning "この警告メッセージはエラーとして扱われ、スクリプトの実行が停止します"
# 注意: 上記のコードを実行するとスクリプトが停止する

# 警告の表示前に確認を求める
$WarningPreference = 'Inquire'
Write-Warning "この警告メッセージは確認ダイアログを表示します"

# 通常の警告表示に戻す
$WarningPreference = 'Continue'
Write-Warning "通常の警告メッセージに戻りました"
```

出力例（順に）:
```
Continue

# 2番目のWrite-Warningは何も表示されない

Write-Warning : この警告メッセージはエラーとして扱われ、スクリプトの実行が停止します
発生場所 行:1 文字:1
+ Write-Warning "この警告メッセージはエラーとして扱われ、スクリプトの実行が停止 ...
+ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    + CategoryInfo          : OperationStopped: (この警告メッセー...が停止します:String) [Write-Warning], ParentContainsErrorRecordException
    + FullyQualifiedErrorId : WarningPreferenceStop,Microsoft.PowerShell.Commands.WriteWarningCommand

# 4番目のWrite-Warningでは確認ダイアログが表示される
WARNING: この警告メッセージは確認ダイアログを表示します
[Y] はい(Y)  [A] すべて続行(A)  [H] このコマンドを停止(H)  [S] 中断(S)  [?] ヘルプ (既定値は "Y"): 

WARNING: 通常の警告メッセージに戻りました
```

### -WarningAction

個別のコマンド実行時に警告の動作を制御する。`$WarningPreference`よりも優先される。

```powershell
# 警告を非表示にする
Write-Warning "このメッセージは表示されません" -WarningAction SilentlyContinue

# 警告をエラーとして扱う
try {
    Write-Warning "この警告はエラーとして扱われます" -WarningAction Stop
} catch {
    Write-Output "警告がエラーとなり、キャッチしました"
}

# 通常どおり警告を表示する
Write-Warning "通常の警告メッセージです" -WarningAction Continue
```

出力例:
```
# 1番目のWrite-Warningは何も表示されない

警告がエラーとなり、キャッチしました

WARNING: 通常の警告メッセージです
```

### -WarningVariable

警告メッセージを変数に格納する。これにより、警告の発生を追跡し、あとから処理することができる。

```powershell
# 警告を変数に保存する
Write-Warning "この警告メッセージは変数にも保存されます" -WarningVariable warnVar
$warnVar

# 変数への追加
Write-Warning "2つ目の警告メッセージ" -WarningVariable +warnVar
$warnVar.Count
$warnVar[0]
$warnVar[1]
```

出力例:
```
WARNING: この警告メッセージは変数にも保存されます
この警告メッセージは変数にも保存されます

WARNING: 2つ目の警告メッセージ
2
この警告メッセージは変数にも保存されます
2つ目の警告メッセージ
```

## 高度な使用例

### 条件に基づいた警告表示

```powershell
# 条件に基づいて警告を表示
$cpuUsage = 85
if ($cpuUsage -gt 80) {
    Write-Warning "CPU使用率が ${cpuUsage}% と高くなっています。処理を最適化してください。"
}
```

出力例:
```
WARNING: CPU使用率が 85% と高くなっています。処理を最適化してください。
```

### 関数内での使用

```powershell
# 関数内でWrite-Warningを使用
function Test-Configuration {
    param (
        [string]$ConfigPath
    )
    
    if (-not (Test-Path $ConfigPath)) {
        Write-Warning "構成ファイル '$ConfigPath' が存在しません。デフォルト設定を使用します。"
        return $false
    }
    
    return $true
}

# 関数を呼び出す
Test-Configuration -ConfigPath "C:\NonExistentConfig.json"
```

出力例:
```
WARNING: 構成ファイル 'C:\NonExistentConfig.json' が存在しません。デフォルト設定を使用します。
False
```

### 警告メッセージのログ記録

```powershell
# 警告メッセージをログファイルに記録
function Log-Warning {
    param (
        [string]$Message
    )
    
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    Write-Warning $Message
    Add-Content -Path "C:\Temp\warnings.log" -Value "[$timestamp] WARNING: $Message"
}

# 関数を使用
Log-Warning "データベース接続がタイムアウトしました。再試行します。"
```

出力例:
```
WARNING: データベース接続がタイムアウトしました。再試行します。
# また、C:\Temp\warnings.logファイルに以下の行が追加される
# [2025-03-29 00:36:31] WARNING: データベース接続がタイムアウトしました。再試行します。
```

## 出力例の見方

`Write-Warning`コマンドレットの出力は、以下の特徴を持つ：

1. **プレフィックス**:
   - 警告メッセージは常に「WARNING: 」というプレフィックスで始まる
   - このプレフィックスは自動的に追加され、メッセージ本体には含める必要がない

2. **表示色**:
   - 警告メッセージはデフォルトで黄色で表示される
   - コンソールの色設定によっては異なる色で表示される場合がある

3. **出力ストリーム**:
   - 警告メッセージは警告ストリーム（3番）に出力される
   - 通常の出力ストリーム（1番）やエラーストリーム（2番）とは別扱いである

4. **$WarningPreferenceによる動作変更**:
   - `Continue`: 警告を表示して続行（デフォルト）
   - `SilentlyContinue`: 警告を表示せずに続行
   - `Stop`: 警告をエラーとして扱い、実行を停止
   - `Inquire`: 警告表示後、続行するかどうかをユーザーに確認

5. **WarningVariable**:
   - `-WarningVariable`で指定した変数には、プレフィックスなしのメッセージ文字列が格納される
   - 変数名の前に`+`を付けると、既存の変数に追加される

## 使用可能なPowerShellバージョン

`Write-Warning`コマンドレットは以下のPowerShellバージョンで使用可能である：

- PowerShell 1.0以降（基本機能）
- PowerShell 3.0以降（メッセージパラメータが位置パラメータとして使用可能）
- PowerShell 5.1および7.x（すべての機能が一貫して利用可能）

PowerShell 7.x（Core）では、Windows以外のプラットフォーム（Linux, macOS）でも同様に使用できる。

各バージョンでのコマンドの確認方法：

```powershell
# PowerShellバージョンの確認
$PSVersionTable.PSVersion

# コマンドのヘルプ表示
Get-Help Write-Warning -Full
```

出力例:
```
Major  Minor  Build  Revision
-----  -----  -----  --------
7      3      4      500

名前
    Write-Warning

概要
    警告メッセージを書き込みます。
...
```

**参考**: `Write-Warning`は、PowerShellの一連の出力コマンドレット（`Write-Output`, `Write-Host`, `Write-Error`, `Write-Debug`, `Write-Verbose`, `Write-Progress`）の一つである。それぞれ異なる目的と動作を持っているため、適切な状況で適切なコマンドレットを使用することが重要である。警告メッセージは、スクリプトの実行は継続するが、ユーザーに注意を促したい場合に最適である。

## 参考サイト

- [Microsoft PowerShell Documentation: Write-Warning](https://docs.microsoft.com/en-us/powershell/module/microsoft.powershell.utility/write-warning)
- [PowerShell.org - PowerShell Output Streams](https://powershell.org/2019/08/the-powershell-output-streams/)
- [技術評論社 - PowerShellポケットリファレンス](https://gihyo.jp/book/2021/978-4-297-11986-6)
- [PowerShell実践リファレンス - 警告メッセージの制御](https://xtech.nikkei.com/it/article/COLUMN/20060922/248292/)
- [GitHub PowerShell Community Documentation](https://github.com/PowerShell/PowerShell-Docs)