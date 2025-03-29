---
title:  【PowerShell】Write-Verboseコマンドの基本と応用
published: 2025-03-29
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-29

## 概要

PowerShellの`Write-Verbose`コマンドレットは、詳細な情報メッセージ（冗長メッセージ）をユーザーに表示するための機能である。このコマンドレットは主にスクリプトやコマンドの内部動作の詳細を表示し、デバッグやトラブルシューティング時に役立つ。通常、このメッセージはデフォルトでは表示されず、`-Verbose`パラメータが指定された場合や、`$VerbosePreference`変数が調整されている場合にのみ表示される。これにより、スクリプトやコマンドレットの実行過程を詳細に追跡したい場合と、通常の実行時には情報を抑制したい場合の両方に対応できる柔軟性を提供している。

自作のpowershell関数を作成した際に警告を表示するときに役に立つ。


## 基本的な使い方

### 単純な詳細メッセージの定義

```powershell
# 基本的な詳細メッセージの定義（このままでは表示されない）
Write-Verbose "設定ファイルを読み込んでいます..."
```

出力例:
```
# デフォルトでは何も表示されない
```

### -Verboseパラメータを使用して表示

```powershell
# -Verboseパラメータを付けて実行することで表示される
Write-Verbose "設定ファイルを読み込んでいます..." -Verbose
```

出力例:
```
VERBOSE: 設定ファイルを読み込んでいます...
```

### $VerbosePreference変数を使用した制御

```powershell
# 現在の詳細メッセージ表示設定を確認
$VerbosePreference

# 詳細メッセージを常に表示するよう設定
$VerbosePreference = 'Continue'
Write-Verbose "この詳細メッセージは自動的に表示されます"

# 詳細メッセージ表示設定を元に戻す
$VerbosePreference = 'SilentlyContinue'
Write-Verbose "この詳細メッセージは表示されません"
```

出力例:
```
SilentlyContinue

VERBOSE: この詳細メッセージは自動的に表示されます

# 最後の行は表示されない
```

## 全オプションに対するハンズオン

### -Message

表示する詳細メッセージを指定する。これは`Write-Verbose`の主要なパラメータである。

```powershell
# 明示的にMessageパラメータを使用
Write-Verbose -Message "データベース接続を初期化しています..." -Verbose
```

出力例:
```
VERBOSE: データベース接続を初期化しています...
```

### -Verbose

詳細メッセージを表示するよう指示するための共通パラメータ。これは`Write-Verbose`だけでなく、多くのPowerShellコマンドレットで使用できる。

```powershell
# -Verboseパラメータを使用して詳細メッセージを表示
Write-Verbose "ファイルをスキャンしています..." -Verbose
```

出力例:
```
VERBOSE: ファイルをスキャンしています...
```

### -VerboseVariable

詳細メッセージを変数に格納する。これにより、メッセージを画面に表示しながら同時に変数にも保存できる。

```powershell
# 詳細メッセージを変数に保存する
Write-Verbose "システム情報を収集しています..." -Verbose -VerboseVariable verbLog
$verbLog

# 変数への追加
Write-Verbose "ネットワーク設定を確認しています..." -Verbose -VerboseVariable +verbLog
$verbLog.Count
$verbLog[0]
$verbLog[1]
```

出力例:
```
VERBOSE: システム情報を収集しています...
システム情報を収集しています...

VERBOSE: ネットワーク設定を確認しています...
2
システム情報を収集しています...
ネットワーク設定を確認しています...
```

### 関数内での使用（CmdletBinding）

`CmdletBinding`属性を持つ関数内で使用すると、関数呼び出し時に`-Verbose`パラメータを指定できるようになる。

```powershell
# CmdletBinding属性を持つ関数の定義
function Test-Connection {
    [CmdletBinding()]
    param (
        [string]$Server
    )
    
    Write-Verbose "サーバー '$Server' への接続をテストしています..."
    Write-Verbose "パケットを送信しています..."
    Write-Verbose "応答を待機しています..."
    
    return $true
}

# 関数を通常実行（詳細メッセージは表示されない）
Test-Connection -Server "example.com"

# 関数を詳細モードで実行
Test-Connection -Server "example.com" -Verbose
```

出力例:
```
True

VERBOSE: サーバー 'example.com' への接続をテストしています...
VERBOSE: パケットを送信しています...
VERBOSE: 応答を待機しています...
True
```

## 高度な使用例

### 詳細レベルに応じたメッセージ表示

```powershell
function Process-Files {
    [CmdletBinding()]
    param (
        [string]$Path,
        [switch]$DeepAnalysis
    )
    
    Write-Verbose "処理を開始しています..."
    
    # 基本的な処理の詳細
    Get-ChildItem -Path $Path -File | ForEach-Object {
        Write-Verbose "ファイル '$($_.Name)' を処理しています..."
        # 実際の処理...
    }
    
    # より詳細な処理（特定の条件下でのみ実行）
    if ($DeepAnalysis) {
        $VerbosePreferenceBak = $VerbosePreference
        $VerbosePreference = 'Continue'
        
        Write-Verbose "詳細分析モードを有効化しています..."
        # 詳細な分析処理...
        
        $VerbosePreference = $VerbosePreferenceBak
    }
    
    Write-Verbose "処理が完了しました"
}

# 通常モードで実行
Process-Files -Path "C:\Temp" -Verbose

# 詳細分析モードで実行
Process-Files -Path "C:\Temp" -DeepAnalysis -Verbose
```

出力例:
```
VERBOSE: 処理を開始しています...
VERBOSE: ファイル 'test.txt' を処理しています...
VERBOSE: ファイル 'data.csv' を処理しています...
VERBOSE: 処理が完了しました

VERBOSE: 処理を開始しています...
VERBOSE: ファイル 'test.txt' を処理しています...
VERBOSE: ファイル 'data.csv' を処理しています...
VERBOSE: 詳細分析モードを有効化しています...
VERBOSE: 処理が完了しました
```

### 進行状況の報告

```powershell
function Copy-LargeFiles {
    [CmdletBinding()]
    param (
        [string]$Source,
        [string]$Destination
    )
    
    Write-Verbose "ソース '$Source' からファイルを収集しています..."
    $files = Get-ChildItem -Path $Source -File
    Write-Verbose "$($files.Count) ファイルが見つかりました"
    
    foreach ($file in $files) {
        Write-Verbose "ファイル '$($file.Name)' を '$Destination' にコピーしています... ($([math]::Round($file.Length / 1MB, 2)) MB)"
        # 実際のコピー処理...
        Start-Sleep -Milliseconds 500 # 処理をシミュレート
    }
    
    Write-Verbose "すべてのファイルのコピーが完了しました"
}

# 詳細モードで実行
Copy-LargeFiles -Source "C:\Temp" -Destination "D:\Backup" -Verbose
```

出力例:
```
VERBOSE: ソース 'C:\Temp' からファイルを収集しています...
VERBOSE: 3 ファイルが見つかりました
VERBOSE: ファイル 'largefile.zip' を 'D:\Backup' にコピーしています... (25.46 MB)
VERBOSE: ファイル 'document.pdf' を 'D:\Backup' にコピーしています... (2.34 MB)
VERBOSE: ファイル 'image.jpg' を 'D:\Backup' にコピーしています... (5.78 MB)
VERBOSE: すべてのファイルのコピーが完了しました
```

### 詳細メッセージのログ記録

```powershell
function Invoke-ProcessWithLogging {
    [CmdletBinding()]
    param (
        [string]$Task,
        [string]$LogPath = "C:\Logs\process_log.txt"
    )
    
    # 詳細メッセージを記録するための関数
    function Log-Verbose {
        param ([string]$Message)
        
        $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
        $logMessage = "[$timestamp] VERBOSE: $Message"
        
        # 画面に表示
        Write-Verbose $Message
        
        # ログファイルに書き込み
        Add-Content -Path $LogPath -Value $logMessage
    }
    
    Log-Verbose "タスク '$Task' の処理を開始します"
    # 実際の処理...
    Start-Sleep -Seconds 1 # 処理をシミュレート
    Log-Verbose "タスクが完了しました"
}

# 詳細モードで実行
Invoke-ProcessWithLogging -Task "データバックアップ" -Verbose
```

出力例:
```
VERBOSE: タスク 'データバックアップ' の処理を開始します
VERBOSE: タスクが完了しました

# また、C:\Logs\process_log.txtファイルには以下のような行が記録される：
# [2025-03-29 00:38:29] VERBOSE: タスク 'データバックアップ' の処理を開始します
# [2025-03-29 00:38:30] VERBOSE: タスクが完了しました
```

## 出力例の見方

`Write-Verbose`コマンドレットの出力は、以下の特徴を持つ：

1. **プレフィックス**:
   - 詳細メッセージは常に「VERBOSE: 」というプレフィックスで始まる
   - このプレフィックスは自動的に追加され、メッセージ本体には含める必要がない

2. **表示色**:
   - 詳細メッセージはデフォルトで青緑色（シアン）で表示される
   - コンソールの色設定によっては異なる色で表示される場合がある

3. **表示制御**:
   - デフォルトでは表示されない
   - `-Verbose`パラメータを指定するか、`$VerbosePreference`を「Continue」に設定することで表示される

4. **出力ストリーム**:
   - 詳細メッセージは詳細ストリーム（4番）に出力される
   - 通常の出力ストリーム（1番）やエラーストリーム（2番）とは別扱いである

5. **$VerbosePreferenceによる動作変更**:
   - `SilentlyContinue`: 詳細メッセージを表示せずに続行（デフォルト）
   - `Continue`: 詳細メッセージを表示して続行
   - `Stop`: 詳細メッセージを表示し、エラーとして扱い、実行を停止
   - `Inquire`: 詳細メッセージを表示後、続行するかどうかをユーザーに確認

6. **-VerboseVariable**:
   - 指定した変数にはプレフィックスなしのメッセージのみが格納される
   - 変数名の前に`+`を付けると、既存の変数に追加する形で保存される

## 使用可能なPowerShellバージョン

`Write-Verbose`コマンドレットは以下のPowerShellバージョンで使用可能である：

- PowerShell 1.0以降（基本機能）
- PowerShell 3.0以降（メッセージパラメータが位置パラメータとして使用可能）
- PowerShell 5.1および7.x（すべての機能が一貫して利用可能）

PowerShell 7.x（Core）では、Windows以外のプラットフォーム（Linux, macOS）でも同様に使用できる。

各バージョンでのコマンドの確認方法：

```powershell
# PowerShellバージョンの確認
$PSVersionTable.PSVersion

# コマンドのヘルプ表示
Get-Help Write-Verbose -Full
```

出力例:
```
Major  Minor  Build  Revision
-----  -----  -----  --------
7      3      4      500

名前
    Write-Verbose

概要
    詳細出力用のストリームにメッセージを書き込みます。
...
```

**補足**: `Write-Verbose`コマンドレットは、適切なレベルの詳細情報を提供するためのものである。再利用可能なスクリプトやモジュール、特に他の人が使用する可能性のあるものを作成する場合は、`Write-Verbose`を活用して内部動作の透明性を確保することが推奨される。ユーザーは必要に応じて`-Verbose`スイッチを使用して詳細情報を表示できるため、通常の使用を妨げることなく、デバッグやトラブルシューティングの際に役立つ情報を提供できる。

## 参考サイト

- [Microsoft PowerShell Documentation: Write-Verbose](https://docs.microsoft.com/en-us/powershell/module/microsoft.powershell.utility/write-verbose)
- [PowerShell.org - PowerShell Output Streams](https://powershell.org/2019/08/the-powershell-output-streams/)
- [Windows PowerShell クックブック](https://www.amazon.co.jp/dp/4873113822/)
- [PowerShell実践ガイド](https://gihyo.jp/book/2017/978-4-7741-8855-3)
- [GitHub PowerShell Community Documentation](https://github.com/PowerShell/PowerShell-Docs)
- [Microsoft TechNet - Advanced Function Development](https://technet.microsoft.com/ja-jp/library/hh847743.aspx)