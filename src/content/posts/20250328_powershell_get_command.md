---
title:  【PowerShell】Get-Commandコマンドの基本と応用
published: 2025-03-28
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-28

## 概要

PowerShellの`Get-Command`コマンドレット（エイリアス：`gcm`）は、システムにインストールされているすべてのコマンドを検索し、表示するための強力なツールである。コマンドレット、関数、エイリアス、アプリケーションなどの各種コマンドタイプを識別し、特定の条件に基づいてフィルタリングすることができる。PowerShellスクリプト開発や日々の作業において、使用可能なコマンドを見つけたり、特定のコマンドの詳細情報を取得したりするのに非常に役立つコマンドである。

## 基本的な使い方

### すべてのコマンドの一覧表示

```powershell
# システム上のすべてのコマンドを一覧表示
Get-Command
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Add-Content                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Add-History                                        7.0.0.0    Microsoft.PowerShell.Core
Cmdlet          Add-Member                                         7.0.0.0    Microsoft.PowerShell.Utility
...
Function        Clear-RecycleBin                                   7.0.0.0    Microsoft.PowerShell.Management
Function        ConvertFrom-Markdown                               7.0.0.0    Microsoft.PowerShell.Utility
...
Application     cmd.exe                                            10.0.19... C:\Windows\system32\cmd.exe
Application     code.cmd                                           0.0.0.0    C:\Program Files\Microsoft VS Code\bin\code.cmd
...
```

### 特定のコマンドの検索

```powershell
# "process"という単語を含むコマンドを検索
Get-Command -Name "*process*"
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Debug-Process                                      7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-Process                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Start-Process                                      7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Stop-Process                                       7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Wait-Process                                       7.0.0.0    Microsoft.PowerShell.Management
```

### コマンドタイプによるフィルタリング

```powershell
# コマンドレット（Cmdlet）のみを表示
Get-Command -CommandType Cmdlet
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Add-Content                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Add-History                                        7.0.0.0    Microsoft.PowerShell.Core
Cmdlet          Add-Member                                         7.0.0.0    Microsoft.PowerShell.Utility
...
```

## 全オプションに対するハンズオン

### -Name

特定の名前パターンに一致するコマンドを検索する。ワイルドカードが使用可能。

```powershell
# "item"を名前に含むコマンドを検索
Get-Command -Name "*item*"
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Clear-Item                                         7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Copy-Item                                          7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-ChildItem                                      7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-Item                                           7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Invoke-Item                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Move-Item                                          7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          New-Item                                           7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Remove-Item                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Rename-Item                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Set-Item                                           7.0.0.0    Microsoft.PowerShell.Management
```

### -CommandType

特定の種類のコマンドのみをフィルタリングする。

```powershell
# 関数のみを表示
Get-Command -CommandType Function

# エイリアスのみを表示
Get-Command -CommandType Alias

# 複数のコマンドタイプを同時に表示
Get-Command -CommandType Cmdlet, Function
```

出力例（エイリアス）:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Alias           % -> ForEach-Object                                           
Alias           ? -> Where-Object                                             
Alias           ac -> Add-Content                                             
Alias           cat -> Get-Content                                            
Alias           cd -> Set-Location                                            
Alias           chdir -> Set-Location                                         
Alias           cls -> Clear-Host                                             
...
```

### -Module

特定のモジュール内のコマンドのみを表示する。

```powershell
# Microsoft.PowerShell.Securityモジュールのコマンドを表示
Get-Command -Module Microsoft.PowerShell.Security
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          ConvertFrom-SecureString                           7.0.0.0    Microsoft.PowerShell.Security
Cmdlet          ConvertTo-SecureString                             7.0.0.0    Microsoft.PowerShell.Security
Cmdlet          Get-Acl                                            7.0.0.0    Microsoft.PowerShell.Security
Cmdlet          Get-AuthenticodeSignature                          7.0.0.0    Microsoft.PowerShell.Security
Cmdlet          Get-CmsMessage                                     7.0.0.0    Microsoft.PowerShell.Security
Cmdlet          Get-Credential                                     7.0.0.0    Microsoft.PowerShell.Security
Cmdlet          Get-ExecutionPolicy                                7.0.0.0    Microsoft.PowerShell.Security
...
```

### -Verb

特定の動詞（操作）を持つコマンドのみを表示する。

```powershell
# "Get"動詞を持つコマンドを表示
Get-Command -Verb Get
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Get-Acl                                            7.0.0.0    Microsoft.PowerShell.Security
Cmdlet          Get-Alias                                          7.0.0.0    Microsoft.PowerShell.Utility
Cmdlet          Get-AuthenticodeSignature                          7.0.0.0    Microsoft.PowerShell.Security
Cmdlet          Get-ChildItem                                      7.0.0.0    Microsoft.PowerShell.Management
...
```

### -Noun

特定の名詞（対象）を持つコマンドのみを表示する。

```powershell
# "Service"名詞を持つコマンドを表示
Get-Command -Noun Service
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Get-Service                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          New-Service                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Restart-Service                                    7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Resume-Service                                     7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Set-Service                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Start-Service                                      7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Stop-Service                                       7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Suspend-Service                                    7.0.0.0    Microsoft.PowerShell.Management
```

### -ParameterName

特定のパラメータを持つコマンドを検索する。

```powershell
# "ComputerName"パラメータを持つコマンドを検索
Get-Command -ParameterName ComputerName
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Connect-PSSession                                  7.0.0.0    Microsoft.PowerShell.Core
Cmdlet          Enter-PSSession                                    7.0.0.0    Microsoft.PowerShell.Core
Cmdlet          Get-Process                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-Service                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-WinEvent                                       7.0.0.0    Microsoft.PowerShell.Diagnostics
Cmdlet          Invoke-Command                                     7.0.0.0    Microsoft.PowerShell.Core
...
```

### -ParameterType

特定の型のパラメータを持つコマンドを検索する。

```powershell
# System.IO.FileInfo型のパラメータを持つコマンドを検索
Get-Command -ParameterType System.IO.FileInfo
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Add-Content                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Clear-Content                                      7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-Content                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-FileHash                                       7.0.0.0    Microsoft.PowerShell.Utility
Cmdlet          Set-Content                                        7.0.0.0    Microsoft.PowerShell.Management
...
```

### -TotalCount

表示するコマンドの最大数を指定する。

```powershell
# 最初の5つのコマンドのみ表示
Get-Command -TotalCount 5
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Add-Content                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Add-History                                        7.0.0.0    Microsoft.PowerShell.Core
Cmdlet          Add-Member                                         7.0.0.0    Microsoft.PowerShell.Utility
Cmdlet          Add-PSSnapin                                       7.0.0.0    Microsoft.PowerShell.Core
Cmdlet          Add-Type                                           7.0.0.0    Microsoft.PowerShell.Utility
```

### -ListImported

読み込まれているモジュールからのコマンドのみを表示する。

```powershell
# 現在のセッションに読み込まれているコマンドのみ表示
Get-Command -ListImported
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Add-Content                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Add-History                                        7.0.0.0    Microsoft.PowerShell.Core
...
Function        Clear-RecycleBin                                   7.0.0.0    Microsoft.PowerShell.Management
...
```

### -ShowCommandInfo

コマンドの情報オブジェクトを表示する。

```powershell
# Get-Processコマンドの詳細情報を表示
Get-Command -Name Get-Process -ShowCommandInfo
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Get-Process                                        7.0.0.0    Microsoft.PowerShell.Management

Syntax
   Get-Process [[-Name] <System.String[]>] [-ComputerName <System.String[]>] [-FileVersionInfo] [-Module] [<CommonParameters>]
   
   Get-Process [-ComputerName <System.String[]>] [-FileVersionInfo] [-Module] -Id <System.Int32[]> [<CommonParameters>]
   
   Get-Process [-ComputerName <System.String[]>] [-FileVersionInfo] [-Module] -InputObject <System.Diagnostics.Process[]> [<CommonParameters>]
...
```

### -All

すべてのコマンドを表示する（通常は隠れたコマンドも含む）。

```powershell
# すべてのコマンド（隠れたものも含む）を表示
Get-Command -All
```

出力例:
```
# 通常のGet-Commandよりもさらに多くのコマンドが表示される
```

### -Syntax

コマンドの構文を表示する。

```powershell
# Get-Contentコマンドの構文を表示
Get-Command -Name Get-Content -Syntax
```

出力例:
```
Get-Content [-Path] <System.String[]> [-ReadCount <System.Int64>] [-TotalCount <System.Int64>] [-Tail <System.Int32>] [-Encoding {ascii | bigendianunicode | bigendianutf32 | byte | default | oem | string | unicode | unknown | utf7 | utf8 | utf8BOM | utf8NoBOM | utf32}] [-Delimiter <System.String>] [-Wait] [-Raw] [-Force] [-Stream <System.String>] [<CommonParameters>]

Get-Content [-ReadCount <System.Int64>] [-TotalCount <System.Int64>] [-Tail <System.Int32>] [-Encoding {ascii | bigendianunicode | bigendianutf32 | byte | default | oem | string | unicode | unknown | utf7 | utf8 | utf8BOM | utf8NoBOM | utf32}] [-Delimiter <System.String>] [-Wait] [-Raw] [-Force] -LiteralPath <System.String[]> [-Stream <System.String>] [<CommonParameters>]
```

## 高度な使用例

### モジュールとコマンドタイプの組み合わせ

```powershell
# Microsoft.PowerShell.Utilityモジュール内の関数のみを検索
Get-Command -Module Microsoft.PowerShell.Utility -CommandType Function
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Function        ConvertFrom-Markdown                               7.0.0.0    Microsoft.PowerShell.Utility
Function        ConvertFrom-ScriptExtent                           7.0.0.0    Microsoft.PowerShell.Utility
Function        ConvertTo-Markdown                                 7.0.0.0    Microsoft.PowerShell.Utility
...
```

### 特定のパターンに一致するパラメータを持つコマンド

```powershell
# "Path"という名前を含むパラメータを持つ"Get"動詞のコマンドを検索
Get-Command -Verb Get -ParameterName *Path*
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Get-ChildItem                                      7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-Content                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-Item                                           7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Get-ItemProperty                                   7.0.0.0    Microsoft.PowerShell.Management
...
```

### コマンドソースによるフィルタリング

```powershell
# 特定のモジュール内でさらに詳細な条件でコマンドを検索
Get-Command -Name "*service*" -Module Microsoft.PowerShell.Management
```

出力例:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Get-Service                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          New-Service                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Restart-Service                                    7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Resume-Service                                     7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Set-Service                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Start-Service                                      7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Stop-Service                                       7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Suspend-Service                                    7.0.0.0    Microsoft.PowerShell.Management
```

### コマンドのエクスポート

```powershell
# 検索結果をCSVファイルにエクスポートする
Get-Command -CommandType Cmdlet | Export-Csv -Path "C:\Temp\PowerShell_Cmdlets.csv"
```

出力例:
```
# 出力なし（ファイルに書き込まれる）
```

### コマンド情報の詳細表示

```powershell
# 特定のコマンドの詳細情報を表示
Get-Command -Name Get-Process | Format-List *
```

出力例:
```
HelpUri                       : https://go.microsoft.com/fwlink/?LinkID=2096820
CommandType                   : Cmdlet
Name                          : Get-Process
ModuleName                    : Microsoft.PowerShell.Management
Module                        : Microsoft.PowerShell.Management
RemotingCapability            : SupportedByCommand
Parameters                    : {[Id, System.Management.Automation.ParameterMetadata], [InputObject, System.Management.Automation.ParameterMetadata], [Name, System.Management.Automation.ParameterMetadata], [ComputerName, System.Management.Automation.ParameterMetadata]...}
ParameterSets                 : {[[-Name] <String[]>] [-ComputerName <String[]>] [-FileVersionInfo] [-Module] [<CommonParameters>], [-ComputerName <String[]>] [-FileVersionInfo] [-Module] -Id <Int32[]> [<CommonParameters>], [-ComputerName <String[]>] [-FileVersionInfo] [-Module] -InputObject <Process[]> [<CommonParameters>]}
Definition                    : Get-Process [[-Name] <String[]>] [-ComputerName <String[]>] [-FileVersionInfo] [-Module] [<CommonParameters>]
                                Get-Process [-ComputerName <String[]>] [-FileVersionInfo] [-Module] -Id <Int32[]> [<CommonParameters>]
                                Get-Process [-ComputerName <String[]>] [-FileVersionInfo] [-Module] -InputObject <Process[]> [<CommonParameters>]
...
```

## 出力例の見方

`Get-Command`コマンドレットの標準出力は表形式で、以下の情報が含まれる：

1. **CommandType**: コマンドの種類を示す
   - Cmdlet: PowerShellのネイティブコマンド
   - Function: PowerShellスクリプトで定義された関数
   - Alias: コマンドの別名
   - Application: 外部アプリケーション（実行可能ファイル）
   - Script: PowerShellスクリプト
   - ExternalScript: .ps1ファイルとして保存されたスクリプト
   - Filter: 特殊な種類の関数
   - Configuration: DSC構成

2. **Name**: コマンドの名前
   - エイリアスの場合は、「エイリアス名 -> 実際のコマンド」という形式

3. **Version**: コマンドが属するモジュールのバージョン
   - 例: 7.0.0.0 は PowerShell 7.0のモジュールを示す

4. **Source**: コマンドが含まれるモジュールや場所
   - PowerShellコマンドレットの場合は、Microsoft.PowerShell.XXX などのモジュール名
   - アプリケーションの場合は、実行可能ファイルのパス

`-Syntax`パラメータを使用した場合は、コマンドの使用法が表示される：

- 角括弧 `[]` で囲まれたパラメータはオプション
- 波括弧 `{}` はパラメータの取り得る値の選択肢
- `<>` で囲まれた部分はパラメータの型
- `|` は選択肢を示す
- `[<CommonParameters>]` は共通パラメータを示す（-Verbose, -Debug, -ErrorAction など）

`Format-List *`と組み合わせた場合は、より詳細な情報が表示される：

- HelpUri: ヘルプドキュメントのURI
- Definition: コマンドの定義（構文）
- ParameterSets: パラメータの組み合わせパターン
- Parameters: サポートされるパラメータの詳細情報



**注意点**: `Get-Command`の結果は、現在のPowerShellセッションで利用できるコマンドに基づいており、まだ読み込まれていないモジュール内のコマンドは表示されない場合がある。また、Windowsの場合とLinuxの場合で利用可能なコマンドセットが若干異なる点にも注意が必要である。