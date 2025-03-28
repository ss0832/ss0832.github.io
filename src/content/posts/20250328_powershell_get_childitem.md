---
title:  【PowerShell】Get-ChildItemコマンドの基本と応用
published: 2025-03-28
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-28

## 概要

PowerShellの`Get-ChildItem`コマンドレット（エイリアス：`dir`, `ls`, `gci`）は、指定した場所に含まれるアイテム（ファイル、フォルダ、レジストリのキーなど）を取得して表示するためのコマンドである。Windowsのコマンドプロンプトにおける`dir`コマンドやUNIX系OSの`ls`コマンドに相当する基本的な機能を持ちながら、より高度なフィルタリングやパイプラインとの連携が可能である。ファイルシステム以外にもレジストリやPowerShellドライブなど、様々なプロバイダー上で使用できる汎用性の高いコマンドレットである。

## 基本的な使い方

### 現在のディレクトリの内容を表示する

```powershell
# 最も基本的な使い方
Get-ChildItem
```

出力例:
```
    ディレクトリ: C:\Users\ss0832\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/15    14:22                PowerShell
d----          2025/03/10    09:15                Projects
-a---          2025/03/20    18:45           2391 memo.txt
-a---          2025/03/25    12:30          15420 report.docx
```

### 特定のパスの内容を表示する

```powershell
# 指定したパスの内容を取得
Get-ChildItem -Path C:\Windows\System32
```

出力例:
```
    ディレクトリ: C:\Windows\System32

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2024/12/25    10:15                drivers
d----          2025/01/15    08:30                wbem
-a---          2024/11/30    15:20         123456 ntdll.dll
-a---          2025/02/10    09:45          56789 kernel32.dll
```

## 全オプションに対するハンズオン

### -Path

指定したパスの内容を表示する。複数のパスを指定することも可能。

```powershell
# 複数のパスを指定
Get-ChildItem -Path C:\Windows\System32\drivers, C:\Windows\Fonts
```

出力例:
```
    ディレクトリ: C:\Windows\System32\drivers

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/01/15    13:40          45678 disk.sys
-a---          2025/02/05    09:22          23456 ntfs.sys

    ディレクトリ: C:\Windows\Fonts

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2024/11/20    14:30         189543 arial.ttf
-a---          2024/11/20    14:30         205678 times.ttf
```

### -Filter

指定したパターンに一致するアイテムのみを表示する。ワイルドカード（*）を使用可能。

```powershell
# .txtファイルのみを表示
Get-ChildItem -Filter *.txt
```

出力例:
```
    ディレクトリ: C:\Users\ss0832\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/20    18:45           2391 memo.txt
-a---          2025/03/15    10:22           5432 notes.txt
```
### -Recurse

指定したパスとそのすべてのサブディレクトリの内容を再帰的に表示する。

```powershell
# サブフォルダも含めて検索
Get-ChildItem -Path C:\Scripts -Recurse
```

出力例:
```
    ディレクトリ: C:\Scripts

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/15    14:22                Utils
-a---          2025/03/20    18:45           2391 main.ps1

    ディレクトリ: C:\Scripts\Utils

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/18    10:30           1234 helper.ps1
-a---          2025/03/19    14:15           5678 functions.ps1
```


### -Include

指定したパターンに一致するアイテムだけを含める。複数のパターンを指定可能。-Recurseと組み合わせるとより効果的。

```powershell
# .txt と .log ファイルのみを取得
Get-ChildItem -Path C:\Logs -Include *.txt,*.log -Recurse
```

出力例:
```
    ディレクトリ: C:\Logs

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/26    14:22          45678 system.log
-a---          2025/03/27    09:15          12345 error.log

    ディレクトリ: C:\Logs\Archive

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/25    18:30           3456 archive.txt
-a---          2025/03/20    11:45          78901 old.log
```

### -Exclude

指定したパターンに一致するアイテムを除外する。複数のパターンを指定可能。

```powershell
# .tmpファイルを除外して表示
Get-ChildItem -Exclude *.tmp
```

出力例:
```
    ディレクトリ: C:\Users\ss0832\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/15    14:22                PowerShell
-a---          2025/03/20    18:45           2391 memo.txt
-a---          2025/03/25    12:30          15420 report.docx
```


### -Depth

-Recurseと併用し、再帰的に探索する深さを制限する。

```powershell
# 最大2階層まで探索
Get-ChildItem -Path C:\Projects -Recurse -Depth 2
```

出力例:
```
    ディレクトリ: C:\Projects

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/20    10:00                ProjectA
d----          2025/03/22    14:30                ProjectB

    ディレクトリ: C:\Projects\ProjectA

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/20    10:15                src
-a---          2025/03/20    10:10            345 README.md

    ディレクトリ: C:\Projects\ProjectB

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/22    14:35                src
-a---          2025/03/22    14:40            456 README.md
```

### -Directory

ディレクトリのみを表示する。

```powershell
# フォルダのみを取得
Get-ChildItem -Directory
```

出力例:
```
    ディレクトリ: C:\Users\ss0832\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/15    14:22                PowerShell
d----          2025/03/10    09:15                Projects
```

### -File

ファイルのみを表示する。

```powershell
# ファイルのみを取得
Get-ChildItem -File
```

出力例:
```
    ディレクトリ: C:\Users\ss0832\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/20    18:45           2391 memo.txt
-a---          2025/03/25    12:30          15420 report.docx
```

### -Hidden

隠しファイルとフォルダを表示する。

```powershell
# 隠しファイルを表示
Get-ChildItem -Hidden
```

出力例:
```
    ディレクトリ: C:\Users\ss0832\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d-h--          2025/03/15    15:30                .git
-ah--          2025/03/18    11:25             45 .gitignore
```

### -System

システムファイルを表示する。

```powershell
# システムファイルを表示
Get-ChildItem -System
```

出力例:
```
    ディレクトリ: C:\Windows

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d--s-          2025/01/15    14:35                System32
-a-s-          2025/02/10    09:20          12345 system.dat
```

### -Force

通常表示されない隠しファイルやシステムファイルをすべて表示する。

```powershell
# すべての隠しファイル、システムファイルを含める
Get-ChildItem -Force
```

出力例:
```
    ディレクトリ: C:\Users\ss0832\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/15    14:22                PowerShell
d----          2025/03/10    09:15                Projects
d-h--          2025/03/15    15:30                .git
-a---          2025/03/20    18:45           2391 memo.txt
-a---          2025/03/25    12:30          15420 report.docx
-ah--          2025/03/18    11:25             45 .gitignore
```

### -Name

名前のみを表示する（その他のプロパティは表示しない）。

```powershell
# 名前のみ表示
Get-ChildItem -Name
```

出力例:
```
PowerShell
Projects
memo.txt
report.docx
```

### -ReadOnly

読み取り専用ファイルのみを表示する。

```powershell
# 読み取り専用ファイルを取得
Get-ChildItem -ReadOnly
```

出力例:
```
    ディレクトリ: C:\Users\ss0832\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-ar--          2025/03/15    10:40          23456 template.docx
```

### -LiteralPath

特殊文字を含むパスをそのままのパスとして扱う（ワイルドカード展開を行わない）。

```powershell
# 特殊文字を含むファイル名を指定
Get-ChildItem -LiteralPath 'C:\Reports\Report[2025].xlsx'
```

出力例:
```
    ディレクトリ: C:\Reports

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/25    11:30          34567 Report[2025].xlsx
```

### -Attributes

指定した属性を持つアイテムのみを表示する。

```powershell
# 隠しファイルと読み取り専用属性を持つファイル
Get-ChildItem -Attributes Hidden, ReadOnly
```

出力例:
```
    ディレクトリ: C:\Users\ss0832\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-ahr-          2025/03/20    16:20          12345 config.ini
```

### 複数オプションの組み合わせ

複数のオプションを組み合わせて高度なフィルタリングを実現できる。

```powershell
# 再帰的に.logファイルを探し、"Error"という文字列を含む行を表示
Get-ChildItem -Path C:\Logs -Include *.log -Recurse | Select-String "Error"
```

出力例:
```
C:\Logs\system.log:15:Error: Failed to connect to database
C:\Logs\error.log:3:Error: Invalid configuration file
C:\Logs\Archive\old.log:120:Error: Service stopped unexpectedly
```

## 出力例の見方

`Get-ChildItem`コマンドの標準出力には以下の情報が含まれる：

1. **Mode**：アイテムの属性を示す記号
   - `d`：ディレクトリ（フォルダ）
   - `a`：アーカイブ（通常のファイル）
   - `h`：隠しファイル/フォルダ
   - `s`：システムファイル/フォルダ
   - `r`：読み取り専用
   - `l`：シンボリックリンク
   - `-`：属性なし（該当位置の属性がない）

2. **LastWriteTime**：最終更新日時

3. **Length**：ファイルサイズ（バイト単位）
   - ディレクトリの場合は表示されない

4. **Name**：ファイルまたはディレクトリの名前

より詳細な情報を取得するには、以下のようにFormat-Listコマンドレットとパイプラインを使用する：

```powershell
Get-ChildItem | Format-List *
```

出力例:
```
PSPath            : Microsoft.PowerShell.Core\FileSystem::C:\Users\ss0832\Documents\memo.txt
PSParentPath      : Microsoft.PowerShell.Core\FileSystem::C:\Users\ss0832\Documents
PSChildName       : memo.txt
PSDrive           : C
PSProvider        : Microsoft.PowerShell.Core\FileSystem
PSIsContainer     : False
Mode              : -a---
VersionInfo       : File:             C:\Users\ss0832\Documents\memo.txt
                    InternalName:     
                    OriginalFilename: 
                    FileVersion:      
                    FileDescription:  
                    Product:          
                    ProductVersion:   
                    Debug:            False
                    Patched:          False
                    PreRelease:       False
                    PrivateBuild:     False
                    SpecialBuild:     False
BaseName          : memo
Target            : {}
LinkType          : 
Name              : memo.txt
Length            : 2391
DirectoryName     : C:\Users\ss0832\Documents
Directory         : C:\Users\ss0832\Documents
IsReadOnly        : False
Exists            : True
FullName          : C:\Users\ss0832\Documents\memo.txt
Extension         : .txt
CreationTime      : 2025/03/15 08:30:45
CreationTimeUtc   : 2025/03/14 23:30:45
LastAccessTime    : 2025/03/26 16:40:22
LastAccessTimeUtc : 2025/03/26 07:40:22
LastWriteTime     : 2025/03/20 18:45:30
LastWriteTimeUtc  : 2025/03/20 09:45:30
Attributes        : Archive
```

