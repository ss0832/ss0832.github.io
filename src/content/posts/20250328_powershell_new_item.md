---
title:  【PowerShell】New-Itemコマンドの基本と応用
published: 2025-03-28
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-28

## 概要

PowerShellの`New-Item`コマンドレット（エイリアス：`ni`）は、新しいアイテムを作成するためのコマンドである。ファイルやディレクトリの作成だけでなく、レジストリキー、シンボリックリンク、ハードリンク、さらには各種PowerShellプロバイダーで使用できる様々なアイテムの作成に対応している。このコマンドレットを使いこなすことで、スクリプト内でのファイル操作や環境設定の自動化が可能になる。

## 基本的な使い方

### 新しいファイルの作成

```powershell
# 基本的なファイル作成
New-Item -Path "C:\Temp\test.txt" -ItemType File
```

出力例:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:15              0 test.txt
```

### 新しいディレクトリの作成

```powershell
# ディレクトリの作成
New-Item -Path "C:\Temp\NewFolder" -ItemType Directory
```

出力例:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/28    12:16                NewFolder
```

### コンテンツを含むファイルの作成

```powershell
# 内容を指定してファイルを作成
New-Item -Path "C:\Temp\log.txt" -ItemType File -Value "This is a log file."
```

出力例:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:17             18 log.txt
```

## 全オプションに対するハンズオン

### -Path

アイテムを作成する場所のパスを指定する。

```powershell
# パスを指定してファイル作成
New-Item -Path "C:\Temp\Documents\report.docx" -ItemType File
```

出力例:
```
    ディレクトリ: C:\Temp\Documents

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:18              0 report.docx
```

### -Name

アイテムの名前を指定する。Pathと併用する場合は、Pathで親ディレクトリを指定する。

```powershell
# 名前を指定してディレクトリ作成
New-Item -Path "C:\Temp" -Name "Scripts" -ItemType Directory
```

出力例:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/28    12:19                Scripts
```

### -ItemType

作成するアイテムの種類を指定する。主な種類は以下の通り。

```powershell
# ファイル作成
New-Item -Path "C:\Temp\config.ini" -ItemType File

# ディレクトリ作成
New-Item -Path "C:\Temp\Data" -ItemType Directory

# シンボリックリンク作成（管理者権限が必要）
New-Item -Path "C:\Temp\link" -ItemType SymbolicLink -Value "C:\Windows"

# ジャンクション作成
New-Item -Path "C:\Temp\junction" -ItemType Junction -Value "C:\Windows"

# レジストリキー作成
New-Item -Path "HKCU:\Software\MyApp" -ItemType Key

# ハードリンク作成
New-Item -Path "C:\Temp\hardlink.txt" -ItemType HardLink -Value "C:\Temp\original.txt"
```

出力例（ファイル）:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:20              0 config.ini
```

出力例（シンボリックリンク）:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
l----          2025/03/28    12:21       <SYMLINK>    link [C:\Windows]
```

### -Value

アイテムに設定する値やコンテンツを指定する。ファイルの場合はその内容、リンクの場合はリンク先となる。

```powershell
# 内容を指定してファイル作成
New-Item -Path "C:\Temp\readme.md" -ItemType File -Value "# README\nThis is a markdown file."
```

出力例:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:22             32 readme.md
```

### -Force

既存のアイテムを上書きするか、存在しない親ディレクトリを作成する。

```powershell
# 既存のファイルを上書き
New-Item -Path "C:\Temp\log.txt" -ItemType File -Value "Overwritten content." -Force

# 存在しない親ディレクトリも含めて作成
New-Item -Path "C:\Temp\SubDir1\SubDir2\file.txt" -ItemType File -Force
```

出力例（上書き）:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:23             20 log.txt
```

出力例（親ディレクトリ自動作成）:
```
    ディレクトリ: C:\Temp\SubDir1\SubDir2

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:24              0 file.txt
```

### -Credential

別のユーザーとして操作するための認証情報を指定する。

```powershell
# 別のユーザーとしてファイル作成
$cred = Get-Credential
New-Item -Path "\\Server\Share\file.txt" -ItemType File -Credential $cred
```

出力例:
```
# 認証情報入力ダイアログが表示される
# 認証後、以下のようにファイルが作成される

    ディレクトリ: \\Server\Share

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:25              0 file.txt
```

### -WhatIf

コマンドを実際に実行せず、実行された場合の動作を表示する。

```powershell
# 実行されるアクションを確認
New-Item -Path "C:\Temp\test-whatif.txt" -ItemType File -WhatIf
```

出力例:
```
WhatIf: C:\Temp\test-whatif.txt に対してアイテムを作成します
```

### -Confirm

コマンド実行前に確認を求める。

```powershell
# 実行前に確認を求める
New-Item -Path "C:\Temp\test-confirm.txt" -ItemType File -Confirm
```

出力例:
```
New-Item を続行しますか?
C:\Temp\test-confirm.txt に対してアイテムを作成します
[Y] はい(Y)  [A] すべて続行(A)  [N] いいえ(N)  [L] すべて無視(L)  [S] 中断(S)  [?] ヘルプ (既定値は "Y"): y

    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:26              0 test-confirm.txt
```

## 高度な使用例

### 複数のアイテムを一度に作成

```powershell
# 複数のファイルをまとめて作成
1..5 | ForEach-Object { New-Item -Path "C:\Temp\file$_.txt" -ItemType File }
```

出力例:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:27              0 file1.txt
-a---          2025/03/28    12:27              0 file2.txt
-a---          2025/03/28    12:27              0 file3.txt
-a---          2025/03/28    12:27              0 file4.txt
-a---          2025/03/28    12:27              0 file5.txt
```

### パイプラインを使って複数の異なるパスにアイテムを作成

```powershell
# 異なるパスに同じタイプのアイテムを作成
"C:\Temp\Project1", "C:\Temp\Project2", "C:\Temp\Project3" | ForEach-Object { New-Item -Path $_ -ItemType Directory }
```

出力例:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/28    12:28                Project1

    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/28    12:28                Project2

    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/28    12:28                Project3
```

### ファイルに複数行の内容を設定

```powershell
# 複数行の内容を持つファイル作成
$content = @"
# Configuration File
AppName = MyApp
Version = 1.0
Environment = Production
"@
New-Item -Path "C:\Temp\config.cfg" -ItemType File -Value $content
```

出力例:
```
    ディレクトリ: C:\Temp

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
-a---          2025/03/28    12:29             69 config.cfg
```

## 出力例の見方

`New-Item`コマンドレットの出力には以下の情報が含まれる：

1. **ディレクトリ行**:
   - 作成されたアイテムの親ディレクトリのパスを表示

2. **Mode**:
   - アイテムの属性を示す記号のセット
   - `d`: ディレクトリ
   - `a`: アーカイブファイル（通常のファイル）
   - `l`: シンボリックリンク
   - `h`: 隠しファイル
   - `r`: 読み取り専用ファイル
   - `-`: 該当する属性なし

3. **LastWriteTime**:
   - アイテムの最終更新日時

4. **Length**:
   - ファイルのサイズ（バイト単位）
   - ディレクトリには表示されない
   - シンボリックリンクには特別な表示（`<SYMLINK>`）がされる場合がある

5. **Name**:
   - アイテムの名前
   - シンボリックリンクやジャンクションの場合は、`[リンク先のパス]`が表示される

アイテムの詳細を確認するには、以下のコマンドを使用する：

```powershell
# 作成したアイテムの詳細情報を表示
Get-Item -Path "C:\Temp\test.txt" | Format-List *
```

出力例:
```
PSPath            : Microsoft.PowerShell.Core\FileSystem::C:\Temp\test.txt
PSParentPath      : Microsoft.PowerShell.Core\FileSystem::C:\Temp
PSChildName       : test.txt
PSDrive           : C
PSProvider        : Microsoft.PowerShell.Core\FileSystem
PSIsContainer     : False
Mode              : -a---
VersionInfo       : File:             C:\Temp\test.txt
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
BaseName          : test
Target            : {}
LinkType          : 
Name              : test.txt
Length            : 0
DirectoryName     : C:\Temp
Directory         : C:\Temp
IsReadOnly        : False
Exists            : True
FullName          : C:\Temp\test.txt
Extension         : .txt
CreationTime      : 2025/03/28 12:15:00
CreationTimeUtc   : 2025/03/28 03:15:00
LastAccessTime    : 2025/03/28 12:15:00
LastAccessTimeUtc : 2025/03/28 03:15:00
LastWriteTime     : 2025/03/28 12:15:00
LastWriteTimeUtc  : 2025/03/28 03:15:00
Attributes        : Archive
```



**注：** 一部の操作（特にシンボリックリンクの作成など）は管理者権限が必要な場合がある。また、ファイルシステム以外のプロバイダー（レジストリなど）での操作では、使用できるアイテムタイプが制限される。