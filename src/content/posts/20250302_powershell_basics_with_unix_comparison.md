---
title:  【PowerShell】基本コマンドとUnixコマンドとの対比
published: 2025-03-02
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-02

# 概要

PowerShellはWindowsのシェル環境およびスクリプト言語であり、コマンドレット（Cmdlet）と呼ばれる「動詞-名詞」形式のコマンドで構成される。システム管理やタスクの自動化に優れており、オブジェクト指向のパイプラインを特徴とする。

現在主に使われているPowerShellのバージョンは以下のとおりである：
- Windows PowerShell 5.1 (Windowsに標準搭載)
- PowerShell 7.x (クロスプラットフォーム対応の最新版)

## PowerShellの基本構造

PowerShellのコマンドレットは基本的に「動詞-名詞」の命名規則に従っている。例えば、`Get-Process`は「プロセスを取得する」という意味である。この一貫した命名規則により、コマンドの役割が把握しやすくなっている。


## 基本コマンド一覧

### ヘルプと情報表示

- `Get-Command`: 利用可能なすべてのコマンドを表示する
- `Get-Help`: 指定したコマンドのヘルプ情報を表示する
- `Get-Alias`: エイリアス（コマンドの別名）一覧を表示する

### ファイルシステム操作

- `Get-ChildItem` (エイリアス: `dir`, `ls`): ディレクトリの内容を表示する
- `Set-Location` (エイリアス: `cd`): カレントディレクトリを変更する
- `New-Item`: 新しいファイルやディレクトリを作成する
- `Remove-Item` (エイリアス: `rm`, `del`): ファイルやディレクトリを削除する
- `Copy-Item` (エイリアス: `cp`, `copy`): ファイルやディレクトリをコピーする
- `Move-Item` (エイリアス: `mv`, `move`): ファイルやディレクトリを移動またはリネームする
- `Get-Content` (エイリアス: `cat`, `type`): ファイルの内容を表示する
- `Set-Content`: ファイルに内容を書き込む

### オブジェクト操作

- `Where-Object` (エイリアス: `where`, `?`): 条件でフィルタリングする
- `Select-Object` (エイリアス: `select`): プロパティを選択して表示する
- `Sort-Object` (エイリアス: `sort`): オブジェクトをソートする
- `ForEach-Object` (エイリアス: `foreach`, `%`): コレクションの各要素に対して処理を行う

### システム管理

- `Get-Process` (エイリアス: `ps`): 実行中のプロセス一覧を表示する
- `Stop-Process` (エイリアス: `kill`): プロセスを終了する
- `Get-Service`: Windowsサービス一覧を表示する
- `Start-Service`: サービスを開始する
- `Stop-Service`: サービスを停止する
- `Restart-Service`: サービスを再起動する

---

## UnixコマンドとPowerShellの対応表

| Unix コマンド | PowerShell コマンド | 説明 |
|--------------|-------------------|-----|
| ls | Get-ChildItem (alias: dir, ls) | ディレクトリ内容の表示 |
| cd | Set-Location (alias: cd) | ディレクトリ変更 |
| mkdir | New-Item -ItemType Directory | ディレクトリ作成 |
| rm | Remove-Item (alias: rm, del) | ファイル/ディレクトリ削除 |
| cp | Copy-Item (alias: cp, copy) | ファイルコピー |
| mv | Move-Item (alias: mv, move) | ファイル移動/リネーム |
| cat | Get-Content (alias: cat, type) | ファイル内容の表示 |
| grep | Select-String | テキスト検索 |
| ps | Get-Process (alias: ps) | プロセス一覧表示 |
| kill | Stop-Process (alias: kill) | プロセス終了 |
| touch | New-Item -ItemType File | 空ファイル作成 |
| man | Get-Help | コマンドのヘルプ表示 |
| pwd | Get-Location (alias: pwd) | 現在のディレクトリを表示 |
| chmod | Set-Acl | アクセス権変更 |
| chown | Set-Acl | 所有者変更 |
| df | Get-PSDrive | ディスク使用状況表示 |
| find | Get-ChildItem -Recurse | ファイル検索 |
| ifconfig | Get-NetIPAddress | ネットワークインターフェース情報 |
| top | Get-Process \| Sort-Object CPU -Descending | リソース使用状況表示 |
| wc | (Get-Content file).Length | 行数/単語数カウント |
| date | Get-Date | 日付・時刻表示 |

---

## ハンズオン実例とコマンド出力例

以下に、基本コマンドの使用例と実際の出力例を示す。

### ヘルプと情報表示

**コマンド一覧を取得**
```powershell
Get-Command
```
出力例（一部抜粋）:
```
CommandType     Name                                               Version    Source
-----------     ----                                               -------    ------
Cmdlet          Add-Content                                        7.0.0.0    Microsoft.PowerShell.Management
Cmdlet          Add-History                                        7.0.0.0    Microsoft.PowerShell.Core
...
```

**特定のコマンドのヘルプを表示**
```powershell
Get-Help Get-Process -Examples
```
出力例:
```
NAME
    Get-Process

EXAMPLES
    Example 1: Get a list of all active processes on the local computer
    Get-Process

    Example 2: Get all available data about one or more processes
    Get-Process -Name winword, explorer | Format-List *
...
```

**PowerShellのバージョンを確認**
```powershell
$PSVersionTable
```
出力例:
```
Name                           Value
----                           -----
PSVersion                      7.3.0
PSEdition                      Core
GitCommitId                    7.3.0
OS                             Microsoft Windows 10.0.19044
Platform                       Win32NT
PSCompatibleVersions           {1.0, 2.0, 3.0, 4.0…}
...
```

### ファイルシステム操作

**ディレクトリの内容を表示**
```powershell
Get-ChildItem
```
または
```powershell
dir
```
出力例:
```
    Directory: C:\Users\ss0832\Documents

Mode                LastWriteTime         Length Name
----                -------------         ------ ----
d-----        2/15/2025   9:30 AM                Projects
d-----        2/28/2025   3:45 PM                Reports
-a----        3/01/2025  10:23 AM           2458 notes.txt
-a----        3/02/2025   2:15 PM          15230 report.docx
```

**ディレクトリを移動**
```powershell
Set-Location C:\Windows
```
出力例（通常は何も表示されない）:
```
PS C:\Windows>
```

**現在のディレクトリを表示**
```powershell
Get-Location
```
出力例:
```
Path
----
C:\Windows
```

**新しいディレクトリを作成**
```powershell
New-Item -ItemType Directory -Name "TestFolder"
```
出力例:
```
    Directory: C:\Windows

Mode                LastWriteTime         Length Name
----                -------------         ------ ----
d-----        3/02/2025   3:47 PM                TestFolder
```

**新しいファイルを作成**
```powershell
New-Item -ItemType File -Name "test.txt" -Value "Hello, PowerShell!"
```
出力例:
```
    Directory: C:\Windows

Mode                LastWriteTime         Length Name
----                -------------         ------ ----
-a----        3/02/2025   3:48 PM             18 test.txt
```

**ファイルの内容を表示**
```powershell
Get-Content test.txt
```
出力例:
```
Hello, PowerShell!
```

**ファイルをコピー**
```powershell
Copy-Item test.txt test_copy.txt
```

**ファイルを削除**
```powershell
Remove-Item test_copy.txt
```

### オブジェクト操作

**フィルタリング：メモリ使用量が100MB以上のプロセスを表示**
```powershell
Get-Process | Where-Object { $_.WorkingSet -gt 100MB }
```
出力例:
```
 NPM(K)    PM(M)      WS(M)     CPU(s)      Id  SI ProcessName
 ------    -----      -----     ------      --  -- -----------
     88    201.75     218.06      69.92    1280   1 chrome
    152    329.47     366.90     127.25    4680   1 firefox
     78    143.32     176.13      59.71    7852   1 explorer
```

**特定のプロパティのみを選択して表示**
```powershell
Get-Process | Select-Object ProcessName, Id, CPU
```
出力例:
```
ProcessName      Id       CPU
-----------      --       ---
chrome         1280     69.92
csrss           548      0.58
dwm            1040     10.20
explorer       7852     59.71
...
```

**テキストファイル内の特定の文字列を検索（Unixのgrep相当）**
```powershell
Select-String -Path "*.txt" -Pattern "PowerShell"
```
出力例:
```
test.txt:1:Hello, PowerShell!
```

**プロセスをCPU使用率でソート（降順）**
```powershell
Get-Process | Sort-Object CPU -Descending | Select-Object ProcessName, CPU -First 5
```
出力例:
```
ProcessName       CPU
-----------       ---
firefox        127.25
chrome          69.92
explorer        59.71
dwm             10.20
powershell       5.43
```

### システム管理

**実行中のプロセス一覧**
```powershell
Get-Process
```
出力例（一部抜粋）:
```
 NPM(K)    PM(M)      WS(M)     CPU(s)      Id  SI ProcessName
 ------    -----      -----     ------      --  -- -----------
     14     2.39      11.12       0.06    9476   0 audiodg
     88    201.75     218.06      69.92    1280   1 chrome
     12     2.23       8.91       0.58     548   0 csrss
     24     3.90      21.40       0.59     644   1 csrss
     21    11.63      24.63      10.20    1040   1 dwm
```

**特定のサービスの状態を確認**
```powershell
Get-Service -Name wuauserv
```
出力例:
```
Status   Name               DisplayName
------   ----               -----------
Running  wuauserv           Windows Update
```

**ネットワーク設定を表示 (Unixのifconfig相当)**
```powershell
Get-NetIPAddress
```
出力例（抜粋）:
```
IPAddress         : 192.168.1.100
InterfaceIndex    : 12
InterfaceAlias    : Ethernet
AddressFamily     : IPv4
...
```

**日付と時刻を表示**
```powershell
Get-Date
```
出力例:
```
Sunday, March 2, 2025 3:48:40 AM
```

**日付と時刻をフォーマット指定で表示**
```powershell
Get-Date -Format "yyyy-MM-dd HH:mm:ss"
```
出力例:
```
2025-03-02 03:48:40
```


PowerShellの特徴は、単なるテキスト処理ではなく、.NETオブジェクトを扱えることである。これにより、複雑なデータ操作も直感的に行える。また、伝統的なUnixコマンドのエイリアスも多く実装されているため、Unixコマンドに慣れた人でも違和感なく使い始めることができる。

バージョン5.1はWindows標準搭載だが、最新のPowerShell 7.xではより多くの機能と改善が提供されており、クロスプラットフォーム（Windows、macOS、Linux）で使用可能である。