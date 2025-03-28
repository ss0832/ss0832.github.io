---
title:  【PowerShell】Copy-Itemコマンドの基本と応用
published: 2025-03-28
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-28

## 概要

PowerShellの`Copy-Item`コマンドレット（エイリアス：`copy`, `cp`, `cpi`）は、ファイル、ディレクトリ、レジストリ項目などのアイテムを、ある場所から別の場所にコピーするために使用するコマンドである。単なるファイルコピーだけでなく、PowerShellのプロバイダーシステムと連携して、レジストリやその他の名前空間での操作も可能である。Windowsのコマンドプロンプトにおける`copy`コマンドやUNIX系OSの`cp`コマンドに相当する機能を提供するが、PowerShellのパイプライン処理やオブジェクト指向の特性を活かした高度な操作が可能である。

## 基本的な使い方

### ファイルのコピー

```powershell
# 基本的なファイルコピー
Copy-Item -Path "C:\Temp\source.txt" -Destination "C:\Temp\Backup\source.txt"
```

出力例:
```
# 成功時は出力なし
```

### ディレクトリのコピー

```powershell
# ディレクトリのコピー
Copy-Item -Path "C:\Temp\SourceFolder" -Destination "C:\Temp\Backup\" -Recurse
```

出力例:
```
# 成功時は出力なし
```

### ワイルドカードを使った複数ファイルのコピー

```powershell
# すべてのテキストファイルをコピー
Copy-Item -Path "C:\Temp\*.txt" -Destination "C:\Temp\Backup\"
```

出力例:
```
# 成功時は出力なし
```

## 全オプションに対するハンズオン

### -Path

コピー元となるアイテムのパスを指定する。ワイルドカードが使用可能。

```powershell
# 特定のパターンに一致するファイルをコピー
Copy-Item -Path "C:\Temp\*report*.xlsx" -Destination "C:\Reports\"
```

出力例:
```
# 成功時は出力なし
```

### -LiteralPath

特殊文字（ワイルドカードなど）をそのままの文字として扱い、コピー元のパスを指定する。

```powershell
# 特殊文字を含むファイル名をそのままコピー
Copy-Item -LiteralPath "C:\Temp\Report[2025].xlsx" -Destination "C:\Reports\"
```

出力例:
```
# 成功時は出力なし
```

### -Destination

コピー先となるパスを指定する。

```powershell
# 別のドライブにファイルをコピー
Copy-Item -Path "C:\Temp\data.db" -Destination "D:\Backup\data.db"
```

出力例:
```
# 成功時は出力なし
```

### -Container

コンテナオブジェクト（ディレクトリなど）をコピーする際に、中身を含めずに空のコンテナだけをコピーする。

```powershell
# 空のディレクトリ構造のみをコピー
Copy-Item -Path "C:\Projects\Structure" -Destination "C:\Temp\Structure" -Container
```

出力例:
```
# 成功時は出力なし
```

### -Force

既存のファイルを上書きするなど、通常はブロックされる操作を強制的に実行する。

```powershell
# 既存のファイルを強制的に上書き
Copy-Item -Path "C:\Temp\config.ini" -Destination "C:\Config\config.ini" -Force
```

出力例:
```
# 成功時は出力なし
```

### -Filter

指定した場所にあるアイテムのうち、指定したパターンに一致するもののみをコピーする。

```powershell
# logで始まるテキストファイルのみをコピー
Copy-Item -Path "C:\Logs\" -Filter "log*.txt" -Destination "C:\Temp\Logs\" -Recurse
```

出力例:
```
# 成功時は出力なし
```

### -Include

指定したパターンに一致するアイテムのみをコピー対象に含める。

```powershell
# テキストファイルとCSVファイルのみをコピー
Copy-Item -Path "C:\Data\*" -Include "*.txt", "*.csv" -Destination "C:\Exports\"
```

出力例:
```
# 成功時は出力なし
```

### -Exclude

指定したパターンに一致するアイテムをコピー対象から除外する。

```powershell
# 一時ファイル以外をすべてコピー
Copy-Item -Path "C:\Project\*" -Exclude "*.tmp", "*.temp" -Destination "C:\Backup\Project\"
```

出力例:
```
# 成功時は出力なし
```

### -Recurse

指定したパスとそのすべてのサブディレクトリ内のアイテムを再帰的にコピーする。

```powershell
# ディレクトリとそのすべての内容を再帰的にコピー
Copy-Item -Path "C:\Projects\App" -Destination "D:\Backup\App" -Recurse
```

出力例:
```
# 成功時は出力なし
```

### -PassThru

コピーしたアイテムを表すオブジェクトをパイプラインに渡す。

```powershell
# コピーしたファイルの情報を取得
$copiedFile = Copy-Item -Path "C:\Temp\data.xml" -Destination "C:\Data\" -PassThru
$copiedFile | Format-List *
```

出力例:
```
PSPath            : Microsoft.PowerShell.Core\FileSystem::C:\Data\data.xml
PSParentPath      : Microsoft.PowerShell.Core\FileSystem::C:\Data
PSChildName       : data.xml
PSDrive           : C
PSProvider        : Microsoft.PowerShell.Core\FileSystem
PSIsContainer     : False
Mode              : -a---
VersionInfo       : File:             C:\Data\data.xml
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
BaseName          : data
Target            : {}
LinkType          : 
Name              : data.xml
Length            : 12345
DirectoryName     : C:\Data
Directory         : C:\Data
IsReadOnly        : False
Exists            : True
FullName          : C:\Data\data.xml
Extension         : .xml
CreationTime      : 2025/03/28 12:10:11
CreationTimeUtc   : 2025/03/28 12:10:11
LastAccessTime    : 2025/03/28 12:10:11
LastAccessTimeUtc : 2025/03/28 12:10:11
LastWriteTime     : 2025/03/28 12:10:11
LastWriteTimeUtc  : 2025/03/28 12:10:11
Attributes        : Archive
```

### -Credential

指定した認証情報を使用して、別のユーザーとしてコマンドを実行する。

```powershell
# 別のユーザー権限でファイルをコピー
$cred = Get-Credential
Copy-Item -Path "C:\Reports\quarterly.xlsx" -Destination "\\Server\Share\" -Credential $cred
```

出力例:
```
# 認証情報入力ダイアログが表示される
# 成功時は出力なし
```

### -WhatIf

コマンドを実際に実行せず、実行された場合の動作を表示する。

```powershell
# 実行されるコピー操作を確認
Copy-Item -Path "C:\Important\*.docx" -Destination "D:\Backup\" -WhatIf
```

出力例:
```
WhatIf: C:\Important\report.docx を D:\Backup\report.docx にコピーします
WhatIf: C:\Important\contract.docx を D:\Backup\contract.docx にコピーします
WhatIf: C:\Important\proposal.docx を D:\Backup\proposal.docx にコピーします
```

### -Confirm

コマンドを実行する前に確認を求める。

```powershell
# コピー前に確認を求める
Copy-Item -Path "C:\Temp\largefile.iso" -Destination "D:\Backup\" -Confirm
```

出力例:
```
Confirm
C:\Temp\largefile.iso を D:\Backup\largefile.iso にコピーしますか？
[Y] はい(Y)  [A] すべて続行(A)  [N] いいえ(N)  [L] すべて無視(L)  [S] 中断(S)  [?] ヘルプ (既定値は "Y"): y
```

## 高度な使用例

### パイプラインを使った条件付きコピー

```powershell
# 特定の日付以降に作成されたファイルのみをコピー
Get-ChildItem -Path "C:\Data\*.xlsx" | Where-Object { $_.LastWriteTime -gt (Get-Date).AddDays(-7) } | Copy-Item -Destination "C:\Reports\Recent\"
```

出力例:
```
# 成功時は出力なし
```

### フォルダ構造の維持

```powershell
# フォルダ構造を維持したまま特定の拡張子のファイルだけをコピー
Get-ChildItem -Path "C:\Projects" -Include "*.cs" -Recurse | ForEach-Object {
    $destPath = $_.FullName.Replace("C:\Projects", "C:\Backup\CodeOnly")
    $destDir = Split-Path -Path $destPath -Parent
    if (!(Test-Path -Path $destDir)) {
        New-Item -Path $destDir -ItemType Directory -Force
    }
    Copy-Item -Path $_.FullName -Destination $destPath
}
```

出力例:
```
# ディレクトリ作成時に出力がある場合があります
    ディレクトリ: C:\Backup\CodeOnly\SubProject1

Mode                 LastWriteTime         Length Name
----                 -------------         ------ ----
d----          2025/03/28    12:10                

# その他は成功時に出力なし
```

### レジストリキーのコピー

```powershell
# レジストリキーのコピー
Copy-Item -Path "HKCU:\Software\MyApp\Settings" -Destination "HKCU:\Software\MyApp\SettingsBackup" -Recurse
```

出力例:
```
# 成功時は出力なし
```

### 進行状況の表示

```powershell
# 大きなファイルをコピーする際の進行状況表示
$sourceFile = "C:\Temp\largefile.iso"
$destFile = "D:\Backup\largefile.iso"
$source = [System.IO.File]::OpenRead($sourceFile)
$dest = [System.IO.File]::Create($destFile)
$buffer = New-Object byte[] 1MB
$totalSize = $source.Length
$bytesTransferred = 0

try {
    while (($bytesRead = $source.Read($buffer, 0, $buffer.Length)) -gt 0) {
        $dest.Write($buffer, 0, $bytesRead)
        $bytesTransferred += $bytesRead
        $percentComplete = [math]::Round(($bytesTransferred / $totalSize) * 100, 2)
        Write-Progress -Activity "ファイルをコピーしています" -Status "$percentComplete% 完了" -PercentComplete $percentComplete
    }
} finally {
    $source.Close()
    $dest.Close()
}
```

出力例:
```
# プログレスバーが表示される（コンソール上でグラフィカルに表示）
```

## 出力例の見方

`Copy-Item`コマンドレットは、基本的に成功した場合は何も出力しない。これはPowerShellの多くのコマンドレットに共通する「沈黙は金」の原則に基づいている。ただし、以下の場合には出力がある：

1. **-PassThruパラメータ使用時**:
   - コピーされたアイテムを表すオブジェクトが返される
   - このオブジェクトには、アイテムのパス、サイズ、作成日時などの属性が含まれる
   - パイプラインを通じて次のコマンドに渡すことができる

2. **-WhatIfパラメータ使用時**:
   - `WhatIf: [元のパス] を [コピー先のパス] にコピーします` という形式で、コピーされる予定のアイテムが表示される
   - 実際のコピー操作は行われない

3. **-Confirmパラメータ使用時**:
   - コピー前に確認ダイアログが表示される
   - 選択肢: [Y]はい、[A]すべて続行、[N]いいえ、[L]すべて無視、[S]中断、[?]ヘルプ

4. **-Verboseパラメータ使用時**:
   - `VERBOSE: [元のパス] を [コピー先のパス] にコピーしています` という形式で、実行中の操作が表示される
   - 処理の詳細な流れを把握できる

5. **エラー時**:
   - アイテムが存在しない場合: `Copy-Item : [パス] が見つかりません`
   - アクセス拒否の場合: `Copy-Item : アクセスが拒否されました`
   - 既存のファイルがあり-Forceが指定されていない場合: `Copy-Item : [パス] が既に存在します`

6. **実行結果の確認方法**:
   - コピーの成功は、通常、戻り値やエラーメッセージがないことで確認できる
   - コピーが成功したかを明示的に確認するには、`Test-Path`コマンドレットやスクリプト内での`$?`変数の確認が有効



**注意点**: `Copy-Item`コマンドは存在するファイルを上書きする可能性があるため、特に重要なファイルを扱う場合は、`-WhatIf`パラメータを使用して実行前に操作内容を確認することをお勧めする。また、複雑なコピー操作では、スクリプトを作成して段階的に実行・検証することで、予期せぬデータ損失を防ぐことができる。