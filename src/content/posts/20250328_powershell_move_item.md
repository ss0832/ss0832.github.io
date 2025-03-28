---
title:  【PowerShell】Move-Itemコマンドの基本と応用
published: 2025-03-28
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-28

## 概要

PowerShellの`Move-Item`コマンドレット（エイリアス：`move`, `mv`, `mi`）は、ファイル、ディレクトリ、レジストリキーなどのアイテムを、ある場所から別の場所に移動するために使用するコマンドである。同一ドライブ内での移動はファイルシステム上で実際に移動が行われるが、異なるドライブ間では内部的にはコピー後に元のアイテムを削除する処理となる。また、同一ディレクトリ内で移動先に異なる名前を指定した場合は、アイテムの名前変更として機能する。WindowsのコマンドプロンプトにおけるMOVEコマンドやUNIX系OSのmvコマンドに相当する機能だが、PowerShellのプロバイダーモデルにより、ファイルシステム以外でも一貫した操作が可能である。

## 基本的な使い方

### ファイルの移動

```powershell
# ファイルを別のディレクトリに移動
Move-Item -Path "C:\Temp\report.docx" -Destination "C:\Reports\"
```

出力例:
```
# 成功時は出力なし
```

### ディレクトリの移動

```powershell
# ディレクトリを別の場所に移動
Move-Item -Path "C:\Temp\ProjectFiles" -Destination "C:\Projects\"
```

出力例:
```
# 成功時は出力なし
```

### ファイル名の変更

```powershell
# ファイル名を変更
Move-Item -Path "C:\Reports\old-report.xlsx" -Destination "C:\Reports\Q1-2025-report.xlsx"
```

出力例:
```
# 成功時は出力なし
```

## 全オプションに対するハンズオン

### -Path

移動するアイテムのパスを指定する。ワイルドカードが使用可能。

```powershell
# 複数のテキストファイルを移動
Move-Item -Path "C:\Temp\*.txt" -Destination "C:\TextFiles\"
```

出力例:
```
# 成功時は出力なし
```

### -LiteralPath

特殊文字（ワイルドカードなど）をそのままの文字として扱い、移動元のパスを指定する。

```powershell
# 角括弧を含むファイル名をそのまま移動
Move-Item -LiteralPath "C:\Temp\Data[2025].csv" -Destination "C:\Reports\"
```

出力例:
```
# 成功時は出力なし
```

### -Destination

移動先のパスを指定する。ディレクトリか、移動先の完全なパス名である必要がある。

```powershell
# 別のドライブに移動
Move-Item -Path "C:\Temp\backup.zip" -Destination "D:\Archives\"
```

出力例:
```
# 成功時は出力なし
```

### -Force

読み取り専用ファイルを移動したり、既存のファイルを上書きするなど、通常はブロックされる操作を強制的に実行する。

```powershell
# 読み取り専用ファイルを強制的に移動
Move-Item -Path "C:\ReadOnly\config.ini" -Destination "C:\Config\" -Force
```

出力例:
```
# 成功時は出力なし
```

### -Filter

指定した場所にあるアイテムのうち、指定したパターンに一致するもののみを移動する。

```powershell
# logで始まるテキストファイルのみを移動
Move-Item -Path "C:\Logs\" -Filter "log*.txt" -Destination "C:\Archives\Logs\"
```

出力例:
```
# 成功時は出力なし
```

### -Include

指定したパターンに一致するアイテムのみを移動対象に含める。

```powershell
# CSVファイルとExcelファイルのみを移動
Move-Item -Path "C:\Data\*" -Include "*.csv", "*.xlsx" -Destination "C:\Reports\"
```

出力例:
```
# 成功時は出力なし
```

### -Exclude

指定したパターンに一致するアイテムを移動対象から除外する。

```powershell
# 一時ファイル以外をすべて移動
Move-Item -Path "C:\Project\*" -Exclude "*.tmp", "*.bak" -Destination "C:\Archive\Project\"
```

出力例:
```
# 成功時は出力なし
```

### -PassThru

移動したアイテムを表すオブジェクトをパイプラインに渡す。

```powershell
# 移動したファイルの情報を取得
$movedFile = Move-Item -Path "C:\Temp\data.xml" -Destination "C:\Data\" -PassThru
$movedFile | Format-List *
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
CreationTime      : 2025/03/28 11:30:00
CreationTimeUtc   : 2025/03/28 02:30:00
LastAccessTime    : 2025/03/28 12:05:15
LastAccessTimeUtc : 2025/03/28 03:05:15
LastWriteTime     : 2025/03/28 12:05:15
LastWriteTimeUtc  : 2025/03/28 03:05:15
Attributes        : Archive
```

### -Credential

指定した認証情報を使用して、別のユーザーとしてコマンドを実行する。

```powershell
# 別のユーザー権限でファイルを移動
$cred = Get-Credential
Move-Item -Path "C:\Reports\quarterly.xlsx" -Destination "\\Server\Share\" -Credential $cred
```

出力例:
```
# 認証情報入力ダイアログが表示される
# 成功時は出力なし
```

### -WhatIf

コマンドを実際に実行せず、実行された場合の動作を表示する。

```powershell
# 実行される移動操作を確認
Move-Item -Path "C:\Important\*.docx" -Destination "D:\Archive\" -WhatIf
```

出力例:
```
WhatIf: C:\Important\report.docx を D:\Archive\report.docx に移動します
WhatIf: C:\Important\contract.docx を D:\Archive\contract.docx に移動します
WhatIf: C:\Important\proposal.docx を D:\Archive\proposal.docx に移動します
```

### -Confirm

コマンドを実行する前に確認を求める。

```powershell
# 移動前に確認を求める
Move-Item -Path "C:\Temp\important-data.xlsx" -Destination "D:\Backup\" -Confirm
```

出力例:
```
Confirm
C:\Temp\important-data.xlsx を D:\Backup\important-data.xlsx に移動しますか？
[Y] はい(Y)  [A] すべて続行(A)  [N] いいえ(N)  [L] すべて無視(L)  [S] 中断(S)  [?] ヘルプ (既定値は "Y"): y
```

## 高度な使用例

### パイプラインを使った条件付き移動

```powershell
# 30日以上前に作成されたログファイルを移動
Get-ChildItem -Path "C:\Logs\*.log" | Where-Object { $_.CreationTime -lt (Get-Date).AddDays(-30) } | Move-Item -Destination "C:\Archives\OldLogs\"
```

出力例:
```
# 成功時は出力なし
```

### レジストリキーの移動

```powershell
# レジストリキーの移動（名前変更）
Move-Item -Path "HKCU:\Software\MyApp\OldSettings" -Destination "HKCU:\Software\MyApp\Settings"
```

出力例:
```
# 成功時は出力なし
```

### ファイル名にタイムスタンプを追加して移動

```powershell
# ファイルにタイムスタンプを追加して移動
$timestamp = Get-Date -Format "yyyy-MM-dd-HHmmss"
Get-ChildItem -Path "C:\Temp\*.log" | ForEach-Object {
    $newName = "$($_.BaseName)-$timestamp$($_.Extension)"
    Move-Item -Path $_.FullName -Destination "C:\Logs\$newName"
}
```

出力例:
```
# 成功時は出力なし
```

### 安全な移動（バックアップ作成後に移動）

```powershell
# 重要ファイルの移動前にバックアップを作成
$file = "C:\Important\critical-data.db"
$destination = "D:\Production\critical-data.db"
$backup = "C:\Backup\critical-data-$(Get-Date -Format 'yyyyMMdd').db"

# バックアップ作成
Copy-Item -Path $file -Destination $backup
if (Test-Path $backup) {
    # バックアップ確認後に移動
    Move-Item -Path $file -Destination $destination -Verbose
    Write-Host "ファイル移動が完了しました。バックアップ: $backup"
} else {
    Write-Error "バックアップの作成に失敗しました。ファイル移動は実行されません。"
}
```

出力例:
```
VERBOSE: C:\Important\critical-data.db を D:\Production\critical-data.db に移動しています
ファイル移動が完了しました。バックアップ: C:\Backup\critical-data-20250328.db
```

## 出力例の見方

`Move-Item`コマンドレットは、基本的に成功した場合は何も出力しない。これはPowerShellの多くのコマンドレットに共通する「沈黙は金」の原則に基づいている。ただし、以下の場合には出力がある：

1. **-PassThruパラメータ使用時**:
   - 移動されたアイテムを表すオブジェクトが返される
   - このオブジェクトには、アイテムの新しいパス、サイズ、作成日時などの属性が含まれる
   - パイプラインを通じて次のコマンドに渡すことができる

2. **-WhatIfパラメータ使用時**:
   - `WhatIf: [元のパス] を [移動先のパス] に移動します` という形式で、移動される予定のアイテムが表示される
   - 実際の移動操作は行われない

3. **-Confirmパラメータ使用時**:
   - 移動前に確認ダイアログが表示される
   - 選択肢: [Y]はい、[A]すべて続行、[N]いいえ、[L]すべて無視、[S]中断、[?]ヘルプ

4. **-Verboseパラメータ使用時**:
   - `VERBOSE: [元のパス] を [移動先のパス] に移動しています` という形式で、実行中の操作が表示される
   - 処理の詳細な流れを把握できる

5. **エラー時**:
   - アイテムが存在しない場合: `Move-Item : [パス] が見つかりません`
   - アクセス拒否の場合: `Move-Item : アクセスが拒否されました`
   - 既存のファイルがあり-Forceが指定されていない場合: `Move-Item : [パス] が既に存在します`
   - 移動先フォルダが存在しない場合: `Move-Item : [パス] が見つかりません`

6. **実行結果の確認方法**:
   - 移動の成功は、通常、戻り値やエラーメッセージがないことで確認できる
   - 移動が成功したかを明示的に確認するには、`Test-Path`コマンドレットや`$?`変数の確認が有効



**注意点**: `Move-Item`コマンドは元のアイテムを削除するため、特に重要なファイルを扱う場合は、`-WhatIf`パラメータを使用して実行前に操作内容を確認するか、バックアップを取っておくことをお勧めする。異なるドライブ間での移動は内部的にコピー＋削除となるため、大きなファイルの場合は時間がかかることがある。また、移動先に既に同名のアイテムが存在する場合、`-Force`パラメータを指定しない限り、エラーとなる点に注意が必要である。