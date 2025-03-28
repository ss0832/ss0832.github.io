---
title:  【PowerShell】Remove-Itemコマンドの基本と応用
published: 2025-03-28
description: ""
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-28

## 概要

PowerShellの`Remove-Item`コマンドレット（エイリアス：`rm`, `del`, `erase`, `ri`）は、ファイルやディレクトリ、レジストリキー、変数など、PowerShellプロバイダーがサポートする様々なアイテムを削除するためのコマンドである。このコマンドレットは、Windowsのコマンドプロンプトにおける`del`コマンドやUNIX系OSの`rm`コマンドに相当する機能を提供するが、PowerShellのパイプラインやプロバイダーモデルと統合されており、より柔軟で強力な削除操作が可能である。

## 基本的な使い方

### ファイルの削除

```powershell
# 基本的なファイル削除
Remove-Item -Path "C:\Temp\test.txt"
```

出力例:
```
# 成功時は出力なし
```

### ディレクトリの削除

```powershell
# 空のディレクトリを削除
Remove-Item -Path "C:\Temp\EmptyFolder"
```

出力例:
```
# 成功時は出力なし
```

### ワイルドカードを使った複数ファイルの削除

```powershell
# 拡張子.tmpのファイルをすべて削除
Remove-Item -Path "C:\Temp\*.tmp"
```

出力例:
```
# 成功時は出力なし
```

## 全オプションに対するハンズオン

### -Path

削除するアイテムのパスを指定する。ワイルドカード（`*`, `?`）も使用可能。

```powershell
# 特定のディレクトリ内のテキストファイルをすべて削除
Remove-Item -Path "C:\Temp\*.txt"

# 複数のパスを指定（カンマ区切り）
Remove-Item -Path "C:\Temp\file1.txt", "C:\Temp\file2.txt"
```

出力例:
```
# 成功時は出力なし
```

### -LiteralPath

特殊文字（ワイルドカードなど）をそのままの文字として扱い、パスを指定する。

```powershell
# 名前に角括弧を含むファイルを削除
Remove-Item -LiteralPath "C:\Temp\Report[2025].xlsx"
```

出力例:
```
# 成功時は出力なし
```

### -Filter

指定した場所にあるアイテムのうち、指定したパターンに一致するもののみを削除する。

```powershell
# テキストファイルのうち、"log"を含むファイルのみ削除
Remove-Item -Path "C:\Temp\" -Filter "*log*.txt"
```

出力例:
```
# 成功時は出力なし
```

### -Include

指定したパスに含まれるアイテムのうち、特定のパターンに一致するもののみを削除する。

```powershell
# テキストファイルと設定ファイルのみ削除
Remove-Item -Path "C:\Temp\*" -Include "*.txt", "*.ini", "*.cfg"
```

出力例:
```
# 成功時は出力なし
```

### -Exclude

指定したパターンに一致するアイテムを削除対象から除外する。

```powershell
# テキストファイル以外のファイルを削除
Remove-Item -Path "C:\Temp\*" -Exclude "*.txt"
```

出力例:
```
# 成功時は出力なし
```

### -Recurse

指定したパスとそのすべてのサブディレクトリ内のアイテムを再帰的に削除する。

```powershell
# フォルダとその中身をすべて削除
Remove-Item -Path "C:\Temp\Project" -Recurse
```

出力例:
```
# 成功時は出力なし
# フォルダ内のアイテムが読み取り専用などの場合は確認メッセージが表示されることがある
```

### -Force

読み取り専用属性が設定されているファイルなど、通常は削除できないアイテムを強制的に削除する。

```powershell
# 読み取り専用ファイルを強制削除
Remove-Item -Path "C:\Temp\readonly.txt" -Force
```

出力例:
```
# 成功時は出力なし
```

### -Credential

別のユーザーとしてコマンドを実行するための認証情報を指定する。

```powershell
# 別ユーザーとしてファイルを削除
$credential = Get-Credential
Remove-Item -Path "\\Server\Share\file.txt" -Credential $credential
```

出力例:
```
# 認証情報入力ダイアログが表示される
# 成功時は出力なし
```

### -WhatIf

コマンドを実際に実行せず、実行された場合の動作を表示する。

```powershell
# 削除されるファイルを事前確認
Remove-Item -Path "C:\Temp\*.log" -WhatIf
```

出力例:
```
WhatIf: C:\Temp\app.log のアイテムを削除します
WhatIf: C:\Temp\system.log のアイテムを削除します
WhatIf: C:\Temp\error.log のアイテムを削除します
```

### -Confirm

コマンドを実行する前に確認を求める。

```powershell
# 削除前に確認を求める
Remove-Item -Path "C:\Temp\important.docx" -Confirm
```

出力例:
```
Delete item C:\Temp\important.docx?
[Y] はい(Y)  [A] すべて続行(A)  [N] いいえ(N)  [L] すべて無視(L)  [S] 中断(S)  [?] ヘルプ (既定値は "Y"): n
```

### -Stream

NTFS代替データストリームを削除する（Windows専用）。

```powershell
# ファイルに添付された代替データストリームを削除
Remove-Item -Path "C:\Temp\document.docx" -Stream "Zone.Identifier"
```

出力例:
```
# 成功時は出力なし
```

## 高度な使用例

### パイプラインを使った条件付き削除

```powershell
# 30日以上前に作成されたログファイルを削除
Get-ChildItem -Path "C:\Logs\*.log" | Where-Object { $_.CreationTime -lt (Get-Date).AddDays(-30) } | Remove-Item
```

出力例:
```
# 成功時は出力なし
```

### レジストリキーの削除

```powershell
# レジストリキーの削除
Remove-Item -Path "HKCU:\Software\TempApp" -Recurse
```

出力例:
```
# 成功時は出力なし
```

### 複数の条件を組み合わせた削除

```powershell
# 複数条件による削除：読み取り専用属性以外の.tmpと.bakファイルを再帰的に削除
Get-ChildItem -Path "C:\Projects" -Recurse -Include "*.tmp", "*.bak" | Where-Object { !$_.IsReadOnly } | Remove-Item
```

出力例:
```
# 成功時は出力なし
```

## 安全な削除方法

### -WhatIfと-Verboseを組み合わせた確認

```powershell
# 実行前に詳細確認
Remove-Item -Path "C:\Important\*" -Include "*.bak" -Recurse -WhatIf -Verbose
```

出力例:
```
VERBOSE: パス C:\Important\*.bak のマッチング中です...
VERBOSE: サブディレクトリ C:\Important\Project1 を処理中です...
WhatIf: C:\Important\Project1\backup.bak のアイテムを削除します
VERBOSE: サブディレクトリ C:\Important\Project2 を処理中です...
WhatIf: C:\Important\Project2\old.bak のアイテムを削除します
```

### 段階的な削除プロセス

```powershell
# 段階的に確認しながら削除
$files = Get-ChildItem -Path "C:\Temp" -Include "*.tmp" -Recurse
$files | Format-Table Name, Length, LastWriteTime -AutoSize
$confirmation = Read-Host "上記のファイルを削除しますか？ (Y/N)"
if ($confirmation -eq 'Y') {
    $files | Remove-Item -Verbose
}
```

出力例:
```
Name           Length LastWriteTime
----           ------ -------------
temp1.tmp        1024 2025/03/25 14:30:22
temp2.tmp        2048 2025/03/26 09:15:45
cache.tmp       10240 2025/03/27 16:45:12

上記のファイルを削除しますか？ (Y/N): Y
VERBOSE: C:\Temp\temp1.tmp のアイテムを削除します
VERBOSE: C:\Temp\temp2.tmp のアイテムを削除します
VERBOSE: C:\Temp\cache.tmp のアイテムを削除します
```

## 出力例の見方

`Remove-Item`コマンドレットは、基本的に成功した場合は何も出力しない。これはUnixの「沈黙は金」の哲学に従っている。ただし、以下の場合には出力がある：

1. **-WhatIfパラメータ使用時**:
   - `WhatIf: [パス] のアイテムを削除します` という形式で、削除される予定のアイテムが表示される
   - 実際の削除は行われない

2. **-Confirmパラメータ使用時**:
   - 削除前に確認ダイアログが表示される
   - 選択肢: [Y]はい、[A]すべて続行、[N]いいえ、[L]すべて無視、[S]中断、[?]ヘルプ

3. **-Verboseパラメータ使用時**:
   - `VERBOSE: [パス] のアイテムを削除します` という形式で、実行中の操作が表示される
   - 処理の詳細な流れを把握できる

4. **エラー時**:
   - アイテムが存在しない場合: `Remove-Item : [パス] にアイテムが存在しません`
   - アクセス拒否の場合: `Remove-Item : アクセスが拒否されました`
   - ディレクトリが空でない場合（-Recurseなしの場合）: `Remove-Item : ディレクトリが空ではありません`

5. **後続の処理用**:
   - 通常は出力がないため、成功したかどうかは `$?` 変数や `$LASTEXITCODE` で確認できる



**注意点**: `Remove-Item`コマンドは削除操作を実行するため、特に`-Recurse`や`-Force`パラメータを使用する場合は十分な注意が必要である。重要なデータを誤って削除するリスクがあるため、実行前に`-WhatIf`パラメータで確認するか、バックアップを取ることを強く推奨する。