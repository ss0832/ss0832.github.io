---
title: 【PowerShell】フォルダ内の複数ファイルを一括処理する方法
published: 2025-11-29
description: "Windowsバッチファイルのforループ処理をPowerShellに移行し、Get-ChildItemとForEach-Objectを用いて複数ファイルに対する外部コマンド実行を自動化する方法を解説する。"
tags: [PowerShell, Automation, Batch]
category: Windows
draft: false
---

最終更新：2025-11-29

※備忘録として、Geminiが出力した内容を交えながらここに記す。

## **概要**


PowerShellでは、Get-ChildItemとパイプライン、そしてオブジェクトのプロパティを利用することで、より柔軟で可読性の高い一括処理が可能である。

本記事では、ワイルドカードで特定した複数の入力ファイル（例：.xyz形式）に対し、ファイル名の一部を加工して設定ファイル名（例：.json形式）を生成し、外部プログラムを実行する自動化スクリプトの作成方法を解説する。

数値計算やデータ処理の自動化において、パラメータを変えながら連続実行する際に役立つ。

## **基本的な使い方**

### **バッチファイルでの実装例（移行元）**

まず、移行元となるバッチファイルの処理を確認する。このスクリプトは、特定のパターンに一致するファイルを検索し、そのファイル名を引数としてPythonプログラムを実行している。
```
@echo off  
echo Starting AutoTS runs...

REM sn2_Et*_*.xyz に一致するファイルを検索しループ処理  
for %%f in (sn2_Et*_*.xyz) do (  
    echo Processing %%f  
    # REM %%~nf は拡張子を除いたファイル名に展開される  
    python run_autots.py %%f -cfg config_%%~nf.json  
)

echo All AutoTS runs finished.  
pause
```
### **PowerShellでの実装例**

上記と同じ動作をPowerShellで行う場合の基本コードである。PowerShellでは、文字列操作ではなくオブジェクトのプロパティとしてファイル名を扱うことができる。

```
Write-Host "Starting AutoTS runs..."

# Get-ChildItemでファイルを取得し、ForEach-Objectでループ処理  
Get-ChildItem "sn2_Et*_*.xyz" | ForEach-Object {  
    Write-Host "Processing $($_.Name)"  
      
    # 外部プログラムの実行  
    # $_.Name はファイル名（拡張子あり）、$_.BaseName はファイル名（拡張子なし）  
    python run_autots.py $_.Name -cfg "config_$($_.BaseName).json"  
}

Write-Host "All AutoTS runs finished."  
Read-Host "Press Enter to continue..."
```
出力例:
```
Starting AutoTS runs...  
Processing sn2_Et_01.xyz  
Processing sn2_Et_02.xyz  
All AutoTS runs finished.  
Press Enter to continue...:
```
## **全オプションに対するハンズオン**

PowerShellでこの処理を構築する各要素について詳しく解説する。

### **Get-ChildItem（ファイル検索）**

バッチの dir や for ... in (...) に相当する。ワイルドカードを使用して処理対象を絞り込む。
```
# 特定のパターンにマッチするファイルのみオブジェクトとして取得  
$files = Get-ChildItem -Path ".\data" -Filter "sn2_Et*_*.xyz"

# 取得した数を確認  
Write-Host "対象ファイル数: $($files.Count)"
```
### **ForEach-Object（ループ処理）**

パイプライン | を通して渡されたファイルオブジェクトを1つずつ処理する。$_ は現在処理中のファイルオブジェクト（この場合は System.IO.FileInfo オブジェクト）を表す。
```
Get-ChildItem "*.xyz" | ForEach-Object {  
    # ここに処理を書く  
    Write-Host "現在処理中: $_"  
}
```
### **ファイル名プロパティの利用（重要）**

バッチファイルの %%~nf のような特殊な記号による修飾子は不要である。PowerShellではファイルオブジェクトが持つプロパティを利用する。
```
Get-ChildItem "sn2_Et_01.xyz" | ForEach-Object {  
    Write-Host "完全パス ($_.FullName):" $_.FullName  
    Write-Host "ファイル名 ($_.Name):"     $_.Name  
    Write-Host "ベース名   ($_.BaseName):" $_.BaseName  
    Write-Host "拡張子     ($_.Extension):" $_.Extension  
}
```
出力例:
```
完全パス ($_.FullName): C:\Work\sn2_Et_01.xyz  
ファイル名 ($_.Name):     sn2_Et_01.xyz  
ベース名   ($_.BaseName): sn2_Et_01  
拡張子     ($_.Extension): .xyz
```
これにより、`config_$($_.BaseName).json` のように記述することで、sn2_Et_01.xyz から config_sn2_Et_01.json という文字列を動的に生成できる。文字列展開を行うために $() で囲んでいる点に注意する。

### **外部コマンドの呼び出し演算子 &**

外部プログラムのパスにスペースが含まれる場合や、コマンド名を変数に入れて動的に変更したい場合は、呼び出し演算子 & を使用すると安全である。
```
$prog = "python"  
$script = "run_autots.py"

Get-ChildItem "*.xyz" | ForEach-Object {  
    # コマンドとして実行  
    & $prog $script $_.Name  
}
```
## **高度な使用例**

### **実行結果（ExitCode）の確認とエラーハンドリング**

外部プログラムが正常終了したかどうかを確認し、エラーがあった場合に警告を表示する。計算が収束しなかった場合などの検知に役立つ。
```
Get-ChildItem "sn2_Et*_*.xyz" | ForEach-Object {  
    Write-Host "Running for $($_.Name)..." -ForegroundColor Cyan  
      
    # コマンドを実行  
    python run_autots.py $_.Name -cfg "config_$($_.BaseName).json"  
      
    # $LASTEXITCODE に直前のコマンドの終了コードが格納される  
    if ($LASTEXITCODE -ne 0) {  
        Write-Warning "エラーが発生しました: $($_.Name) (ExitCode: $LASTEXITCODE)"  
    }  
}
```
### **Start-Process を使用した制御**

より厳密にプロセス管理を行いたい場合（例えば、新しいウィンドウで立ち上げたい、明示的に待機したいなど）は Start-Process を使用する。
```
Get-ChildItem "*.xyz" | ForEach-Object {  
    $argsList = @(  
        "run_autots.py",  
        $_.Name,  
        "-cfg",  
        "config_$($_.BaseName).json"  
    )

    # プロセスを開始し、終了まで待機(-Wait)する  
    # -NoNewWindowを指定しない場合、別のコンソール窓が開く  
    Start-Process -FilePath "python" -ArgumentList $argsList -Wait -NoNewWindow  
}
```
### **ログファイルへの出力**

処理内容を画面だけでなくファイルにも保存する。Tee-Object を使うと画面表示とファイル保存を同時に行える。長時間かかる計算を放置する場合に推奨される。
```
$logFile = "execution_log.txt"

Get-ChildItem "*.xyz" | ForEach-Object {  
    $msg = "$(Get-Date): Processing $($_.Name)"  
    Write-Output $msg  
      
    # 実際のコマンド実行（出力もキャプチャしたい場合）  
    python run_autots.py $_.Name | Out-String | Write-Output  
} | Tee-Object -FilePath $logFile
```
## **出力例の見方**

本スクリプトを実行した際の挙動と、各ストリームの関係を以下に示す。

1. **Write-Host**:  
   * コンソール画面に情報を表示する。  
   * 色付け（-ForegroundColor）が可能で、進捗状況の視認性を高める。ログファイル等にはリダイレクトされないため、あくまでユーザーへの通知用である。  
2. **外部プログラム（Pythonなど）の出力**:  
   * PowerShellコンソールにそのまま標準出力として表示される。  
   * バッチファイルと同様に振る舞う。  
3. **Read-Host**:  
   * バッチファイルの pause に相当。  
   * ユーザーがEnterキーを押すまでスクリプトの終了を待機する。

## **使用可能なPowerShellバージョン**

本記事のコードは、以下の環境で動作確認済みである。

* PowerShell 5.1（Windows 10/11 デフォルト）  
* PowerShell 7.x（Core）

バージョンによる主な違い：

* **PowerShell 5.1**: ForEach-Object はシーケンシャル（順次）処理のみ。  
* **PowerShell 7.0以降**: ForEach-Object \-Parallel オプションが使用可能になり、並列処理（マルチスレッド処理）が簡単に実装できる。計算資源に余裕がある場合は7.xへのアップグレードを推奨する。

各バージョンでの確認方法：

$PSVersionTable

出力例（PowerShell 5.1の場合）:
```
Name                           Value  
----                           -----  
PSVersion                      5.1.26100.7019  
PSEdition                      Desktop  
PSCompatibleVersions           {1.0, 2.0, 3.0, 4.0...}  
BuildVersion                   10.0.26100.7019  
CLRVersion                     4.0.30319.42000  
WSManStackVersion              3.0  
PSRemotingProtocolVersion      2.3  
SerializationVersion           1.1.0.1
```
**参考**: 外部コマンドへの引数の渡し方でトラブルが起きた場合（スペースを含むパスなど）、cmd /c を経由するか、Start-Process の \-ArgumentList を配列として渡す方法を検討すると解決することが多い。

## **参考サイト**

* [Microsoft Docs: Get-ChildItem](https://learn.microsoft.com/ja-jp/powershell/module/microsoft.powershell.management/get-childitem)  
* [Microsoft Docs: ForEach-Object](https://learn.microsoft.com/ja-jp/powershell/module/microsoft.powershell.core/foreach-object)  
* [Microsoft Docs: Start-Process](https://learn.microsoft.com/ja-jp/powershell/module/microsoft.powershell.management/start-process)  
* [Qiita \- PowerShellで外部コマンドを呼び出すいくつかの方法](https://www.google.com/search?q=https://qiita.com/rawr/items/325eb2230198150499e7)