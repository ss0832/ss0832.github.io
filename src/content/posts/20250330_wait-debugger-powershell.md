---
title:  【PowerShell】Wait-Debuggerコマンドとデバッグモードの活用法
published: 2025-03-30
description: "PowerShellのWait-Debuggerコマンドレットとデバッグモードの使い方を解説"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Wait-Debugger`はPowerShell 5.0以降で導入されたコマンドレットで、スクリプトの実行を一時停止してデバッガーを起動するためのものである。このコマンドを使用すると、スクリプト内の特定のポイントで実行を停止し、変数の値の検査や、ステップ実行などのデバッグ作業を行うことができる。従来のブレークポイントと異なり、コード内に直接組み込むことができるため、条件付きデバッグや動的な状況での問題解決に特に有用である。また、起動されるデバッグモードでは、変数の検査や変更、実行フローの制御などの操作が可能となる。

## 基本的な使い方

### シンプルなデバッグの一時停止

最も基本的な使い方は、スクリプト内にWait-Debuggerを配置して実行を一時停止することである。

```powershell
function Test-WaitDebugger {
    $i = 1
    $j = 2
    
    # ここでデバッガーが起動する
    Wait-Debugger
    
    $result = $i + $j
    return $result
}

Test-WaitDebugger
```

### 条件付きデバッグ

特定の条件が満たされた場合にのみデバッガーを起動することもできる。

```powershell
function Process-Items {
    param (
        [int[]]$Items
    )
    
    foreach ($item in $Items) {
        # 値が10より大きい場合にのみデバッガーを起動
        if ($item -gt 10) {
            Write-Host "大きな値を検出: $item"
            Wait-Debugger
        }
        
        $processedItem = $item * 2
        Write-Host "処理結果: $processedItem"
    }
}

Process-Items -Items @(5, 12, 3, 15)
```

### PSBreakpointと組み合わせて使用

PSBreakpointと組み合わせることで、より柔軟なデバッグが可能になる。

```powershell
function Complex-Operation {
    param($data)
    
    # デバッグモードが有効な場合のみデバッガーを起動
    if ($PSDebugContext) {
        Wait-Debugger
    }
    
    # 処理ロジック
    $result = $data | ForEach-Object { $_ * 2 }
    return $result
}

# デバッグモードを有効にして関数を呼び出す
Set-PSBreakpoint -Command Complex-Operation
Complex-Operation -data @(1, 2, 3)
```

## デバッグモードの基本操作

デバッガーが起動すると、`[DBG]:`というプロンプトが表示され、以下のコマンドを使用できる。

| コマンド | 別名 | 説明 |
|---------|------|------|
| step | s | 1行実行（関数やスクリプトの中に入る） |
| stepover | v | 1行実行（関数やスクリプトの中には入らない） |
| stepout | o | 現在の関数やスクリプトから抜ける |
| continue | c | 実行を再開（次のブレークポイントまで） |
| quit | q | デバッグを終了しスクリプトを停止 |
| Get-PSCallStack | k | コールスタックを表示 |
| list | l | 現在のスクリプトのソースを表示 |
| help | ? または h | ヘルプを表示 |

### 変数の検査と変更

デバッグモードでは、変数の値を検査したり、変更したりすることができる。

```powershell
[DBG]: PS C:\> $i  # 変数の値を表示
1

[DBG]: PS C:\> $i = 10  # 変数の値を変更
[DBG]: PS C:\> $i
10
```

### デバッグモードでソースを表示

`list`コマンドを使用して、現在のスクリプトソースを表示できる。

```powershell
[DBG]: PS C:\> list

    1: function Test-WaitDebugger {
    2:     $i = 1
    3:     $j = 2
    4:     
    5:     # ここでデバッガーが起動する
    6:*    Wait-Debugger
    7:     
    8:     $result = $i + $j
    9:     return $result
   10: }
```

## 応用的な使い方

### エラーハンドリングでのデバッグ

エラーが発生した場合にのみデバッガーを起動する例である。

```powershell
function Divide-Numbers {
    param (
        [int]$Numerator,
        [int]$Denominator
    )
    
    try {
        $result = $Numerator / $Denominator
        return $result
    } catch {
        Write-Host "エラーが発生しました: $_"
        Wait-Debugger
        # エラー後の処理
        return $null
    }
}

Divide-Numbers -Numerator 10 -Denominator 0
```

### 複雑なステップ実行の制御

複雑な入れ子構造の関数でのデバッグ制御例。

```powershell
function Main-Function {
    param($input)
    
    Write-Host "メイン関数開始"
    Wait-Debugger
    
    # 副関数を呼び出す
    $result = Sub-Function -data $input
    
    Write-Host "メイン関数終了"
    return $result
}

function Sub-Function {
    param($data)
    
    Write-Host "副関数開始"
    
    $processed = $data * 2
    
    Write-Host "副関数終了"
    return $processed
}

# 関数を実行
Main-Function -input 5
```

このスクリプトでは、デバッグモードで以下のような操作が可能：

1. `stepover` (v) - Sub-Function内部に入らずに実行
2. `step` (s) - Sub-Function内部に入って詳細にデバッグ
3. `stepout` (o) - 現在の関数から抜ける
4. `list` (l) - 現在位置のコードを表示

## ハンズオン：Wait-Debuggerとデバッグモードの実践

以下はWait-Debuggerコマンドとデバッグモードの実践例である。

1. **基本的なデバッグ停止とデバッグモード操作**

```powershell
function Test-BasicDebugger {
    $name = "PowerShell"
    $version = 7.3
    
    # デバッガーを起動
    Wait-Debugger
    
    $message = "Hello from $name $version!"
    return $message
}

# 関数を実行
Test-BasicDebugger
```

デバッガーが起動すると、以下のような出力が表示される：

```
Entering debug mode. Use h or ? for help.

Hit Debug breakpoint
At line:6 char:5
+     Wait-Debugger
+     ~~~~~~~~~~~~

[DBG]: PS C:\>
```

デバッグモードの操作例：

```
# 変数の値を確認
[DBG]: PS C:\> $name
PowerShell

[DBG]: PS C:\> $version
7.3

# 変数の値を変更
[DBG]: PS C:\> $version = 7.4
[DBG]: PS C:\> $version
7.4

# ヘルプを表示
[DBG]: PS C:\> h
s, step                  Single step (step into functions, scripts, etc.)
v, stepover              Step to next statement (step over functions, scripts, etc.)
o, stepout               Step out of the current function, script, etc.
c, continue              Continue operation
q, quit                  Stop operation and exit the debugger
k, Get-PSCallStack       Display call stack
l, list                  List source code for the current script
<enter>                  Repeat last command if it was step (s), stepover (v) or list (l)
?, h                     Displays this help

# ソースを表示
[DBG]: PS C:\> l
    1: function Test-BasicDebugger {
    2:     $name = "PowerShell"
    3:     $version = 7.3
    4:     
    5:     # デバッガーを起動
    6:*    Wait-Debugger
    7:     
    8:     $message = "Hello from $name $version!"
    9:     return $message
   10: }

# ステップ実行
[DBG]: PS C:\> s
At line:8 char:5
+     $message = "Hello from $name $version!"
+     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 次のステップを実行
[DBG]: PS C:\> s
At line:9 char:5
+     return $message
+     ~~~~~~~~~~~~~~~

# 変数の内容を確認
[DBG]: PS C:\> $message
Hello from PowerShell 7.4

# 実行を継続
[DBG]: PS C:\> c

Hello from PowerShell 7.4
```

2. **コールスタックの確認と複数のデバッグポイント**

```powershell
function Outer-Function {
    param($value)
    Write-Host "外部関数開始"
    $doubledValue = Inner-Function -data $value
    Write-Host "外部関数終了"
    return $doubledValue
}

function Inner-Function {
    param($data)
    Write-Host "内部関数開始"
    Wait-Debugger  # 最初のデバッグポイント
    $result = $data * 2
    Wait-Debugger  # 2番目のデバッグポイント
    Write-Host "内部関数終了"
    return $result
}

# 関数を実行
Outer-Function -value 10
```

デバッグセッションの操作例：

```
外部関数開始
内部関数開始

Entering debug mode. Use h or ? for help.

Hit Debug breakpoint
At line:11 char:5
+     Wait-Debugger  # 最初のデバッグポイント
+     ~~~~~~~~~~~~

# コールスタックを確認
[DBG]: PS C:\> k

Command                                 Arguments
-------                                 ---------
<ScriptBlock>                           {}
Inner-Function                          {-data 10}
Outer-Function                          {-value 10}
<ScriptBlock>                           {}

# 内部関数の変数を確認
[DBG]: PS C:\> $data
10

# 外部関数のスコープにアクセス
[DBG]: PS C:\> $value
# 別のスコープのため直接アクセスできない

# ステップオーバーで次の行へ
[DBG]: PS C:\> v
At line:12 char:5
+     $result = $data * 2
+     ~~~~~~~~~~~~~~~~~~

# もう一度ステップオーバー
[DBG]: PS C:\> v
At line:13 char:5
+     Wait-Debugger  # 2番目のデバッグポイント
+     ~~~~~~~~~~~~

# 計算結果を確認
[DBG]: PS C:\> $result
20

# 実行を継続
[DBG]: PS C:\> c

Entering debug mode. Use h or ? for help.

Hit Debug breakpoint
At line:13 char:5
+     Wait-Debugger  # 2番目のデバッグポイント
+     ~~~~~~~~~~~~

# 結果を変更してみる
[DBG]: PS C:\> $result = 50

# 実行を継続
[DBG]: PS C:\> c
内部関数終了
外部関数終了
50
```

3. **条件付きデバッグとエラーハンドリング**

```powershell
function Process-DataItem {
    param(
        [string]$Item,
        [int]$Index
    )
    
    Write-Host "処理中: [$Index] $Item"
    
    try {
        # 特定の条件でデバッガーを起動
        if ($Item -eq "error") {
            throw "エラーアイテムを検出"
        }
        
        if ($Index -eq 2) {
            Write-Host "特定のインデックスを検出"
            Wait-Debugger
        }
        
        return $Item.ToUpper()
    }
    catch {
        Write-Host "エラーが発生: $_" -ForegroundColor Red
        Wait-Debugger
        return $null
    }
}

# テストデータで関数を呼び出す
$items = @("apple", "banana", "cherry", "error", "grape")
$results = @()

for ($i = 0; $i -lt $items.Count; $i++) {
    $result = Process-DataItem -Item $items[$i] -Index $i
    if ($result) {
        $results += $result
    }
}

Write-Host "結果: $results"
```

このようなスクリプトでは、インデックス2（cherry）とエラーアイテムでデバッガーが起動する。デバッグモードでは、さまざまなコマンドを試すことができる：

```
処理中: [0] apple
処理中: [1] banana
処理中: [2] cherry
特定のインデックスを検出

Entering debug mode. Use h or ? for help.

Hit Debug breakpoint
At line:15 char:13
+             Wait-Debugger
+             ~~~~~~~~~~~~

# 現在の変数を検査
[DBG]: PS C:\> $Item
cherry
[DBG]: PS C:\> $Index
2

# ループカウンターにアクセス
[DBG]: PS C:\> $i
2

# 全アイテムリストを確認
[DBG]: PS C:\> $items
apple
banana
cherry
error
grape

# 収集した結果を確認
[DBG]: PS C:\> $results
APPLE
BANANA

# 実行を継続
[DBG]: PS C:\> c
処理中: [3] error
エラーが発生: エラーアイテムを検出

Entering debug mode. Use h or ? for help.

Hit Debug breakpoint
At line:20 char:9
+         Wait-Debugger
+         ~~~~~~~~~~~~

# エラー情報を確認
[DBG]: PS C:\> $_
エラーアイテムを検出

# エラー処理用の修正を試す
[DBG]: PS C:\> return "ERROR_FIXED"
# 直接returnはできないが、変数は変更可能

# 実行を継続
[DBG]: PS C:\> c
処理中: [4] grape
結果: APPLE BANANA CHERRY GRAPE
```

## 対応PowerShellバージョン

Wait-Debuggerコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 5.0以降
- PowerShell Core 6.0以降
- PowerShell 7.0以降

PowerShell 4.0以前では、このコマンドレットは利用できないため、代わりに`Set-PSBreakpoint`コマンドレットを使用する必要がある。

## 参考サイト

- [Microsoft公式ドキュメント: Wait-Debugger](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.utility/wait-debugger)
- [Microsoft Learn - PowerShellのデバッグ機能](https://docs.microsoft.com/ja-jp/powershell/scripting/learn/deep-dives/everything-about-debugging)
- [PowerShell.org - Advanced Debugging Techniques](https://powershell.org/2020/05/advanced-debugging-techniques-in-powershell/)
- [Microsoft Learn - PowerShellデバッガーの使用方法](https://docs.microsoft.com/ja-jp/powershell/scripting/windows-powershell/ise/how-to-debug-scripts-in-windows-powershell-ise)
- [SS64.com PowerShell Commands - Wait-Debugger](https://ss64.com/ps/wait-debugger.html)