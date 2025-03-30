---
title:  【PowerShell】Get-WmiObjectコマンドでシステム情報を取得する方法
published: 2025-03-30
description: "PowerShellのGet-WmiObjectコマンドレットを使用してシステム情報を取得する基本的な方法"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Get-WmiObject`はPowerShellにおいて、Windows Management Instrumentation（WMI）を使用してローカルおよびリモートのコンピュータからシステム情報を取得するためのコマンドレットである。このコマンドレットによりオペレーティングシステム、ハードウェア構成、ソフトウェア設定、およびその他のシステム情報に関するデータにアクセスできる。ユーザーはWMIクラスを指定することで、特定のハードウェアコンポーネントやシステム設定に関する詳細情報を取得し、システムの監視、管理、トラブルシューティングを効率的に行うことが可能となる。なお、PowerShell 3.0以降では`Get-CimInstance`コマンドレットも提供されており、より新しい環境では代替として使用されることが多い。

## 基本的な使い方

### 特定のWMIクラスの情報を取得

特定のWMIクラスに関する情報を取得するには、`-Class`パラメータを使用する。

```powershell
Get-WmiObject -Class Win32_OperatingSystem
```

### リモートコンピュータからの情報取得

リモートコンピュータから情報を取得するには、`-ComputerName`パラメータを使用する。

```powershell
Get-WmiObject -Class Win32_BIOS -ComputerName "Server01"
```

### 特定プロパティの選択

取得したオブジェクトから特定のプロパティのみを選択するには、`Select-Object`と組み合わせる。

```powershell
Get-WmiObject -Class Win32_LogicalDisk -Filter "DriveType=3" | Select-Object DeviceID, Size, FreeSpace
```

### フィルターの使用

特定の条件に一致するオブジェクトのみを取得するには、`-Filter`パラメータを使用する。

```powershell
Get-WmiObject -Class Win32_Process -Filter "Name='chrome.exe'"
```

## 応用的な使い方

### WMI名前空間の探索

デフォルト以外の名前空間内のWMIクラスを取得するには、`-Namespace`パラメータを使用する。

```powershell
Get-WmiObject -Class __NAMESPACE -Namespace "root"
Get-WmiObject -Class Win32_PerfFormattedData_PerfOS_Memory -Namespace "root\cimv2"
```

### 認証情報を使用したリモート接続

別のユーザーアカウントの認証情報を使用してリモートコンピュータに接続するには、`-Credential`パラメータを使用する。

```powershell
$cred = Get-Credential
Get-WmiObject -Class Win32_ComputerSystem -ComputerName "Server01" -Credential $cred
```

### WQLクエリの実行

より複雑なクエリを実行するには、WMI Query Language（WQL）を使用する。

```powershell
Get-WmiObject -Query "SELECT * FROM Win32_Process WHERE ExecutablePath LIKE '%system32%'"
```

### インスタンスメソッドの呼び出し

WMIオブジェクトのメソッドを呼び出して、システムに対してアクションを実行する。

```powershell
$proc = Get-WmiObject -Class Win32_Process -Filter "Name='notepad.exe'"
$proc.Terminate()
```

## ハンズオン：Get-WmiObjectの実践

以下はGet-WmiObjectコマンドの実践例とその出力である。

1. **オペレーティングシステムの情報を取得**

```powershell
# OSの基本情報を取得
Get-WmiObject -Class Win32_OperatingSystem | Select-Object Caption, Version, BuildNumber, OSArchitecture
```

出力例：
```
Caption                   Version    BuildNumber OSArchitecture
-------                   -------    ----------- --------------
Microsoft Windows 11 Pro  10.0.22621 22621       64-bit
```

**出力の見方**: この出力では、Windowsのエディション名（Caption）、バージョン番号（Version）、ビルド番号（BuildNumber）、OSのアーキテクチャ（OSArchitecture）が表示されている。この例ではWindows 11 Pro、バージョン10.0.22621、64ビット版がインストールされていることがわかる。

2. **コンピュータのシステム情報を取得**

```powershell
# コンピュータ名、メーカー、モデルの情報を取得
Get-WmiObject -Class Win32_ComputerSystem | Select-Object Name, Manufacturer, Model
```

出力例：
```
Name        Manufacturer      Model
----        ------------      -----
DESKTOP-ABC Dell Inc.         Latitude 5420
```

**出力の見方**: この出力からはコンピュータのホスト名（Name）、製造メーカー名（Manufacturer）、モデル名（Model）がわかる。この例では、「DESKTOP-ABC」というホスト名のDell製Latitude 5420というモデルのコンピュータであることを示している。

3. **BIOSの情報を取得**

```powershell
# BIOSの情報を取得
Get-WmiObject -Class Win32_BIOS | Select-Object Manufacturer, Version, SerialNumber
```

出力例：
```
Manufacturer Version      SerialNumber
------------ -------      ------------
Dell Inc.    1.8.2        ABC123XYZ
```

**出力の見方**: この出力ではBIOSの製造元（Manufacturer）、バージョン（Version）、シリアル番号（SerialNumber）が表示されている。デバイスの識別やBIOSのアップデートが必要かどうかを判断する際に役立つ情報である。

4. **物理メモリの情報を取得**

```powershell
# 物理メモリの容量を計算（バイト単位からGB単位に変換）
$totalRAM = [Math]::Round((Get-WmiObject -Class Win32_ComputerSystem).TotalPhysicalMemory / 1GB, 2)
Write-Host "搭載物理メモリ容量: $totalRAM GB"

# メモリモジュールの詳細情報を取得
Get-WmiObject -Class Win32_PhysicalMemory | 
    Select-Object @{Name="容量(GB)";Expression={[Math]::Round($_.Capacity / 1GB, 2)}}, 
                 @{Name="メモリタイプ";Expression={
                    switch($_.MemoryType) {
                        21 {"DDR2"}
                        24 {"DDR3"}
                        26 {"DDR4"}
                        default {"Unknown"}
                    }
                 }},
                 Speed, DeviceLocator
```

出力例：
```
搭載物理メモリ容量: 16 GB

容量(GB) メモリタイプ Speed DeviceLocator
-------- ---------- ----- -------------
      8 DDR4        3200 DIMM 1
      8 DDR4        3200 DIMM 2
```

**出力の見方**: この出力では2つの情報が表示されている。まず、コンピュータに搭載されている物理メモリ（RAM）の総容量が16GBであることが示されている。次に、個々のメモリモジュールの詳細が表示され、このコンピュータには8GB容量のDDR4メモリが2枚あり、それぞれのメモリ速度は3200MHzで、「DIMM 1」と「DIMM 2」というスロットに装着されていることがわかる。

5. **ディスクドライブの情報を取得**

```powershell
# 論理ディスクの情報を取得
Get-WmiObject -Class Win32_LogicalDisk -Filter "DriveType=3" |
    Select-Object DeviceID, 
                 @{Name="容量(GB)";Expression={[Math]::Round($_.Size / 1GB, 2)}}, 
                 @{Name="空き容量(GB)";Expression={[Math]::Round($_.FreeSpace / 1GB, 2)}},
                 @{Name="使用率(%)";Expression={[Math]::Round(($_.Size - $_.FreeSpace) / $_.Size * 100, 1)}}
```

出力例：
```
DeviceID 容量(GB) 空き容量(GB) 使用率(%)
-------- -------- ------------ ---------
C:         476.92       322.58      32.4
D:         523.08       453.29      13.3
```

**出力の見方**: この出力は固定ディスク（DriveType=3はローカルディスクを意味する）についての情報を表示している。C:ドライブは総容量約477GBで、そのうち約323GB（67.6%）が空き容量であり、使用率は32.4%である。D:ドライブは総容量約523GBで、そのうち約453GB（86.7%）が空き容量であり、使用率は13.3%である。この情報からディスク容量の使用状況を一目で把握できる。

6. **実行中のプロセスの情報を取得**

```powershell
# メモリ使用量上位5件のプロセスを取得
Get-WmiObject -Class Win32_Process |
    Select-Object Name, ProcessId, 
                 @{Name="メモリ使用量(MB)";Expression={[Math]::Round($_.WorkingSetSize / 1MB, 2)}} |
    Sort-Object -Property "メモリ使用量(MB)" -Descending |
    Select-Object -First 5
```

出力例：
```
Name                  ProcessId メモリ使用量(MB)
----                  --------- ---------------
chrome.exe                 4568          285.72
explorer.exe               2340          156.31
MicrosoftEdge.exe          3928          125.68
powershell.exe             5824           98.45
Teams.exe                  3568           85.23
```

**出力の見方**: この出力はメモリ使用量の多い上位5つのプロセスを表示している。各行にはプロセス名（Name）、プロセスID（ProcessId）、メモリ使用量（MB単位）が表示されている。この例では、Google Chrome（chrome.exe）がもっとも多くのメモリを使用しており約286MBを消費している。次いでexplorer.exe（Windowsエクスプローラー）が約156MB、Microsoft Edgeブラウザが約126MB使用していることがわかる。システムのメモリ使用状況を監視する際に役立つ情報である。

7. **ネットワークアダプタの情報を取得**

```powershell
# 物理ネットワークアダプタの情報を取得
Get-WmiObject -Class Win32_NetworkAdapter | 
    Where-Object { $_.PhysicalAdapter -eq $true -and $_.NetEnabled -eq $true } |
    Select-Object Name, MACAddress, Speed, AdapterType
```

出力例：
```
Name                        MACAddress        Speed       AdapterType
----                        ----------        -----       -----------
Intel(R) Wi-Fi 6 AX201      00:11:22:33:44:55 1200000000  イーサネット802.3
Realtek PCIe GBE Family Cont56:67:78:89:9A:AB 1000000000  イーサネット802.3
```

**出力の見方**: この出力は有効な物理ネットワークアダプタの情報を表示している。コンピュータには2つのネットワークアダプタがあり、1つはIntel製のWi-Fi 6アダプタで、MACアドレスは「00:11:22:33:44:55」、通信速度は1.2Gbps（1200000000ビット/秒）である。もう1つはRealtek製の有線LANアダプタで、MACアドレスは「56:67:78:89:9A:AB」、通信速度は1Gbps（1000000000ビット/秒）である。どちらも「イーサネット802.3」規格に準拠している。

8. **サービスの情報を取得**

```powershell
# 自動起動設定かつ実行中のサービスを取得
Get-WmiObject -Class Win32_Service -Filter "StartMode='Auto' AND State='Running'" |
    Select-Object Name, DisplayName, State, StartMode |
    Sort-Object -Property Name |
    Select-Object -First 5
```

出力例：
```
Name       DisplayName                              State   StartMode
----       -----------                              -----   ---------
AudioEndpo Windowsオーディオエンドポイントビルダー   Running Auto
BFE        Base Filtering Engine                    Running Auto
BrokerInfr Background Tasks Infrastructure Service   Running Auto
CDPSvc     Connected Devices Platform Service       Running Auto
CoreMessa... CoreMessaging                          Running Auto
```

**出力の見方**: この出力はシステム上で「自動起動」（Auto）に設定されており、かつ現在「実行中」（Running）のサービスを、名前のアルファベット順に並べた最初の5つを表示している。各行にはサービスの内部名（Name）、表示名（DisplayName）、現在の状態（State）、起動モード（StartMode）が表示されている。この情報はシステムの起動時に自動的に開始されるサービスを確認する際に役立つ。

9. **WMIイベントを監視**

```powershell
# 新しいプロセスの作成を監視する（Ctrl+Cで終了）
Write-Host "新しいプロセスを監視しています。(Ctrl+Cで停止)..." -ForegroundColor Yellow
try {
    Register-WmiEvent -Query "SELECT * FROM __InstanceCreationEvent WITHIN 1 WHERE TargetInstance ISA 'Win32_Process'" -SourceIdentifier "ProcessStarted" -Action {
        $process = $Event.SourceEventArgs.NewEvent.TargetInstance
        Write-Host "新しいプロセスが開始されました: $($process.Name) (PID: $($process.ProcessId))" -ForegroundColor Cyan
    }
    
    # イベントを待機（10秒間のみ、実際の使用ではCtrl+Cで終了）
    Start-Sleep -Seconds 10
    
} finally {
    # イベント監視を解除
    Unregister-Event -SourceIdentifier "ProcessStarted" -ErrorAction SilentlyContinue
    Write-Host "プロセス監視を終了しました。" -ForegroundColor Yellow
}
```

出力例：
```
新しいプロセスを監視しています。(Ctrl+Cで停止)...
新しいプロセスが開始されました: notepad.exe (PID: 6784)
新しいプロセスが開始されました: calc.exe (PID: 7120)
プロセス監視を終了しました。
```

**出力の見方**: この出力はWMIイベント監視の結果を表示している。最初に「新しいプロセスを監視しています」というメッセージが表示され、その後システム上で新しく起動されたプロセスが検出されるたびに通知が表示される。この例では監視中にメモ帳（notepad.exe）と電卓（calc.exe）が起動され、それぞれのプロセスIDとともに表示されている。最後に監視が終了したことを示すメッセージが表示されている。この機能はリアルタイムでのシステム活動監視に役立つ。

10. **インストールされているアプリケーションの一覧取得**

```powershell
# インストールされているアプリケーションの一覧を取得
Get-WmiObject -Class Win32_Product | 
    Select-Object Name, Version, Vendor |
    Sort-Object -Property Name |
    Select-Object -First 5
```

出力例：
```
Name                      Version    Vendor
----                      -------    ------
7-Zip                     21.07      Igor Pavlov
Adobe Acrobat Reader DC   23.003.20244 Adobe Systems Incorporated
Brave                     120.1.60.4 Brave Software Inc.
Git                       2.42.0     The Git Development Community
Google Chrome             123.0.6312.58 Google LLC
```

**出力の見方**: この出力はWindows Installerによってインストールされたアプリケーションの一覧を表示している。アルファベット順に並べられた最初の5つのアプリケーションが表示され、各行にはアプリケーション名（Name）、バージョン（Version）、開発元（Vendor）の情報が含まれている。この例では7-Zip、Adobe Acrobat Reader DC、Brave、Git、Google Chromeがインストールされており、それぞれのバージョン番号も確認できる。この情報はインストールされているソフトウェアの管理やバージョン確認に役立つ。

## 対応PowerShellバージョン

Get-WmiObjectコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 1.0～5.1

注意点：
- PowerShell Core 6.0以降では、Get-WmiObjectは非推奨となり、Get-CimInstanceに置き換えられている
- Windows PowerShell 3.0以降ではGet-CimInstanceの使用が推奨されている
- Windowsでのみサポートされており、Linux/macOSなどのクロスプラットフォーム環境では使用できない

PowerShell 6.0以降で同等の機能を使用するには以下のように置き換える：
- `Get-WmiObject -Class Win32_Process` ➡ `Get-CimInstance -ClassName Win32_Process`

## 参考サイト

- [Microsoft公式ドキュメント: Get-WmiObject](https://docs.microsoft.com/ja-jp/powershell/module/microsoft.powershell.management/get-wmiobject)
- [Windows Management Instrumentation](https://docs.microsoft.com/ja-jp/windows/win32/wmisdk/wmi-start-page)
- [PowerShell.org - WMI入門ガイド](https://powershell.org/2013/08/wmi-query-basics/)
- [Microsoft Tech Community - PowerShellでのWMI/CIM使用方法](https://techcommunity.microsoft.com/t5/itops-talk-blog/powershell-basics-how-to-use-wmi-with-powershell/ba-p/885715)
- [SS64.com PowerShell Commands - Get-WmiObject](https://ss64.com/ps/get-wmiobject.html)