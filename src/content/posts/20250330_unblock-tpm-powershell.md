---
title:  【PowerShell】Unblock-TpmコマンドでロックされたTPMをリセットする方法
published: 2025-03-30
description: "PowerShellのUnblock-Tpmコマンドレットを使ったTPM(Trusted Platform Module)のロック解除方法"
tags: [PowerShell]
category: Windows
draft: false
---
最終更新：2025-03-30

## 概要

`Unblock-Tpm`はWindows PowerShellのTrustedPlatformModule(TpmProvider)モジュールに含まれるコマンドレットである。このコマンドレットは、ロックアウト状態になったTPM（Trusted Platform Module）チップをリセットするために使用される。TPMはハードウェアセキュリティチップであり、BitLockerドライブ暗号化やWindows Helloなどの機能で使用されるが、誤ったPINの連続入力などによりロックされることがある。そのような場合に、このコマンドレットを使用して回復することが可能である。TPM所有者認証情報（オーナーパスワード）が必要となる場合があり、管理者権限でのみ実行可能である。

## 基本的な使い方

### TPM状態の確認

TPMの状態を確認するには、まず`Get-Tpm`コマンドレットを使用する。

```powershell
Get-Tpm
```

### 標準的なTPMロック解除

TPMがロックされている場合、以下のコマンドで解除を試みる。

```powershell
Unblock-Tpm
```

### 所有者認証情報を使用したロック解除

TPMの所有者認証情報（オーナーパスワード）がファイルに保存されている場合は、それを使用してロック解除できる。

```powershell
Unblock-Tpm -OwnerAuth (Get-Content C:\TPM\ownerAuth.bin)
```

### TPM所有者認証をバイト配列で指定

TPMの所有者認証情報をバイト配列として直接指定することも可能である。

```powershell
$ownerAuth = [byte[]](0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08)
Unblock-Tpm -OwnerAuth $ownerAuth
```

## 応用的な使い方

### TPMステータス確認と条件付きロック解除

TPMの状態を確認し、ロックされている場合のみロック解除を試みる。

```powershell
$tpmStatus = Get-Tpm
if ($tpmStatus.LockedOut) {
    Write-Host "TPMはロックされています。ロック解除を試みます..."
    Unblock-Tpm
} else {
    Write-Host "TPMはロックされていません。"
}
```

### OwnerAuthファイルの検索とロック解除

システム上の所有者認証ファイルを探して、それを使用したロック解除を試みる例。

```powershell
$ownerAuthPath = "C:\Windows\System32\TPM\TPMOwnerAuth.bin"
if (Test-Path $ownerAuthPath) {
    Write-Host "TPM所有者認証ファイルが見つかりました。ロック解除を試みます..."
    Unblock-Tpm -OwnerAuth (Get-Content $ownerAuthPath -Encoding Byte)
} else {
    Write-Host "TPM所有者認証ファイルが見つかりません。標準的なロック解除を試みます..."
    Unblock-Tpm
}
```

### エラー処理を含んだTPMロック解除

エラーハンドリングを含めたTPMロック解除の例。

```powershell
try {
    Write-Host "TPMのロック解除を試みています..."
    Unblock-Tpm -ErrorAction Stop
    Write-Host "TPMのロック解除に成功しました。"
} catch {
    Write-Host "TPMのロック解除中にエラーが発生しました: $_"
    Write-Host "所有者認証が必要かもしれません。"
}
```

## ハンズオン：Unblock-Tpmの実践

以下はUnblock-Tpmコマンドの実践例とその出力である。これらのコマンドは管理者権限のPowerShellで実行する必要がある。

1. **TPM状態の確認**

```powershell
# TPMの状態を確認
Get-Tpm
```

出力例：
```
TpmPresent                : True
TpmReady                  : True
TpmEnabled                : True
TpmActivated              : True
TpmOwned                  : True
RestartPending            : False
ManufacturerId            : 1398033696
ManufacturerVersion       : 7.63
ManagedAuthLevel          : Full
OwnerAuth                 : 
OwnerClearDisabled        : False
AutoProvisioning          : Enabled
LockedOut                 : True
LockoutHealTime           : 10
LockoutCount              : 3
LockoutMax                : 10
SelfTest                  : {}
```

上記の例では、`LockedOut: True`が表示されており、TPMがロックアウト状態であることを示している。

2. **標準的なロック解除の試行**

```powershell
# TPMのロック解除を試みる
Unblock-Tpm
```

出力例（成功時）:
```
# 成功時は出力がなく、コマンドプロンプトに戻る
```

出力例（失敗時）:
```
Unblock-Tpm : TPMがロックアウト状態です。所有者認証が必要です。
発生場所 行:1 文字:1
+ Unblock-Tpm
+ ~~~~~~~~~~~
    + CategoryInfo          : NotSpecified: (:) [Unblock-Tpm], Exception
    + FullyQualifiedErrorId : System.Exception,Microsoft.Tpm.Commands.UnblockTpmCommand
```

3. **ロック解除後のTPM状態の確認**

```powershell
# ロック解除後のTPM状態を確認
Get-Tpm
```

出力例（ロック解除成功後）:
```
TpmPresent                : True
TpmReady                  : True
TpmEnabled                : True
TpmActivated              : True
TpmOwned                  : True
RestartPending            : False
ManufacturerId            : 1398033696
ManufacturerVersion       : 7.63
ManagedAuthLevel          : Full
OwnerAuth                 : 
OwnerClearDisabled        : False
AutoProvisioning          : Enabled
LockedOut                 : False
LockoutHealTime           : 10
LockoutCount              : 0
LockoutMax                : 10
SelfTest                  : {}
```

4. **所有者認証情報を使用したロック解除（ファイルから）**

```powershell
# 所有者認証ファイルが存在するか確認
$ownerAuthPath = "$env:SystemRoot\System32\Tpm\TPMOwnerAuth.bin"
if (Test-Path $ownerAuthPath) {
    # 所有者認証ファイルからロック解除
    $ownerAuth = Get-Content -Path $ownerAuthPath -Encoding Byte
    Unblock-Tpm -OwnerAuth $ownerAuth
    Write-Host "所有者認証情報を使用してTPMのロック解除を試みました。"
} else {
    Write-Host "所有者認証ファイルが見つかりません: $ownerAuthPath"
}
```

出力例：
```
所有者認証情報を使用してTPMのロック解除を試みました。
```

または：
```
所有者認証ファイルが見つかりません: C:\Windows\System32\Tpm\TPMOwnerAuth.bin
```

5. **TPM状態の監視と自動ロック解除（スクリプト例）**

```powershell
# TPM状態監視と自動ロック解除スクリプト
$maxRetries = 3
$retryCount = 0
$success = $false

while (-not $success -and $retryCount -lt $maxRetries) {
    $tpm = Get-Tpm
    if ($tpm.LockedOut) {
        Write-Host "TPMがロックされています。ロック解除を試みます... (試行 $($retryCount + 1)/$maxRetries)"
        try {
            Unblock-Tpm -ErrorAction Stop
            $success = $true
            Write-Host "TPMのロック解除に成功しました。"
        } catch {
            Write-Host "TPMのロック解除に失敗しました: $_"
            $retryCount++
            if ($retryCount -lt $maxRetries) {
                Write-Host "5秒後に再試行します..."
                Start-Sleep -Seconds 5
            }
        }
    } else {
        Write-Host "TPMはロックされていません。処理は不要です。"
        $success = $true
    }
}

if (-not $success) {
    Write-Host "最大試行回数に達しました。手動での対応が必要かもしれません。" -ForegroundColor Red
}
```

出力例：
```
TPMがロックされています。ロック解除を試みます... (試行 1/3)
TPMのロック解除に失敗しました: TPMがロックアウト状態です。所有者認証が必要です。
5秒後に再試行します...
TPMがロックされています。ロック解除を試みます... (試行 2/3)
TPMのロック解除に失敗しました: TPMがロックアウト状態です。所有者認証が必要です。
5秒後に再試行します...
TPMがロックされています。ロック解除を試みます... (試行 3/3)
TPMのロック解除に失敗しました: TPMがロックアウト状態です。所有者認証が必要です。
最大試行回数に達しました。手動での対応が必要かもしれません。
```

または（TPMがロックされていない場合）：
```
TPMはロックされていません。処理は不要です。
```

6. **TPM仕様バージョンの確認（情報収集の例）**

```powershell
# TPMのバージョン情報を取得
$tpm = Get-Tpm
$wmiTpm = Get-WmiObject -Namespace root\cimv2\Security\MicrosoftTpm -Class Win32_Tpm

Write-Host "TPM情報:"
Write-Host "ロックアウト状態: $($tpm.LockedOut)"
Write-Host "メーカーID: $($tpm.ManufacturerId)"
Write-Host "メーカーバージョン: $($tpm.ManufacturerVersion)"
if ($wmiTpm) {
    Write-Host "TPM仕様バージョン: $($wmiTpm.SpecVersion)"
    Write-Host "物理的な存在: $($wmiTpm.PhysicalPresenceVersionInfo)"
}
```

出力例：
```
TPM情報:
ロックアウト状態: False
メーカーID: 1398033696
メーカーバージョン: 7.63
TPM仕様バージョン: 2.0
物理的な存在: 1.3
```

## 対応PowerShellバージョン

Unblock-Tpmコマンドレットは以下のPowerShellバージョンで利用可能である：
- Windows PowerShell 4.0以降（Windows 8.1/Server 2012 R2以降）
- Windows PowerShell 5.0および5.1（Windows 10/Server 2016以降で標準搭載）
- PowerShell Core 6.0以降は非対応（TpmProviderモジュールはWindows PowerShellのみ）
- PowerShell 7.0以降も非対応

このコマンドレットを使用するには、以下が必要：
- Windows 8.1/Server 2012 R2以降のオペレーティングシステム
- TPMモジュール搭載のコンピュータ
- 管理者権限でのPowerShell実行
- TrustedPlatformModuleモジュールのインポート（必要に応じて）

PowerShellでTPMプロバイダモジュールを明示的にインポートするには：
```powershell
Import-Module TrustedPlatformModule
```

## 参考サイト

- [Microsoft公式ドキュメント: Unblock-Tpm](https://docs.microsoft.com/ja-jp/powershell/module/trustedplatformmodule/unblock-tpm)
- [Microsoft Docs - TPM基礎](https://docs.microsoft.com/ja-jp/windows/security/information-protection/tpm/trusted-platform-module-overview)
- [Microsoft Docs - TPMのトラブルシューティング](https://docs.microsoft.com/ja-jp/windows/security/information-protection/tpm/tpm-recommendations)
- [TechNet - TPMのロックアウトからの回復](https://social.technet.microsoft.com/wiki/contents/articles/50002.tpm-lockout-and-recovery.aspx)