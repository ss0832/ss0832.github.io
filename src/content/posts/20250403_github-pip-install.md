---
title:  【Python】GitHubリポジトリから自作モジュールをpipでインストールする方法
published: 2025-04-03
description: "自作Pythonモジュールをpipコマンドでお手軽にGitHubリポジトリからインストールする方法"
tags: [python]
category: Programming
draft: false
---
最終更新：2025-04-03

## 概要

Pythonのパッケージ管理ツール「pip」は、PyPI（Python Package Index）だけでなく、GitHubなどのバージョン管理システムに公開されているコードもインストールできる。本記事では、PyPIに登録せずに、GitHubリポジトリから直接自作モジュールをインストールする方法を解説する。


## GitHubリポジトリからインストールする基本形式

GitHubリポジトリからパッケージをインストールするには、以下の形式を使用する：

```bash
pip install git+https://github.com/ユーザー名/リポジトリ名.git
```

例えば：
```bash
pip install git+https://github.com/ss0832/ase_afir_for_pip.git
```

これだけで、リポジトリに含まれるパッケージが標準のPythonパッケージと同じようにインストールされる。

## リポジトリの構成要件

GitHubからインストール可能なパッケージには、最低限以下のファイル構成が必要である：

```
リポジトリ/
├── パッケージ名/
│   ├── __init__.py
│   └── (その他のPythonファイル)
└── setup.py
```

特に`setup.py`は必須で、以下のように記述する：

```python
from setuptools import setup, find_packages

setup(
    name="パッケージ名",
    version="0.1.0",
    description="簡単な説明",
    author="作者名",
    packages=find_packages(),
    install_requires=[
        "必要な依存パッケージ",
    ],
)
```

## インストールの応用例

### 特定のブランチからインストール

```bash
pip install git+https://github.com/ユーザー名/リポジトリ名.git@ブランチ名
```

例：
```bash
pip install git+https://github.com/ss0832/ase_afir_for_pip.git@develop
```

### 特定のタグからインストール

```bash
pip install git+https://github.com/ユーザー名/リポジトリ名.git@v1.0.0
```

### 特定のコミットからインストール

```bash
pip install git+https://github.com/ユーザー名/リポジトリ名.git@コミットハッシュ
```

### SSHプロトコルを使用する場合

```bash
pip install git+ssh://git@github.com/ユーザー名/リポジトリ名.git
```

### サブディレクトリにあるパッケージをインストール

リポジトリのサブディレクトリにパッケージがある場合：

```bash
pip install git+https://github.com/ユーザー名/リポジトリ名.git#subdirectory=サブディレクトリ名
```

### 開発モードでインストール

コードを編集しながら開発する場合は、`-e`オプションを使用する：

```bash
pip install -e git+https://github.com/ユーザー名/リポジトリ名.git
```

## requirements.txtでの指定

`requirements.txt`ファイル内でも同様に指定できる：

```
numpy>=1.19.0
git+https://github.com/ユーザー名/リポジトリ名.git@ブランチ名
```

## ハンズオン：ASE_AFIRパッケージのセットアップとGitHubからのインストール

ここでは、実際にASE_AFIRパッケージをGitHubで公開し、そこからインストールする手順を示す。

### 1. プロジェクト構成の作成

```
ase_afir_for_pip/
├── ase_afir/
│   ├── __init__.py
|   ├── __version__.py
│   ├── (プログラム).py
├── setup.py
└── README.md

```

### 2. setup.pyの作成

```python
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="ase_afir",
    version="0.1.0",
    author="ss0832",
    author_email="your.email@example.com",
    description="ASE implementation of AFIR method for chemical reaction path search",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ss0832/ase_afir_for_pip",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.19.0",
        "ase>=3.20.0",
    ],
)
```

### 3. __init__.py及び__version__.pyの作成

`ase_afir/__init__.py`に以下を記述する：

```python
from .__version__ import __version__

from ase_afir import afirpot4ase
```

`ase_afir/__version__.py`に以下を記述する：
```python
VERSION = (0, 1)
__version__ = '.'.join(map(str, VERSION))
```

### 4. GitHubにプッシュする


```bash
# リポジトリを初期化
git init
git add .
git commit -m "Initial commit"

# GitHubリポジトリを作成し、リモートとして追加
git remote add origin https://github.com/ss0832/ase_afir_for_pip.git

# リポジトリにプッシュ
git push -u origin main
```

### 5. GitHubからインストール

リポジトリを公開した後、以下のコマンドでインストールできる：

```bash
# 最新の開発版をインストール
pip install git+https://github.com/ss0832/ase_afir_for_pip.git

# または特定のブランチ/タグをインストール
pip install git+https://github.com/ss0832/ase_afir_for_pip.git@develop
```

### 6. インストールの確認とパッケージの使用

```python
# Pythonを起動
python

# パッケージのインポートを確認
>>> import ase_afir
>>> print(ase_afir.__version__)
'0.1'

# 機能を使用
>>> import ase_afir
>>> test_class = ase_afir.afirpot4ase.ASE_AFIRCalculator()
```


## プライベートリポジトリからのインストール

プライベートリポジトリからインストールする場合は、認証が必要になる：

### HTTPSを使用する場合

```bash
pip install git+https://ユーザー名:パスワード@github.com/ユーザー名/プライベートリポジトリ.git
```

またはPersonal Access Tokenを使用：

```bash
pip install git+https://ユーザー名:個人アクセストークン@github.com/ユーザー名/プライベートリポジトリ.git
```

### SSHを使用する場合（推奨）

```bash
pip install git+ssh://git@github.com/ユーザー名/プライベートリポジトリ.git
```

この場合は事前にSSHキーの設定が必要。

## 依存関係の指定方法

パッケージが他のGitHubリポジトリに依存している場合、`setup.py`の`install_requires`または`dependency_links`で指定できる：

```python
setup(
    # ...
    install_requires=[
        "numpy>=1.19.0",
        "ase>=3.20.0",
        "other-package",
    ],
    dependency_links=[
        "git+https://github.com/user/other-package.git@v1.0.0#egg=other-package-1.0.0",
    ],
)
```

## まとめ

GitHubリポジトリから直接パッケージをインストールする方法は、以下のような場合に特に有用である：

1. パッケージをPyPIに公開する前のテスト
2. 組織内限定または特定ユーザー向けのパッケージ配布
3. 開発中の最新版を試したい場合
4. フォークした修正版を使いたい場合

この方法により、PyPIに登録しなくても、Pythonモジュールを簡単に配布・共有できる。GitHubのバージョン管理機能を活用することで、様々なバージョンやブランチからインストールすることも可能である。

## 参考サイト

- [Pythonオリジナルライブラリをpipで公開する方法](https://irukanobox.blogspot.com/2021/05/pythonpip.html)
- [pip公式ドキュメント: VCSサポート](https://pip.pypa.io/en/stable/topics/vcs-support/)
- [GitHubからpipでインストールする方法](https://pip.pypa.io/en/stable/reference/pip_install/#git)
- [GitHub: ss0832/ase_afir_for_pip](https://github.com/ss0832/ase_afir_for_pip)