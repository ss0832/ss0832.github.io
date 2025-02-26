---
title: 【Javascript】 Puppeteerの使い方：ウェブスクレイピングによる画像ファイルとPDFファイルの自動取得ガイド
published: 2025-02-26
description: ""
tags: [Puppeteer]
category: JavaScript
draft: false
---
最終更新：2025-02-26

## 目次

1. [Puppeteerとは](#puppeteerとは)
2. [環境構築](#環境構築)
3. [プロジェクトのセットアップ](#プロジェクトのセットアップ)
4. [Puppeteerの基本的な使い方](#puppeteerの基本的な使い方)
5. [画像とPDFファイルの抽出](#画像とpdfファイルの抽出)
6. [完全なダウンローダープログラム](#完全なダウンローダープログラム)

## Puppeteerとは

Puppeteerは、GoogleのChromeチームによって開発されたNode.jsのライブラリである。このライブラリを使うと、ヘッドレスChrome（画面表示なし）やChromiumを自動的に操作可能である。主な用途は以下の通りである：

- ウェブページのスクリーンショット取得
- SPAやJavaScriptを使用したウェブサイトのテスト自動化
- キーボード入力やフォーム送信などのブラウザ操作の自動化
- ウェブスクレイピング
- PDFの生成

## 環境構築

Puppeteerを使うには、Node.js環境が必要である。

### 1. Node.jsのインストール

[Node.js公式サイト](https://nodejs.org/)から最新のLTS版をダウンロードしてインストールする。

インストールが完了したら、バージョンを確認してNode.jsが使用可能か確認する。：

```bash
node -v
npm -v
```

### 2. プロジェクト理ディレクトリの作成
```bash
mkdir puppeteer-workshop
cd puppeteer-workshop
```

## プロジェクトのセットアップ

### 1. プロジェクトの初期化
```bash
npm init -y
```

### 2. Puppeteerのインストール
```bash
npm install puppeteer
```
### 3. Linux依存関係のインストール
```bash
# Debian/Ubuntu
sudo apt-get update
sudo apt-get install -y ca-certificates fonts-liberation libasound2 libatk-bridge2.0-0 libatk1.0-0 libc6 libcairo2 libcups2 libdbus-1-3 libexpat1 libfontconfig1 libgbm1 libgcc1 libglib2.0-0 libgtk-3-0 libnspr4 libnss3 libpango-1.0-0 libpangocairo-1.0-0 libstdc++6 libx11-6 libx11-xcb1 libxcb1 libxcomposite1 libxcursor1 libxdamage1 libxext6 libxfixes3 libxi6 libxrandr2 libxrender1 libxss1 libxtst6 lsb-release wget xdg-utils

# Fedora/CentOS
sudo dnf install alsa-lib atk cups-libs gtk3 libXcomposite libXcursor libXdamage libXext libXi libXrandr libXScrnSaver libXtst pango at-spi2-atk libXt xorg-x11-server-Xvfb
```


## Puppeteerの基本的な使い方

まずは、基本的なブラウザ操作とスクリーンショットの取得から行う。

### 1. Webサイトのスクリーンショット取得
以下のコードをscreenshot.jsとして保存する。：
```JavaScript
const puppeteer = require('puppeteer');

 (async () => {
  // サンドボックスを無効化してブラウザを起動
  const browser = await puppeteer.launch({
    headless: "new",
    args: ['--no-sandbox', '--disable-setuid-sandbox']
  });
  // 新しいページを開く
  const page = await browser.newPage();

  // URLにアクセス
  await page.goto('https://example.com');

  // スクリーンショットを撮る
  await page.screenshot({ path: 'example.png' });

  // ブラウザを閉じる
  await browser.close();

  console.log('スクリーンショットを保存しました: example.png');
})();
```
以下のコマンドで実行する。：
```bash
node screenshot.js
```

### 2.ブラウザを可視化する

デバッグ時には、ブラウザを実際に表示したい場合がある：
```JavaScript
const puppeteer = require('puppeteer');

(async () => {
  // ヘッドレスモードをオフにしてブラウザを起動
  const browser = await puppeteer.launch({
    headless: false,  // ブラウザウィンドウを表示
    slowMo: 100,      // 各操作間に100msの遅延を追加
    args: ['--no-sandbox', '--disable-setuid-sandbox'],
  });
  
  const page = await browser.newPage();
  await page.goto('https://example.com');
  
  // 3秒間待機
  await new Promise(resolve => setTimeout(resolve, 3000));
  
  await browser.close();
})();
```

### 3.ページ情報の取得
ページのタイトルやテキストを取得する例：

```JavaScript
const puppeteer = require('puppeteer');

(async () => {
  const browser = await puppeteer.launch({
    args: ['--no-sandbox', '--disable-setuid-sandbox'],
  });
  const page = await browser.newPage();
  await page.goto('https://example.com');
  
  // ページタイトルを取得
  const title = await page.title();
  console.log('ページタイトル:', title);
  
  // ページのテキストを取得
  const bodyText = await page.evaluate(() => document.body.innerText);
  console.log('ページ本文:', bodyText);
  
  await browser.close();
})();
```

## 画像とPDFファイルの抽出

### 1. ページ内の画像URLを抽出
```JavaScript
 const puppeteer = require('puppeteer');

(async () => {
  const browser = await puppeteer.launch({
    args: ['--no-sandbox', '--disable-setuid-sandbox'],
  });
  const page = await browser.newPage();
  await page.goto('https://www.yahoo.com/');
  
  // ページ内のすべての画像URLを取得
  const imageUrls = await page.evaluate(() => {
    const images = Array.from(document.querySelectorAll('img'));
    return images.map(img => {
      return {
        src: img.src,
        alt: img.alt,
        width: img.width,
        height: img.height
      };
    });
  });
  
  console.log('見つかった画像:', imageUrls);
  await browser.close();
})();
```
### 2. ページ内のリンクからPDFを見つける
```JavaScript
const puppeteer = require('puppeteer');

(async () => {
  const browser = await puppeteer.launch({
    args: ['--no-sandbox', '--disable-setuid-sandbox'],
  });
  const page = await browser.newPage();
  await page.goto('https://example.com');
  
  // ページ内のPDFへのリンクを検索
  const pdfLinks = await page.evaluate(() => {
    const links = Array.from(document.querySelectorAll('a'));
    return links
      .filter(link => link.href.toLowerCase().endsWith('.pdf'))
      .map(link => {
        return {
          url: link.href,
          text: link.textContent.trim()
        };
      });
  });
  
  console.log('見つかったPDFリンク:', pdfLinks);
  await browser.close();
})();
```

## 完全なダウンロードダウンロード

以下のプログラムは、指定したURLにアクセスして、そのページに含まれるすべての画像とPDFファイルを重複なくダウンロードする。

以下のスクリプトをfile-downloader.jsとして保存する。

```JavaScript
const puppeteer = require('puppeteer');
const fs = require('fs');
const path = require('path');
const https = require('https');
const http = require('http');
const { URL } = require('url');

// ダウンロードディレクトリの設定
const DOWNLOAD_DIR = path.join(__dirname, 'downloads');
const IMAGES_DIR = path.join(DOWNLOAD_DIR, 'images');
const PDFS_DIR = path.join(DOWNLOAD_DIR, 'pdfs');

// ディレクトリが存在しなければ作成
function ensureDirectoryExistence(dirPath) {
  if (!fs.existsSync(dirPath)) {
    fs.mkdirSync(dirPath, { recursive: true });
    console.log(`ディレクトリを作成しました: ${dirPath}`);
  }
}

// ファイル名から拡張子を取得
function getFileExtension(filename) {
  return path.extname(filename).toLowerCase();
}

// URLからファイル名を抽出
function getFilenameFromUrl(urlString) {
  try {
    const url = new URL(urlString);
    const pathname = decodeURIComponent(url.pathname);
    let filename = path.basename(pathname);
    
    // クエリパラメータがある場合はそれも含める
    if (url.search) {
      // 不正なファイル名文字を置換
      const safeQuery = url.search.replace(/[/?<>\\:*|"]/g, '_');
      filename += safeQuery;
    }
    
    return filename;
  } catch (error) {
    console.error('URLの解析に失敗しました:', error);
    return `file_${Date.now()}`;
  }
}

// ファイルのダウンロード関数
async function downloadFile(url, destPath) {
  return new Promise((resolve, reject) => {
    // 既に存在するかチェック
    if (fs.existsSync(destPath)) {
      console.log(`ファイルは既に存在します: ${destPath}`);
      resolve(false);
      return;
    }
    
    // URLがhttpまたはhttpsで始まるか確認
    if (!url.startsWith('http://') && !url.startsWith('https://')) {
      console.error(`サポートされていないURL形式です: ${url}`);
      resolve(false);
      return;
    }
    
    const protocol = url.startsWith('https') ? https : http;
    
    protocol.get(url, (response) => {
      // リダイレクトの処理
      if (response.statusCode === 301 || response.statusCode === 302) {
        const redirectUrl = response.headers.location;
        console.log(`リダイレクト: ${url} -> ${redirectUrl}`);
        downloadFile(redirectUrl, destPath).then(resolve).catch(reject);
        return;
      }
      
      // エラーレスポンスの処理
      if (response.statusCode !== 200) {
        console.error(`ダウンロード失敗: ${url} (Status: ${response.statusCode})`);
        resolve(false);
        return;
      }
      
      // ファイルをディスクに保存
      const fileStream = fs.createWriteStream(destPath);
      response.pipe(fileStream);
      
      fileStream.on('finish', () => {
        fileStream.close();
        console.log(`ダウンロード完了: ${destPath}`);
        resolve(true);
      });
      
      fileStream.on('error', (err) => {
        fs.unlink(destPath, () => {}); // エラー時にファイルを削除
        console.error(`ファイル書き込みエラー: ${err.message}`);
        reject(err);
      });
    }).on('error', (err) => {
      console.error(`リクエストエラー: ${err.message}`);
      reject(err);
    });
  });
}

// URLの正規化
function normalizeUrl(url, baseUrl) {
  try {
    return new URL(url, baseUrl).href;
  } catch (error) {
    console.error(`URLの正規化に失敗しました: ${url}`, error);
    return null;
  }
}

// メイン関数
async function downloadMediaFromUrl(url) {
  ensureDirectoryExistence(DOWNLOAD_DIR);
  ensureDirectoryExistence(IMAGES_DIR);
  ensureDirectoryExistence(PDFS_DIR);
  
  console.log(`処理開始: ${url}`);
  
  const browser = await puppeteer.launch({
    args: ['--no-sandbox', '--disable-setuid-sandbox'],
    headless: "new"  // 新しいヘッドレスモード
  });
  
  try {
    const page = await browser.newPage();
    
    // ユーザーエージェントを設定
    await page.setUserAgent('Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36');
    
    // リクエストを監視してPDFを検出
    const directPdfUrls = new Set();
    page.on('response', async (response) => {
      const url = response.url();
      const contentType = response.headers()['content-type'] || '';
      
      if (contentType.includes('application/pdf') || url.toLowerCase().endsWith('.pdf')) {
        const normalizedUrl = normalizeUrl(url, page.url());
        if (normalizedUrl) directPdfUrls.add(normalizedUrl);
      }
    });
    
    // ページに移動
    await page.goto(url, { waitUntil: 'networkidle0', timeout: 60000 });
    console.log(`ページ読み込み完了: ${page.url()}`);
    
    // ページ内の画像とPDFリンクを収集
    const mediaUrls = await page.evaluate((baseUrl) => {
      const result = {
        images: new Set(),
        pdfs: new Set()
      };
      
      // 画像を収集
      const images = document.querySelectorAll('img');
      images.forEach(img => {
        if (img.src) {
          try {
            const fullUrl = new URL(img.src, baseUrl).href;
            result.images.add(fullUrl);
          } catch (e) {
            console.error('Invalid image URL:', img.src);
          }
        }
      });
      
      // 背景画像を持つ要素を検索
      const allElements = document.querySelectorAll('*');
      allElements.forEach(el => {
        const style = window.getComputedStyle(el);
        const bgImage = style.backgroundImage;
        if (bgImage && bgImage !== 'none') {
          const urlMatch = bgImage.match(/url\(['"]?(.*?)['"]?\)/);
          if (urlMatch && urlMatch[1]) {
            try {
              const fullUrl = new URL(urlMatch[1], baseUrl).href;
              result.images.add(fullUrl);
            } catch (e) {
              console.error('Invalid background image URL:', urlMatch[1]);
            }
          }
        }
      });
      
      // PDFリンクを収集
      const links = document.querySelectorAll('a');
      links.forEach(link => {
        if (link.href) {
          try {
            const href = link.href.toLowerCase();
            if (href.endsWith('.pdf')) {
              const fullUrl = new URL(link.href, baseUrl).href;
              result.pdfs.add(fullUrl);
            }
          } catch (e) {
            console.error('Invalid link URL:', link.href);
          }
        }
      });
      
      return {
        images: [...result.images],
        pdfs: [...result.pdfs]
      };
    }, page.url());
    
    // 画像のダウンロード
    console.log(`${mediaUrls.images.length}個の画像を見つけました`);
    let downloadedImages = 0;
    
    for (const imageUrl of mediaUrls.images) {
      try {
        const filename = getFilenameFromUrl(imageUrl);
        // ファイル名に拡張子がない場合はデフォルトで.jpgを付与
        const ext = getFileExtension(filename) || '.jpg';
        const imageExt = ext.match(/\.(jpg|jpeg|png|gif|webp|svg|avif)/i) ? ext : '.jpg';
        
        // ファイル名の衝突を防ぐためにタイムスタンプを追加
        const safeFilename = `${path.parse(filename).name.substring(0, 100)}_${Date.now()}${imageExt}`;
        const imagePath = path.join(IMAGES_DIR, safeFilename);
        
        const success = await downloadFile(imageUrl, imagePath);
        if (success) downloadedImages++;
      } catch (error) {
        console.error(`画像ダウンロードエラー (${imageUrl}):`, error);
      }
    }
    
    // PDFのダウンロード (ページから見つけたものとレスポンスから直接検出したもの)
    const allPdfUrls = [...new Set([...mediaUrls.pdfs, ...directPdfUrls])];
    console.log(`${allPdfUrls.length}個のPDFを見つけました`);
    let downloadedPdfs = 0;
    
    for (const pdfUrl of allPdfUrls) {
      try {
        const filename = getFilenameFromUrl(pdfUrl);
        const ext = getFileExtension(filename);
        // PDFでない場合はスキップ
        if (ext !== '.pdf') continue;
        
        // ファイル名の衝突を防ぐためにタイムスタンプを追加
        const safeFilename = `${path.parse(filename).name.substring(0, 100)}_${Date.now()}.pdf`;
        const pdfPath = path.join(PDFS_DIR, safeFilename);
        
        const success = await downloadFile(pdfUrl, pdfPath);
        if (success) downloadedPdfs++;
      } catch (error) {
        console.error(`PDFダウンロードエラー (${pdfUrl}):`, error);
      }
    }
    
    console.log('ダウンロード統計:');
    console.log(`- 画像: ${downloadedImages}/${mediaUrls.images.length} 個ダウンロード完了`);
    console.log(`- PDF: ${downloadedPdfs}/${allPdfUrls.length} 個ダウンロード完了`);
    
  } catch (error) {
    console.error('実行中にエラーが発生しました:', error);
  } finally {
    await browser.close();
    console.log('ブラウザを閉じました');
  }
}

// コマンドライン引数からURLを取得
const targetUrl = process.argv[2];

if (!targetUrl) {
  console.error('使用方法: node file-downloader.js <URL>');
  process.exit(1);
}

downloadMediaFromUrl(targetUrl)
  .then(() => console.log('処理が完了しました'))
  .catch(err => console.error('エラーが発生しました:', err));
```

### 使用方法

```
node file-downloader.js https://example.com
```

このプログラムは次のことを行う：

1. 指定されたURLにアクセス
2. ページ内のすべての画像を検索（img要素や背景画像など）
3. ページ内のPDFリンクを検索
4. ネットワークレスポンスからPDFを直接検出
5. 見つかったメディアを重複なくダウンロード
6. 「downloads」ディレクトリ内に「images」と「pdfs」フォルダを作成して保存

