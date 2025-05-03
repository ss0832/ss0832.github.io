---
title:  【Markdown】Markdownで折り畳み表示を実装する方法
published: 2025-05-03
description: "Markdownで折り畳み表示（アコーディオン）を実装する方法について解説する"
tags: [Markdown]
category: text
draft: false
---
最終更新：2025-05-03

# Markdownにおける折り畳み表示の実装方法

## 概要

Markdownで折り畳み表示（アコーディオン）を実装するには、HTML記法の`<details>`と`<summary>`タグを使用する。これらのタグはMarkdown内でもHTML要素として認識され、折り畳み機能を提供する。GitHubやQiitaなど多くのMarkdownをサポートするプラットフォームで使用可能である。

## 実装方法

### 基本構文

```
<details>
<summary>ここをクリックして開く</summary>

ここに折り畳まれる内容を記述する。
Markdown記法も使用可能である。

</details>
```

### 表示例

以下は実際の表示例である：

<details>
<summary>ここをクリックして開く</summary>

折り畳まれた内容がここに表示される。
- リスト表示も可能
- 画像やコードブロックも挿入できる

```javascript
console.log('コードブロックも使用可能');
```

</details>

### 注意点

1. `<summary>`タグと`</summary>`タグの間にクリック時に表示されるタイトルを記述する
2. `<summary>`タグ閉じた後に必ず改行を入れること
3. 折り畳まれる内容を記述した後、`</details>`タグの前にも改行を入れること
4. 折り畳み内部でもMarkdown記法を使用できるが、環境によっては対応していない場合がある
5. GitHubでは折り畳み内のコード構文ハイライトも正常に機能する

## 応用例

### デフォルトで開いた状態にする

最初から開いた状態にする場合は以下のように`open`属性を追加する：

```
<details open>
<summary>最初から開いた状態</summary>

この部分は最初から表示されている。

</details>
```

### ネスト（入れ子）構造

折り畳み要素は入れ子にすることも可能である：

```
<details>
<summary>親の折り畳み</summary>

親の内容

<details>
<summary>子の折り畳み</summary>

子の内容

</details>
</details>
```

### スタイルのカスタマイズ

`style`属性を使用してカスタマイズすることも可能である：

```
<details>
<summary style="color: red; font-weight: bold;">カスタマイズされたタイトル</summary>

内容をここに記述する。

</details>
```

### マークダウン内での見出し

`summary`タグ内でMarkdownの見出し記法を使うこともできる：

```
<details>
<summary><h3>見出しスタイルのサマリー</h3></summary>

内容をここに記述する。

</details>
```

ただし、プラットフォームによっては正しく表示されない場合がある。

## まとめ

Markdownでの折り畳み機能は、`<details>`と`<summary>`タグを使用することで簡単に実装できる。長文の資料や補足情報をすっきりと整理したい場合に有効な手法である。サポートされている環境であれば、記事の可読性を大きく向上させることができる。

折り畳み内でもMarkdown記法が使用できることが大きな利点だが、プラットフォームによって対応状況が異なるため、実装時には対象環境での動作確認が重要である。

## 参考サイト
- [Markdownで折り畳み表示（アコーディオン）を実装する方法](https://qiita.com/matagawa/items/31e26e9cd53c3e61ae07)
- [Markdownの折り畳み(details)記法を紹介するよ](https://qiita.com/P-man_Brown/items/067bfa132eb3c4b49bc4)
- [Markdown as supported in GitHub, including extensions](https://gist.github.com/Phroneris/e7e6c869640b95bd42434bdc995cd4f6)