import rss from '@astrojs/rss';
import {siteConfig} from '@/config';
import { getCollection } from 'astro:content';
import sanitizeHtml from 'sanitize-html';
import MarkdownIt from 'markdown-it';
const parser = new MarkdownIt();

export async function GET(context: any) {
  const blog = await getCollection('posts');

  // Sort posts by date (descending: newest first)
  const sortedPosts = blog.sort((a, b) => 
    new Date(b.data.published).valueOf() - new Date(a.data.published).valueOf()
  );

  return rss({
    title: siteConfig.title,
    description: siteConfig.subtitle || 'No description',
    site: context.site,
    // Limit to the latest 30 posts to avoid the 25 MiB file size limit.
    // Adjust the number (30) if you need more or less history.
    items: sortedPosts.slice(0, 30).map((post) => ({
        title: post.data.title,
        pubDate: post.data.published,
        description: post.data.description,
        link: `/posts/${post.slug}/`,
        content: sanitizeHtml(parser.render(post.body), {
          allowedTags: sanitizeHtml.defaults.allowedTags.concat(['img'])
        }),
      })),
    customData: `<language>${siteConfig.lang}</language>`,
  });
}
