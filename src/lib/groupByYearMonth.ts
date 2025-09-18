import type { CollectionEntry } from 'astro:content'

export interface YearMonthCount {
  key: string // YYYY-MM
  year: string
  month: string // 01..12
  count: number
}

export function groupPostsByYearMonth(
  posts: CollectionEntry<'posts'>[],
): YearMonthCount[] {
  const map = new Map<string, number>()

  for (const post of posts) {
    if (post.data.draft) continue

    const date = new Date(post.data.published)
    const year = date.getFullYear().toString()
    const month = (date.getMonth() + 1).toString().padStart(2, '0')
    const key = `${year}-${month}`
    map.set(key, (map.get(key) || 0) + 1)
  }

  const arr: YearMonthCount[] = Array.from(map.entries()).map(
    ([key, count]) => {
      const [year, month] = key.split('-')
      return { key, year, month, count }
    },
  )

  // Sort descending by year-month (lexicographic works for YYYY-MM)
  arr.sort((a, b) => b.key.localeCompare(a.key))
  return arr
}
