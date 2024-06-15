// LeetCode Contest
//
// NO:      3161
// LEVEL:   HARD
// TITLE:   Block Placement Queries
// URL:     https://leetcode.com/problems/block-placement-queries/
//
// STATUS: Accepted, Jun 14, 2024 22:44
// RUNTIME: 718ms | Beats 99.20%
// MEMORY: 294.48 MB | Beats 80.46%

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

constexpr size_t tree_height(size_t n)
{
    size_t h = 0;
    for (; (1 << h) < n; ++h);
    return 0 == n ? 0 : 1 + h; // = 0 or 1 + ceil(log2(n))
}

constexpr size_t tree_size(size_t n)
{
    return (1 << tree_height(n)) - 1; // = 2^{1 + ceil(log2(n))} - 1
}

constexpr size_t NMAX = 5E4;


/**
 * Segment tree specific to this contest.
 *
 * Following contest hints, the tree holds values of (d[x]-1) for each x, where
 * d[x] is distance to an obstacle next to x. We use (d[x]-1) instead of d[x],
 * as this allows us to use uint16_t type (if suitable) to represent tree
 * values for all possible values of x \in [0, 5E4]. This is possible, because
 * the nearest power of 2 next to 5E4 is 2^16 = 65536, and uint16_t can hold
 * values 0 .. 65565 inclusive.
 */
template <typename T = int> class Tree {
public:
    T tree[tree_size(NMAX)];
    size_t height;
    size_t length;

    // Initialize tree using vertices [2^h - 1, ..., 1, 0]
    void build(size_t nmax) noexcept
    {
        height = tree_height(nmax);
        length = 1 << (height - 1);

        for (size_t i = (1 << height) - 1, x = 0, s = 1; i != 0;) {
            tree[--i] = x;
            if (length - 1 < (x += s)) {
                s *= 2;
                x = s - 1;
            }
        }
    }

    // Put obstacle at x
    void obstacle(int x) noexcept
    {
        size_t r = length - 1; // element [0] of the last tree level
        size_t i = r + x - 1;  // start from "x - 1" (x is at obstacle, it must remain untouched)
        int y = 0;             // "x - 1" is "1"-far from "x", but we use 0 to denote 1

        for (; i >= r && y < tree[i];) {
            tree[i--] = y++;
        }

        update(x - y, x - 1, 0, 0, r);
    }

    // Update non-vertex tree nodes
    void update(size_t xl, size_t xr, size_t i, size_t l, size_t r) noexcept
    {
        if (l == r) {
            return;
        } else {
            const size_t li = 2 * i + 1, ri = 2 * i + 2;
            const size_t m = (l + r) / 2ul;

            if (xl <= m) {
                update(xl, xr, li, l, m);
            }

            if (xr >= m + 1) {
                update(xl, xr, ri, m + 1, r);
            }

            tree[i] = std::max(tree[li], tree[ri]);
        }
    }

    // Returns max(d[0] - 1, d[1] - 1, ..., d[x] - 1)
    int query(int x, size_t i, size_t l, size_t r) noexcept
    {
        if (x < 0) {
            return 0;
        }

        if (r <= x) {
            return tree[i];
        }

        size_t li = 2 * i + 1, ri = 2 * i + 2;
        size_t m  = (l + r) / 2ul;
        if (m + 1 <= x) {
            return std::max(query(x, li, l, m), query(x, ri, m + 1, r));
        } else {
            return query(x, li, l, m);
        }
    }
};

class Solution {
public:
    Tree<int> tree;

    std::vector<bool> getResults(std::vector<std::vector<int>> const& queries) {
        std::vector<bool> results;

        size_t xmax = std::min(NMAX, 3 * queries.size());

        tree.build(1 + xmax);

        for (auto const& q :queries) {
            if (1 == q[0]) {
                tree.obstacle(q[1]);
            } else {
                int max_gap = 1 + tree.query(q[1] - q[2], 0, 0, tree.length - 1);
                results.push_back(q[2] <= max_gap);
            }
        }

        return results;
    }
};
