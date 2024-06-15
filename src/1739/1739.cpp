// LeetCode Contest
//
// NO:      1739
// LEVEL:   HARD
// TITLE:   Building Boxes
// URL:     https://leetcode.com/problems/building-boxes/description/
//
// STATUS:  Accepted, Jun 15, 2024 12:12
// RUNTIME: 0ms | Beats 100%
// MEMORY:  7.18 MB | Beats 48.81%
//
// Note that the % performance (Beats) is not reliable here. The solution
// has O(1) complexity and takes virtually no time and no memory. The execution
// time varies between 0ms .. 5ms changing Beats % between 100% and 20..40%.

#include <cmath>
#include <iostream>


/**
 * Returns k'th the tetrahedral number.
 */
constexpr unsigned long int th(unsigned long int k) noexcept
{
    return k * (k + 1) * (k + 2) / 6;
}

/**
 * For a given number n, solves the equation
 *
 *  x * (x + 1) * (x + 2) / 6 = n
 *
 * for "x". The formula on the left-hand side
 * is the value of x'th tetrahedral number.
 * See: https://en.wikipedia.org/wiki/Tetrahedral_number
 */
constexpr long double thr(long double n) noexcept
{
    long double const u = std::sqrt(9.0 * n * n - 1/27.0);
    return std::cbrt(3.0 * n + u) + std::cbrt(3.0 * n - u) - 1.0;
}

struct thr_t {
    unsigned int k;
    unsigned long m;
};

/**
 * Integer version of the thr().
 */
constexpr thr_t thr(unsigned int n) noexcept
{
    unsigned int k = static_cast<unsigned int>(thr(static_cast<long double>(n)));
    unsigned long m = th(k);
    for (; m < n; m = th(++k));
    return { k, m };
}

class Solution {
public:
    unsigned int minimumBoxes(unsigned int n)
    {
        auto const [ k, m ] = thr(n);
        unsigned const int r = m - n;
        long double const b = (2.0 * k + 1);
        long double const d = b * b - 8.0 * r;
        long double const l = (b - std::sqrt(d)) / 2.0;

        return k * (k + 1) / 2 - static_cast<unsigned int>(l);
    }
};
