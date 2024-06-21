// LeetCode Contest
//
// NO:      1675
// LEVEL:   HARD
// TITLE:   Minimize Deviation in Array
// URL:     https://leetcode.com/problems/minimize-deviation-in-array/description/
//
// STATUS:      Accepted
// SUBMISSION:  https://leetcode.com/problems/minimize-deviation-in-array/submissions/1294000780/
// RUNTIME:     100 ms | Beats 98.33%
// MEMORY:      58.07 MB | Beats 95.32%

#include <climits>
#include <cmath>
#include <functional>
#include <queue>
#include <ranges>
#include <set>
#include <vector>

constexpr size_t NMAX = 5E4;

class Solution {
public:
    int minimumDeviation(std::vector<int>& nums) const {
        int marr[NMAX];

        int mmin = INT_MAX, mmax = 0, omin = INT_MAX;
        for (size_t i = 0; i != nums.size(); ++i) {
            auto m = nums[i];
            omin = std::min(omin, m << (m & 1));    // multiply odd numbers once
            for (; m && !(m & 1); m >>= 1);         // divide m until it gets odd
            marr[i] = m;
            mmin = std::min(m, mmin);
            mmax = std::max(m, mmax);
        }
        omin = std::min(omin, mmax + 1);

        if (mmin == mmax) {
            return 0; // solution found in O(n)
        }

        if (3 * mmin > 2 * mmax) {
            return mmax - mmin; // solution found in O(n)
        }

        if (2 * mmin > mmax) {
            std::set<int> mset{marr, marr + nums.size()};
            return mindev(mset, mmax);
        }

        std::set<int> cset;
        for (size_t i = 0; i < nums.size(); ++i) {
            int m = marr[i];
            for (; 2 * m < omin; m *= 2);
            if (mmax < 2 * m && m < omin) {
                cset.insert(m);
            }
        }
        cset.insert(omin);

        return mindev(cset, mmax);
    }

    int mindev(std::set<int> const& cset, int bmax) const noexcept
    {
        auto citr = cset.cbegin();
        auto const cend = cset.cend();

        int bmin = *citr++, dmin = bmax - bmin;
        for (; citr != cend; ++citr) {
            bmax = 2 * bmin;
            bmin = *citr;
            dmin = std::min(dmin, bmax - bmin);
        }

        return dmin;
    }
};
