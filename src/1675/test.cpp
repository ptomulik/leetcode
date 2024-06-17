#include <vector>
#include <iostream>
#include <iterator>

#include "solution.cpp"

int main()
{
    Solution s;
    std::vector<std::vector<int>> cases({
        { 4, 1, 5, 20, 3 },
        { 3, 5 },
        { 10, 4, 3 },
        { 2, 9, 12, 4 },
        { 4, 9, 4, 5 },
        { 136, 465, 87 },
        { 399, 908, 648, 357, 693, 502, 331, 649, 596, 698 },
        { 165, 319, 305 },
        { 610, 778, 846, 733, 395 },
        { 200000, 199998, 199996, 199994, 199992 }
    });

    for (auto& nums: cases) {
        std::cout << "============" << std::endl;
        std::cout << "a[i]: ";
        std::set<int>aset(nums.begin(), nums.end());
        std::copy(aset.begin(), aset.end(), std::ostream_iterator<int>(std::cout, ", "));
        std::cout << std::endl;
        auto d = s.minimumDeviation(nums);
        std::cout << "d: " << d << std::endl;
    }
    return 0;
}
