// LeetCode Contest
//
// NO:      1601
// LEVEL:   HARD
// TITLE:   Maximum Number of Achievable Transfer Requests
// URL:     https://leetcode.com/problems/maximum-number-of-achievable-transfer-requests/description/
//
// STATUS:
// SUBMISSION:
// RUNTIME:     ??? ms | Beats ??.??%
// MEMORY:      ??.?? MB | Beats ??.??%

#include <climits>
//#include <cmath>
//#include <functional>
//#include <queue>
//#include <ranges>
#include <set>
#include <vector>

// Storage index in adjacency matrix
constexpr static int index(int i, int j) {
    return (i & 0x1f) | ((j & 0x1f) << 5);
}

constexpr int WMAX = 1 + index(31, 31);

// Adjacency matrix, with weights.
class Weights {
    int storage[WMAX] = { 0 };
public:
    constexpr int& operator () (int i, int j) noexcept {
        return storage[index(i, j)];
    }

    constexpr int operator () (int i, int j) const noexcept {
        return storage[index(i, j)];
    }

    constexpr void clear() noexcept {
        for (auto &s: storage) { s = 0; }
    }
};

class Graph {
    Weights weights;
public:
    void insert(int from, int to) {
        ++(weights(from, to));
    }

    constexpr int& operator () (int i, int j) noexcept {
        return weights(i, j);
    }

    constexpr int operator () (int i, int j) const noexcept {
        return weights(i, j);
    }

    void clear() noexcept {
        weights.clear();
    }
};

class Solution {
public:
    int maximumRequests(int n, std::vector<std::vector<int>>& requests) {
        Graph g;

        int total = 0;

        return total;
    }
};
