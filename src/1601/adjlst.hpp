#ifndef LC_1601_ADJLST_HPP
#define LC_1601_ADJLST_HPP

#include "vector.hpp"

/**
 * T - data type stored
 * N - max number of nodes in graph
 * M - max number of directed arcs from i to j (multigraphs).
 */
template<typename T, size_t N, size_t M>
class AdjLstStar {
    Vector<size_t, N>     _point;
    Vector<T, M * N * N>  _edges;

    constexpr size_t _edge_index(size_t i, size_t j, size_t m) const noexcept {
        return  M * (N * i  + j) + m;
    }
public:


};

#endif /* LC_1601_ADJLST_HPP */
