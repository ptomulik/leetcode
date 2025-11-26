#ifndef LC_1601_EDGE_HPP
#define LC_1601_EDGE_HPP

#include <cstddef>

namespace Templates {

/**
 * Edge endpoint.
 *
 * Encapsulates node index i of node at endpoint and number m (multi-graph).
 */
template<typename T, size_t M>
struct Endpoint {
    typedef T value_type;

    T i;
    T m;

    constexpr Endpoint(T encoded) noexcept : i(encoded / M), m(encoded % M) {}
    constexpr Endpoint(T i, T m) noexcept : i(i), m(m) {}

    constexpr T encoded() const noexcept {
        return M * i + m;
    }
};

template <size_t M, typename T>
constexpr auto make_endpoint(T i, T m) noexcept {
    return Endpoint<T, M>(i, m);
}

template <size_t M, typename T>
constexpr auto make_endpoint(T encoded) noexcept {
    return Endpoint<T, M>(encoded);
}

template<typename T, size_t M>
struct Edge {
    constexpr static const size_t max_m = M;

    typedef T value_type;
    typedef Endpoint<T, M> endpoint_type;
    typedef endpoint_type tail_type;
    typedef endpoint_type head_type;

    T i;
    T j;
    T m;

    constexpr Edge(T i, T j, T m) noexcept : i(i), j(j), m(m) {}
    constexpr Edge(T i, head_type h) noexcept : i(i), j(h.i), m(h.m) {}
    constexpr Edge(tail_type t, T j) noexcept : i(t.i), j(j), m(t.m) {}

    constexpr auto tail() const noexcept {
        return tail_type(i, m);
    }

    constexpr auto head() const noexcept {
        return head_type(j, m);
    }

    constexpr auto reversed() const noexcept {
        return Edge<T, M>(j, i, m);
    }
};

template <size_t M, typename T>
constexpr auto make_edge(T i, T j, T m) noexcept {
    return Edge<T, M>(i, j, m);
}

template <size_t M, typename T>
constexpr auto make_edge(T i, Endpoint<T, M> head) noexcept {
    return Endpoint<T, M>(i, head);
}

template <size_t M, typename T>
constexpr auto make_edge(Endpoint<T, M> tail, T j) noexcept {
    return Endpoint<T, M>(tail, j);
}

} /* namespace Structures */

#endif /* LC_1601_EDGE_HPP */
