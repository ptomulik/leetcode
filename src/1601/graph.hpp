#ifndef LC_1601_GRAPH_HPP
#define LC_1601_GRAPH_HPP

#include "matrix.hpp"
#include "adjlst.hpp"

#include <type_traits>
#include <limits>

namespace Templates {

template <typename T, size_t N, size_t M>
class Graph {
public:
    typedef T node_type;

    typedef typename std::make_unsigned<T>::type unsigned_int;
    typedef typename std::make_signed<T>::type signed_int;

    typedef signed_int flow_t;
    typedef signed_int cost_t;
    typedef signed_int caps_t;

    typedef Matrix<bool, N>   adjmat_type;
    typedef Matrix<flow_t, N> flow_matrix;
    typedef Matrix<cost_t, N> cost_matrix;
    typedef Matrix<caps_t, N> caps_matrix;

    typedef Edge<T, M> edge_type;
    typedef AdjLst<edge_type, N> adjlst_type;
    typedef SortSeq<node_type, N> nodeseq_type;
private:

    nodeseq_type _nodes;
    adjlst_type  _outbound;
    adjlst_type  _inbound;
    adjmat_type  _adjm[M];
    flow_matrix  _flow[M];
    cost_matrix  _cost[M];
    caps_matrix  _caps[M];

public:

    constexpr Graph() noexcept {}

    constexpr Graph(size_t n) noexcept {
        reset(n);
    }

    constexpr auto const& nodes() const noexcept {
        return _nodes;
    }

    constexpr auto const& outbound() const noexcept {
        return _outbound;
    }

    constexpr auto const& outbound(node_type i) const noexcept {
        return _outbound(i);
    }

    constexpr auto const& inbound() const noexcept {
        return _inbound;
    }

    constexpr auto const& inbound(node_type i) const noexcept {
        return _inbound(i);
    }

    /**
     * Returns the value provided to the constructor or to reset().
     *
     * The returned value determines maximum possible node index. It's not
     * the actual number of nodes existing in the graph.
     */
    constexpr auto size() const noexcept {
        return _adjm[0].size();
    }

    constexpr void reset(size_t n, flow_t x = 0, cost_t c = std::numeric_limits<cost_t>::max(), caps_t u = 0) noexcept {
        for (size_t m = 0; m < M; ++m) {
            _adjm[m].reset(n, 0);
            _flow[m].reset(n, x);
            _cost[m].reset(n, c);
            _caps[m].reset(n, u);
        }
        _outbound.reset();
        _inbound.reset();
    }

    constexpr void connect(edge_type e, cost_t c = 1) noexcept {
        adjm(e) = true;
        cost(e) = c;
        _outbound.connect(e);
        _inbound.connect(e.reversed());
        _nodes.insert(e.i);
        _nodes.insert(e.j);
    }

    constexpr bool connected(edge_type e) const noexcept {
        return adjm(e);
    }

    constexpr void disconnect(edge_type e) noexcept {
        _adjm[e.m](e.i, e.j) = false;
        _outbound.disconnect(e);
        _inbound.disconnect(e.reversed());
        if (0 == _outbound(e.i).size() && 0 == _inbound(e.i).size()) {
            _nodes.remove(e.i);
        }
        if (0 == _outbound(e.j).size() && 0 == _inbound(e.j).size()) {
            _nodes.remove(e.j);
        }
    }

    constexpr void remove(node_type i) noexcept {
        for (size_t m = 0; m < M; ++m) {
            _adjm[m].set(i, false);
            _flow[m].set(i, 0);
            _cost[m].set(i, std::numeric_limits<cost_t>::max());
            _caps[m].set(i, std::numeric_limits<caps_t>::max());
        }
        _outbound.remove(i);
        _inbound.remove(i);
        _nodes.remove(i);
    }

    // Add capacity to edge
    constexpr void add(edge_type e, caps_t du = 1) noexcept {
        if (!connected(e)) {
            connect(e);
        }

        if (caps(e) < std::numeric_limits<caps_t>::max() - du) {
            caps(e) += du;
        } else {
            caps(e) = std::numeric_limits<caps_t>::max();
        }
    }

    // Remve capacity from edge
    constexpr void sub(edge_type e, caps_t du = 1) noexcept {
        if (du < caps(e)) {
            caps(e) -= du;
        } else {
            caps(e) = 0;
        }
    }

    constexpr bool adjm(edge_type e) const noexcept {
        return _adjm[e.m](e.i, e.j);
    }

    constexpr bool& adjm(edge_type e) noexcept {
        return _adjm[e.m](e.i, e.j);
    }

    constexpr flow_t flow(edge_type e) const noexcept {
        return _flow[e.m](e.i, e.j);
    }

    constexpr flow_t& flow(edge_type e) noexcept {
        return _flow[e.m](e.i, e.j);
    }

    constexpr cost_t cost(edge_type e) const noexcept {
        return _cost[e.m](e.i, e.j);
    }

    constexpr cost_t& cost(edge_type e) noexcept {
        return _cost[e.m](e.i, e.j);
    }

    constexpr caps_t caps(edge_type e) const noexcept {
        return _caps[e.m](e.i, e.j);
    }

    constexpr caps_t& caps(edge_type e) noexcept {
        return _caps[e.m](e.i, e.j);
    }

    constexpr auto outdegree(node_type i) const noexcept {
        return _outbound.degree(i);
    }

    constexpr auto indegree(node_type i) const noexcept {
        return _inbound.degree(i);
    }

    constexpr auto outflow(node_type i) const noexcept {
        flow_t x = 0;
        for (auto head: _outbound(i)) {
            x += flow(edge_type(i, head));
        }
        return x;
    }

    constexpr auto inflow(node_type i) const noexcept {
        flow_t x = 0;
        for (auto head: _inbound(i)) {
            x += flow(edge_type(head, i));
        }
        return x;
    }

    constexpr auto balance(node_type i) const noexcept {
        return inflow(i) - outflow(i);
    }
};

}

#endif /* LC_1601_GRAPH_HPP */
