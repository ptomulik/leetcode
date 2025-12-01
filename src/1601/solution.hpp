// LeetCode Contest
//
// NO:      1601
// LEVEL:   HARD
// TITLE:   Maximum Number of Achievable Transfer Requests
// URL:     https://leetcode.com/problems/maximum-number-of-achievable-transfer-requests/description/
//
// STATUS:      Accepted
// SUBMISSION:  https://leetcode.com/problems/maximum-number-of-achievable-transfer-requests/submissions/1836052186/
// RUNTIME:     0ms | Beats 100.00%
// MEMORY:      11.35 MB | Beats 99.39%

#ifndef LC_1601_SOLUTION_HPP
#define LC_1601_SOLUTION_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

typedef int8_t int_t; //! Base integer type.
typedef typename std::make_unsigned<int_t>::type uint_t;   //! Unsigned integer value.
typedef typename std::make_signed<int_t>::type sint_t;     //! Signed integer value.
typedef uint_t node_t; //! Node identifier.
typedef uint_t edge_t; //! Arc number for multiarc.
typedef sint_t flow_t; //! Arc flow.
typedef sint_t cost_t; //! Arc flow cost.
typedef sint_t caps_t; //! Arc capacity.
typedef cost_t dist_t; //! Node distance on path.
typedef node_t pred_t; //! Node predecessor on path.
typedef sint_t potn_t; //! Node potential.

constexpr const size_t NMAX = 20; //! Max supported no. nodes in graph.
constexpr const size_t MMAX = 2;  //! Residual network is a multigraph with up to MMAX parallel arcs.
// constexpr uint_t FMAX = 16; //! Max supported total capacity (sum of arc capacities).

constexpr cost_t cost_max() noexcept {
    return std::numeric_limits<cost_t>::max();
}

constexpr caps_t caps_max() noexcept {
    return std::numeric_limits<caps_t>::max();
}

constexpr caps_t caps_min() noexcept {
    return std::numeric_limits<caps_t>::min();
}

/**
 * Vector of size N with inplace storage.
 */
template <typename T, size_t N>
class Vector {
public:
    typedef T value_type;
    typedef value_type* pointer;
    typedef value_type const* const_pointer;
    typedef value_type* iterator;
    typedef value_type const* const_iterator;
    typedef value_type& reference;
    typedef value_type const& const_reference;
private:
    size_t     _n;
    value_type _vec[N];
public:

    constexpr Vector() noexcept: _n(0u) {}

    constexpr size_t size() const noexcept {
        return _n;
    }

    constexpr auto&& operator()(this auto&& self, size_t i) noexcept {
        return self._vec[i];
    }

    constexpr void reset(size_t n, value_type const& v = value_type()) noexcept {
        resize(n);
        fill(v);
    }

    constexpr void resize(size_t n) noexcept {
        _n = n;
    }

    constexpr iterator erase(iterator pos) noexcept {
        auto const* const e = end();
        if (pos >= begin() && pos < e) {
            for (iterator ptr = pos; ptr < e; ++ptr) {
                *ptr = *(ptr + 1);
            }
            --_n;
        }
        return pos;
    }

    constexpr void fill(value_type const& v) noexcept {
        std::fill_n(&_vec[0], _n, v);
    }

    constexpr void clear() noexcept {
        _n = 0;
    }

    constexpr void push_back(value_type value) noexcept {
        _vec[_n++] = value;
    }

    constexpr value_type pop_back() noexcept {
        return _vec[--_n];
    }

    constexpr auto* begin(this auto&& self) noexcept {
        return &self._vec[0];
    }

    constexpr auto* end(this auto&& self) noexcept {
        return &self._vec[self._n];
    }
};

/**
 * Square NxN matrix with inplace storage.
 */
template <typename T, size_t N>
class Matrix {
public:
    typedef T value_type;
    typedef value_type* pointer;
    typedef value_type const* const_pointer;
private:
    size_t     _n;
    value_type _mat[N][N];
public:
    constexpr Matrix() noexcept: _n(0u) { }
    constexpr Matrix(size_t n) noexcept : _n(n) { }
    constexpr Matrix(size_t n, value_type const& x) noexcept : _n(n), _mat{x} { }

    constexpr void reset(size_t n, value_type const& v = value_type()) noexcept {
        resize(n);
        fill(v);
    }

    constexpr void resize(size_t n) noexcept {
        _n = n;
    }

    constexpr void fill(value_type const& v) noexcept {
        std::fill_n(&_mat[0][0], N * N, v);
    }

    constexpr auto size() const noexcept {
        return _n;
    }

    constexpr auto&& operator () (this auto&& self, size_t i, size_t j) noexcept {
        return self._mat[i][j];
    }

    // Set value to all elements of row and column i.
    constexpr void set(size_t i, value_type const& v) noexcept {
        for (size_t j = 0; j < _n; ++j) {
            _mat[i][j] = v;
            _mat[j][i] = v;
        }
    }
};

/**
 * Sorted sequence of numbers.
 */
template <typename T, size_t N>
class SortSeq {
public:
    typedef T value_type;
    typedef value_type* iterator;
    typedef value_type const* const_iterator;
    typedef value_type& reference;
    typedef value_type const& const_reference;
private:
    size_t _n;
    value_type _seq[N];

    constexpr ptrdiff_t _insertion_index(value_type v) const noexcept {
        if (0 == _n) {
            return 0;
        }

        if (v < _seq[0]) {
            return 0;
        }

        if (v > _seq[_n-1]) {
            return _n;
        }

        value_type lo = 0, hi = _n-1;
        // Find the index using bisection
        for (auto mid = (hi + lo)/2; (hi - lo) > 1; mid = (lo + hi)/2) {
            if (v < _seq[mid]) {
                hi = mid;
            } else if (v > _seq[mid]) {
                lo = mid;
            } else {
                return -1; // already exists
            }
        }

        if (v == _seq[lo] || v == _seq[hi]) {
            return -1;
        }

        return hi;
    }

    constexpr ptrdiff_t _removal_index(value_type v) const noexcept {
        if (0 == _n) {
            return -1;
        }

        if (v < _seq[0] || v > _seq[_n-1]) {
            return -1;
        }

        value_type lo = 0, hi = _n-1;
        // Find the index using bisection
        for (auto mid = (hi + lo)/2; (hi - lo) > 1; mid = (lo + hi)/2) {
            if (v < _seq[mid]) {
                hi = mid;
            } else if (v > _seq[mid]){
                lo = mid;
            } else {
                return mid; // found
            }
        }

        if (v == _seq[lo]) {
            return lo;
        }

        if (v == _seq[hi]) {
            return hi;
        }

        return -1;
    }

    constexpr void _unshift(value_type v, size_t at = 0) noexcept {
        for (auto j = _n; j > at; --j) {
            _seq[j] = _seq[j-1];
        }
        _seq[at] = v;
        ++_n;
    }

    constexpr auto _shift(size_t at = 0) noexcept {
        auto v = _seq[at];

        for (value_type j = at + 1; j < _n; ++j) {
            _seq[j-1] = _seq[j];
        }
        --_n;

        return v;
    }

public:
    constexpr SortSeq() noexcept: _n(0) {}

    constexpr auto size() const noexcept {
        return _n;
    }

    constexpr void reset() noexcept {
        _n = 0;
    }

    constexpr auto operator() (size_t k) const noexcept {
        return _seq[k];
    }

    constexpr void insert(value_type v) noexcept {
        auto const at = _insertion_index(v);
        if (-1 != at) {
            _unshift(v, at);
        }
    }

    constexpr void remove(value_type v) noexcept {
        auto const at = _removal_index(v);
        if (-1 != at) {
            _shift(at);
        }
    }

    constexpr auto* begin(this auto&& self) noexcept {
        return &self._seq[0];
    }

    constexpr auto* end(this auto&& self) noexcept {
        return &self._seq[self._n];
    }
};

/**
 * Reference to an "encoded" object (forgive me naming).
 *
 * Encapsulates a reference to an "encoded" value. Provides methods for
 * converting between the encoded value type and TObject.
 */
template <typename TObject>
class EncodedObjectRef {
public:
    typedef typename std::remove_cv<TObject>::type object_type;
    typedef object_type value_type;
    typedef typename TObject::encoded_type encoded_type;
    typedef typename std::conditional<
        std::is_const<TObject>::value,
        typename std::add_const<encoded_type>::type,
        encoded_type
    >::type& encoded_reference;
private:
    encoded_reference _ref;
public:
    constexpr EncodedObjectRef(encoded_reference ref) noexcept : _ref(ref) {}

    constexpr auto const& operator= (object_type obj) const noexcept {
        _ref = obj.encoded();
        return *this;
    }

    constexpr operator object_type() const noexcept {
        return object_type(_ref);
    }

    constexpr encoded_type encoded() const noexcept {
        return _ref;
    }
};

template <typename TObject>
class EncodedObjectIter {
public:
    typedef typename std::remove_cv<TObject>::type object_type;
    typedef object_type value_type;
    typedef typename TObject::encoded_type encoded_type;
    typedef typename std::conditional<
        std::is_const<TObject>::value,
        typename std::add_const<encoded_type>::type,
        typename std::remove_const<encoded_type>::type
    >::type* encoded_pointer;
    typedef EncodedObjectIter<TObject> iterator;
    typedef EncodedObjectRef<TObject> reference;
private:
    encoded_pointer _ptr;
public:
    constexpr EncodedObjectIter(encoded_pointer ptr) noexcept : _ptr(ptr) {}

    constexpr reference operator*() const noexcept {
        return reference(*_ptr);
    }

    constexpr void operator++() noexcept { ++_ptr; }
    constexpr void operator--() noexcept { --_ptr; }

    constexpr reference operator[](size_t n) const noexcept {
        return _ptr[n];
    }

    template<class TOther>
    constexpr bool operator == (const TOther& other) const noexcept {
        return _ptr == other._ptr;
    }

    template<class TOther>
    constexpr bool operator != (const TOther& other) const noexcept {
        return !(*this == other);
    }
};

template <typename TObject, size_t N>
class EncodedObjectSortSeq
{
public:
    typedef TObject object_type;
    typedef TObject value_type;
    typedef typename TObject::encoded_type encoded_type;
    typedef SortSeq<encoded_type, N> sequence_type;
    typedef EncodedObjectIter<object_type> iterator;
    typedef EncodedObjectIter<const object_type> const_iterator;
    typedef EncodedObjectRef<object_type> reference;
    typedef EncodedObjectRef<const object_type> const_reference;
private:
    sequence_type _seq;
public:
    constexpr auto size() const noexcept {
        return _seq.size();
    }

    constexpr void reset() noexcept {
        _seq.reset();
    }

    constexpr auto operator() (size_t k) const noexcept {
        return _seq(k);
    }

    constexpr void insert(object_type obj) noexcept {
        _seq.insert(obj.encoded());
    }

    constexpr void remove(object_type obj) noexcept {
        _seq.remove(obj.encoded());
    }

    constexpr auto begin() noexcept {
        return iterator(_seq.begin());
    }

    constexpr auto begin() const noexcept {
        return const_iterator(_seq.begin());
    }

    constexpr auto end() noexcept {
        return iterator(_seq.end());
    }

    constexpr auto end() const noexcept {
        return const_iterator(_seq.end());
    }
};

template <typename HeadLstT>
class AdjLstRef {
public:
    typedef HeadLstT head_list;
private:
    node_t _i;
    head_list& _list;
public:
    constexpr AdjLstRef(node_t i, head_list& list) noexcept
        : _i(i), _list(list) {
    }

    constexpr auto i() const noexcept {
        return _i;
    }

    constexpr head_list& list() const noexcept {
        return _list;
    }
};

template <typename NodeLstT, typename HeadLstT>
class AdjLstIter {
public:
    typedef NodeLstT node_list;
    typedef HeadLstT head_list;
    typedef AdjLstRef<head_list> reference;
private:
    node_list& _nodes;
    head_list* _lists;
    size_t _k;
public:
    constexpr AdjLstIter(node_list& nodes, head_list* lists, size_t k) noexcept
        : _nodes(nodes), _lists(lists), _k(k) {
    }

    constexpr void operator++() noexcept { ++_k; }
    constexpr void operator--() noexcept { --_k; }

    template<class TOther>
    constexpr bool operator == (const TOther& other) const noexcept {
        return _lists == other._lists && &_nodes == &other._nodes && _k == other._k;
    }

    template<class TOther>
    constexpr bool operator != (const TOther& other) const noexcept {
        return !(*this == other);
    }

    constexpr auto operator* () const noexcept {
        auto i = _nodes(_k);
        return reference(i, _lists[i]);
    }
};

/**
 * Edge endpoint.
 *
 * Encapsulates node index i of node at endpoint and number m (multi-graph).
 */
struct Endpoint {
    typedef uint_t encoded_type;

    node_t node;
    edge_t edge;

    constexpr Endpoint(encoded_type encoded) noexcept : node(encoded / MMAX), edge(encoded % MMAX) {}
    constexpr Endpoint(encoded_type node, encoded_type edge) noexcept : node(node), edge(edge) {}

    constexpr encoded_type encoded() const noexcept {
        return MMAX * node + edge;
    }
};

typedef Endpoint Head;
typedef Endpoint Tail;

/**
 * Graph edge.
 *
 * e.i - source node
 * e.j - target node
 * e.m - edge number for multi-graph.
 */
struct Edge {
    node_t i;
    node_t j;
    edge_t m;

    constexpr Edge(node_t i, node_t j, edge_t m) noexcept : i(i), j(j), m(m) {}
    constexpr Edge(node_t i, Head h) noexcept : i(i), j(h.node), m(h.edge) {}
    constexpr Edge(Tail t, node_t j) noexcept : i(t.node), j(j), m(t.edge) {}

    constexpr auto tail() const noexcept {
        return Tail(i, m);
    }

    constexpr auto head() const noexcept {
        return Head(j, m);
    }

    constexpr auto reversed() const noexcept {
        return Edge(j, i, m);
    }

    constexpr auto reversed(edge_t m) const noexcept {
        return Edge(j, i, m);
    }
};

/**
 * Ajdacency lists.
 */
class AdjLst {
public:
    typedef SortSeq<node_t, NMAX> node_list;
    typedef EncodedObjectSortSeq<Head, NMAX> head_list;

    typedef AdjLstIter<node_list, head_list> iterator;
    typedef AdjLstIter<const node_list, const head_list> const_iterator;
private:
    node_list _nodes;
    head_list _lists[NMAX];
public:
    constexpr void reset() noexcept {
        for(auto& list: _lists) {
            list.reset();
        }
    }

    constexpr void connect(Edge edge) noexcept {
        if (0 == _lists[edge.i].size()) {
            _nodes.insert(edge.i);
        }
        _lists[edge.i].insert(edge.head());
    }

    constexpr void disconnect(Edge edge) noexcept {
        _lists[edge.i].remove(edge.head());
        if (0 == _lists[edge.i].size()) {
            _nodes.remove(edge.i);
        }
    }

    constexpr void remove(node_t i) noexcept {
        node_t rj[NMAX];
        node_t rn = 0;

        _lists[i].reset();
        _nodes.remove(i);
        for (auto j: _nodes) {
            _lists[j].remove(i);
            if (0 == _lists[j].size()) {
                rj[rn++] = j;
            }
        }

        for (; rn > 0;) {
            _nodes.remove(rj[--rn]);
        }
    }

    constexpr size_t size() const noexcept {
        return _nodes.size();
    }

    constexpr auto& operator() (this auto&& self, size_t i) noexcept {
        return self._lists[i];
    }

    constexpr size_t degree(size_t i) const noexcept {
        return _lists[i].size();
    }

    constexpr auto begin() noexcept {
        return iterator(_nodes, _lists, 0);
    }

    constexpr auto end() noexcept {
        return iterator(_nodes, _lists, _nodes.size());
    }

    constexpr auto begin() const noexcept {
        return const_iterator(_nodes, _lists, 0);
    }

    constexpr auto end() const noexcept {
        return const_iterator(_nodes, _lists, _nodes.size());
    }
};

class Graph {
public:
    typedef SortSeq<node_t, NMAX> node_list;
    typedef std::array<Matrix<bool, NMAX>, MMAX>   adjc_matrix;
    typedef std::array<Matrix<flow_t, NMAX>, MMAX> flow_matrix;
    typedef std::array<Matrix<cost_t, NMAX>, MMAX> cost_matrix;
    typedef std::array<Matrix<caps_t, NMAX>, MMAX> caps_matrix;
private:

    node_list   _nodes;
    AdjLst      _outbound;
    AdjLst      _inbound;
    adjc_matrix _adjc; //! adjacency matrix
    flow_matrix _flow; //! arc flows x_{ij}
    cost_matrix _cost; //! arc costs c_{ij}
    caps_matrix _caps; //! (residual) capacities r_{ij}

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

    constexpr auto const& outbound(node_t i) const noexcept {
        return _outbound(i);
    }

    constexpr auto const& inbound() const noexcept {
        return _inbound;
    }

    constexpr auto const& inbound(node_t i) const noexcept {
        return _inbound(i);
    }

    /**
     * Returns the value provided to the constructor or to reset().
     *
     * The returned value determines maximum possible node index. It's not
     * the actual number of nodes existing in the graph.
     */
    constexpr auto size() const noexcept {
        return _adjc[0].size();
    }

    constexpr void reset(
        size_t n,
        flow_t x = 0,
        cost_t c = cost_max(),
        caps_t u = 0
    ) noexcept {
        for (size_t m = 0; m < MMAX; ++m) {
            _adjc[m].reset(n, 0);
            _flow[m].reset(n, x);
            _cost[m].reset(n, c);
            _caps[m].reset(n, u);
        }
        _outbound.reset();
        _inbound.reset();
    }

    constexpr void connect(Edge e, cost_t c = 1) noexcept {
        adjc(e) = true;
        cost(e) = c;
        _outbound.connect(e);
        _inbound.connect(e.reversed());
        _nodes.insert(e.i);
        _nodes.insert(e.j);
    }

    constexpr void disconnect(Edge e) noexcept {
        adjc(e) = false;
        _outbound.disconnect(e);
        _inbound.disconnect(e.reversed());
        if (0 == _outbound(e.i).size() && 0 == _inbound(e.i).size()) {
            _nodes.remove(e.i);
        }
        if (0 == _outbound(e.j).size() && 0 == _inbound(e.j).size()) {
            _nodes.remove(e.j);
        }
    }

    constexpr bool connected(Edge e) const noexcept {
        return adjc(e);
    }

    constexpr void remove(node_t i) noexcept {
        for (size_t m = 0; m < MMAX; ++m) {
            _adjc[m].set(i, false);
            _flow[m].set(i, 0);
            _cost[m].set(i, cost_max());
            _caps[m].set(i, caps_max());
        }
        _outbound.remove(i);
        _inbound.remove(i);
        _nodes.remove(i);
    }

    // Add capacity to edge
    constexpr void add(Edge e, caps_t du = 1) noexcept {
        if (!connected(e)) {
            connect(e);
        }

        if (caps(e) < caps_max() - du) {
            caps(e) += du;
        } else {
            caps(e) = caps_max();
        }
    }

    // Remove capacity from edge
    constexpr void sub(Edge e, caps_t du = 1) noexcept {
        if (du < caps(e)) {
            caps(e) -= du;
        } else {
            caps(e) = 0;
        }
    }

    constexpr bool adjc(Edge e) const noexcept {
        return _adjc[e.m](e.i, e.j);
    }

    constexpr bool& adjc(Edge e) noexcept {
        return _adjc[e.m](e.i, e.j);
    }

    constexpr flow_t flow(Edge e) const noexcept {
        return _flow[e.m](e.i, e.j);
    }

    constexpr flow_t& flow(Edge e) noexcept {
        return _flow[e.m](e.i, e.j);
    }

    constexpr cost_t cost(Edge e) const noexcept {
        return _cost[e.m](e.i, e.j);
    }

    constexpr cost_t& cost(Edge e) noexcept {
        return _cost[e.m](e.i, e.j);
    }

    constexpr caps_t caps(Edge e) const noexcept {
        return _caps[e.m](e.i, e.j);
    }

    constexpr caps_t& caps(Edge e) noexcept {
        return _caps[e.m](e.i, e.j);
    }

    constexpr auto outdegree(node_t i) const noexcept {
        return _outbound.degree(i);
    }

    constexpr auto indegree(node_t i) const noexcept {
        return _inbound.degree(i);
    }

    constexpr auto outflow(node_t i) const noexcept {
        flow_t x = 0;
        for (auto head: _outbound(i)) {
            x += flow(Edge(i, head));
        }
        return x;
    }

    constexpr auto inflow(node_t i) const noexcept {
        flow_t x = 0;
        for (auto head: _inbound(i)) {
            x += flow(Edge(head, i));
        }
        return x;
    }

    constexpr auto balance(node_t i) const noexcept {
        return inflow(i) - outflow(i);
    }
};

/**
 * Encapsulates shortest paths solution.
 */
class Paths {
public:
    typedef Vector<dist_t, NMAX> dist_vector;
    typedef Vector<pred_t, NMAX> pred_vector;

private:
    size_t      _i;
    dist_vector _dists; //! Vector of distances.
    pred_vector _preds; //! Vector of predecessors.
public:
    constexpr Paths() noexcept : _i(), _dists(), _preds() { }

    constexpr size_t i() const noexcept  {
        return _i;
    }

    constexpr size_t size() const noexcept {
        return _dists.size();
    }

    constexpr void reset(size_t n, size_t i) noexcept {
        _i = i;
        _dists.reset(n, NMAX);
        _preds.reset(n, -1);
    }

    constexpr auto&& dists(this auto&& self) noexcept {
        return self._dists;
    }

    constexpr auto&& dist(this auto&& self, size_t j) noexcept {
        return self._dists(j);
    }

//    constexpr auto preds() const noexcept {
//        return EncodedObjectIter<const Tail>(_preds.begin());
//    }
//
//    constexpr auto preds() noexcept {
//        return EncodedObjectIter<Tail>(_preds.begin());
//    }

    constexpr auto pred(size_t j) const noexcept {
        return EncodedObjectRef<const Tail>(_preds(j));
    }

    constexpr auto pred(size_t j) noexcept {
        return EncodedObjectRef<Tail>(_preds(j));
    }

    constexpr bool exists(size_t j) const noexcept {
        return NMAX != _dists(j);
    }

    template<typename Function, typename U>
    constexpr U reduce(size_t j, Function func, U value) const noexcept {
        if (exists(j)) {
            while(j != _i) {
                Tail t = pred(j);
                value = func(Edge(t, j), value);
                j = t.node;
            }
        }
        return value;
    }

    template<typename Function>
    constexpr void walk(size_t j, Function func) const noexcept {
        if (exists(j)) {
            while (j != _i) {
                Tail t = pred(j);
                func(Edge(t, j));
                j = t.node;
            }
        }
    }

    constexpr sint_t flow(size_t j, Graph const& graph) const noexcept {
        if(_i == j || !exists(j)) {
            return 0;
        }

        const sint_t fmax = std::numeric_limits<sint_t>::max();

        return reduce(j, [&graph](Edge e, sint_t x) { return std::min(x, graph.flow(e)); }, fmax);
    }
};

/**
 * Dijkstra shortest path with Dial modification (bucket).
 */
class Dijkstra {
    mutable Vector<Vector<dist_t, NMAX>, NMAX> _bucket;

public:
    constexpr void shortest_paths(Graph const& graph, size_t s, Paths& paths) const noexcept {
        paths.reset(graph.size(), s);

        _bucket.reset(graph.size());

        paths.dist(s) = 0;
        _bucket(0).push_back(s);

        auto first_nonempty = [](auto& container) {
            return std::find_if(container.begin(), container.end(), [](auto const& d) {
                return d.size() > 0;
            });
        };

        for (auto stack = first_nonempty(_bucket); stack != _bucket.end(); stack = first_nonempty(_bucket)) {
            dist_t i = stack->pop_back();

            for (Head h: graph.outbound(i)) {
                const node_t j = h.node;
                dist_t dist = paths.dist(i) + graph.cost(Edge(i, j, h.edge));
                if (dist < paths.dist(j)) {
                    _bucket(dist).push_back(j);
                    paths.dist(j) = dist;
                    paths.pred(j) = Tail(i, h.edge);
                }
            }
        }
    }
};

/**
 * The workhorse class.
 */
class Optimizer {
public:
    typedef Vector<node_t, NMAX> node_vector;
    typedef Vector<sint_t, NMAX> sint_vector;
//    typedef Vector<uint_t, NMAX> uint_vector;
private:
    mutable node_vector _nodes1; //! Vector 1 holding nodes.
    mutable node_vector _nodes2; //! Vector 2 holding nodes.
    mutable sint_vector _sints1; //! Vector holding signed integers
    mutable Paths _paths1;       //! Structure of shortest paths

    Dijkstra _dijkstra;
public:

    constexpr void find_deficit_and_excess_nodes(
        Graph const& graph,
        node_vector& deficit,
        node_vector& excess,
        sint_vector& balances
    ) const noexcept {
        deficit.clear();
        excess.clear();
        balances.reset(graph.size(), 0);
        for (auto i: graph.nodes()) {
            auto e = graph.balance(i);
            balances(i) = e;
            if (e < 0) {
                deficit.push_back(i);
            } else if (e > 0) {
                excess.push_back(i);
            }
        }
    }

//    constexpr void remove_bridges(Graph& graph) const noexcept {
//        auto& bridges = _pairs1;
//        bridges.clear();
//        _tarjan.bridges(graph, [&bridges](Pair const& bridge) {
//            bridges.push_back(bridge);
//        });
//        graph.disconnect(bridges);
//    }

    constexpr uint_t max_circulation(Graph& graph) const noexcept {
        auto& deficit = _nodes1;
        auto& excess = _nodes2;
        auto& balances = _sints1;
        auto& paths = _paths1;

        // Find and remove minimal path flow.
        find_deficit_and_excess_nodes(graph, deficit, excess, balances);
//        std::cout << "deficit: [" << deficit << "]; excess: [" << excess << "]" << std::endl;
        while (deficit.size() != 0) {
            for (auto ki = deficit.begin(); ki != deficit.end();) {
                auto k = *ki;
                _dijkstra.shortest_paths(graph, k, paths);
//                std::cout << "sortest: " << paths.i() << " -> " << paths.dists() << ";" << paths.preds() << std::endl;
                for (auto li = excess.begin(); li != excess.end();) {
                    auto l = *li;
                    if (paths.exists(l)) {
                        // Flow f(P_{kl}) along the path P_{kl}
                        uint_t f = std::min(std::min((sint_t)-balances(k), balances(l)), paths.flow(l, graph));
                        paths.walk(l, [&graph, f](Edge e) {
                            const auto r = e.reversed((e.m + 1) % MMAX);
                            graph.sub(e, f);
                            graph.add(r, f);
                            graph.flow(e) += f;
                        });
                        balances(k) += f;
                        balances(l) -= f;
                        if (0 == balances(k)) {
                            break;
                        }
                    }
                    li = (0 == balances(l)) ? excess.erase(li) : li + 1;
                }
                ki = (0 == balances(k)) ? deficit.erase(ki) : ki + 1;
            }
        }

        return flow_cost(graph);
    }

    constexpr flow_t flow_cost(Graph const& graph) const noexcept {
        flow_t total = 0;
        for (auto const& ref: graph.outbound()) {
            auto i = ref.i();
            for (Head h: ref.list()) {
                Edge e(i, h);
                total += graph.cost(e) * graph.flow(e);
            }
        }
        return total;
    }
};

/**
 * Tarjan's algorithm -- identifies strongly connected components in graph
 * (actually, we identify only bridges between them).
 */
//class Tarjan {
//public:
//    typedef Vector<bool, NMAX> BoolVec;
//    typedef Vector<sint_t, NMAX> SnumVec;
//    typedef Vector<uint_t, NMAX> NodeStack;
//private:
//    mutable BoolVec _onstack;
//    mutable SnumVec _low;
//    mutable SnumVec _tin;
//    mutable uint_t  _timer;
//    mutable NodeStack _stack;
//
//    template<typename BridgeFunc>
//    constexpr void _bridges_dfs(
//        Graph const& graph,
//        BridgeFunc bridge,
//        sint_t i,
//        sint_t p = -1
//    ) const noexcept {
//        _tin(i) = _low(i) = _timer++;
//
//        _stack.push_back(i);
//        _onstack(i) = true;
//
//        for (auto j: graph.outbound(i)) {
//            if (_tin(j) == -1) {
//                // Successor j has not yet been visited. Recurse on it.
//                _bridges_dfs(graph, bridge, j, i);
//                _low(i) = std::min(_low(i), _low(j));
//            } else if (_onstack(j)) {
//                // Successor j is on stack and hence in the current SCC.
//                _low(i) = std::min(_low(i), _tin(j));
//            } else {
//                // If j is not on stack, then (i, j) is an arc pointing to an
//                // SCC already found, and must be ignored. See below, regarding
//                // the next line.
//                bridge(Pair(i, j));
//            }
//        }
//
//        // If i is a root node, pop the stack and generate an SCC.
//        if (_low(i) == _tin(i)) {
//            if (-1 != p) {
//                bridge(Pair(p, i));
//            }
//            sint_t j;
//            do {
//                j = _stack.pop_back();
//                _onstack(j) = false;
//            } while(i != j);
//        }
//    }
//
//    constexpr void _dfs_init(size_t n) const noexcept {
//        _stack.clear();
//
//        _onstack.reset(n, false);
//        _tin.reset(n, -1);
//        _low.reset(n, -1);
//
//        _timer = 0;
//    }
//
//public:
//    //! Identifies bridges between strongly connected components in graph.
//    template<typename BridgeFunc>
//    constexpr void bridges(Graph const& graph, BridgeFunc bridge) const noexcept {
//        _dfs_init(graph.size());
//
//        for (auto i: graph.nodes()) {
//            if (-1 == _tin(i)) {
//                _bridges_dfs(graph, bridge, i);
//            }
//        }
//    }
//};
//#include "io.cpp"

class Solution {
    auto _setup_graph(int n, std::vector<std::vector<int>> const& requests) {
        graph.reset(n);

        uint_t loops = 0;

        for(auto const& r: requests) {
            if (r[0] == r[1]) {
                // Loops can be handled immediately.
                ++loops;
            } else {
                const auto e = Edge(r[0], r[1], 0);
                const auto r = e.reversed(1);
                graph.add(e);
                if (!graph.connected(r)) {
                    // reversed arc in residual network
                    graph.connect(r, -1);
                }
            }
        }

//        // Bridges do not contribute, and fool our opimizer.
//        optimizer.remove_bridges(graph);

        return loops;
    }

public:
    Graph graph;
    Optimizer optimizer;

    int maximumRequests(int n, std::vector<std::vector<int>> const& requests) {
        uint_t loops = _setup_graph(n, requests);
        // return loops + optimizer.max_circulation(graph);
        return loops;
    }
};

#endif /* LC_1601_SOLUTION_HPP */
