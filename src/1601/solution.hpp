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

namespace Templates {

/**
 * Vector of size N with inplace storage.
 */
template <typename T, size_t N>
class Vector {
    size_t _n;
    T      _vec[N];
public:

    typedef T value_type;

    constexpr Vector() noexcept: _n(0u) {}

    constexpr size_t size() const noexcept {
        return _n;
    }

    constexpr auto&& operator()(this auto&& self, size_t i) noexcept {
        return self._vec[i];
    }

    constexpr void reset(size_t n, T const& v = T()) noexcept {
        resize(n);
        fill(v);
    }

    constexpr void resize(size_t n) noexcept {
        _n = n;
    }

    constexpr T* erase(T* pos) noexcept {
        auto const* const e = end();
        if (pos >= begin() && pos < e) {
            for (T* ptr = pos; ptr < e; ++ptr) {
                *ptr = *(ptr + 1);
            }
            --_n;
        }
        return pos;
    }

    constexpr void fill(T const& v) noexcept {
        std::fill_n(&_vec[0], _n, v);
    }

    constexpr void clear() noexcept {
        _n = 0;
    }

    constexpr void push_back(T value) noexcept {
        _vec[_n++] = value;
    }

    constexpr T pop_back() noexcept {
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
    size_t _n;
    T      _mat[N][N];
public:
    constexpr Matrix() noexcept: _n(0u) { }
    constexpr Matrix(size_t n) noexcept : _n(n) { }
    constexpr Matrix(size_t n, T const& x) noexcept : _n(n), _mat{x} { }

    constexpr void reset(size_t n, T const& v = T()) noexcept {
        resize(n);
        fill(v);
    }

    constexpr void resize(size_t n) noexcept {
        _n = n;
    }

    constexpr void fill(T const& v) noexcept {
        std::fill_n(&_mat[0][0], N * N, v);
    }

    constexpr auto size() const noexcept {
        return _n;
    }

    constexpr auto&& operator () (this auto&& self, size_t i, size_t j) noexcept {
        return self._mat[i][j];
    }

    // Set value to all elements of row and column i.
    constexpr void set(size_t i, T const& v) noexcept {
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
 * Edge endpoint.
 *
 * Encapsulates node index i of node at endpoint and number m (multi-graph).
 */
template<typename T, size_t M>
struct Endpoint {
    typedef T value_type;

    value_type node;
    value_type edge;

    constexpr Endpoint(value_type encoded) noexcept : node(encoded / M), edge(encoded % M) {}
    constexpr Endpoint(value_type node, value_type edge) noexcept : node(node), edge(edge) {}

    constexpr value_type encoded() const noexcept {
        return M * node + edge;
    }
};

template <typename T>
struct TypeTraits {
    typedef typename std::make_unsigned<T>::type unsigned_int;
    typedef typename std::make_signed<T>::type signed_int;
    typedef unsigned_int node_type;
};

template <typename T, size_t M>
struct EdgeTraits {
    typedef typename TypeTraits<T>::node_type node_type;
    typedef Endpoint<T, M> endpoint_type;
    typedef endpoint_type tail_type;
    typedef endpoint_type head_type;
};

template<typename T, size_t M>
struct Edge {
    constexpr static const size_t max_m = M;

    typedef T value_type;
    typedef typename EdgeTraits<T, M>::node_type node_type;
    typedef typename EdgeTraits<T, M>::endpoint_type endpoint_type;
    typedef typename EdgeTraits<T, M>::tail_type tail_type;
    typedef typename EdgeTraits<T, M>::head_type head_type;

    node_type i;
    node_type j;
    T m;

    constexpr Edge(node_type i, node_type j, T m) noexcept : i(i), j(j), m(m) {}
    constexpr Edge(node_type i, head_type h) noexcept : i(i), j(h.node), m(h.edge) {}
    constexpr Edge(tail_type t, node_type j) noexcept : i(t.node), j(j), m(t.edge) {}

    constexpr auto tail() const noexcept {
        return tail_type(i, m);
    }

    constexpr auto head() const noexcept {
        return head_type(j, m);
    }

    constexpr auto reversed() const noexcept {
        return Edge<T, M>(j, i, m);
    }

    constexpr auto reversed(T m) const noexcept {
        return Edge<T, M>(j, i, m);
    }
};

/**
 * Reference to an "encoded" object (forgive me naming).
 *
 * Encapsulates a reference to an "encoded" value. Provides methods for
 * converting between the encoded value type and ObjectT.
 */
template <typename ObjectT>
class EncodedObjectRef {
public:
    typedef typename std::remove_cv<ObjectT>::type object_type;
    typedef typename std::remove_cv<ObjectT>::type value_type;
    typedef typename ObjectT::value_type encoded_type;
    typedef typename std::conditional<
        std::is_const<ObjectT>::value,
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

    constexpr operator value_type() const noexcept {
        return value_type(_ref);
    }

    constexpr encoded_type encoded() const noexcept {
        return _ref;
    }
};

template <typename ObjectT>
class EncodedObjectIter {
public:
    typedef typename std::remove_cv<ObjectT>::type object_type;
    typedef typename std::remove_cv<ObjectT>::type value_type;
    typedef typename ObjectT::value_type encoded_type;
    typedef typename std::conditional<
        std::is_const<ObjectT>::value,
        typename std::add_const<encoded_type>::type,
        typename std::remove_const<encoded_type>::type
    >::type* encoded_pointer;
    typedef EncodedObjectIter<ObjectT> iterator;
    typedef EncodedObjectRef<ObjectT> reference;
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

template <typename ObjectT, size_t N>
class EncodedObjectSortSeq
{
public:
    typedef ObjectT object_type;
    typedef ObjectT value_type;
    typedef typename ObjectT::value_type encoded_type;
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

template <typename T, typename HeadLstT>
class AdjLstRef {
public:
    typedef T node_type;
    typedef HeadLstT head_list_type;
private:
    node_type _i;
    head_list_type& _list;
public:
    constexpr AdjLstRef(node_type i, head_list_type& list) noexcept
        : _i(i), _list(list) {
    }

    constexpr auto i() const noexcept {
        return _i;
    }

    constexpr head_list_type& list() const noexcept {
        return _list;
    }
};

template <typename NodeLstT, typename HeadLstT>
class AdjLstIter {
public:
    typedef NodeLstT node_list_type;
    typedef HeadLstT head_list_type;
    typedef typename NodeLstT::value_type node_type;
    typedef typename HeadLstT::value_type head_type;
    typedef AdjLstRef<node_type, head_list_type> reference;
private:
    node_list_type& _nodes;
    head_list_type* _lists;
    size_t _k;
public:
    constexpr AdjLstIter(node_list_type& nodes, head_list_type* lists, size_t k) noexcept
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
 * N - max number of nodes in graph.
 */
template <typename EdgeT, size_t N>
class AdjLst {
public:
    typedef EdgeT edge_type;
    typedef typename EdgeT::value_type value_type;
    typedef typename EdgeT::tail_type tail_type;
    typedef typename EdgeT::head_type head_type;

    typedef SortSeq<value_type, N> node_list_type;
    typedef EncodedObjectSortSeq<head_type, N> head_list_type;

    typedef AdjLstIter<node_list_type, head_list_type> iterator;
    typedef AdjLstIter<const node_list_type, const head_list_type> const_iterator;
private:
    node_list_type _nodes;
    head_list_type _lists[N];
public:
    constexpr void reset() noexcept {
        for(auto& list: _lists) {
            list.reset();
        }
    }

    constexpr void connect(edge_type edge) noexcept {
        if (0 == _lists[edge.i].size()) {
            _nodes.insert(edge.i);
        }
        _lists[edge.i].insert(edge.head());
    }

    constexpr void disconnect(edge_type edge) noexcept {
        _lists[edge.i].remove(edge.head());
        if (0 == _lists[edge.i].size()) {
            _nodes.remove(edge.i);
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

template <typename T, size_t N, size_t M>
struct GraphTraits {
    typedef typename TypeTraits<T>::node_type node_type;
    typedef typename TypeTraits<T>::unsigned_int unsigned_int;
    typedef typename TypeTraits<T>::signed_int signed_int;

    typedef signed_int flow_type;
    typedef signed_int cost_type;
    typedef signed_int caps_type;
    // typedef signed_int pots_type;

    typedef std::array<Matrix<bool, N>, M>   adjc_matrix;
    typedef std::array<Matrix<flow_type, N>, M> flow_matrix;
    typedef std::array<Matrix<cost_type, N>, M> cost_matrix;
    typedef std::array<Matrix<caps_type, N>, M> caps_matrix;
    // typedef std::array<Vector<pots_type, N>, M> pots_vector;

    typedef Edge<node_type, M> edge_type;
    typedef typename EdgeTraits<T, M>::endpoint_type endpoint_type;
    typedef typename EdgeTraits<T, M>::head_type head_type;
    typedef typename EdgeTraits<T, M>::tail_type tail_type;
    typedef AdjLst<edge_type, N> adjc_list;
    typedef SortSeq<node_type, N> node_list;
};

template <typename T, size_t N, size_t M>
class Graph {
public:
    typedef typename GraphTraits<T, N, M>::node_type node_type;

    typedef typename TypeTraits<T>::unsigned_int unsigned_int;
    typedef typename TypeTraits<T>::signed_int signed_int;
    typedef typename GraphTraits<T, N, M>::flow_type flow_type;
    typedef typename GraphTraits<T, N, M>::cost_type cost_type;
    typedef typename GraphTraits<T, N, M>::caps_type caps_type;
    // typedef typename GraphTraits<T, N, M>::pots_type pots_type;

    typedef typename GraphTraits<T, N, M>::adjc_matrix adjc_matrix;
    typedef typename GraphTraits<T, N, M>::flow_matrix flow_matrix;
    typedef typename GraphTraits<T, N, M>::cost_matrix cost_matrix;
    typedef typename GraphTraits<T, N, M>::caps_matrix caps_matrix;
    // typedef typename GraphTraits<T, N, M>::pots_vector pots_vector;

    typedef typename GraphTraits<T, N, M>::edge_type edge_type;
    typedef typename GraphTraits<T, N, M>::endpoint_type endpoint_type;
    typedef typename GraphTraits<T, N, M>::head_type head_type;
    typedef typename GraphTraits<T, N, M>::tail_type tail_type;

    typedef typename GraphTraits<T, N, M>::adjc_list adjc_list;
    typedef typename GraphTraits<T, N, M>::node_list node_list;
private:

    node_list   _nodes;
    adjc_list   _outbound;
    adjc_list   _inbound;
    adjc_matrix _adjc; //! adjacency matrix
    flow_matrix _flow; //! arc flows x_{ij}
    cost_matrix _cost; //! arc costs c_{ij}
    caps_matrix _caps; //! (residual) capacities r_{ij}
    // pots_vector _pots; //! potentials pi_{i}

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
        return _adjc[0].size();
    }

    constexpr static cost_type cost_max() noexcept {
        return std::numeric_limits<cost_type>::max();
    }

    constexpr static caps_type caps_max() noexcept {
        return std::numeric_limits<caps_type>::max();
    }

    constexpr void reset(
        size_t n,
        flow_type x = 0,
        cost_type c = cost_max(),
        caps_type u = 0/*,
        pots_type pi = 0*/
    ) noexcept {
        for (size_t m = 0; m < M; ++m) {
            _adjc[m].reset(n, 0);
            _flow[m].reset(n, x);
            _cost[m].reset(n, c);
            _caps[m].reset(n, u);
            // _pots[m].reset(n, pi);
        }
        _outbound.reset();
        _inbound.reset();
    }

    constexpr void connect(edge_type e, cost_type c = 1) noexcept {
        adjc(e) = true;
        cost(e) = c;
        _outbound.connect(e);
        _inbound.connect(e.reversed());
        _nodes.insert(e.i);
        _nodes.insert(e.j);
    }

    constexpr void disconnect(edge_type e, cost_type c = cost_max()) noexcept {
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

    constexpr bool connected(edge_type e) const noexcept {
        return adjc(e);
    }

    constexpr void remove(node_type i) noexcept {
        for (size_t m = 0; m < M; ++m) {
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
    constexpr void add(edge_type e, caps_type du = 1) noexcept {
        if (!connected(e)) {
            connect(e);
        }

        if (caps(e) < caps_max() - du) {
            caps(e) += du;
        } else {
            caps(e) = caps_max();
        }
    }

    // Remve capacity from edge
    constexpr void sub(edge_type e, caps_type du = 1) noexcept {
        if (du < caps(e)) {
            caps(e) -= du;
        } else {
            caps(e) = 0;
        }
    }

    constexpr bool adjc(edge_type e) const noexcept {
        return _adjc[e.m](e.i, e.j);
    }

    constexpr bool& adjc(edge_type e) noexcept {
        return _adjc[e.m](e.i, e.j);
    }

    constexpr flow_type flow(edge_type e) const noexcept {
        return _flow[e.m](e.i, e.j);
    }

    constexpr flow_type& flow(edge_type e) noexcept {
        return _flow[e.m](e.i, e.j);
    }

    constexpr cost_type cost(edge_type e) const noexcept {
        return _cost[e.m](e.i, e.j);
    }

    constexpr cost_type& cost(edge_type e) noexcept {
        return _cost[e.m](e.i, e.j);
    }

    constexpr caps_type caps(edge_type e) const noexcept {
        return _caps[e.m](e.i, e.j);
    }

    constexpr caps_type& caps(edge_type e) noexcept {
        return _caps[e.m](e.i, e.j);
    }

    constexpr auto outdegree(node_type i) const noexcept {
        return _outbound.degree(i);
    }

    constexpr auto indegree(node_type i) const noexcept {
        return _inbound.degree(i);
    }

    constexpr auto outflow(node_type i) const noexcept {
        flow_type x = 0;
        for (auto head: _outbound(i)) {
            x += flow(edge_type(i, head));
        }
        return x;
    }

    constexpr auto inflow(node_type i) const noexcept {
        flow_type x = 0;
        for (auto head: _inbound(i)) {
            x += flow(edge_type(head, i));
        }
        return x;
    }

    constexpr auto balance(node_type i) const noexcept {
        return inflow(i) - outflow(i);
    }
};

/**
 * Encapsulates shortest paths solution.
 */
template <typename T, size_t N, size_t M>
class Paths {
public:
    typedef TypeTraits<T>::unsigned_int unsigned_int;
    typedef TypeTraits<T>::signed_int signed_int;
    typedef EdgeTraits<T, M>::tail_type tail_type;
    typedef GraphTraits<T, N, M>::edge_type edge_type;

    typedef Vector<unsigned_int, N> dist_vector;
    typedef Vector<unsigned_int, N> pred_vector;
    typedef Graph<T, N, M> graph_type;

private:
    size_t      _i;
    dist_vector _distances;
    pred_vector _predecessors;
public:
    constexpr Paths() noexcept : _i(), _distances(), _predecessors() { }

    constexpr size_t i() const noexcept  {
        return _i;
    }

    constexpr size_t size() const noexcept {
        return _distances.size();
    }

    constexpr void reset(size_t n, size_t i) noexcept {
        _i = i;
        _distances.reset(n, N);
        _predecessors.reset(n, -1);
    }

    constexpr auto&& distances(this auto&& self) noexcept {
        return self._distances;
    }

    constexpr auto&& distance(this auto&& self, size_t j) noexcept {
        return self._distances(j);
    }

//    constexpr auto predecessors() const noexcept {
//        return EncodedObjectIter<const tail_type>(_predecessors.begin());
//    }
//
//    constexpr auto predecessors() noexcept {
//        return EncodedObjectIter<tail_type>(_predecessors.begin());
//    }

    constexpr auto predecessor(size_t j) const noexcept {
        return EncodedObjectRef<const tail_type>(_predecessors(j));
    }

    constexpr auto predecessor(size_t j) noexcept {
        return EncodedObjectRef<tail_type>(_predecessors(j));
    }

    constexpr bool exists(size_t j) const noexcept {
        return N != _distances(j);
    }

    template<typename Function, typename U>
    constexpr U reduce(size_t j, Function func, U value) const noexcept {
        if (exists(j)) {
            while(j != _i) {
                tail_type t = predecessor(j);
                value = func(edge_type(t, j), value);
                j = t.node;
            }
        }
        return value;
    }

    template<typename Function>
    constexpr void walk(size_t j, Function func) const noexcept {
        if (exists(j)) {
            while (j != _i) {
                tail_type t = predecessor(j);
                func(edge_type(t, j));
                j = t.node;
            }
        }
    }

    constexpr signed_int flow(size_t j, graph_type const& graph) const noexcept {
        if(_i == j || !exists(j)) {
            return 0;
        }

        const signed_int fmax = std::numeric_limits<signed_int>::max();

        return reduce(j, [&graph](edge_type e, signed_int x) { return std::min(x, graph.flow(e)); }, fmax);
    }
};

/**
 * Dijkstra shortest path with Dial modification (bucket).
 */
template <typename T, size_t N, size_t M>
class Dijkstra {
public:
    typedef typename GraphTraits<T, N, M>::edge_type edge_type;
    typedef typename GraphTraits<T, N, M>::node_type node_type;
    typedef typename GraphTraits<T, N, M>::cost_type cost_type;
    typedef typename EdgeTraits<T, M>::head_type head_type;
    typedef Graph<T, N, M> graph_type;
    typedef Paths<T, N, M> paths_type;
    typedef cost_type dist_type;

private:
    mutable Vector<Vector<dist_type, N>, N> _bucket;

public:
    constexpr void shortest_paths(graph_type const& graph, size_t s, paths_type& paths) const noexcept {
        paths.reset(graph.size(), s);

        _bucket.reset(graph.size());

        paths.distance(s) = 0;
        _bucket(0).push_back(s);

        auto first_nonempty = [](auto& container) {
            return std::find_if(container.begin(), container.end(), [](auto const& d) {
                return d.size() > 0;
            });
        };

        for (auto stack = first_nonempty(_bucket); stack != _bucket.end(); stack = first_nonempty(_bucket)) {
            dist_type i = stack->pop_back();

            for (head_type h: graph.outbound(i)) {
                const node_type j = h.node;
                dist_type dist = paths.distance(i) + graph.cost(edge_type(i, j, h.edge));
                if (dist < paths.distance(j)) {
                    _bucket(dist).push_back(j);
                    paths.distance(j) = dist;
                    paths.predecessor(j) = Endpoint<node_type, M>(i, h.edge);
                }
            }
        }
    }
};

/**
 * The workhorse class.
 */
template <typename T, size_t N, size_t M>
class Optimizer {
public:
    typedef typename TypeTraits<T>::unsigned_int unsigned_int;
    typedef typename TypeTraits<T>::signed_int signed_int;
    typedef typename GraphTraits<T, N, M>::node_type node_type;
    typedef typename GraphTraits<T, N, M>::endpoint_type endpoint_type;
    typedef typename GraphTraits<T, N, M>::head_type head_type;
    typedef typename GraphTraits<T, N, M>::tail_type tail_type;
    typedef typename GraphTraits<T, N, M>::edge_type edge_type;
    typedef Graph<T, N, M> graph_type;
    typedef Vector<node_type, N> node_vector;
    typedef Vector<signed_int, N> sint_vector;
    typedef Vector<unsigned_int, N> uint_vector;
private:
    mutable node_vector _nodes1;    //! Vector 1 holding nodes.
    mutable node_vector _nodes2;    //! Vector 2 holding nodes.
    mutable sint_vector _snums1;    //! Vector holding signed integers
//    mutable pair_vector _pairs1;  //! Vector holding node pairs.
    mutable Paths<T, N, M> _paths1; //! Structure of shortest paths

    Dijkstra<T, N, M> _dijkstra;
//    Tarjan _tarjan;
public:

    constexpr void find_deficit_and_excess_nodes(
        graph_type const& graph,
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

    constexpr unsigned_int max_circulation(graph_type& graph) const noexcept {
        auto& deficit = _nodes1;
        auto& excess = _nodes2;
        auto& balances = _snums1;
        auto& paths = _paths1;

        // Find and remove minimal path flow.
        find_deficit_and_excess_nodes(graph, deficit, excess, balances);
//        std::cout << "deficit: [" << deficit << "]; excess: [" << excess << "]" << std::endl;
        while (deficit.size() != 0) {
            for (auto ki = deficit.begin(); ki != deficit.end();) {
                auto k = *ki;
                _dijkstra.shortest_paths(graph, k, paths);
//                std::cout << "sortest: " << paths.i() << " -> " << paths.distances() << ";" << paths.predecessors() << std::endl;
                for (auto li = excess.begin(); li != excess.end();) {
                    auto l = *li;
                    if (paths.exists(l)) {
                        // Flow f(P_{kl}) along the path P_{kl}
                        unsigned_int f = std::min(std::min((signed_int)-balances(k), balances(l)), paths.flow(l, graph));
                        paths.walk(l, [&graph, f](edge_type e) {
                            const auto r = e.reversed((e.m + 1) %M);
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

    constexpr unsigned_int flow_cost(graph_type const& graph) const noexcept {
        unsigned_int total = 0;
        for (auto const& ref: graph.outbound()) {
            auto i = ref.i();
            for (head_type h: ref.list()) {
                edge_type e(i, h);
                total += graph.cost(e) * graph.flow(e);
            }
        }
        return total;
    }
};

} /* namespace Templates */

typedef uint8_t unum_t; // unsigned integer value
typedef typename std::make_signed<unum_t>::type snum_t; // signed integer value

constexpr unum_t NMAX = 20; //! Max supported no. nodes in graph.
constexpr unum_t FMAX = 16; //! Max supported total capacity (sum of arc capacities).

typedef Templates::Graph<unum_t, NMAX, 2> Graph;
typedef Templates::Optimizer<unum_t, NMAX, 2> Optimizer;
typedef Templates::Edge<unum_t, 2> Edge;

/**
 * Tarjan's algorithm -- identifies strongly connected components in graph
 * (actually, we identify only bridges between them).
 */
//class Tarjan {
//public:
//    typedef Vector<bool, NMAX> BoolVec;
//    typedef Vector<snum_t, NMAX> SnumVec;
//    typedef Vector<unum_t, NMAX> NodeStack;
//private:
//    mutable BoolVec _onstack;
//    mutable SnumVec _low;
//    mutable SnumVec _tin;
//    mutable unum_t  _timer;
//    mutable NodeStack _stack;
//
//    template<typename BridgeFunc>
//    constexpr void _bridges_dfs(
//        Graph const& graph,
//        BridgeFunc bridge,
//        snum_t i,
//        snum_t p = -1
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
//            snum_t j;
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

        unum_t loops = 0;

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
        unum_t loops = _setup_graph(n, requests);
        // return loops + optimizer.max_circulation(graph);
        return loops;
    }
};

#endif /* LC_1601_SOLUTION_HPP */
