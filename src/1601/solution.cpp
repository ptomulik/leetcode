// LeetCode Contest
//
// NO:      1601
// LEVEL:   HARD
// TITLE:   Maximum Number of Achievable Transfer Requests
// URL:     https://leetcode.com/problems/maximum-number-of-achievable-transfer-requests/description/
//
// STATUS:      Accepted
// SUBMISSION: https://leetcode.com/problems/maximum-number-of-achievable-transfer-requests/submissions/1836052186/
// RUNTIME:     0ms | Beats 100.00%
// MEMORY:      11.35 MB | Beats 99.39%

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <functional>
#include <type_traits>
#include <vector>

typedef int8_t snum_t;     // signed integer value
typedef uint8_t unum_t;    // unsigned integer value

constexpr unum_t NMAX = 20; //! Max supported no. nodes in graph.
constexpr unum_t FMAX = 16; //! Max supported total capacity (sum of arc capacities).

/**
 * Vector with inplace storage.
 */
template<size_t N, typename T> class Vector {
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
 * Square matrix with inplace storage.
 */
template<size_t N, typename T> class Matrix {
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

    // Remove node from graph (node indices and _n remain unchanged!)
    constexpr void remove(size_t i) noexcept {
        for (size_t j = 0; j < _n; ++j) {
            _mat[i][j] = 0;
            _mat[j][i] = 0;
        }
    }
};


/**
 * Flow matrix (a kind of adjacency matrix).
 */
typedef Matrix<NMAX, snum_t> FlowMat;

/**
 * Sorted sequence of numbers.
 */
template <size_t N, typename T>
class SortSeq {
    size_t _n;
    T _seq[N];

    constexpr ptrdiff_t _insertion_index(T v) const noexcept {
        if (0 == _n) {
            return 0;
        }

        if (v < _seq[0]) {
            return 0;
        }

        if (v > _seq[_n-1]) {
            return _n;
        }

        T lo = 0, hi = _n-1;
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

    constexpr ptrdiff_t _removal_index(T v) const noexcept {
        if (0 == _n) {
            return -1;
        }

        if (v < _seq[0] || v > _seq[_n-1]) {
            return -1;
        }

        T lo = 0, hi = _n-1;
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

    constexpr void _unshift(T v, size_t at = 0) noexcept {
        for (auto j = _n; j > at; --j) {
            _seq[j] = _seq[j-1];
        }
        _seq[at] = v;
        ++_n;
    }

    constexpr auto _shift(size_t at = 0) noexcept {
        auto v = _seq[at];

        for (T j = at + 1; j < _n; ++j) {
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

    constexpr void insert(T v) noexcept {
        auto const at = _insertion_index(v);
        if (-1 != at) {
            _unshift(v, at);
        }
    }

    constexpr void remove(T v) noexcept {
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

template<size_t N = NMAX>
using NodeSeq = SortSeq<N, unum_t>;

template <class TNodeLst> class AdjLstRef {
    unum_t    _i;
    TNodeLst& _list;
public:
    constexpr AdjLstRef(unum_t i, TNodeLst& list) noexcept: _i(i), _list(list) {}

    constexpr auto i() const noexcept {
        return _i;
    }

    constexpr TNodeLst& list() const noexcept {
        return _list;
    }

    constexpr operator TNodeLst& () const noexcept {
        return _list;
    }
};

template <class TNodeLst, class TNodeSeq> class AdjLstIterator {
    TNodeLst* _lists;
    TNodeSeq& _nodes;
    unum_t   _k;
public:
    constexpr AdjLstIterator(TNodeLst* lists, TNodeSeq& vtxseq, int pos = 0) noexcept
        : _lists(lists), _nodes(vtxseq), _k(pos)
    { }

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

    constexpr AdjLstRef<TNodeLst> operator* () const noexcept {
        auto _i = _nodes(_k);
        return AdjLstRef<TNodeLst>(_i, _lists[_i]);
    }
};

/**
 * Adjacency list for a graph.
 */
class AdjLst {
    NodeSeq<NMAX>   _nodes;
    SortSeq<NMAX-1,unum_t> _lists[NMAX];
public:
    constexpr void reset() noexcept {
        for(auto& list: _lists) {
            list.reset();
        }
    }

    constexpr auto& operator() (this auto&& self, unum_t i) noexcept {
        return self._lists[i];
    }

    constexpr void connect(unum_t i, unum_t j) noexcept {
        if (0 == _lists[i].size()) {
            _nodes.insert(i);
        }
        _lists[i].insert(j);
    }

    constexpr void disconnect(unum_t i, unum_t j) noexcept {
        _lists[i].remove(j);
        if (0 == _lists[i].size()) {
            _nodes.remove(i);
        }
    }

    constexpr void remove(unum_t i) noexcept {
        unum_t rj[NMAX];
        unum_t rn = 0;

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

    constexpr auto size() const noexcept {
        return _nodes.size();
    }

    constexpr auto degree(unum_t i) const noexcept {
        return _lists[i].size();
    }

    constexpr auto begin() const noexcept {
        return AdjLstIterator<SortSeq<NMAX-1,unum_t> const, NodeSeq<NMAX> const>(_lists, _nodes, 0);
    }

    constexpr auto end() const noexcept {
        return AdjLstIterator<SortSeq<NMAX-1,unum_t> const, NodeSeq<NMAX> const>(_lists, _nodes, _nodes.size());
    }
};

class Pair {
    snum_t _i;
    snum_t _j;
public:
    constexpr Pair(): _i(-1), _j(-1) {}
    constexpr Pair(snum_t i, snum_t j): _i(i), _j(j) {}
    constexpr snum_t i() const noexcept { return _i; }
    constexpr snum_t j() const noexcept { return _j; }
};

class Graph {
    NodeSeq<NMAX> _nodes;
    FlowMat       _flow;
    AdjLst        _outbound;
    AdjLst        _inbound;

    constexpr void _remove_arc(unum_t i, unum_t j) noexcept {
        _outbound.disconnect(i, j);
        _inbound.disconnect(j, i);
        if (0 == _outbound(i).size() && 0 == _inbound(i).size()) {
            _nodes.remove(i);
        }
        if (0 == _outbound(j).size() && 0 == _inbound(j).size()) {
            _nodes.remove(j);
        }
    }

public:

    constexpr Graph() noexcept: _flow(0) {}
    constexpr Graph(int n) noexcept: _flow(n) {}

    constexpr auto const& nodes() const noexcept {
        return _nodes;
    }

    constexpr auto& flow() const noexcept {
        return _flow;
    }

    constexpr auto flow(unum_t i, unum_t j) const noexcept {
        return _flow(i, j);
    }

    constexpr auto const& outbound() const noexcept {
        return _outbound;
    }

    constexpr auto const& outbound(unum_t i) const noexcept {
        return _outbound(i);
    }

    constexpr auto const& inbound() const noexcept {
        return _inbound;
    }

    constexpr auto const& inbound(unum_t i) const noexcept {
        return _inbound(i);
    }

    /**
     * Returns the value provided to the constructor or to reset().
     *
     * The returned value determines maximum possible node index. It's not
     * the actual number of nodes existing in the graph.
     */
    constexpr auto size() const noexcept {
        return _flow.size();
    }

    constexpr void reset(unum_t n) noexcept {
        _flow.reset(n);
        _outbound.reset();
        _inbound.reset();
    }

    constexpr void add(unum_t i, unum_t j) noexcept {
        add(i, j, 1);
    }

    constexpr void add(unum_t i, unum_t j, unum_t f) noexcept {
        _outbound.connect(i, j);
        _inbound.connect(j, i);
        _flow(i, j) += f;
        _nodes.insert(i);
        _nodes.insert(j);
    }

    constexpr void sub(unum_t i, unum_t j) noexcept {
        sub(i, j, 1);
    }

    constexpr void sub(unum_t i, unum_t j, unum_t f) noexcept {
        auto f2 = (_flow(i, j) -= f);
        if (0 == f2) {
            _remove_arc(i, j);
        }
    }

    constexpr void remove(unum_t i) noexcept {
        _flow.remove(i);
        _outbound.remove(i);
        _inbound.remove(i);
        _nodes.remove(i);
    }

    template<size_t N, typename T>
    constexpr void remove(Vector<N, T> const& nodes) noexcept {
        for (auto i: nodes) {
            remove(i);
        }
    }

    constexpr void disconnect(unum_t i, unum_t j) noexcept {
        _flow(i, j) = 0;
        _remove_arc(i, j);
    }

    constexpr void disconnect(Pair const& arc) noexcept {
        disconnect(arc.i(), arc.j());
    }

    template<size_t N, typename T>
    constexpr void disconnect(Vector<N, T> const& arcs) noexcept {
        for (auto arc: arcs) {
            disconnect(arc);
        }
    }

    constexpr unum_t outdegree(unum_t i) const noexcept {
        return _outbound.degree(i);
    }

    constexpr unum_t indegree(unum_t i) const noexcept {
        return _inbound.degree(i);
    }

    constexpr snum_t outflow(unum_t i) const noexcept {
        snum_t x = 0;
        for (auto j: _outbound(i)) {
            x += _flow(i, j);
        }
        return x;
    }

    constexpr snum_t inflow(unum_t i) const noexcept {
        snum_t x = 0;
        for (auto j: _inbound(i)) {
            x += _flow(j, i);
        }
        return x;
    }

    constexpr snum_t balance(unum_t i) const noexcept {
        return inflow(i) - outflow(i);
    }

    constexpr bool is_terminal(unum_t i) const noexcept {
        return 0 == outdegree(i) || 0 == indegree(i);
    }

    constexpr bool is_isolated(unum_t i) const noexcept {
        return 0 == outdegree(i) && 0 == indegree(i);
    }

    constexpr bool is_source_or_sink(unum_t i) const noexcept {
        return 0 == outflow(i) || 0 == inflow(i);
    }
};

/**
 * Encapsulates shortest paths solution.
 */
class ShortestPaths {
public:
    typedef Vector<NMAX, unum_t> Distances;
    typedef Vector<NMAX, snum_t> Predecessors;

private:
    size_t       _i;
    Distances    _distances;
    Predecessors _predecessors;
public:
    constexpr ShortestPaths() noexcept : _i(), _distances(), _predecessors() { }

    constexpr size_t i() const noexcept  {
        return _i;
    }

    constexpr size_t size() const noexcept {
        return _distances.size();
    }

    constexpr void reset(size_t n, size_t i) noexcept {
        _i = i;
        _distances.reset(n, NMAX);
        _predecessors.reset(n, -1);
    }

    constexpr auto&& distances(this auto&& self) noexcept {
        return self._distances;
    }

    constexpr auto&& distance(this auto&& self, size_t j) noexcept {
        return self._distances(j);
    }

    constexpr auto const& predecessors() const noexcept {
        return _predecessors;
    }

    constexpr auto&& predecessor(this auto&& self, size_t j) noexcept {
        return self._predecessors(j);
    }

    constexpr bool exists(size_t j) const noexcept {
        return NMAX != _distances(j);
    }

    template<typename Function, typename T>
    constexpr T reduce(size_t j, Function func, T value) const noexcept {
        if (exists(j)) {
            while(j != _i) {
                size_t k = predecessor(j);
                value = func(k, j, value);
                j = k;
            }
        }
        return value;
    }

    template<typename Function>
    constexpr void walk(size_t j, Function func) const noexcept {
        if (exists(j)) {
            while (j != _i) {
                size_t k = predecessor(j);
                func(k, j);
                j = k;
            }
        }
    }

    constexpr snum_t flow(size_t j, Graph const& graph) const noexcept {
        if(_i == j || !exists(j)) {
            return 0;
        }

        return reduce(j, [&graph](size_t i, size_t j, snum_t x) { return std::min(x, graph.flow(i, j)); }, FMAX);
    }
};

/**
 * Dijkstra shortest path with Dial modification (bucket).
 */
class Dijkstra {
    mutable Vector<NMAX, Vector<NMAX, unum_t>> _bucket;

public:
    constexpr void shortest_paths(Graph const& graph, size_t s, ShortestPaths& paths) const noexcept {
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
            unum_t i = stack->pop_back();

            unum_t dist = paths.distance(i) + 1;

            for (auto j: graph.outbound(i)) {
                if (dist < paths.distance(j)) {
                    _bucket(dist).push_back(j);
                    paths.distance(j) = dist;
                    paths.predecessor(j) = i;
                }
            }
        }
    }
};

/**
 * Tarjan's algorithm -- identifies strongly connected components in graph
 * (actually, we identify only bridges between them).
 */
class Tarjan {
public:
    typedef Vector<NMAX, bool> BoolVec;
    typedef Vector<NMAX, snum_t> SnumVec;
    typedef Vector<NMAX, unum_t> NodeStack;
private:
    mutable BoolVec _onstack;
    mutable SnumVec _low;
    mutable SnumVec _tin;
    mutable unum_t  _timer;
    mutable NodeStack _stack;

    template<typename BridgeFunc>
    constexpr void _bridges_dfs(
        Graph const& graph,
        BridgeFunc bridge,
        snum_t i,
        snum_t p = -1
    ) const noexcept {
        _tin(i) = _low(i) = _timer++;

        _stack.push_back(i);
        _onstack(i) = true;

        for (auto j: graph.outbound(i)) {
            if (_tin(j) == -1) {
                // Successor j has not yet been visited. Recurse on it.
                _bridges_dfs(graph, bridge, j, i);
                _low(i) = std::min(_low(i), _low(j));
            } else if (_onstack(j)) {
                // Successor j is on stack and hence in the current SCC.
                _low(i) = std::min(_low(i), _tin(j));
            } else {
                // If j is not on stack, then (i, j) is an arc pointing to an
                // SCC already found, and must be ignored. See below, regarding
                // the next line.
                bridge(Pair(i, j));
            }
        }

        // If i is a root node, pop the stack and generate an SCC.
        if (_low(i) == _tin(i)) {
            if (-1 != p) {
                bridge(Pair(p, i));
            }
            snum_t j;
            do {
                j = _stack.pop_back();
                _onstack(j) = false;
            } while(i != j);
        }
    }

    constexpr void _dfs_init(size_t n) const noexcept {
        _stack.clear();

        _onstack.reset(n, false);
        _tin.reset(n, -1);
        _low.reset(n, -1);

        _timer = 0;
    }

public:
    //! Identifies bridges between strongly connected components in graph.
    template<typename BridgeFunc>
    constexpr void bridges(Graph const& graph, BridgeFunc bridge) const noexcept {
        _dfs_init(graph.size());

        for (auto i: graph.nodes()) {
            if (-1 == _tin(i)) {
                _bridges_dfs(graph, bridge, i);
            }
        }
    }
};

/**
 * The workhorse class.
 */
class Optimizer {
public:
    typedef Vector<NMAX, unum_t> NodeVec;
    typedef Vector<NMAX, snum_t> SnumVec;
    typedef Vector<NMAX, Pair>   PairVec;
private:
    mutable NodeVec _nodes1;            //! Vector 1 holding nodes.
    mutable NodeVec _nodes2;            //! Vector 2 holding nodes.
    mutable SnumVec _snums1;            //! Vector holding signed integers
    mutable PairVec _pairs1;            //! Vector holding node pairs.
    mutable ShortestPaths _shortest1;   //! Structure of shortest paths
    Dijkstra _dijkstra;
    Tarjan _tarjan;
public:

    constexpr void find_deficit_and_excess_nodes(
        Graph const& graph,
        NodeVec& deficit,
        NodeVec& excess,
        SnumVec& balances
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

    constexpr void remove_bridges(Graph& graph) const noexcept {
        auto& bridges = _pairs1;
        bridges.clear();
        _tarjan.bridges(graph, [&bridges](Pair const& bridge) {
            bridges.push_back(bridge);
        });
        graph.disconnect(bridges);
    }

    constexpr unum_t max_circulation(Graph& graph) const noexcept {
        NodeVec& deficit = _nodes1;
        NodeVec& excess = _nodes2;
        SnumVec& balances = _snums1;
        ShortestPaths& shortest = _shortest1;

        // Find and remove minimal path flow.
        find_deficit_and_excess_nodes(graph, deficit, excess, balances);
        while (deficit.size() != 0) {
            for (auto ki = deficit.begin(); ki != deficit.end();) {
                auto k = *ki;
                _dijkstra.shortest_paths(graph, k, shortest);
                for (auto li = excess.begin(); li != excess.end();) {
                    auto l = *li;
                    if (shortest.exists(l)) {
                        // Flow f(P_{kl}) along the path P_{kl}
                        unum_t f = std::min(std::min((snum_t)-balances(k), balances(l)), shortest.flow(l, graph));
                        shortest.walk(l, [&graph, f](size_t i, size_t j) { graph.sub(i, j, f); });
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

    constexpr unum_t flow_cost(Graph const& graph) const noexcept {
        unum_t total = 0;
        for (auto const& ref: graph.outbound()) {
            auto i = ref.i();
            for (auto j: ref.list()) {
                total += graph.flow(i, j);
            }
        }
        return total;
    }
};

class Solution {
    auto _setup_graph(int n, std::vector<std::vector<int>> const& requests) {
        graph.reset(n);

        unum_t loops = 0;

        for(auto const& r: requests) {
            if (r[0] == r[1]) {
                // Loops can be handled immediately.
                ++loops;
            } else {
                graph.add(r[0], r[1]);
            }
        }

        // Bridges do not contribute, and fool our opimizer.
        optimizer.remove_bridges(graph);

        return loops;
    }

public:
    Graph graph;
    Optimizer optimizer;

    int maximumRequests(int n, std::vector<std::vector<int>> const& requests) {
        unum_t loops = _setup_graph(n, requests);
        return loops + optimizer.max_circulation(graph);
    }
};
