#ifndef LC_1601_ADJLST_HPP
#define LC_1601_ADJLST_HPP

#include "sortseq.hpp"
#include "edge.hpp"
#include "encoded.hpp"

namespace Templates {

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

} /* Templates */

#endif /* LC_1601_ADJLST_HPP */
