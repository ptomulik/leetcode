#ifndef LC_1601_SORTSEQ_HPP
#define LC_1601_SORTSEQ_HPP

#include <cstddef>

namespace Templates {

/**
 * Sorted sequence of numbers.
 */
template <typename T, size_t N>
class SortSeq {
public:
    typedef T value_type;
    typedef T* iterator;
    typedef T const* const_iterator;
    typedef T& reference;
    typedef T const& const_reference;
private:
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

} /* namespace Templates */

#endif /* LC_1601_SORTSEQ_HPP */
