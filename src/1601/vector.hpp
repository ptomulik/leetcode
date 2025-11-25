#ifndef 
#include <algorithm>

/**
 * Vector with inplace storage.
 */
template<typename T, size_t N> class Vector {
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
