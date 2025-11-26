#ifndef LC_1601_MATRIX_HPP
#define LC_1601_MATRIX_HPP

#include <algorithm>

namespace Templates {

/**
 * Square matrix with inplace storage.
 */
template<typename T, size_t N> class Matrix {
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

} /* namespace Templates */

#endif /* LC_1601_MATRIX_HPP */
