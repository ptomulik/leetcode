#ifndef LC_1601_ENCODED_HPP
#define LC_1601_ENCODED_HPP

#include "sortseq.hpp"
#include <type_traits>

namespace Templates {

/**
 * Endpoint value reference.
 */
template<typename ObjectT>
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

template<typename ObjectT>
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

} /* namespace Templates */

#endif /* LC_1601_ENCODED_HPP */
