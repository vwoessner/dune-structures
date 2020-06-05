#ifndef DUNE_STRUCTURES_ENUMERATE_HH
#define DUNE_STRUCTURES_ENUMERATE_HH

/** This is a C++ equivalent of Python's enumerate to be used with
 *  structured bindings. It is taken from http://reedbeta.com/blog/python-like-enumerate-in-cpp17/
 *  Copyright is with the author Nathan Reed which has granted the public rights
 *  to use it in their code.
 *
 *  Usage example: https://godbolt.org/z/XyxYDS
 */

#include <tuple>


template <typename T,
          typename TIter = decltype(std::begin(std::declval<T>())),
          typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T && iterable)
{
    struct iterator
    {
        size_t i;
        TIter iter;
        bool operator != (const iterator & other) const { return iter != other.iter; }
        void operator ++ () { ++i; ++iter; }
        auto operator * () const { return std::tie(i, *iter); }
    };
    struct iterable_wrapper
    {
        T iterable;
        auto begin() { return iterator{ 0, std::begin(iterable) }; }
        auto end() { return iterator{ 0, std::end(iterable) }; }
    };
    return iterable_wrapper{ std::forward<T>(iterable) };
}

#endif
