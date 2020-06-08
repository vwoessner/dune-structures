#ifndef DUNE_BLOCKLAB_UTILITIES_UNIQUEVARIANT_HH
#define DUNE_BLOCKLAB_UTILITIES_UNIQUEVARIANT_HH

/** This header introduces unique_variant<T...>.
 *  It is an alternative to std::variant that allows duplication of entries
 *  in the variadic template parameter list. Forwards to std::variant
 *  after removing the duplicates in a meta-programming nightmare.
 */

#include<variant>
#include<tuple>
#include<type_traits>


namespace Dune::BlockLab {

  namespace impl {
    // This is available from std:: starting with C++2a
    template<typename T>
    struct type_identity
    {
      using type = T;
    };

    template <typename T, typename... Ts>
    struct unique : type_identity<T> {};

    template <typename... Ts, typename U, typename... Us>
    struct unique<std::tuple<Ts...>, U, Us...>
	: std::conditional_t<(std::is_same_v<U, Ts> || ...)
			   , unique<std::tuple<Ts...>, Us...>
			   , unique<std::tuple<Ts..., U>, Us...>> {};

    template <typename... Ts>
    using unique_tuple = typename unique<std::tuple<>, Ts...>::type;


    template<typename T>
    struct unique_variant_impl
    {};

    template<typename... Ts>
    struct unique_variant_impl<std::tuple<Ts...>>
    {
      using type = std::variant<Ts...>;
    };
  }

  template<typename... A>
  using unique_variant = typename impl::unique_variant_impl<impl::unique_tuple<A...>>::type;

} // namespace Dune::BlockLab

#endif
