#ifndef DUNE_BLOCKLAB_UTILITIES_TUPLECAT_HH
#define DUNE_BLOCKLAB_UTILITIES_TUPLECAT_HH

/** The C++ standard library has std::tuple_cat to
 *  concatenate tuple, but it does not have std::tuple_cat_t
 *  which just gives you the type of the concatenated tuple.
 *  Implementation from https://stackoverflow.com/a/53398815
 */

#include<tuple>

namespace Dune::BlockLab {

  template<typename... T>
  using tuple_cat_t = decltype(std::tuple_cat(std::declval<T>()...));

} // namespace Dune::BlockLab

#endif
