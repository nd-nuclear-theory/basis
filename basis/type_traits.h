/****************************************************************

  type_traits.h

  Language: C++17

  Patrick J. Fasano
  Argonne National Laboratory

  + 09/25/24: Created.

****************************************************************/

#ifndef BASIS_TYPE_TRAITS_H_
#define BASIS_TYPE_TRAITS_H_

#include <limits>
#include <type_traits>

namespace basis {

  static constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();
  // Flag value for missing target in index lookups.

  ////////////////////////////////////////////////////////////////
  // type trait helpers
  //
  // for use with SFINAE
  ////////////////////////////////////////////////////////////////

  template<typename, typename, typename = void>
  struct is_subspace
      : std::false_type
  {};

  template<typename Subspace, typename Space>
  struct is_subspace<
      Subspace, Space,
      std::enable_if_t<std::is_same_v<Subspace, typename Space::SubspaceType>>
    >
      : std::true_type
  {};

  template<typename Subspace, typename Space>
  struct is_subspace<
      Subspace, Space,
      std::enable_if_t<
          is_subspace<Subspace, typename Space::SubspaceType>::value
        >
    >
      : std::true_type
  {};

  template<typename T, typename U>
  constexpr bool is_subspace_v = is_subspace<T, U>::value;

}  // namespace basis

#endif  // BASIS_TYPE_TRAITS_H_
