/****************************************************************

  map.h

  Language: C++17

  Patrick J. Fasano
  Argonne National Laboratory

  + 10/01/24: Created, defining customization points for equal_to, less,
    and hash.

****************************************************************/

#ifndef BASIS_MAP_H_
#define BASIS_MAP_H_

#ifdef BASIS_HASH

#include <unordered_map>

#if defined(BASIS_BOOST_HASH)
#include <boost/container_hash/hash.hpp>
#elif defined(BASIS_STD_HASH)
#endif

#else  // BASIS_HASH

#include <map>

#endif  // BASIS_HASH

namespace basis {

// forward declarations
template<typename T> class equal_to;
template<typename T> class less;
template<typename T> class hash;


template<typename T> class equal_to : public std::equal_to<T> {};

template<typename T> class less : public std::less<T> {};

#ifdef BASIS_HASH

#if defined(BASIS_BOOST_HASH)
template<typename T> class hash : public boost::hash<T> {};
#elif defined(BASIS_STD_HASH)
template<typename T> class hash : public std::hash<T> {};
#endif

template<typename Key, typename T>
using map = std::unordered_map<Key, T, basis::hash<Key>, basis::equal_to<Key>>;

#else  // BASIS_HASH
template<typename Key, typename T>
using map = std::map<Key, T, basis::less<Key>>;

#endif // BASIS_HASH

}  // namespace basis

#endif  // BASIS_MAP_H_
