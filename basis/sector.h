/****************************************************************

  sectors.h

  Language: C++17

  Mark A. Caprio
  University of Notre Dame
  Patrick J. Fasano
  Argonne National Laboratory

  + 09/25/24: Split basis.h into separate headers.
  + 09/27/24: Move implementations into impl namespace to hide template
    parameters which should be invisible to library users.

****************************************************************/

#ifndef BASIS_SECTORS_H_
#define BASIS_SECTORS_H_

#include <cassert>
#include <cstddef>
#include <iterator>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// for DebugStr
#include <iosfwd>
#include <sstream>
#include <string>

#include <mcutils/meta.h>

#ifdef BASIS_ALLOW_DEPRECATED
// emit warnings on deprecated
#include "mcutils/deprecated.h"
#endif

#include "basis/type_traits.h"
#include "basis/map.h"

namespace basis {

  // SectorDirection -- sector direction specifier
  //
  //   kCanonical : specifies bra_sector_index <= ket_sector_index
  //   kBoth : specifies that both directions are allowed
  //
  // Note: It is up to the derived class constructor to accept a
  // sector direction argument and to implement this symmetry
  // constraint, if so decided.

  enum class SectorDirection {kCanonical,kBoth};

namespace impl {

  ////////////////////////////////////////////////////////////////
  // sector indexing
  ////////////////////////////////////////////////////////////////

  // Note: A "sector" is a pair of subspaces (which thus define a
  // "sector" or block in the matrix representation of an operator on
  // the space).
  //
  // The sector may also be labeled with a multiplicity index.  The
  // concept of a multiplicity index on a sector applies when the
  // symmetry group has outer multiplicites, so reduced matrix
  // elements are labeled not only by the bra and ket but also by a
  // multiplicity index.


  // declare BaseSectors template
  template <
      typename tBraSubspaceType, typename tKetSubspaceType = tBraSubspaceType,
      bool = std::is_same_v<tBraSubspaceType, tKetSubspaceType>
    >
    class BaseSector;

  /// BaseSector -- provide storage of information for a single sector

  // Here we specialize the template for the case where the bra and ket
  // subspaces are of different types. For the case where they are the same
  // (below), we will inherit from this specialization.
  template <typename tBraSubspaceType, typename tKetSubspaceType>
    class BaseSector<tBraSubspaceType, tKetSubspaceType, false>
    // Store indexing and subspace reference information for a sector.
    //
    // Allows for multiplicity index on sector, for symmetry sectors
    // of groups with outer multiplicities.
    {

      ////////////////////////////////////////////////////////////////
      // typedefs
      ////////////////////////////////////////////////////////////////

      public:
      using BraSubspaceType = tBraSubspaceType;
      using KetSubspaceType = tKetSubspaceType;
      typedef std::tuple<std::size_t,std::size_t,std::size_t> KeyType;

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSector()
        : bra_subspace_index_{}, ket_subspace_index_{}, multiplicity_index_{}
      {}
      // default constructor

      BaseSector(
          std::size_t bra_subspace_index, std::size_t ket_subspace_index,
          std::shared_ptr<const BraSubspaceType> bra_subspace_ptr,
          std::shared_ptr<const KetSubspaceType> ket_subspace_ptr,
          std::size_t multiplicity_index=1
        )
        : bra_subspace_index_{bra_subspace_index},
          ket_subspace_index_{ket_subspace_index},
          multiplicity_index_{multiplicity_index},
          bra_subspace_ptr_{std::move(bra_subspace_ptr)},
          ket_subspace_ptr_{std::move(ket_subspace_ptr)}
      {}

      ////////////////////////////////////////////////////////////////
      // accessors
      ////////////////////////////////////////////////////////////////

      inline KeyType Key() const
      // Return tuple key identifying sector for sorting/lookup
      // purposes.
      {
        return KeyType(bra_subspace_index(),ket_subspace_index(),multiplicity_index());
      }

      std::size_t bra_subspace_index() const {return bra_subspace_index_;}
      std::size_t ket_subspace_index() const {return ket_subspace_index_;}
      // Return integer index of bra/ket subspace.

      const BraSubspaceType& bra_subspace() const {return *bra_subspace_ptr_;}
      const KetSubspaceType& ket_subspace() const {return *ket_subspace_ptr_;}
      // Return reference to bra/ket subspace.

      std::shared_ptr<const BraSubspaceType> bra_subspace_ptr() const
      {
        return bra_subspace_ptr_;
      }
      std::shared_ptr<const KetSubspaceType> ket_subspace_ptr() const
      {
        return ket_subspace_ptr_;
      }
      // Return shared pointer to bra/ket subspace.

      std::size_t multiplicity_index() const {return multiplicity_index_;}
      // Return multiplicity index of this sector.

      inline bool IsDiagonal() const
      // Test if sector is diagonal (i.e., within a single subspace).
      {
        return (bra_subspace_index()==ket_subspace_index());
      }

      inline bool IsUpperTriangle() const
      // Test if sector is in upper triangle (includes diagonal).
      {
        return (bra_subspace_index()<=ket_subspace_index());
      }

      inline std::size_t num_elements() const
      {
        return (bra_subspace().dimension() * ket_subspace().dimension());
      }

      private:
      std::size_t bra_subspace_index_, ket_subspace_index_;
      std::size_t multiplicity_index_;
      std::shared_ptr<const BraSubspaceType> bra_subspace_ptr_;
      std::shared_ptr<const KetSubspaceType> ket_subspace_ptr_;
    };

  // Here we specialize to the case where the bra and ket subspaces have the
  // same type. We inherit all the functionality from the case where they
  // differ, and add an additional type alias.
  template<typename tSubspaceType>
    class BaseSector<tSubspaceType, tSubspaceType, true>
      : public BaseSector<tSubspaceType, tSubspaceType, false>
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SubspaceType = tSubspaceType;

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSector() = default;

      inline BaseSector(
          std::size_t bra_subspace_index, std::size_t ket_subspace_index,
          std::shared_ptr<const SubspaceType> bra_subspace_ptr,
          std::shared_ptr<const SubspaceType> ket_subspace_ptr,
          std::size_t multiplicity_index=1
        )
        : BaseSector<tSubspaceType, tSubspaceType, false>{
            bra_subspace_index, ket_subspace_index,
            std::move(bra_subspace_ptr), std::move(ket_subspace_ptr),
            multiplicity_index
        }
        {}
    };

  // BaseSectors -- container to hold a set of sectors with reverse
  // lookup by sector labels
  //
  // Template arguments:
  //   tBraSpaceType (typename) : type for bra space
  //   tKetSpaceType (typename) : type for ket space
  //   tSectorType (typename) : type for sector
  //   (bool): true if tBraSpaceType and tKetSpaceType are the same

  // declare BaseSectors template
  template <
      typename tBraSpaceType, typename tKetSpaceType,
      typename tSectorType,
      bool = std::is_same_v<tBraSpaceType, tKetSpaceType>
    >
    class BaseSectors;

  // Here we specialize the template for the case where the bra and ket spaces
  // are of different types. For the case where they are the same (below), we
  // will inherit from this specialization.
  template<typename tBraSpaceType, typename tKetSpaceType, typename tSectorType>
    class BaseSectors<tBraSpaceType, tKetSpaceType, tSectorType, false>
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using BraSpaceType = tBraSpaceType;
      using KetSpaceType = tKetSpaceType;

      using SectorType = tSectorType;
      using BraSubspaceType = typename SectorType::BraSubspaceType;
      static_assert(
          is_subspace_v<BraSubspaceType,BraSpaceType>,
          "BraSubspaceType must be subspace of BraSpaceType"
        );
      using KetSubspaceType = typename SectorType::KetSubspaceType;
      static_assert(
          is_subspace_v<KetSubspaceType,KetSpaceType>,
          "KetSubspaceType must be subspace of KetSpaceType"
        );


      protected:

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSectors()
        : num_elements_{}
      {}
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      template<
          typename T, typename U,
          typename std::enable_if_t<std::is_constructible_v<
              std::shared_ptr<const BraSpaceType>,
              T
            >>* = nullptr,
          typename std::enable_if_t<std::is_constructible_v<
              std::shared_ptr<const KetSpaceType>,
              U
            >>* = nullptr
        >
      BaseSectors(
            T&& bra_space_ptr,
            U&& ket_space_ptr
          )
          : bra_space_ptr_{std::forward<T>(bra_space_ptr)},
            ket_space_ptr_{std::forward<U>(ket_space_ptr)},
            num_elements_{}
      {}

      template<
          typename T, typename U,
          typename std::enable_if_t<std::is_same_v<std::decay_t<T>, BraSpaceType>>* = nullptr,
          typename std::enable_if_t<std::is_same_v<std::decay_t<U>, KetSpaceType>>* = nullptr
        >
      BaseSectors(T&& bra_space, U&& ket_space)
        // : bra_space_ptr_(bra_space.weak_from_this().lock()),
        //   ket_space_ptr_(ket_space.weak_from_this().lock())
        : num_elements_{}
      {
        // if bra_space or ket_space is not managed by a shared_ptr, make copy
        // or move it to a locally-created shared_ptr.
        //
        // Note: this is probably safe, since the spaces themselves store
        // vectors of shared_ptrs to subspaces, so a copy only involves copying
        // those vectors and maps
        // if (!bra_space_ptr_)
          bra_space_ptr_
            = std::make_shared<const BraSpaceType>(std::forward<T>(bra_space));
        // if (!ket_space_ptr_)
        {
          // if bra and ket spaces are identical, we can share the new copy we
          // just made for the bra_space_ptr_
          //
          // note: this also ensures that if bra_space and ket_space are the
          // same object, then if we just moved from bra_space we don't try to
          // move from it again (which would be UB)
          if (std::addressof(ket_space) == std::addressof(bra_space))
            ket_space_ptr_ = bra_space_ptr_;
          else
            ket_space_ptr_
              = std::make_shared<const KetSpaceType>(std::forward<U>(ket_space));
        }
      }

      public:

      ////////////////////////////////////////////////////////////////
      // space accessors
      ////////////////////////////////////////////////////////////////

      const BraSpaceType& bra_space() const
      {
        return *bra_space_ptr_;
      }

      const KetSpaceType& ket_space() const
      {
        return *ket_space_ptr_;
      }

      std::shared_ptr<const BraSpaceType> bra_space_ptr() const
      {
        return bra_space_ptr_;
      }

      std::shared_ptr<const KetSpaceType> ket_space_ptr() const
      {
        return ket_space_ptr_;
      }

      ////////////////////////////////////////////////////////////////
      // iterators
      ////////////////////////////////////////////////////////////////

      using value_type = SectorType;
      using const_reference = const SectorType&;
      using iterator = decltype(std::declval<const std::vector<SectorType>>().begin());
      using const_iterator = decltype(std::declval<const std::vector<SectorType>>().cbegin());
      using difference_type = typename iterator::difference_type;
      using size_type = std::make_unsigned_t<difference_type>;

      iterator begin()  const { return sectors_.begin(); }
      iterator end()    const { return sectors_.end(); }
      const_iterator cbegin() const { return sectors_.cbegin(); }
      const_iterator cend()   const { return sectors_.cend(); }

      ////////////////////////////////////////////////////////////////
      // sector lookup and retrieval
      ////////////////////////////////////////////////////////////////

      const SectorType& GetSector(std::size_t sector_index) const
      // Given sector index, return reference to sector itself.
      {
        return sectors_.at(sector_index);
      };

      bool ContainsSector(typename SectorType::KeyType key) const
      // Given the labels for a sector, returns whether or not the sector
      // is found within the the sector set.
      {
        return lookup_.count(key);
      };

      bool ContainsSector(
          std::size_t bra_subspace_index,
          std::size_t ket_subspace_index,
          std::size_t multiplicity_index=1
        ) const
      // Given the labels for a sector, returns whether or not the sector
      // is found within the the sector set.
      {
        return ContainsSector(
            typename SectorType::KeyType{bra_subspace_index,ket_subspace_index,multiplicity_index}
          );
      };

      std::size_t LookUpSectorIndex(const typename SectorType::KeyType& key) const
      // Given the key for a sector, look up its index within the
      // sector set.
      //
      // If no such labels are found, basis::kNone is returned.
      {
        auto pos = lookup_.find(key);
        if (pos==lookup_.end())
          return kNone;
        else
          return pos->second;
      };

      std::size_t LookUpSectorIndex(
          std::size_t bra_subspace_index,
          std::size_t ket_subspace_index,
          std::size_t multiplicity_index=1
        ) const
      // Given the labels for a sector, look up its index within the
      // sector set.
      //
      // If no such labels are found, basis::kNone is returned.
      {
        return LookUpSectorIndex(
            typename SectorType::KeyType{bra_subspace_index,ket_subspace_index,multiplicity_index}
          );
      }

      std::size_t GetSectorOffset(std::size_t i) const
      /// Given the index for a sector, return the offset of the sector in
      /// the full matrix indexing.
      {
        return sector_offsets_.at(i);
      }

      ////////////////////////////////////////////////////////////////
      // size retrieval
      ////////////////////////////////////////////////////////////////

      std::size_t size() const
      // Return number of sectors within sector set.
      {
        return sectors_.size();
      };

      std::size_t num_elements() const
      // Return number of matrix elements within sector set
      {
        return num_elements_;
      };

      ////////////////////////////////////////////////////////////////
      // diagnostic strings
      ////////////////////////////////////////////////////////////////

      std::string DebugStr() const;
      // Generate string dump of contents, for debugging purposes.
      //
      // Requires subspace to have a LabelStr() member function.

      protected:

      ////////////////////////////////////////////////////////////////
      // sector push (for initial construction)
      ////////////////////////////////////////////////////////////////

      template<
          typename T,
          typename std::enable_if_t<std::is_same_v<std::decay_t<T>,SectorType>>* = nullptr
        >
      void PushSector(T&& sector)
      // Create indexing information (in both directions, index <->
      // labels) for a sector.
      {
        sector_offsets_.push_back(num_elements());
        lookup_[sector.Key()] = size(); // index for lookup
        num_elements_ += sector.num_elements();
        sectors_.push_back(std::forward<T>(sector));  // save sector
      };

      template<typename... Args>
      void EmplaceSector(Args&&... args)
      // Create indexing information (in both directions, index <->
      // labels) for a sector.
      {
        const std::size_t index = sectors_.size();
        sector_offsets_.push_back(num_elements());
        sectors_.emplace_back(std::forward<Args>(args)...);  // save sector
        const SectorType& sector = sectors_.back();
        lookup_[sector.Key()] = index; // index for lookup
        num_elements_ += sector.num_elements();
      };

      template<
          typename T = SectorType,
          std::enable_if_t<std::is_same_v<typename T::BraSubspaceType,typename BraSpaceType::SubspaceType>>* = nullptr,
          std::enable_if_t<std::is_same_v<typename T::KetSubspaceType,typename KetSpaceType::SubspaceType>>* = nullptr
        >
      inline void PushSector(
          std::size_t bra_subspace_index,
          std::size_t ket_subspace_index,
          std::size_t multiplicity_index=1
        )
      // Create indexing information (in both directions, index <->
      // labels) for a sector, given indices.
      {
        EmplaceSector(
            bra_subspace_index, ket_subspace_index,
            bra_space_ptr_->GetSubspacePtr(bra_subspace_index),
            ket_space_ptr_->GetSubspacePtr(ket_subspace_index),
            multiplicity_index
          );  // save sector
      }

      template<
          typename T = SectorType,
          std::enable_if_t<std::is_same_v<typename T::BraSubspaceType,typename BraSpaceType::SubspaceType>>* = nullptr,
          std::enable_if_t<std::is_same_v<typename T::KetSubspaceType,typename KetSpaceType::SubspaceType>>* = nullptr
        >
      inline void PushSector(const typename SectorType::KeyType& key)
      // Create indexing information (in both directions, index <->
      // labels) for a sector, given key.
      {
        const auto& [bra_subspace_index,ket_subspace_index,multiplicity_index] = key;
        PushSector(bra_subspace_index, ket_subspace_index, multiplicity_index);
      };

      private:
      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      // spaces
      std::shared_ptr<const BraSpaceType> bra_space_ptr_;
      std::shared_ptr<const KetSpaceType> ket_space_ptr_;

      // sectors (accessible by index)
      std::vector<SectorType> sectors_;

      // sector index lookup by subspace indices
      basis::map<typename SectorType::KeyType,std::size_t> lookup_;

      // indexing
      std::vector<std::size_t> sector_offsets_;
      std::size_t num_elements_;

    };

  // Here we specialize to the case where the bra and ket spaces have the same
  // type. We inherit all the functionality from the case where they differ, and
  // add some additional type aliases and constructors.
  template<typename tSpaceType, typename tSectorType>
    class BaseSectors<tSpaceType, tSpaceType, tSectorType, true>
      : public BaseSectors<tSpaceType, tSpaceType, tSectorType, false>
    {
      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SpaceType = tSpaceType;
      // TODO (pjf): Fix up template type here so that SubspaceType only
      // gets defined if unambiguous
      // template<
      //     typename T = tSectorType,
      //     std::enable_if_t<std::is_same_v<typename T::BraSubspaceType, typename T::KetSubspaceType>>* = nullptr
      //   >
      static_assert(
          std::is_same_v<typename tSectorType::BraSubspaceType, typename tSectorType::KetSubspaceType>,
          "same-Space BaseSectors currently allowed with same-Subspace BaseSector"
        );
      using SubspaceType = typename tSectorType::BraSubspaceType;

      protected:

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseSectors() = default;

      template<
          typename T,
          typename std::enable_if_t<mcutils::is_derived_constructible_v<
              BaseSectors<tSpaceType, tSpaceType, tSectorType, false>,
              T, T
            >>* = nullptr
        >
      inline explicit BaseSectors(T&& space_ptr)
        : BaseSectors<tSpaceType, tSpaceType, tSectorType, false>{
              space_ptr, std::forward<T>(space_ptr)
          }
      {}
      // construct from shared_ptr to space

      template<
          typename... Args,
          typename std::enable_if_t<mcutils::is_derived_constructible_v<
              BaseSectors<tSpaceType, tSpaceType, tSectorType, false>,
              Args...
            >>* = nullptr
        >
      inline explicit BaseSectors(Args&&... args)
        : BaseSectors<tSpaceType, tSpaceType, tSectorType, false>{
              std::forward<Args>(args)...
            }
      {}

    };

  template <typename tBraSpaceType, typename tKetSpaceType, typename tSectorType>
    std::string BaseSectors<tBraSpaceType, tKetSpaceType, tSectorType, false>::DebugStr() const
    {
      std::ostringstream os;
      for (std::size_t sector_index=0; sector_index<size(); ++sector_index)
        {
          const SectorType& sector = GetSector(sector_index);

          os << "  sector " << sector_index
             << "  bra index " << sector.bra_subspace_index()
             << " labels " << sector.bra_subspace().LabelStr()
             << " size " << sector.bra_subspace().size()
             << " dim " << sector.bra_subspace().dimension()
             << "  ket index " << sector.ket_subspace_index()
             << " labels " << sector.ket_subspace().LabelStr()
             << " size " << sector.ket_subspace().size()
             << " dim " << sector.ket_subspace().dimension()
             << "  multiplicity index " << sector.multiplicity_index()
             << "  elements " << sector.num_elements()
             << std::endl;
        }
      return os.str();
    }

}  // namespace impl

template <typename tBraSubspaceType, typename tKetSubspaceType = tBraSubspaceType>
using BaseSector = impl::BaseSector<
    tBraSubspaceType, tKetSubspaceType,
    std::is_same_v<tBraSubspaceType, tKetSubspaceType>
  >;

template <
    typename tBraSpaceType, typename tKetSpaceType = tBraSpaceType,
    typename tSectorType
      = BaseSector<typename tBraSpaceType::SubspaceType, typename tKetSpaceType::SubspaceType>
  >
  using BaseSectors = impl::BaseSectors<
      tBraSpaceType, tKetSpaceType,
      tSectorType,
      std::is_same_v<tBraSpaceType, tKetSpaceType>
    >;

}  // namespace basis

#endif  // BASIS_SECTORS_H_
