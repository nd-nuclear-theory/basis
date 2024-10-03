/************************************************************//**
  @file hypersector.h

  Enumeration of sectors for multiple operators.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 04/04/17 (mac): Created, building on code from basis.h.
  + 04/11/17 (mac): Add basic hyperoperator support --
    OperatorHypersectors and SetHyperoperatorToZero.
  + 04/22/17 (aem): Fix error in hypersector constructor.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
  + 08/04/21 (pjf): Fix use of size() vs. dimension().
  + 09/27/24: Move implementations into impl namespace to hide template
    parameters which should be invisible to library users.
****************************************************************/

#ifndef BASIS_HYPERSECTOR_H_
#define BASIS_HYPERSECTOR_H_

#include <cstddef>
#include <iosfwd>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "basis.h"
#include "operator.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // operator storage by hyperblock
  ////////////////////////////////////////////////////////////////


  // typedef for operator hyperblocks
  //
  // EX:
  //   basis::OperatorHyperblocks<double> matrices;  // matrices for all operator planes of all hypersectors
  //   ...
  //   matrix_element = matrices[hypersector_index][operator_index](bra_index,ket_index);  // operator_index is within given hypersector's operator subspace

  template <typename tFloat>
    using OperatorHyperblocks = std::vector<std::vector<OperatorBlock<tFloat>>>;

namespace impl {

  ////////////////////////////////////////////////////////////////
  // hypersector indexing
  ////////////////////////////////////////////////////////////////

  // declare BaseSectors template
  template <
      typename tBraSubspaceType,
      typename tOperatorSubspaceType,
      typename tKetSubspaceType = tBraSubspaceType,
      bool = std::is_same_v<tBraSubspaceType, tKetSubspaceType>
    >
    class BaseHypersector;


  // Note: A "hypersector" represents a triplet consisting of a pair
  // of "state" subspaces (which thus define a "sector" or block in
  // the matrix representation of an operator on the space) and an
  // "operator" subspace.  The hypersector may also be labeled with an
  // optional multiplicity index.


  /// BaseHypersector -- provide storage of information for a single hypersector

  // Here we specialize the template for the case where the bra and ket
  // subspaces are of different types. For the case where they are the same
  // (below), we will inherit from this specialization.
  template <typename tBraSubspaceType, typename tOperatorSubspaceType, typename tKetSubspaceType>
    class BaseHypersector<tBraSubspaceType, tOperatorSubspaceType, tKetSubspaceType, false>
    : public BaseSector<tBraSubspaceType, tKetSubspaceType>
    // Store indexing and subspace reference information for a
    // hypersector.
    {
      private:

      ////////////////////////////////////////////////////////////////
      // private (convenience) typedefs
      ////////////////////////////////////////////////////////////////
      using BaseSectorType = BaseSector<tBraSubspaceType, tKetSubspaceType>;

      ////////////////////////////////////////////////////////////////
      // typedefs
      ////////////////////////////////////////////////////////////////

      public:
      using BraSubspaceType = typename BaseSectorType::BraSubspaceType;
      using KetSubspaceType = typename BaseSectorType::KetSubspaceType;
      using OperatorSubspaceType = tOperatorSubspaceType;
      using KeyType = std::tuple<std::size_t,std::size_t,std::size_t,std::size_t>;

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////
      BaseHypersector() = default;

      inline BaseHypersector(
          std::size_t bra_subspace_index,
          std::size_t ket_subspace_index,
          std::size_t operator_subspace_index,
          std::shared_ptr<const BraSubspaceType> bra_subspace_ptr,
          std::shared_ptr<const KetSubspaceType> ket_subspace_ptr,
          std::shared_ptr<const OperatorSubspaceType> operator_subspace_ptr,
          std::size_t multiplicity_index=1
        )
        : BaseSectorType{
              bra_subspace_index, ket_subspace_index,
              std::move(bra_subspace_ptr), std::move(ket_subspace_ptr),
              multiplicity_index
            },
          operator_subspace_index_{operator_subspace_index},
          operator_subspace_ptr_{std::move(operator_subspace_ptr)}
      {}

      ////////////////////////////////////////////////////////////////
      // accessors
      ////////////////////////////////////////////////////////////////

      inline KeyType Key() const
      // Return tuple key identifying sector for sorting/lookup
      // purposes.
      {
        return KeyType(
            BaseSectorType::bra_subspace_index(), BaseSectorType::ket_subspace_index(),
            operator_subspace_index(), BaseSectorType::multiplicity_index()
          );
      }

      std::size_t operator_subspace_index() const
      {
        return operator_subspace_index_;
      }
      // Return integer index of operator subspace.
      const OperatorSubspaceType& operator_subspace() const
      {
        return *operator_subspace_ptr_;
      }
      // Return reference to operator subspace.
      std::shared_ptr<const OperatorSubspaceType> operator_subspace_ptr() const
      {
        return operator_subspace_ptr_;
      }
      // Return shared pointer to bra/ket subspace.

      private:
      std::size_t operator_subspace_index_;
      std::shared_ptr<const OperatorSubspaceType> operator_subspace_ptr_;
    };



  // Here we specialize to the case where the bra and ket subspaces have the
  // same type. We inherit all the functionality from the case where they
  // differ, and add an additional type alias.
  template<typename tSubspaceType,typename tOperatorSubspaceType>
    class BaseHypersector<tSubspaceType, tOperatorSubspaceType, tSubspaceType, true>
      : public BaseHypersector<tSubspaceType, tOperatorSubspaceType, tSubspaceType, false>
    {
      private:

      ////////////////////////////////////////////////////////////////
      // private (convenience) typedefs
      ////////////////////////////////////////////////////////////////
      using BaseHypersectorType
        = BaseHypersector<tSubspaceType, tOperatorSubspaceType, tSubspaceType, false>;

      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SubspaceType = tSubspaceType;
      using OperatorSubspaceType = tOperatorSubspaceType;

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseHypersector() = default;

      inline BaseHypersector(
          std::size_t bra_subspace_index,
          std::size_t ket_subspace_index,
          std::size_t operator_subspace_index,
          std::shared_ptr<const SubspaceType> bra_subspace_ptr,
          std::shared_ptr<const SubspaceType> ket_subspace_ptr,
          std::shared_ptr<const OperatorSubspaceType> operator_subspace_ptr,
          std::size_t multiplicity_index=1
        )
        : BaseHypersectorType{
              bra_subspace_index, ket_subspace_index,
              operator_subspace_index,
              std::move(bra_subspace_ptr), std::move(ket_subspace_ptr),
              std::move(operator_subspace_ptr),
              multiplicity_index
            }
      {}

    };


  // BaseHypersectors -- container to hold a set of hypersectors with
  // reverse lookup by hypersector labels
  //
  // Template arguments:
  //   tSpaceType (typename): type for bra/ket space
  //   tOperatorSpaceType (typename): type for operator space

  // declare BaseSectors template
  template <
      typename tBraSpaceType,
      typename tOperatorSpaceType,
      typename tKetSpaceType = tBraSpaceType,
      typename tHypersectorType
        = BaseHypersector<
          typename tBraSpaceType::SubspaceType,
          typename tOperatorSpaceType::SubspaceType,
          typename tKetSpaceType::SubspaceType
        >,
      bool = std::is_same_v<tBraSpaceType, tKetSpaceType>
    >
    class BaseHypersectors;

  template <typename tBraSpaceType, typename tOperatorSpaceType, typename tKetSpaceType, typename tHypersectorType>
    class BaseHypersectors<tBraSpaceType, tOperatorSpaceType, tKetSpaceType, tHypersectorType, false>
    : public BaseSectors<tBraSpaceType, tKetSpaceType, tHypersectorType>
    {
      private:

    ////////////////////////////////////////////////////////////////
    // private (convenience) typedefs
    ////////////////////////////////////////////////////////////////
      using BaseSectorsType
        = BaseSectors<tBraSpaceType, tKetSpaceType, tHypersectorType>;

      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////


      using BraSpaceType = typename BaseSectorsType::BraSpaceType;
      using KetSpaceType = typename BaseSectorsType::KetSpaceType;
      using OperatorSpaceType = tOperatorSpaceType;
      using OperatorSubspaceType = typename tOperatorSpaceType::SubspaceType;
      using HypersectorType = tHypersectorType;


      BaseHypersectors() = default;


      template<
          typename T, typename U, typename V,
          typename std::enable_if_t<mcutils::is_derived_constructible_v<
              BaseSectorsType, T, U
            >>* = nullptr,
          typename std::enable_if_t<std::is_constructible_v<
              std::shared_ptr<const OperatorSpaceType>, V
            >>* = nullptr
        >
      BaseHypersectors(
          T&& bra_space_or_ptr,
          U&& ket_space_or_ptr,
          V&& operator_space_ptr
        )
        : BaseSectorsType{
              std::forward<T>(bra_space_or_ptr), std::forward<U>(ket_space_or_ptr)
            },
          operator_space_ptr_{std::forward<V>(operator_space_ptr)}
      {}

      template<
          typename T, typename U, typename V,
          typename std::enable_if_t<mcutils::is_derived_constructible_v<
              BaseSectorsType,
              T, U
            >>* = nullptr,
          typename std::enable_if_t<std::is_same_v<std::decay_t<V>, OperatorSpaceType>>* = nullptr
        >
      BaseHypersectors(
          T&& bra_space_or_ptr,
          U&& ket_space_or_ptr,
          V&& operator_space
        )
        : BaseSectorsType{
              std::forward<T>(bra_space_or_ptr), std::forward<U>(ket_space_or_ptr)
            }
      {
          operator_space_ptr_
            = std::make_shared<const OperatorSpaceType>(std::forward<V>(operator_space));
      }

      ////////////////////////////////////////////////////////////////
      // space accessors
      ////////////////////////////////////////////////////////////////

      const OperatorSpaceType& operator_space() const
      {
        return *operator_space_ptr_;
      }

      std::shared_ptr<const BraSpaceType> operator_space_ptr() const
      {
        return operator_space_ptr_;
      }

      ////////////////////////////////////////////////////////////////
      // sector lookup and retrieval
      ////////////////////////////////////////////////////////////////

      const HypersectorType& GetHypersector(std::size_t hypersector_index) const
      // Given sector index, return reference to sector itself.
      {
        return BaseSectorsType::GetSector(hypersector_index);
      };

      bool ContainsHypersector(typename HypersectorType::KeyType key) const
      // Given the labels for a sector, returns whether or not the sector
      // is found within the the sector set.
      {
        return BaseSectorsType::ContainsSector(key);
      };

      bool ContainsHypersector(
          std::size_t bra_subspace_index, std::size_t ket_subspace_index, std::size_t operator_subspace_index,
          std::size_t multiplicity_index=1
        ) const
      // Given the labels for a sector, returns whether or not the sector
      // is found within the the sector set.
      {
        return ContainsHypersector(
            typename HypersectorType::KeyType{
                bra_subspace_index,ket_subspace_index,operator_subspace_index,multiplicity_index
              }
          );
      };

      std::size_t LookUpHypersectorIndex(const typename HypersectorType::KeyType& key) const
      // Given the key for a sector, look up its index within the
      // sector set.
      //
      // If no such labels are found, basis::kNone is returned.
      {
        return BaseSectorsType::LookUpSectorIndex(key);
      }
      std::size_t LookUpHypersectorIndex(
          std::size_t bra_subspace_index, std::size_t ket_subspace_index, std::size_t operator_subspace_index,
          std::size_t multiplicity_index=1
        ) const
      // Given the labels for a sector, look up its index within the
      // sector set.
      //
      // If no such labels are found, basis::kNone is returned.
      {
        return LookUpHypersectorIndex(
            typename HypersectorType::KeyType{
                bra_subspace_index,ket_subspace_index,operator_subspace_index,multiplicity_index
              }
          );
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
      // hypersector push (for initial construction)
      ////////////////////////////////////////////////////////////////

      template<
          typename T,
          typename std::enable_if_t<std::is_same_v<std::decay_t<T>,HypersectorType>>* = nullptr
        >
      void PushHypersector(T&& hypersector)
      // Create indexing information (in both directions, index <->
      // labels) for a sector.
      {
        BaseSectorsType::PushSector(std::forward<T>(hypersector));
      };

      template<typename... Args>
      void EmplaceHypersector(Args&&... args)
      // Create indexing information (in both directions, index <->
      // labels) for a sector.
      {
        BaseSectorsType::EmplaceSector(std::forward<Args>(args)...);
      }

      template<
          typename T = HypersectorType,
          std::enable_if_t<std::is_same_v<typename T::BraSubspaceType,typename BraSpaceType::SubspaceType>>* = nullptr,
          std::enable_if_t<std::is_same_v<typename T::KetSubspaceType,typename KetSpaceType::SubspaceType>>* = nullptr
        >
      inline void PushHypersector(
          std::size_t bra_subspace_index,
          std::size_t ket_subspace_index,
          std::size_t operator_subspace_index,
          std::size_t multiplicity_index=1
        )
      // Create indexing information (in both directions, index <->
      // labels) for a sector, given indices.
      {
        EmplaceHypersector(
            bra_subspace_index,
            ket_subspace_index,
            operator_subspace_index,
            BaseSectorsType::bra_space().GetSubspacePtr(bra_subspace_index),
            BaseSectorsType::ket_space().GetSubspacePtr(ket_subspace_index),
            operator_space().GetSubspacePtr(operator_subspace_index),
            multiplicity_index
          );  // save sector
      }

      template<
          typename T = HypersectorType,
          std::enable_if_t<std::is_same_v<typename T::BraSubspaceType, typename BraSpaceType::SubspaceType>>* = nullptr,
          std::enable_if_t<std::is_same_v<typename T::KetSubspaceType, typename KetSpaceType::SubspaceType>>* = nullptr
        >
      inline void PushHypersector(const typename HypersectorType::KeyType& key)
      // Create indexing information (in both directions, index <->
      // labels) for a sector, given key.
      {
        const auto& [bra_subspace_index,ket_subspace_index,operator_subspace_index,multiplicity_index] = key;
        PushHypersector(bra_subspace_index, ket_subspace_index, operator_subspace_index, multiplicity_index);
      };


      private:
      ////////////////////////////////////////////////////////////////
      // internal storage
      ////////////////////////////////////////////////////////////////

      // spaces
      std::shared_ptr<const OperatorSpaceType> operator_space_ptr_;

    };

  // Here we specialize to the case where the bra and ket spaces have the same
  // type. We inherit all the functionality from the case where they differ, and
  // add some additional type aliases and constructors.
  template<typename tSpaceType, typename tOperatorSpaceType, typename tHypersectorType>
    class BaseHypersectors<tSpaceType, tOperatorSpaceType, tSpaceType,  tHypersectorType, true>
      : public BaseHypersectors<tSpaceType, tOperatorSpaceType, tSpaceType, tHypersectorType, false>
    {
      private:

      ////////////////////////////////////////////////////////////////
      // private (convenience) typedefs
      ////////////////////////////////////////////////////////////////
      using BaseHypersectorsType
        = BaseHypersectors<tSpaceType, tOperatorSpaceType, tSpaceType, tHypersectorType, false>;

      public:

      ////////////////////////////////////////////////////////////////
      // common typedefs
      ////////////////////////////////////////////////////////////////

      using SpaceType = tSpaceType;
      using OperatorSpaceType = tOperatorSpaceType;
      // TODO (pjf): Fix up template type here so that SubspaceType only
      // gets defined if unambiguous
      // template<
      //     typename T = tSectorType,
      //     std::enable_if_t<std::is_same_v<typename T::BraSubspaceType, typename T::KetSubspaceType>>* = nullptr
      //   >
      static_assert(
          std::is_same_v<typename tHypersectorType::BraSubspaceType, typename tHypersectorType::KetSubspaceType>,
          "same-Space BaseHypersectors currently allowed with same-Subspace BaseHypersector"
        );
      using SubspaceType = typename tHypersectorType::BraSubspaceType;
      using OperatorSubspaceType = typename tHypersectorType::OperatorSubspaceType;

      protected:

      ////////////////////////////////////////////////////////////////
      // constructors
      ////////////////////////////////////////////////////////////////

      BaseHypersectors() = default;



      template<
          typename T, typename U,
          typename std::enable_if_t<mcutils::is_derived_constructible_v<
            BaseHypersectorsType, T, T, U
          >>* = nullptr
        >
      inline BaseHypersectors(T&& space_or_ptr, U&& operator_space_or_ptr)
      : BaseHypersectorsType{
            space_or_ptr,
            std::forward<T>(space_or_ptr),
            std::forward<U>(operator_space_or_ptr)
          }
      {}

      template<
          typename T, typename U, typename V,
          typename std::enable_if_t<mcutils::is_derived_constructible_v<
            BaseHypersectorsType, T, U, V
          >>* = nullptr
        >
      inline BaseHypersectors(
          T&& bra_space_or_ptr,
          U&& ket_space_or_ptr,
          V&& operator_space_or_ptr
        )
        : BaseHypersectorsType{
              std::forward<T>(bra_space_or_ptr),
              std::forward<U>(ket_space_or_ptr),
              std::forward<V>(operator_space_or_ptr)
            }
      {}

    };

  template <typename tBraSpaceType, typename tOperatorSpaceType, typename tKetSpaceType, typename tHypersectorType>
    std::string BaseHypersectors<tBraSpaceType, tOperatorSpaceType, tKetSpaceType, tHypersectorType, false>::DebugStr() const
    {
      std::ostringstream os;
      for (std::size_t hypersector_index=0; hypersector_index<BaseSectorsType::size(); ++hypersector_index)
        {
          const HypersectorType& hypersector = GetHypersector(hypersector_index);

          os << "  hypersector " << hypersector_index
             << "  bra index " << hypersector.bra_subspace_index()
             << " labels " << hypersector.bra_subspace().LabelStr()
             << " size " << hypersector.bra_subspace().size()
             << " dim " << hypersector.bra_subspace().dimension()
             << "  ket index " << hypersector.ket_subspace_index()
             << " labels " << hypersector.ket_subspace().LabelStr()
             << " size " << hypersector.ket_subspace().size()
             << " dim " << hypersector.ket_subspace().dimension()
             << "  operator index " << hypersector.operator_subspace_index()
             << " labels " << hypersector.operator_subspace().LabelStr()
             << " size " << hypersector.operator_subspace().size()
             << " dim " << hypersector.operator_subspace().dimension()
             << "  multiplicity index " << hypersector.multiplicity_index()
             << std::endl;
        }
      return os.str();
    }

}  // namespace impl

  ////////////////////////////////////////////////////////////////
  // template alias definitions (in basis:: namespace)
  ////////////////////////////////////////////////////////////////

  template <
      typename tBraSubspaceType,
      typename tOperatorSubspaceType,
      typename tKetSubspaceType = tBraSubspaceType
    >
    using BaseHypersector = impl::BaseHypersector<
        tBraSubspaceType, tOperatorSubspaceType, tKetSubspaceType,
        std::is_same_v<tBraSubspaceType, tKetSubspaceType>
      >;

  template <
      typename tBraSpaceType,
      typename tOperatorSpaceType,
      typename tKetSpaceType = tBraSpaceType,
      typename tHypersectorType
        = BaseHypersector<
          typename tBraSpaceType::SubspaceType,
          typename tOperatorSpaceType::SubspaceType,
          typename tKetSpaceType::SubspaceType
        >
    >
    using BaseHypersectors = impl::BaseHypersectors<
        tBraSpaceType, tOperatorSpaceType, tKetSpaceType, tHypersectorType,
        std::is_same_v<tBraSpaceType, tKetSpaceType>
      >;


  ////////////////////////////////////////////////////////////////
  // zero operator
  ////////////////////////////////////////////////////////////////

  template <typename tHypersectorsType, typename tFloat>
  void SetHyperoperatorToZero(
      const tHypersectorsType& hypersectors,
      basis::OperatorHyperblocks<tFloat>& matrices
    )
    // Set operator hyperblocks to zero.
    //
    // Arguments:
    //   hypersectors (input): the set of hypersectors on which the
    //     operators are defined
    //   matrices (output): matrices to hold hyperblocks
  {

    // clear vector of matrices
    matrices.clear();
    matrices.resize(hypersectors.size());

    // iterate over hypersectors
    for (std::size_t hypersector_index = 0; hypersector_index < hypersectors.size(); ++hypersector_index)
      {

        // extract hypersector
        const typename tHypersectorsType::HypersectorType& hypersector = hypersectors.GetHypersector(hypersector_index);

        // extract hypersector subspaces
        const typename tHypersectorsType::SubspaceType& bra_subspace = hypersector.bra_subspace();
        const typename tHypersectorsType::SubspaceType& ket_subspace = hypersector.ket_subspace();
        const typename tHypersectorsType::OperatorSubspaceType& operator_subspace = hypersector.operator_subspace();

        // generate matrices for hypersector (by operator)
        matrices[hypersector_index].resize(operator_subspace.size());
        for (std::size_t operator_index = 0; operator_index < operator_subspace.size(); ++operator_index)
          matrices[hypersector_index][operator_index]
            = basis::OperatorBlock<tFloat>::Zero(bra_subspace.dimension(),ket_subspace.dimension());
      }
  }

  ////////////////////////////////////////////////////////////////
  // counting operator
  ////////////////////////////////////////////////////////////////

  template <typename tHypersectorsType>
  std::size_t GetNumHyperoperatorME(
      const tHypersectorsType& hypersectors
    )
    // Set operator hyperblocks to zero.
    //
    // Arguments:
    //   hypersectors (input): the set of hypersectors on which the
    //     operators are defined
    //   matrices (output): matrices to hold hyperblocks
  {

    // iterate over hypersectors
    std::size_t num_rmes=0;
    for (std::size_t hypersector_index = 0; hypersector_index < hypersectors.size(); ++hypersector_index)
      {

        // extract hypersector
        const typename tHypersectorsType::HypersectorType& hypersector = hypersectors.GetHypersector(hypersector_index);

        // extract hypersector subspaces
        const typename tHypersectorsType::SubspaceType& bra_subspace = hypersector.bra_subspace();
        const typename tHypersectorsType::SubspaceType& ket_subspace = hypersector.ket_subspace();
        const typename tHypersectorsType::OperatorSubspaceType& operator_subspace = hypersector.operator_subspace();

        // generate matrices for hypersector (by operator)
        for (std::size_t operator_index = 0; operator_index < operator_subspace.size(); ++operator_index)
          num_rmes+=bra_subspace.dimension()*ket_subspace.dimension();
      }
    return num_rmes;
  }

  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
}  // namespace basis

#endif
