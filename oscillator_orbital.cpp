/****************************************************************
  oscillator_orbital.cpp

  Mark A. Caprio, Patrick J. Fasano, Jakub Herko, Shwetha L. Vittal
  University of Notre Dame

****************************************************************/


// #include <cstddef>
// #include <algorithm>
#include <iomanip>
// #include <set>
// #include <istream>
#include <sstream>

// #include "am/am.h"

#include "oscillator_orbital.h"
#include "basis.h"

namespace basis {


  ////////////////////////////////////////////////////////////////
  // single-particle orbitals
  ////////////////////////////////////////////////////////////////

  /**
   * Construct an Nmax-truncated subpace of harmonic oscillator spatial
   * orbitals.
   *
   * @param[in] l orbital angular momentum for subspace
   * @param[in] Nmax maximum number of oscillator quanta in truncation
   */
  OscillatorOrbitalSubspace::OscillatorOrbitalSubspace(int l, int Nmax)
    : BaseSubspace{{l}}, Nmax_{Nmax}
  {
    
    // validate subspace labels
    assert(ValidLabels());

    // iterate over total oscillator quanta
    for (int n = 0; n <= (Nmax-l)/2; ++n)
      // DEBUGGING: Must use arguments Nmax and l (or data members Nmax_ and l_)
      // here, instead of (hidden) accessors Nmax() and l().
      PushStateLabels(StateLabelsType(n));
  }

  /**
   * Verify validity of proposed subspace labels for subpace of harmonic
   * oscillator spatial orbitals.
   */
  bool OscillatorOrbitalSubspace::ValidLabels() const
  {

    bool valid = true;

    // nonnegativity
    valid &= l()>=0;

    // nonemptiness in given Nmax truncation
    valid &= l()<=Nmax();
    
    return valid;
  }

  /**
   * Generate a string representation of the subspace labels.
   * @return subspace labels as a string
   */
  // TODO (mac): reimplement in terms of fmt library
  std::string OscillatorOrbitalSubspace::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << int(l())
       << " " << "]";

    return os.str();
  }

  /**
   * Generate a string representation, useful for debugging.
   * @return debug string
   */
  // TODO (mac): reimplement in terms of fmt library
  std::string OscillatorOrbitalSubspace::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        StateType state(*this,state_index);

        os
          << " " << "index"
          << " " << std::setw(width) << state_index
          << " " << "n"
          << " " << std::setw(width) << state.n()
          << " " << "=> " << "N"
          << " " << std::setw(width) << state.N()
          << std::endl;
      }

    return os.str();

  }

  /**
   * Generate a string representation of the orbital labels.
   * @return orbital labels as a string
   */
  // TODO (mac): reimplement in terms of fmt library
  std::string OscillatorOrbitalState::LabelStr() const
  {
    std::ostringstream os;
  
    const int width = 0;  // for now, no fixed width
  
    os << "["
       << " " << "l"
       << " " << std::setw(width) << l()
       << " " << "index"
       << " " << std::setw(width) << index()
       << " :"
       << " " << "n"
       << " " << std::setw(width) << n()
       << " " << "]";
  
    return os.str();
  }
  
  // /**
  //  * Construct an Nmax-truncated single-particle space with species subspaces.
  //  *
  //  * @param[in] Nmax number of oscillator quanta
  //  */
  // OrbitalSpacePN::OrbitalSpacePN(int Nmax)
  // {
  // 
  //   // save truncation
  //   weight_max_ = double(Nmax);
  //   is_oscillator_like_ = true;
  //   Nmax_ = Nmax;
  // 
  //   // iterate over species
  //   for (OrbitalSpeciesPN orbital_species : {OrbitalSpeciesPN::kP,OrbitalSpeciesPN::kN})
  //     {
  //       OrbitalSubspacePN subspace(orbital_species,Nmax);
  //       PushSubspace(subspace);
  //     }
  // }
  // 
  // /**
  //  * Construct a space with species subspaces from a list of orbitals.
  //  *
  //  * @param[in] states vector of orbitals
  //  */
  // OrbitalSpacePN::OrbitalSpacePN(const OrbitalPNList& states)
  // {
  //   weight_max_ = 0.0;
  //   // collect orbital_species subspace labels sorted in canonical order
  //   std::set<OrbitalSubspacePNLabels> subspace_labels_set;
  //   for (std::size_t state_index=0; state_index<states.size(); ++state_index)
  //     {
  //       OrbitalPNInfo state = states[state_index];
  //       OrbitalSubspacePNLabels labels(state.orbital_species);
  //       subspace_labels_set.insert(labels);
  //     }
  // 
  //   // construct subspaces
  //   for (const OrbitalSubspacePNLabels& labels : subspace_labels_set)
  //     {
  //       OrbitalSpeciesPN orbital_species;
  //       std::tie(orbital_species) = labels;
  //       OrbitalSubspacePN subspace(orbital_species,states);
  //       PushSubspace(subspace);
  //       weight_max_ = std::max(weight_max_,subspace.weight_max());
  //     }
  // 
  //   // check if oscillator-like and set Nmax
  //   if (IsOscillatorLike_()) {
  //     is_oscillator_like_ = true;
  //     Nmax_ = int(weight_max_);
  //   } else {
  //     is_oscillator_like_ = false;
  //     Nmax_ = -1;
  //   }
  // 
  // }
  // 
  // /**
  //  * Check if space is truncated like an Nmax oscillator truncation.
  //  *
  //  * @return true if all subspaces are Nmax truncated with the same Nmax.
  //  */
  // bool OrbitalSpacePN::IsOscillatorLike_() const
  // {
  // 
  //   // first check that the space contains subspaces
  //   if (size() < 1)
  //     return false;
  // 
  //   // check if subspaces are individually Nmax truncated
  //   for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
  //     {
  //       const SubspaceType& subspace = GetSubspace(subspace_index);
  //       if (!subspace.is_oscillator_like())
  //         return false;
  //     }
  // 
  //   // see if can extract viable Nmax
  //   int Nmax = GetSubspace(0).Nmax();
  //   for (std::size_t subspace_index=1; subspace_index<size(); ++subspace_index)
  //     {
  //       const SubspaceType& subspace = GetSubspace(subspace_index);
  //       if (subspace.Nmax() != Nmax)
  //         return false;
  //     }
  // 
  //   // orbitals are oscillator like!
  //   return true;
  // }
  // 
  // /**
  //  * Generate a string representation, useful for debugging.
  //  * @return debug string
  //  */
  // std::string OrbitalSpacePN::DebugStr() const
  // {
  // 
  //   std::ostringstream os;
  // 
  //   const int width = 3;
  // 
  //   std::string oscillator_like_indicator = (is_oscillator_like() ? "true" : "false");
  //   os << " weight_max " << weight_max();
  //   if (is_oscillator_like())
  //     os << " Nmax " << Nmax();
  //   os << " (oscillator-like: " << oscillator_like_indicator << ")"
  //      << std::endl;
  // 
  //   for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
  //     {
  //       const SubspaceType& subspace = GetSubspace(subspace_index);
  // 
  //       os
  //         << " " << "index"
  //         << " " << std::setw(width) << subspace_index
  //         << " " << "species"
  //         << " " << std::setw(width) << int(subspace.orbital_species())
  //         << " " << "dim"
  //         << " " << std::setw(width) << subspace.size()
  //         << " " << std::endl;
  //     }
  // 
  //   return os.str();
  // 
  // }
  // 
  // /**
  //  * Flatten space into a vector of OrbitalPNInfo objects.
  //  *
  //  * @return vector representation of space
  //  */
  // OrbitalPNList OrbitalSpacePN::OrbitalInfo() const
  // {
  //   OrbitalPNList orbitals;
  // 
  //   for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
  //     {
  //       const SubspaceType& subspace = GetSubspace(subspace_index);
  //       OrbitalPNList subspace_orbitals;
  // 
  //       // get orbitals for subspace and append to vector
  //       subspace_orbitals = subspace.OrbitalInfo();
  //       orbitals.insert(orbitals.end(),
  //                       subspace_orbitals.begin(),
  //                       subspace_orbitals.end()
  //                      );
  // 
  //     }
  // 
  //   return orbitals;
  // 
  // }
  // 
  // 
  // ////////////////////////////////////////////////////////////////
  // // single-particle orbitals - lj subspaces
  // ////////////////////////////////////////////////////////////////
  // 
  // /**
  //  * Construct an Nmax-truncated single-particle subspace with a particular
  //  * species.
  //  *
  //  * @param[in] orbital_species species type for subspace
  //  * @param[in] l orbital angular momentum quantum number
  //  * @param[in] j total angular momentum quantum number
  //  * @param[in] Nmax number of oscillator quanta
  //  */
  // OrbitalSubspaceLJPN::OrbitalSubspaceLJPN(
  //     OrbitalSpeciesPN orbital_species, int l, HalfInt j, int Nmax
  //   )
  //   : BaseSubspace{{orbital_species,l,j}}, weight_max_{double(Nmax)},
  //     is_oscillator_like_{true}, Nmax_{Nmax}
  // {
  //   // iterate over radial quantum number
  //   for (int n = 0; (2*n+l) <= Nmax; ++n) {
  //     // save state
  //     PushStateLabels(StateLabelsType(n));
  // 
  //     // save oscillator quantum number as weight
  //     weights_.push_back(double(2*n+l));
  //   }
  // }
  // 
  // /**
  //  * Construct a subspace with a particular l, j, and species from a list
  //  * of orbitals.
  //  *
  //  * @param[in] orbital_species species type for subspace
  //  * @param[in] l orbital angular momentum quantum number
  //  * @param[in] j total angular momentum quantum number
  //  * @param[in] states vector of orbitals
  //  */
  // OrbitalSubspaceLJPN::OrbitalSubspaceLJPN(
  //     OrbitalSpeciesPN orbital_species, int l, HalfInt j,
  //     const OrbitalPNList& states
  //   )
  //   : BaseSubspace{{orbital_species,l,j}}, weight_max_{0.0}
  // {
  //   for (auto&& state : states) {
  //     if (state.orbital_species == orbital_species
  //         && state.l == l && state.j == j) {
  //       PushStateLabels(StateLabelsType(state.n));
  //       weights_.push_back(state.weight);
  //       weight_max_ = std::max(weight_max_,state.weight);
  //     }
  //   }
  // 
  //   // check if oscillator-like and set Nmax
  //   if (IsOscillatorLike_()) {
  //     is_oscillator_like_ = true;
  //     Nmax_ = static_cast<int>(weight_max_);
  //   } else {
  //     is_oscillator_like_ = false;
  //     Nmax_ = -1;
  //   }
  // }
  // 
  // /**
  //  * Do a deep comparison to oscillator-truncated basis.
  //  * @return true if truncation is like Nmax-truncated oscillator
  //  */
  // bool OrbitalSubspaceLJPN::IsOscillatorLike_() const
  // {
  //   // oscillators have states
  //   if (size()==0)
  //     return false;
  // 
  //   // maximum weight should be an integer
  //   int Nmax = int(weight_max());
  //   if (Nmax != weight_max())
  //     return false;
  // 
  //   // only positive Nmax is allowed
  //   if (Nmax < 0)
  //     return false;
  // 
  //   // compare to equivalent Nmax-truncated subspace
  //   OrbitalSubspacePN reference_subspace(orbital_species(),Nmax);
  //   if (reference_subspace.OrbitalInfo() != OrbitalInfo())
  //     return false;
  // 
  //   // if it looks like an oscillator, swims like an oscillator, and quacks
  //   // like an oscillator, it's probably like an oscillator
  //   return true;
  // }
  // 
  // /**
  //  * Generate a string representation of the subspace labels.
  //  * @return subspace labels as a string
  //  */
  // std::string OrbitalSubspaceLJPN::LabelStr() const {
  //   std::ostringstream os;
  // 
  //   const int width = 3;
  // 
  //   os << "["
  //      << " " << std::setw(width) << int(orbital_species())
  //      << " " << std::setw(width) << l()
  //      << " " << std::setw(width+2) << j().Str()
  //      << " " << "]";
  // 
  //   return os.str();
  // }
  // 
  // /**
  //  * Generate a string representation, useful for debugging.
  //  * @return debug string
  //  */
  // std::string OrbitalSubspaceLJPN::DebugStr() const {
  //   std::ostringstream os;
  // 
  //   const int width = 3;
  // 
  //   for (std::size_t state_index=0; state_index<size(); ++state_index) {
  //     OrbitalStateLJPN state(*this,state_index);
  // 
  //     os
  //       << " " << "index"
  //       << " " << std::setw(width) << state_index
  //       << " " << "nlj"
  //       << " " << std::setw(width) << state.n()
  //       << " " << std::setw(width) << state.l()
  //       << " " << std::setw(width+2) << state.j().Str()
  //       << " " << "weight"
  //       << " " << state.weight()
  //       << std::endl;
  //   }
  // 
  //   return os.str();
  // 
  // }
  // 
  // /**
  //  * Flatten subspace into a vector of OrbitalPNInfo objects.
  //  *
  //  * @return vector representation of subspace
  //  */
  // OrbitalPNList OrbitalSubspaceLJPN::OrbitalInfo() const
  // {
  //   OrbitalPNList orbitals;
  // 
  //   for (std::size_t state_index=0; state_index<size(); ++state_index)
  //     {
  //       OrbitalStateLJPN state(*this,state_index);
  //       orbitals.push_back(state.OrbitalInfo());
  //     }
  // 
  //   return orbitals;
  // 
  // }
  // 
  // /**
  //  * Generate a string representation of the orbital labels.
  //  * @return orbital labels as a string
  //  */
  // std::string OrbitalStateLJPN::LabelStr() const
  // {
  //   std::ostringstream os;
  // 
  //   const int width = 0;  // for now, no fixed width
  // 
  //   os << "["
  //      << " " << std::setw(width) << int(orbital_species())
  //      << " " << std::setw(width) << index()
  //      << " :"
  //      << " " << std::setw(width) << n()
  //      << " " << std::setw(width) << l()
  //      << " " << std::setw(width) << j()
  //      << " :"
  //      << " " << std::setw(width) << weight()
  //      << " " << "]";
  // 
  //   return os.str();
  // }
  // 
  // /**
  //  * Flatten state into an OrbitalPNInfo object.
  //  *
  //  * @return OrbitalPNInfo representation of state
  //  */
  // OrbitalPNInfo OrbitalStateLJPN::OrbitalInfo() const
  // {
  //   OrbitalPNInfo orbital;
  // 
  //   orbital.orbital_species = orbital_species();
  //   orbital.n = n();
  //   orbital.l = l();
  //   orbital.j = j();
  //   orbital.weight = weight();
  // 
  //   return orbital;
  // }
  // 
  // /**
  //  * Construct an Nmax-truncated single-particle space divided into LJPN
  //  * subspaces.
  //  *
  //  * @param[in] Nmax number of oscillator quanta
  //  */
  // OrbitalSpaceLJPN::OrbitalSpaceLJPN(int Nmax)
  //   : weight_max_{double(Nmax)}, is_oscillator_like_{true}, Nmax_{Nmax}
  // {
  //   // iterate over species
  //   for (OrbitalSpeciesPN orbital_species :
  //         {OrbitalSpeciesPN::kP,OrbitalSpeciesPN::kN}) {
  //     for (int l=0; l<=Nmax; ++l) {
  //       for (HalfInt j = l-HalfInt(1,2); j<=(l+HalfInt(1,2)); ++j) {
  //         if (j<0) continue;
  //         OrbitalSubspaceLJPN subspace(orbital_species,l,j,Nmax);
  //         PushSubspace(subspace);
  //       }
  //     }
  //   }
  // }
  // 
  // /**
  //  * Construct a space with LJPN subspaces from a list of orbitals.
  //  *
  //  * @param[in] states vector of orbitals
  //  */
  // OrbitalSpaceLJPN::OrbitalSpaceLJPN(const OrbitalPNList& states)
  //   : weight_max_{0.0}
  // {
  //   // collect (l,j) subspace labels sorted in canonical order
  //   std::set<OrbitalSubspaceLJPNLabels> subspace_labels_set;
  //   for (std::size_t state_index=0; state_index<states.size(); ++state_index)
  //     {
  //       OrbitalPNInfo state = states[state_index];
  //       OrbitalSubspaceLJPNLabels labels(state.orbital_species,state.l,state.j);
  //       subspace_labels_set.insert(labels);
  //     }
  // 
  //   // construct subspaces
  //   for (const OrbitalSubspaceLJPNLabels& labels : subspace_labels_set)
  //     {
  //       OrbitalSpeciesPN orbital_species;
  //       int l;
  //       HalfInt j;
  //       std::tie(orbital_species,l,j) = labels;
  //       OrbitalSubspaceLJPN subspace(orbital_species,l,j,states);
  //       PushSubspace(subspace);
  //       weight_max_ = std::max(weight_max_,subspace.weight_max());
  //     }
  // 
  //   // check if oscillator-like and set Nmax
  //   if (IsOscillatorLike_()) {
  //     is_oscillator_like_ = true;
  //     Nmax_ = int(weight_max_);
  //   } else {
  //     is_oscillator_like_ = false;
  //     Nmax_ = -1;
  //   }
  // }
  // 
  // /**
  //  * Check if space is truncated like an Nmax oscillator truncation.
  //  *
  //  * @return true if all subspaces are Nmax truncated with the same Nmax.
  //  */
  // bool OrbitalSpaceLJPN::IsOscillatorLike_() const
  // {
  // 
  //   // first check that the space contains subspaces
  //   if (size() < 1)
  //     return false;
  // 
  //   // check if subspaces are individually Nmax truncated
  //   for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
  //     {
  //       const SubspaceType& subspace = GetSubspace(subspace_index);
  //       if (!subspace.is_oscillator_like())
  //         return false;
  //     }
  // 
  //   // see if can extract viable Nmax
  //   int Nmax = GetSubspace(0).Nmax();
  //   for (std::size_t subspace_index=1; subspace_index<size(); ++subspace_index)
  //     {
  //       const SubspaceType& subspace = GetSubspace(subspace_index);
  //       if (subspace.Nmax() != Nmax)
  //         return false;
  //     }
  // 
  //   // orbitals are oscillator like!
  //   return true;
  // }
  // 
  // /**
  //  * Generate a string representation, useful for debugging.
  //  * @return debug string
  //  */
  // std::string OrbitalSpaceLJPN::DebugStr() const {
  // 
  //   std::ostringstream os;
  // 
  //   const int width = 3;
  // 
  //   for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index) {
  //     const SubspaceType& subspace = GetSubspace(subspace_index);
  // 
  //     os
  //       << " " << "index"
  //       << " " << std::setw(width) << subspace_index
  //       << " " << "species"
  //       << " " << std::setw(width) << int(subspace.orbital_species())
  //       << " " << "dim"
  //       << " " << std::setw(width) << subspace.size()
  //       << " " << std::endl;
  //   }
  // 
  //   return os.str();
  // 
  // }
  // 
  // /**
  //  * Flatten space into a vector of OrbitalPNInfo objects.
  //  *
  //  * @return vector representation of space
  //  */
  // OrbitalPNList OrbitalSpaceLJPN::OrbitalInfo() const
  // {
  //   OrbitalPNList orbitals;
  // 
  //   for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
  //     {
  //       const SubspaceType& subspace = GetSubspace(subspace_index);
  //       OrbitalPNList subspace_orbitals;
  // 
  //       // get orbitals for subspace and append to vector
  //       subspace_orbitals = subspace.OrbitalInfo();
  //       orbitals.insert(orbitals.end(),
  //                       subspace_orbitals.begin(),
  //                       subspace_orbitals.end()
  //                      );
  // 
  //     }
  // 
  //   return orbitals;
  // 
  // }
  // 
  // #ifdef BASIS_ALLOW_DEPRECATED
  // OrbitalSectorsLJPN::OrbitalSectorsLJPN(
  //     const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space
  //   )
  //   : BaseSectors(bra_space, ket_space), j0_(-1), g0_(-1), Tz0_(1)
  // {
  //   for (std::size_t bra_subspace_index=0; bra_subspace_index<bra_space.size(); ++bra_subspace_index) {
  //     for (std::size_t ket_subspace_index=0; ket_subspace_index<ket_space.size(); ++ket_subspace_index) {
  //       // retrieve subspaces
  //       const SubspaceType& bra_subspace = bra_space.GetSubspace(bra_subspace_index);
  //       const SubspaceType& ket_subspace = ket_space.GetSubspace(ket_subspace_index);
  // 
  //       // push sector
  //       PushSector(bra_subspace_index,ket_subspace_index);
  //     }
  //   }
  // }
  // #endif  // BASIS_ALLOW_DEPRECATED
  // 
  // OrbitalSectorsLJPN::OrbitalSectorsLJPN(
  //     const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space,
  //     int j0, int g0, int Tz0)
  //   : BaseSectors(bra_space, ket_space), j0_(j0), g0_(g0), Tz0_(Tz0)
  // {
  //   for (std::size_t bra_subspace_index=0; bra_subspace_index<bra_space.size(); ++bra_subspace_index) {
  //     for (std::size_t ket_subspace_index=0; ket_subspace_index<ket_space.size(); ++ket_subspace_index) {
  // 
  //       // retrieve subspaces
  //       const SubspaceType& bra_subspace = bra_space.GetSubspace(bra_subspace_index);
  //       const SubspaceType& ket_subspace = ket_space.GetSubspace(ket_subspace_index);
  // 
  //       bool allowed = true;
  //       allowed &= am::AllowedTriangle(ket_subspace.j(), j0, bra_subspace.j());
  //       allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2 == 0);
  //       allowed &= ((bra_subspace.Tz() - ket_subspace.Tz()) == Tz0);
  // 
  //       // push sector
  //       if (allowed) {
  //         PushSector(bra_subspace_index,ket_subspace_index);
  //       }
  //     }
  //   }
  // }
  // 
  // /**
  //  * Generate a string representation, useful for debugging.
  //  * @return debug string
  //  */
  // std::string OrbitalSectorsLJPN::DebugStr() const
  // {
  //   std::ostringstream os;
  //   int width = 3;
  // 
  //   for (std::size_t sector_index=0; sector_index<size(); ++sector_index) {
  //     // const basis::OrbitalSectorLJPN& sector = GetSector(sector_index);
  //     // os << sector_index+1 << sector.DebugStr();
  //     const SectorType& sector = GetSector(sector_index);
  //     os << std::setw(width) << sector_index
  //        << " bra " << std::setw(width) << sector.bra_subspace_index()
  //        << " (" << int(sector.bra_subspace().orbital_species())
  //        << ", " << sector.bra_subspace().l()
  //        << ", " << sector.bra_subspace().j().Str()  << ")"
  //        << " ket " << std::setw(width) << sector.ket_subspace_index()
  //        << " (" << int(sector.ket_subspace().orbital_species())
  //        << ", " << sector.ket_subspace().l()
  //        << ", " << sector.ket_subspace().j().Str() << ")"
  //        << std::endl;
  //   }
  //   return os.str();
  //
  // }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
