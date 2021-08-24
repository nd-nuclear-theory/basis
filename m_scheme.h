/****************************************************************
  @file m_scheme.h

  Defines indexing for simple M-scheme code.

  Language: C++11

  Patrick J. Fasano
  University of Notre Dame

  + 07/10/16 (pjf): Created (m_scheme), building on code from nlj_orbital.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
  + 07/03/21 (pjf): Call base class constructor for initializing labels.
  + 07/04/21 (pjf): Pass derived subspace class as template argument to
    BaseSubspace.
*/
#ifndef BASIS_M_SCHEME_H_
#define BASIS_M_SCHEME_H_

#include <cassert>
#include <cstddef>
#include <array>
#include <tuple>
#include <vector>
#include <string>

#include "am/halfint.h"

#include "basis.h"
#include "proton_neutron.h"
#include "nlj_orbital.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle state label containers
  ////////////////////////////////////////////////////////////////

  typedef std::tuple<OrbitalSpeciesPN,int,int,HalfInt,HalfInt> FullSingleParticleStateLabels;
  ///< Full (species,n,l,j,m) labels for single-particle state.

  ///////////////////////////////////////////////////////////////
  // PN Subspaces
  // Single particle states divided into subspaces with definite
  // orbital species.
  //
  // ## Labeling ##
  //
  // subspace labels: (species)
  //
  //  * species (enum): species (kP=0 for proton, kN=1 for neutron)
  //
  // state labels within subspace: (n,l,j,m)
  //
  //  * n (int): radial quantum number (0,1,...)
  //  * l (int): orbital angular momentum
  //  * j (HalfInt): total angular momentum
  //  * m (HalfInt): angular momentum projection
  //
  // The parity for each state is deduced from the l quantum number:
  //
  //  * g (int): grade (=0,1) for the parity P, given by g~l
  //
  // Each state also has a floating point "weight"
  // associated with it:
  //
  //  * weight (double): orbital weight
  //
  // This quantity is not considered a "label", since it is not used
  // for lookup purposes.  In a traditional oscillator scheme, the
  // weight is set to the oscillator quantum number w -> N=2n+l.
  //
  // The hard-coded oscillator quantum number is also deduced from the
  // n and l quantum numbers, as an integer, to be used only when the
  // orbitals are known to be oscillator orbitals:
  //
  //  * N (int): oscillator quanta (N=2n+l)
  //
  ///////////////////////////////////////////////////////////////
  //
  // ## Subspaces ##
  //
  // Within the full space, subspaces are ordered by:
  //    * enumerated species
  //
  ///////////////////////////////////////////////////////////////
  //
  // ## States ##
  //
  // Within a subspace, the states are ordered by:
  //   * lexicographically increasing (n,l,j,m)
  //
  ///////////////////////////////////////////////////////////////

  // declarations
  class SingleParticleSubspacePN;
  class SingleParticleStatePN;
  class SingleParticleSpacePN;

  // labels

  typedef std::tuple<OrbitalSpeciesPN> SingleParticleSubspacePNLabels;
  typedef std::tuple<int,int,HalfInt,HalfInt> SingleParticleStatePNLabels;


  // subspace

  class SingleParticleSubspacePN
    : public BaseSubspace<SingleParticleSubspacePN,SingleParticleSubspacePNLabels,SingleParticleStatePN,SingleParticleStatePNLabels>
  {

    public:

      // constructors

      SingleParticleSubspacePN() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      explicit SingleParticleSubspacePN(const OrbitalSubspacePN& orbital_subspace);
      // Set up indexing from an orbital subspace.

      SingleParticleSubspacePN(
        OrbitalSpeciesPN orbital_species, const OrbitalPNList& states
      ) : SingleParticleSubspacePN(OrbitalSubspacePN(orbital_species, states)) {}
      // Set up indexing from a list of orbitals.


      // accessors
      OrbitalSpeciesPN orbital_species() const {return std::get<0>(labels());}
      HalfInt Tz() const {return kOrbitalSpeciesPNCodeTz[int(orbital_species())];}
      double weight_max() const {return weight_max_;}
      bool is_oscillator_like() const {return is_oscillator_like_;}
      int Nmax() const {assert(is_oscillator_like()); return Nmax_;}
      // only meaningful if oscillator scheme constructor used
      const std::vector<double>& weights() const {return weights_;}

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

    private:

      // truncation
      double weight_max_;
      bool is_oscillator_like_;
      int Nmax_;  // only meaningful if oscillator scheme constructor used

      // weights
      std::vector<double> weights_;

  };

  // state

  class SingleParticleStatePN
    : public BaseState<SingleParticleSubspacePN>
  {

    public:

      // pass-through constructors

      SingleParticleStatePN(const SubspaceType& subspace, std::size_t index)
        // Construct state by index.
        : BaseState (subspace, index) {}

      SingleParticleStatePN(const SubspaceType& subspace, const StateLabelsType& state_labels)
        // Construct state by reverse lookup on labels.
        : BaseState (subspace, state_labels) {}

      // pass-through accessors
      OrbitalSpeciesPN orbital_species() const {return subspace().orbital_species();}
      HalfInt Tz() const {return subspace().Tz();}

      // state label accessors
      int n() const {return std::get<0>(labels());}
      int l() const {return std::get<1>(labels());}
      HalfInt j() const {return std::get<2>(labels());}
      HalfInt m() const {return std::get<3>(labels());}
      int g() const {return l()%2;}
      FullSingleParticleStateLabels full_labels() const
      {
        return std::make_tuple(orbital_species(),n(),l(),j(),m());
      }

      // state weight accessors
      double weight() const {return subspace().weights()[index()];}
      // Look up floating-point weight.
      int N() const {return 2*n()+l();}
      // Calculate hard-coded oscillator quantum number.

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.

      // comparison
      friend bool operator == (const basis::SingleParticleStatePN& a1, const basis::SingleParticleStatePN& a2)
      // Equality test based on labels (so permits comparison across different subspace indexings).
      {
        return (a1.labels() == a2.labels()) && (a1.subspace().labels() == a2.subspace().labels());
      }
  };


  // space

  class SingleParticleSpacePN
    : public BaseSpace<SingleParticleSpacePN, SingleParticleSubspacePN>
  {

    public:

    // constructor

    SingleParticleSpacePN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit SingleParticleSpacePN(const OrbitalSpacePN& orbital_space);
    // Set up indexing for a list of states.

    explicit SingleParticleSpacePN(const OrbitalPNList& states)
      : SingleParticleSpacePN(OrbitalSpacePN(states)) {}
    // Set up indexing for a list of states.

    // accessors
    double weight_max() const {return weight_max_;}
    bool is_oscillator_like() const {return is_oscillator_like_;}
    int Nmax() const {assert(is_oscillator_like()); return Nmax_;}
    // only meaningful if oscillator scheme constructor used

    // diagnostic string
    std::string DebugStr() const;

    private:

    // truncation
    double weight_max_;
    bool is_oscillator_like_;
    int Nmax_;  // only meaningful if oscillator scheme constructor used

  };

  // sectors -- not applicable

  // ///////////////////////////////////////////////////////////////
  // // MPH Subspaces
  // // Many-body M-scheme states with definite M and number of
  // // particle-hole excitations.
  // //
  // // ## Labeling ##
  // //
  // // subspace labels: (species,M,num_ex)
  // //
  // //  * species (enum): species (kP=0 for proton, kN=1 for neutron)
  // //  * M (HalfInt): total angular momentum projection
  // //  * num_ex (int): number of particle-hole excitations
  // //
  // // state labels within subspace: vector of s.p. indices
  // //
  // ///////////////////////////////////////////////////////////////
  // //
  // // ## Subspaces ##
  // //
  // // Within the full space, subspaces are ordered by:
  // //    * lexicographically increasing (species,M,num_ex)
  // //
  // ///////////////////////////////////////////////////////////////
  // //
  // // ## States ##
  // //
  // // Within a subspace, the states are ordered by:
  // //   * lexicographically increasing (n,l,j,m)
  // //
  // ///////////////////////////////////////////////////////////////
  //
  // // labels
  //
  // typedef std::tuple<int,int> ManyBodySubspaceMPHLabels;
  // typedef std::vector<int> ManyBodyStateMPHLabels;
  //
  //
  // // subspace
  //
  // class ManyBodySubspaceMPH
  //   : public BaseSubspace<ManyBodySubspaceMPHLabels,ManyBodyStateMPHLabels>
  // {
  //
  //   public:
  //
  //     // constructors
  //
  //     ManyBodySubspaceMPH() = default;
  //     // default constructor -- provided since required for certain
  //     // purposes by STL container classes (e.g., std::vector::resize)
  //
  //     ManyBodySubspaceMPH(
  //       OrbitalSpeciesPN orbital_species, const OrbitalPNList& states
  //     );
  //     // Set up indexing from a list of orbitals.
  //
  //     explicit ManyBodySubspaceMPH(OrbitalSubspacePN& orbital_subspace);
  //     // Set up indexing from an orbital subspace.
  //
  //
  //     // accessors
  //     OrbitalSpeciesPN orbital_species() const {return std::get<0>(labels());}
  //     HalfInt Tz() const {return kOrbitalSpeciesPNCodeTz[int(orbital_species())];}
  //     double weight_max() const {return weight_max_;}
  //     bool is_oscillator_like() const {return is_oscillator_like_;}
  //     int Nmax() const {assert(is_oscillator_like()); return Nmax_;}
  //     // only meaningful if oscillator scheme constructor used
  //     const std::vector<double>& weights() const {return weights_;}
  //
  //     // diagnostic strings
  //     std::string LabelStr() const;
  //     // Provide string representation of subspace labels.
  //     std::string DebugStr() const;
  //     // Dump subspace contents.
  //
  //   private:
  //
  //     // truncation
  //     double weight_max_;
  //     bool is_oscillator_like_;
  //     int Nmax_;  // only meaningful if oscillator scheme constructor used
  //
  //     // test for truncation scheme
  //     bool IsOscillatorLike_() const;
  //     // Test if labeling and weights match oscillator truncation.
  //
  //     // weights
  //     std::vector<double> weights_;
  //
  // };
  //
  // // state
  //
  // class ManyBodyStateMPH
  //   : public BaseState<ManyBodySubspaceMPH>
  // {
  //
  //   public:
  //
  //     // pass-through constructors
  //
  //     ManyBodyStateMPH(const SubspaceType& subspace, int index)
  //       // Construct state by index.
  //       : BaseState (subspace, index) {}
  //
  //     ManyBodyStateMPH(const SubspaceType& subspace, const StateLabelsType& state_labels)
  //       // Construct state by reverse lookup on labels.
  //       : BaseState (subspace, state_labels) {}
  //
  //     // pass-through accessors
  //     OrbitalSpeciesPN orbital_species() const {return subspace().orbital_species();}
  //     HalfInt Tz() const {return subspace().Tz();}
  //
  //     // state label accessors
  //     int n() const {return std::get<0>(labels());}
  //     int l() const {return std::get<1>(labels());}
  //     HalfInt j() const {return std::get<2>(labels());}
  //     HalfInt m() const {return std::get<3>(labels());}
  //     int g() const {return l()%2;}
  //     FullSingleParticleStateLabels full_labels() const
  //     {
  //       return std::make_tuple(orbital_species(),n(),l(),j(),m());
  //     }
  //
  //     // state weight accessors
  //     double weight() const {return subspace().weights()[index()];}
  //     // Look up floating-point weight.
  //     int N() const {return 2*n()+l();}
  //     // Calculate hard-coded oscillator quantum number.
  //
  //     // diagnostic strings
  //     std::string LabelStr() const;
  //     // Provide string representation of subspace labels.
  //
  //     // comparison
  //     friend bool operator == (const basis::ManyBodyStateMPH& a1, const basis::ManyBodyStateMPH& a2)
  //     // Equality test based on labels (so permits comparison across different subspace indexings).
  //     {
  //       return (a1.labels() == a2.labels()) && (a1.subspace().labels() == a2.subspace().labels());
  //     }
  // };
  //
  //
  // // space
  //
  // class ManyBodySpaceMPH
  //   : public BaseSpace<ManyBodySubspaceMPH>
  // {
  //
  //   public:
  //
  //   // constructor
  //
  //   ManyBodySpaceMPH() = default;
  //   // default constructor -- provided since required for certain
  //   // purposes by STL container classes (e.g., std::vector::resize)
  //
  //   explicit ManyBodySpaceMPH(int Nmax);
  //   // Set up indexing and weights in traditional oscillator Nmax
  //   // truncation.
  //
  //   explicit ManyBodySpaceMPH(const OrbitalPNList& states);
  //   // Set up indexing for a list of states.
  //
  //   explicit ManyBodySpaceMPH(const OrbitalSpacePN& orbital_space);
  //   // Set up indexing for a list of states.
  //
  //   // accessors
  //   double weight_max() const {return weight_max_;}
  //   bool is_oscillator_like() const {return is_oscillator_like_;}
  //   int Nmax() const {assert(is_oscillator_like()); return Nmax_;}
  //   // only meaningful if oscillator scheme constructor used
  //
  //   // diagnostic string
  //   std::string DebugStr() const;
  //
  //   private:
  //
  //   // truncation
  //   double weight_max_;
  //   bool is_oscillator_like_;
  //   int Nmax_;  // only meaningful if oscillator scheme constructor used
  //
  //   // test for truncation scheme
  //   bool IsOscillatorLike_() const;
  //   // Test if labeling and weights match oscillator truncation for
  //   // all subspaces.
  //
  // };
  //
  // // sectors -- not applicable

}  // namespace basis

#endif  // BASIS_M_SCHEME_H_
