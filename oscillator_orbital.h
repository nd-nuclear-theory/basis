/****************************************************************
  @file oscillator_orbital.h

  Defines spatial (n,l) orbitals for three-dimensional harmonic oscillator.

  This is meant as an introductory example to indexing with the basis package.

  Language: C++11

  Mark A. Caprio, Patrick J. Fasano, Jakub Herko, Shwetha L. Vittal
  University of Notre Dame

  + 02/02/22 (mac): Created, drawing on nlj_orbital and lsjt_scheme.

****************************************************************/

#ifndef BASIS_OSCILLATOR_ORBITAL_H_
#define BASIS_OSCILLATOR_ORBITAL_H_

// #include <cassert>
// #include <cmath>
// #include <cstddef>
// #include <array>
// #include <iosfwd>
// #include <string>
// #include <tuple>
// #include <vector>

// #include "am/halfint.h"

#include "basis.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  /// @defgroup oscillator-subspaces Oscillator orbitals
  /// Harmonic oscillator spatial orbitals.
  ///
  /// ## Labeling ##
  ///
  /// subspace labels: (l)
  ///
  ///  * l (int): orbital angular momentum
  ///
  /// state labels within subspace: (n)
  ///
  ///  * n (int): radial quantum number (0,1,...)
  ///
  /// The parity for each orbital is deduced from the l quantum number:
  ///
  ///  * g (int): grade (=0,1) for the parity P, given by g~l
  ///
  /// The oscillator quantum number is deduced from the
  /// n and l quantum numbers:
  ///
  ///  * N (int): oscillator quanta (N=2n+l)
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## Subspaces ##
  ///
  /// Within the full space, subspaces are ordered by:
  ///
  ///    * increasing l (l=0,1,...,Nmax)
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## States ##
  ///
  /// Within a subspace, the states are ordered by:
  ///
  ///   * increasing n (n=0,1,...,(Nmax-l)/2)
  ///
  ////////////////////////////////////////////////////////////////
  /// @{

  // declarations

  class OscillatorOrbitalSubspace;
  class OscillatorOrbitalState;
  class OscillatorOrbitalSpace;

  // labels

  typedef std::tuple<int> OscillatorOrbitalSubspaceLabels;
  typedef std::tuple<int> OscillatorOrbitalStateLabels;

  // subspace

  class OscillatorOrbitalSubspace
    : public BaseSubspace<
    OscillatorOrbitalSubspace,OscillatorOrbitalSubspaceLabels,
    OscillatorOrbitalState,OscillatorOrbitalStateLabels
    >
    {

      public:

      // constructor

      OscillatorOrbitalSubspace() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      OscillatorOrbitalSubspace(int l, int Nmax);
      // Set up indexing in with truncation by oscillator quanta.

      // accessors

      int l() const {return std::get<0>(labels());}
      int g() const {return l()%2;}
      int Nmax() const {return Nmax_;}

      // diagnostic strings
      
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      // truncation
      int Nmax_;  // maximum oscillator quanta

    };

  // state

  // WIP
  
  class OrbitalStatePN
    : public BaseState<OrbitalSubspacePN>
  {

    public:

    // pass-through constructors

    OrbitalStatePN(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    OrbitalStatePN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // produce flattened orbital information
    OrbitalPNInfo OrbitalInfo() const;

    // pass-through accessors
    OrbitalSpeciesPN orbital_species() const {return subspace().orbital_species();}
    HalfInt Tz() const {return subspace().Tz();}

    // state label accessors
    int n() const {return std::get<0>(labels());}
    int l() const {return std::get<1>(labels());}
    HalfInt j() const {return std::get<2>(labels());}
    int g() const {return l()%2;}
    FullOrbitalLabels full_labels() const
    {
      return std::make_tuple(orbital_species(),n(),l(),j());
    }

    // state weight accessors
    double weight() const {return subspace().weights()[index()];}
    // Look up floating-point weight.
    int N() const {return 2*n()+l();}
    // Calculate hard-coded oscillator quantum number.

    // diagnostic strings
    std::string LabelStr() const;
    // Provide string representation of state labels.

    // comparison
    friend bool operator == (const basis::OrbitalStatePN& a1, const basis::OrbitalStatePN& a2)
    // Equality test based on labels (so permits comparison across different subspace indexings).
    {
      return (a1.labels() == a2.labels()) && (a1.subspace().labels() == a2.subspace().labels());
    }

  };

  // space

  class OrbitalSpacePN
    : public BaseSpace<OrbitalSubspacePN>
  {

    public:

    // constructor

    OrbitalSpacePN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit OrbitalSpacePN(int Nmax);
    // Set up indexing and weights in traditional oscillator Nmax
    // truncation.

    explicit OrbitalSpacePN(const OrbitalPNList& states);
    // Set up indexing for a list of states.

    // produce flattened orbital information
    OrbitalPNList OrbitalInfo() const;

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

    // test for truncation scheme
    bool IsOscillatorLike_() const;
    // Test if labeling and weights match oscillator truncation for
    // all subspaces.

  };

  // sectors -- not applicable

  ////////////////////////////////////////////////////////////////
  /// @}
  ////////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////////
  /// @defgroup ljpn-subspaces LJPN Subspaces
  /// Single particle orbitals divided into subspaces with definite l, j and
  /// orbital species.
  ///
  /// ## Labeling ##
  ///
  /// subspace labels: (species,l,j)
  ///
  ///  * species (enum): species (kP=0 for proton, kN=1 for neutron)
  ///  * l (int): orbital angular momentum
  ///  * j (HalfInt): total angular momentum
  ///
  /// state labels within subspace: (n,l,j)
  ///
  ///  * n (int): radial quantum number (0,1,...)
  ///
  /// The parity for each orbital is deduced from the l quantum number:
  ///
  ///  * g (int): grade (=0,1) for the parity P, given by g~l
  ///
  /// Each orbital (state) also has a floating point "weight"
  /// associated with it:
  ///
  ///  * weight (double): orbital weight
  ///
  /// This quantity is not considered a "label", since it is not used
  /// for lookup purposes.  In a traditional oscillator scheme, the
  /// weight is set to the oscillator quantum number w -> N=2n+l.
  ///
  /// The hard-coded oscillator quantum number is also deduced from the
  /// n and l quantum numbers, as an integer, to be used only when the
  /// orbitals are known to be oscillator orbitals:
  ///
  ///   *N (int): oscillator quanta (N=2n+l)
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## Subspaces ##
  ///
  /// Within the full space, subspaces are ordered by:
  ///    * lexicographically increasing (species,l,j)
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## States ##
  ///
  /// Within a subspace, the states are ordered by:
  ///   * increasing n
  ///
  ////////////////////////////////////////////////////////////////
  /// @{

  // declarations
  class OrbitalSubspaceLJPN;
  class OrbitalStateLJPN;
  class OrbitalSpaceLJPN;

  // labels

  typedef std::tuple<OrbitalSpeciesPN,int,HalfInt> OrbitalSubspaceLJPNLabels;
  typedef std::tuple<int> OrbitalStateLJPNLabels;

  // subspace

  class OrbitalSubspaceLJPN
    : public BaseSubspace<OrbitalSubspaceLJPN,OrbitalSubspaceLJPNLabels,OrbitalStateLJPN,OrbitalStateLJPNLabels>
    {

      public:

      // constructor

      OrbitalSubspaceLJPN() = default;
      /// default constructor -- provided since required for certain
      /// purposes by STL container classes (e.g., std::vector::resize)

      OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species, int l, HalfInt j, int Nmax);
      /// Set up indexing and weights in traditional oscillator Nmax
      /// truncation.

      OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species, int l, HalfInt j,
        const OrbitalPNList& states);
      /// Set up indexing for a list of states.

      OrbitalPNList OrbitalInfo() const;
      /// produce flattened orbital information

      // accessors

      OrbitalSpeciesPN orbital_species() const {return std::get<0>(labels());}
      HalfInt Tz() const {return kOrbitalSpeciesPNCodeTz[int(orbital_species())];}
      int l() const {return std::get<1>(labels());}
      HalfInt j() const {return std::get<2>(labels());}
      int g() const {return l()%2;}
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

      // test for truncation scheme
      bool IsOscillatorLike_() const;
      // Test if labeling and weights match oscillator truncation.

      // weights
      std::vector<double> weights_;
    };

  // state

  class OrbitalStateLJPN
    : public BaseState<OrbitalSubspaceLJPN>
  {

    public:

    // pass-through constructors

    OrbitalStateLJPN(const SubspaceType& subspace, std::size_t index)
      /// Construct state by index.
      : BaseState (subspace, index) {}

    OrbitalStateLJPN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      /// Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    OrbitalPNInfo OrbitalInfo() const;
    /// produce flattened orbital information

    // pass-through accessors
    OrbitalSpeciesPN orbital_species() const {return subspace().orbital_species();}
    HalfInt Tz() const {return subspace().Tz();}
    int l() const {return subspace().l();}
    HalfInt j() const {return subspace().j();}
    int g() const {return l()%2;}
    FullOrbitalLabels full_labels() const
    {
      return std::make_tuple(orbital_species(),n(),l(),j());
    }

    // state label accessors
    int n() const {return std::get<0>(labels());}

    // state weight accessors
    double weight() const {return subspace().weights()[index()];}
    // Look up floating-point weight.
    int N() const {return 2*n()+l();}
    // Calculate hard-coded oscillator quantum number.

    // diagnostic strings
    std::string LabelStr() const;
    // Provide string representation of state labels.
  };

  // space

  class OrbitalSpaceLJPN
    : public BaseSpace<OrbitalSubspaceLJPN>
  {

    public:

    // constructor

    OrbitalSpaceLJPN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit OrbitalSpaceLJPN(int Nmax);
    // Set up indexing and weights in traditional oscillator Nmax
    // truncation.

    explicit OrbitalSpaceLJPN(const OrbitalPNList& states);
    // Set up indexing for a list of states.

    explicit OrbitalSpaceLJPN(const OrbitalSpacePN& space)
      : OrbitalSpaceLJPN(space.OrbitalInfo())
    {}
    // Set up indexing from space divided into pn-sectors.

    // accessors
    double weight_max() const {return weight_max_;}
    bool is_oscillator_like() const {return is_oscillator_like_;}
    int Nmax() const {assert(is_oscillator_like()); return Nmax_;}
      // only meaningful if oscillator scheme constructor used

    // produce flattened orbital information
    OrbitalPNList OrbitalInfo() const;

    // diagnostic string
    std::string DebugStr() const;

    private:

    // truncation
    double weight_max_;
    bool is_oscillator_like_;
    int Nmax_;  // only meaningful if oscillator scheme constructor used

    // test for truncation scheme
    bool IsOscillatorLike_() const;
    // Test if labeling and weights match oscillator truncation for
    // all subspaces.
  };

  // sectors

  class OrbitalSectorsLJPN
    : public BaseSectors<OrbitalSpaceLJPN>
  {

    public:

    // constructor

    OrbitalSectorsLJPN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    #ifdef BASIS_ALLOW_DEPRECATED
    DEPRECATED("all-to-all enumeration is now deprecated")
    OrbitalSectorsLJPN(
      const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space
    );
    DEPRECATED("all-to-all enumeration is now deprecated")
    explicit OrbitalSectorsLJPN(
        const OrbitalSpaceLJPN& space
    ) : OrbitalSectorsLJPN(space, space) {}
    // Enumerate all sector pairs ("all-to-all" sector enumeration).
    //
    // Sectors are enumerated in lexicographical order by (bra)(ket).
    #endif  // BASIS_ALLOW_DEPRECATED

    OrbitalSectorsLJPN(
      const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space,
      int j0, int g0, int Tz0
    );
    OrbitalSectorsLJPN(
        const OrbitalSpaceLJPN& space,
        int j0, int g0, int Tz0
    ) : OrbitalSectorsLJPN(space, space, j0, g0, Tz0) {}
    // Enumerate sector pairs between two spaces connected by an operator of
    // given tensorial and parity character ("constrained" sector enumeration).

    // diagnostic string
    std::string DebugStr() const;

    // accessors
    int j0()     const {return j0_;}
    int g0()     const {return g0_;}
    int Tz0()    const {return Tz0_;}

   private:
    // operator properties
    int j0_, g0_, Tz0_;
  };

  ////////////////////////////////////////////////////////////////
  /// @}
  ////////////////////////////////////////////////////////////////
}  // namespace basis

#endif  // BASIS_OSCILLATOR_ORBITAL_H_
