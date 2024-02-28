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
  ///
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
  ///  * g (int): grade (=0,1) for the parity P, given by g~l (mod 2)
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
      // Set up indexing with truncation by oscillator quanta.

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

      //validation
      bool ValidLabels() const;

      // truncation
      int Nmax_;  // maximum oscillator quanta

    };

  // state

  class OscillatorOrbitalState
    : public BaseState<OscillatorOrbitalSubspace>
  {

    public:

    // pass-through constructors

    OscillatorOrbitalState(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    OscillatorOrbitalState(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors (from subspace)
    int l() const {return subspace().l();}
    int g() const {return subspace().g();}

    // state label accessors
    int n() const {return std::get<0>(labels());}

    // derived quantum numbers
    int N() const {return 2*n()+l();}
    // Calculate oscillator quantum number.

    // diagnostic strings
    std::string LabelStr() const;
    // Provide string representation of state labels.

    // comparison
    friend bool operator == (const basis::OscillatorOrbitalState& a1, const basis::OscillatorOrbitalState& a2)
    // Equality test based on labels (so permits comparison across different subspace indexings).
    {
      return (a1.labels() == a2.labels()) && (a1.subspace().labels() == a2.subspace().labels());
    }

  };

  // space
  
  class OscillatorOrbitalSpace
    : public BaseSpace<OscillatorOrbitalSpace,OscillatorOrbitalSubspace>
  {

    public:

    // constructor

    OscillatorOrbitalSpace() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    explicit OscillatorOrbitalSpace(int Nmax);
    // Set up indexing and weights in traditional oscillator Nmax
    // truncation.

    // accessors
    int Nmax() const {return Nmax_;}

    // diagnostic string
    std::string DebugStr() const;

    private:

    // truncation
    int Nmax_;

  };

  // sectors

  class OscillatorOrbitalSectors
    : public BaseSectors<OscillatorOrbitalSpace>
  {

    public:

    // constructor

    OscillatorOrbitalSectors() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    OscillatorOrbitalSectors(
        const OscillatorOrbitalSpace& space,
        int L0, int g0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).
    //
    // Arguments:
    //
    //   space (const OscillatorOrbitalSpace& space): underlying space
    //
    //   L0 (int): spherical tensor rank of operator
    //
    //   g0 (int): parity "grade" of operator
    //
    //   sector_direction (basis::SectorDirection): whether to retain only
    //     canonical (upper triangle) sectors or both canonical and
    //     non-canonical (upper and lower triangles)

    // accessors
    int L0() const {return L0_;}
    int g0() const {return g0_;}

   private:
    // operator properties
    int L0_, g0_;
    
  };
  
  ////////////////////////////////////////////////////////////////
  /// @}
  ////////////////////////////////////////////////////////////////
}  // namespace basis

#endif  // BASIS_OSCILLATOR_ORBITAL_H_
