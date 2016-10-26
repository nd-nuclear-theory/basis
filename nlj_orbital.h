/// @file
/****************************************************************
  nlj_orbital.h

  Defines general single-particle orbital sets.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 7/7/16 (mac): Created (jjjpnorb_scheme), building on code from jjjt_scheme.
  + 7/19/16 (mac):
   - Add default constructors.
   - Use enum Rank for truncation rank.
   - Add two-body species code definitions.
   - Add GetOrbital accessors.
  + 7/22/16 (mac):
   - Fix reference error in TwoBodySpaceJJJPN.
   - Add debugging strings.
  + 9/28/16 (mac,pjf): Break out into nlj_orbital.
  + 10/6/16 (pjf): Add LJPN classes
  + 10/7/16 (pjf): Add LJPN sectors and general constructors
  + 10/13/16 (mac): Add default constructors.
  + 10/18/16 (pjf): Add PN general constructors and serializers (OrbitalInfo).
  + 10/18/16 (pjf):
   - Add general constructors and serializers for Orbital*PN (OrbitalPNInfo).
   - Refactor OrbitalDefinitionStr out of spaces.
   - Add oscillator-likeness tests to OrbitalSpacePN and OrbitalSubspacePN.
   - Add == operator for OrbitalPNInfo.
  + 10/19/16 (pjf): Update documentation to Doxygen format.
  + 10/19/16 (mac):
    - Add default constructors yet again.
    - Make orbital I/O switchable between standalone and embedded modes.
  + 10/22/16 (mac): Augment debug strings.
  + 10/25/16 (pjf):
    - Fix constrained sector enumeration to use l0max and Tz0.
    - Add sector enumeration between different bra and ket spaces.
    - Define array mappings for orbital species.
    - Define Tz accessor for states and subspaces.

****************************************************************/

#ifndef NLJ_ORBITAL_H_
#define NLJ_ORBITAL_H_

#include <array>
#include <string>

#include "am/halfint.h"
#include "mcutils/parsing.h"

#include "basis/indexing.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle orbitals
  ////////////////////////////////////////////////////////////////

  /**
   * Enumerated type for orbital species
   *
   * Note: Follows same sequence as MFDn, but MFDn uses 1-based
   * numbering.
   */
  enum class OrbitalSpeciesPN {kP=0,kN=1};

  /**
   * Notational definition for orbital species conversion to Tz.
   * Uses "up quark is positive" convention.
   *
   * Use of this array requires conversion of the OrbitalSpeciesPN to int.
   *
   * Example:
   * \code
   *   basis::OrbitalSpeciesPN orbital_species;
   *   ...
   *   os << basis::kOrbitalSpeciesPNCodeTz[int(orbital_species)];
   * \endcode
   */
  extern const std::array<HalfInt, 2> kOrbitalSpeciesPNCodeTz;
  /**
   * Notational definition for orbital species conversion to MFDn decimal code.
   * "p" -> 1, "n" -> 2
   *
   * Use of this array requires conversion of the OrbitalSpeciesPN to int.
   *
   * Example:
   * \code
   *   basis::OrbitalSpeciesPN orbital_species;
   *   ...
   *   os << basis::kOrbitalSpeciesPNCodeDecimal[int(orbital_species)];
   * \endcode
   */
  extern const std::array<int, 2> kOrbitalSpeciesPNCodeDecimal;
  /**
   * Notational definition for orbital species conversion to character name.
   *
   * Use of this array requires conversion of the OrbitalSpeciesPN to int.
   *
   * Example:
   * \code
   *   basis::OrbitalSpeciesPN orbital_species;
   *   ...
   *   os << basis::kOrbitalSpeciesPNCodeChar[int(orbital_species)];
   * \endcode
   */
  extern const std::array<const char*, 2> kOrbitalSpeciesPNCodeChar;

  /**
   * A flattened orbital container.
   *
   * All quantum numbers used by MFDn are contained in this simple struct.
   */
  struct OrbitalPNInfo
  {
    OrbitalSpeciesPN orbital_species;
    int n;
    int l;
    HalfInt j;
    double weight;

    OrbitalPNInfo() = default;
    // default constructor

    OrbitalPNInfo(OrbitalSpeciesPN os, int n, int l, HalfInt j, double weight)
      : orbital_species(os), n(n), l(l), j(j), weight(weight) {};

    static constexpr double kWeightTolerance = 1e-8;
    inline bool operator==(const OrbitalPNInfo& rhs) const
    {
      bool equiv = true;
      equiv &= (this->orbital_species==rhs.orbital_species);
      equiv &= (this->n==rhs.n);
      equiv &= (this->l==rhs.l);
      equiv &= (this->j==rhs.j);
      equiv &= (abs(this->weight-rhs.weight)<kWeightTolerance);
      return equiv;
    }
  };

  // orbital I/O

  std::vector<OrbitalPNInfo> ParseOrbitalPNStream(std::istream& is, bool standalone);
  // Read orbital definitions from a stream.
  //
  // Arguments:
  //   is (std::istream, input) :
  //     input stream containing MFDn-formatted orbital definitions
  //   standalone (bool): whether or not to expect initial comments and version number
  //     as for standalone orbital file
  //
  // Returns:
  //   (std::vector<OrbitalPNInfo>) : list of flattened orbital parameters

  std::string OrbitalDefinitionStr(const std::vector<OrbitalPNInfo>& orbitals, bool standalone = false);
  // Output orbital info as a string suitable for MFDn version 15.
  //
  // Arguments:
  //   orbitals (const std::vector<OrbitalPNInfo>&, input) :
  //     list of flattened orbital parameters
  //   standalone (bool, optional): whether or not to include initial comments and version number
  //     as for standalone orbital file
  //
  // Returns:
  //   (std::string) output stream containing MFDn-formatted orbital definitions

  ////////////////////////////////////////////////////////////////
  /// @defgroup pn-subspaces PN Subspaces
  /// Single particle orbitals divided into subspaces with definite
  /// orbital species.
  ///
  /// ## Labeling ##
  ///
  /// subspace labels: (species)
  ///
  ///  * species (enum): species (kP=0 for proton, kN=1 for neutron)
  ///
  /// state labels within subspace: (n,l,j)
  ///
  ///  * n (int): radial quantum number (0,1,...)
  ///  * l (int): orbital angular momentum
  ///  * j (HalfInt): total angular momentum
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
  ///  * N (int): oscillator quanta (N=2n+l)
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## Subspaces ##
  ///
  /// Within the full space, subspaces are ordered by:
  ///    * enumerated species
  ///
  ////////////////////////////////////////////////////////////////
  ///
  /// ## States ##
  ///
  /// Within a subspace, the states are ordered by:
  ///   * lexicographically increasing (n,l,j)
  ///
  ////////////////////////////////////////////////////////////////
  /// @{

  // labels

  typedef std::tuple<OrbitalSpeciesPN> OrbitalSubspacePNLabels;
  typedef std::tuple<int,int,HalfInt> OrbitalStatePNLabels;


  // subspace

  class OrbitalSubspacePN
    : public BaseSubspace<OrbitalSubspacePNLabels,OrbitalStatePNLabels>
    {

      public:

      // constructor

      OrbitalSubspacePN() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      OrbitalSubspacePN(OrbitalSpeciesPN orbital_species, int Nmax);
      // Set up indexing and weights in traditional oscillator Nmax
      // truncation.

      OrbitalSubspacePN(OrbitalSpeciesPN orbital_species,
        const std::vector<OrbitalPNInfo>& states);
      // Set up indexing for a list of states.

      // produce flattened orbital information
      std::vector<OrbitalPNInfo> OrbitalInfo() const;

      // accessors

      OrbitalSpeciesPN orbital_species() const {return std::get<0>(labels_);}
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

      // test for truncation scheme
      bool IsOscillatorLike_() const;
      // Test if labeling and weights match oscillator truncation.

      // weights
      std::vector<double> weights_;
    };

  // state

  class OrbitalStatePN
    : public BaseState<OrbitalSubspacePN>
  {

    public:

    // pass-through constructors

    OrbitalStatePN(const SubspaceType& subspace, int index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    OrbitalStatePN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // produce flattened orbital information
    OrbitalPNInfo OrbitalInfo() const;

    // pass-through accessors
    OrbitalSpeciesPN orbital_species() const {return Subspace().orbital_species();}
    HalfInt Tz() const {return Subspace().Tz();}

    // state label accessors
    int n() const {return std::get<0>(GetStateLabels());}
    int l() const {return std::get<1>(GetStateLabels());}
    HalfInt j() const {return std::get<2>(GetStateLabels());}
    int g() const {return l()%2;}

    // state weight accessors
    double weight() const {return Subspace().weights()[index()];}
    // Look up floating-point weight.
    int N() const {return 2*n()+l();}
    // Calculate hard-coded oscillator quantum number.

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

    OrbitalSpacePN(int Nmax);
    // Set up indexing and weights in traditional oscillator Nmax
    // truncation.

    OrbitalSpacePN(const std::vector<OrbitalPNInfo>& states);
    // Set up indexing for a list of states.

    // produce flattened orbital information
    std::vector<OrbitalPNInfo> OrbitalInfo() const;

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

  // labels

  typedef std::tuple<OrbitalSpeciesPN,int,HalfInt> OrbitalSubspaceLJPNLabels;
  typedef std::tuple<int> OrbitalStateLJPNLabels;

  // subspace

  class OrbitalSubspaceLJPN
    : public BaseSubspace<OrbitalSubspaceLJPNLabels,OrbitalStateLJPNLabels>
    {

      public:

      // constructor

      OrbitalSubspaceLJPN() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species, int l, HalfInt j, int Nmax);
      // Set up indexing and weights in traditional oscillator Nmax
      // truncation.

      OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species, int l, HalfInt j,
        const std::vector<OrbitalPNInfo>& states);
      // Set up indexing for a list of states.

      // produce flattened orbital information
      std::vector<OrbitalPNInfo> OrbitalInfo() const;

      // accessors

      OrbitalSpeciesPN orbital_species() const {return std::get<0>(labels_);}
      HalfInt Tz() const {return kOrbitalSpeciesPNCodeTz[int(orbital_species())];}
      int l() const {return std::get<1>(labels_);}
      HalfInt j() const {return std::get<2>(labels_);}
      int g() const {return l()%2;}
      double weight_max() const {return weight_max_;}
      int Nmax() const {return Nmax_;}  // only meaningful if oscillator scheme constructor used
      const std::vector<double>& weights() const {return weights_;}

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      // truncation
      double weight_max_;
      int Nmax_;  // only meaningful if oscillator scheme constructor used

      // weights
      std::vector<double> weights_;
    };

  // state

  class OrbitalStateLJPN
    : public BaseState<OrbitalSubspaceLJPN>
  {

    public:

    // pass-through constructors

    OrbitalStateLJPN(const SubspaceType& subspace, int index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    OrbitalStateLJPN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // produce flattened orbital information
    OrbitalPNInfo OrbitalInfo() const;

    // pass-through accessors
    OrbitalSpeciesPN orbital_species() const {return Subspace().orbital_species();}
    HalfInt Tz() const {return Subspace().Tz();}
    int l() const {return Subspace().l();}
    HalfInt j() const {return Subspace().j();}
    int g() const {return l()%2;}

    // state label accessors
    int n() const {return std::get<0>(GetStateLabels());}

    // state weight accessors
    double weight() const {return Subspace().weights()[index()];}
    // Look up floating-point weight.
    int N() const {return 2*n()+l();}
    // Calculate hard-coded oscillator quantum number.

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

    OrbitalSpaceLJPN(int Nmax);
    // Set up indexing and weights in traditional oscillator Nmax
    // truncation.

    OrbitalSpaceLJPN(const std::vector<OrbitalPNInfo>& states);
    // Set up indexing for a list of states.

    // accessors
    double weight_max() const {return weight_max_;}
    int Nmax() const {return Nmax_;}  // only meaningful if oscillator scheme
                                      // constructor used

    // produce flattened orbital information
    std::vector<OrbitalPNInfo> OrbitalInfo() const;

    // diagnostic string
    std::string DebugStr() const;

    private:

    // truncation
    double weight_max_;
    int Nmax_;  // only meaningful if oscillator scheme constructor used

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

    OrbitalSectorsLJPN(
        const OrbitalSpaceLJPN& space,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate all sector pairs ("all-to-all" sector enumeration).
    //
    // Sectors are enumerated in lexicographical order by (bra)(ket).

    OrbitalSectorsLJPN(
        const OrbitalSpaceLJPN& space,
        int l0max, int Tz0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

    OrbitalSectorsLJPN(
        const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space
      );
    // Enumerate all sector pairs between two spaces ("all-to-all" sector
    // enumeration).
    //
    // Sectors are enumerated in lexicographical order by (bra)(ket).

    OrbitalSectorsLJPN(
        const OrbitalSpaceLJPN& bra_space, const OrbitalSpaceLJPN& ket_space,
        int l0max, int Tz0
      );
    // Enumerate sector pairs between two spaces connected by an operator of
    // given tensorial and parity character ("constrained" sector
    // enumeration).

    // diagnostic string
    std::string DebugStr() const;

    // accessors
    int l0max() const {return l0max_;}
    int Tz0() const {return Tz0_;}

   private:
    // operator properties
    int l0max_, Tz0_;
  };

  ////////////////////////////////////////////////////////////////
  /// @}
  ////////////////////////////////////////////////////////////////
}  // namespace basis

#endif  // NLJ_ORBITAL_H_
