/****************************************************************
  jjjpnorb_scheme.h

  Defines two-body state indexing in jjJpn coupling scheme, based on
  general single-particle orbital sets.

  Nomenclature: The "orb" in the filename flags that these definitions
  are for general orbitals, rather than oscillator orbitals (compare,
  e.g., lsjt_scheme.h or jjjt_scheme.h).  However, we do not propagate
  this "orb" nomenclature to the actual typedefs and function names
  defined in this file.  Naming conflicts would therefore arise if we
  were to wish to define a traditional "jjjpn_scheme" with hard-coded
  (N,j) orbital labels, in the spirit of jjjt_scheme.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/7/16 (mac): Created, building on code from jjjt_scheme.
  7/19/16 (mac):
   - Add default constructors.
   - Use enum Rank for truncation rank.
   - Add two-body species code definitions.

****************************************************************/

#ifndef JJJPNORB_SCHEME_H_
#define JJJPNORB_SCHEME_H_

#include <array>
#include <string>

#include "am/halfint.h"

#include "basis/indexing.h"
#include "basis/many_body.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle orbitals
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (s)
  //
  //   s (enum): species (kP=0 for proton, kN=1 for neutron)
  //
  // state labels within subspace: (n,l,j)
  //
  //   n (int): radial quantum number (0,1,...)
  //   l (int): orbital angular momentum
  //   j (HalfInt): total angular momentum
  //
  // The parity for each orbital is deduced from the l quantum number:
  //
  //   g (int): grade (=0,1) for the parity P, given by g~l
  //
  // Each orbital (state) also has a floating point "weight"
  // associated with it:
  //
  //   weight (double): orbital weight
  //
  // This quantity is not considered a "label", since it is not used
  // for lookup purposes.  In a traditional oscillator scheme, the
  // weight is set to the oscillator quantum number w -> N=2n+l.
  //
  // The hard-coded oscillator quantum number is also deduced from the
  // n and l quantum numbers, as an integer, to be used only when the
  // orbitals are known to be oscillator orbitals:
  //
  //   N (int): oscillator quanta (N=2n+l)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- enumerated species
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- lexicographically increasing (n,l,j)
  //
  ////////////////////////////////////////////////////////////////  

  // enumerated type for orbital species
  //
  // Note: Follows same sequence as MFDn, but MFDn uses 1-based
  // numbering.

  enum class OrbitalSpeciesPN {kP=0,kN=1};

  // labels

  typedef std::tuple<OrbitalSpeciesPN> OrbitalSubspacePNLabels;
  typedef std::tuple<int,int,HalfInt> OrbitalStatePNLabels;

  // subspace

  class OrbitalSubspacePN
    : public BaseSubspace<OrbitalSubspacePNLabels,OrbitalStatePNLabels>
    {
    
      public:

      // constructor

      OrbitalSubspacePN(OrbitalSpeciesPN orbital_species, int Nmax);
      // Set up indexing and weights in traditional oscillator Nmax
      // truncation.

      // accessors

      OrbitalSpeciesPN orbital_species() const {return std::get<0>(labels_);}
      double weight_max() const {return weight_max_;}
      int Nmax() const {return Nmax_;}  // only meaningful if oscillator scheme constructor used
      const std::vector<double>& weights() const {return weights_;}
      
      // diagnostic string
      std::string DebugStr() const;

      private:

      // truncation
      double weight_max_;
      int Nmax_;  // only meaningful if oscillator scheme constructor used

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

    // pass-through accessors
    OrbitalSpeciesPN orbital_species() {return Subspace().orbital_species();}

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

    OrbitalSpacePN(int Nmax);
    // Set up indexing and weights in traditional oscillator Nmax
    // truncation.

    // accessors
    double weight_max() const {return weight_max_;}
    int Nmax() const {return Nmax_;}  // only meaningful if oscillator scheme constructor used

    // diagnostic string
    std::string DebugStr() const;

    // orbital tabulation
    std::string OrbitalDefinitionStr() const;
    // Generate orbital tabulation suitable for output as an MFDn
    // Version 15 orbital file.
    //
    // See Pieter Maris's README_SPorbitals_2016June20.txt.


    private:

    // truncation
    double weight_max_;
    int Nmax_;  // only meaningful if oscillator scheme constructor used

  };

  // sectors -- not applicable

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJpn scheme with general orbitals
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (type,J,g)    P=(-)^g
  //
  //   species (enum): two-body species (equivalent to Tz)
  //   J (int): total angular momentum
  //   g (int): grade (=0,1) for the parity P
  //
  // state labels within subspace: (index1,index2)
  //
  //   index1 (int): index of particle 1 within appropriate
  //     (proton or neutron) orbital set
  //   index2 (int): index of particle 2 within appropriate
  //     (proton or neutron) orbital set
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing type (type=pp,nn,pn)
  //    -- increasing J
  //    -- increasing g (g=0,1)
  // 
  // Subspaces are pruned to those of nonzero dimension.
  //
  // Note that ordering of subspaces is lexicographic by (type,J,g).
  //
  // Truncation of the space is by one-body and two-body weights.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing index1
  //   -- increasing index2
  // and subject to:
  //   -- triangularity constraint on (j1,j2,J)
  //   -- parity constraint g1+g2~g
  //   -- in pp/nn subspaces, antisymmetry constraint J~0
  //      if index1==index2
  //
  // This basis is for *identical* particle states:
  //   -- In the pp/nn subspaces, the labels are subject to 
  //      the antisymmetry constraint (J~0) if the orbitals are identical.
  //   -- In the pp/nn subspaces, a canonical (lexicographic) ordering
  //      constraint is applied to the single-particle quantum numbers.
  //      That is, when enumerating the basis, the states
  //    
  //        |(index1,index2)...>  and  |(index2,index1)...>
  //
  //      would be redundant, and only the first (for index1<=index2) is
  //      retained.
  //
  ////////////////////////////////////////////////////////////////  

  // enumerated type for two-body state species
  //
  // Note: Ordering of pp/nn/pn labels follows the same sequence as
  // used in MFDn.  However, note that MFDn uses 1-based numbering.

  enum class TwoBodySpeciesPN {kPP=0,kNN=1,kPN=2};

  // notational definitions for two-body state species
  //
  // Use of these arrays requires conversion of the TwoBodySpeciesPN to int.
  //
  // Example:
  //   basis::TwoBodySpeciesPN two_body_species;
  //   ...
  //   os << basis::kTwoBodySpeciesPNCodeTz[int(two_body_species)];
  extern const std::array<int,3> kTwoBodySpeciesPNCodeTz;  // {+1,-1,0} -- "up quark is positive" convention
  extern const std::array<int,3> kTwoBodySpeciesPNCodeDecimal;  // {11,22,12} -- used by MFDn
  extern const std::array<const char*,3> kTwoBodySpeciesPNCodeChar;  // {"pp","nn","pn"}

  // maximum weight collection
  struct WeightMax
  {

    // constructor

    WeightMax() {};
    // default constructor

    WeightMax(basis::Rank truncation_rank, int truncation_cutoff)
    // Set conventional oscillator one-body/two-body trunctation.
    {
      // extract one-body and two-body cutoffs
      int N1max, N2max;
      std::tie(N1max,N2max) = basis::TwoBodyCutoffs(truncation_rank,truncation_cutoff);

      // save cutoffs
      one_body[0] = N1max;
      one_body[1] = N1max;
      two_body[0] = N2max;
      two_body[1] = N2max;
      two_body[2] = N2max;
    }

    // maximum weights 
    std::array<double,2> one_body;
    std::array<double,3> two_body;
  };

  // labels

  typedef std::tuple<TwoBodySpeciesPN,int,int> TwoBodySubspaceJJJPNLabels;
  typedef std::tuple<int,int> TwoBodyStateJJJPNLabels;

  // subspace

  class TwoBodySubspaceJJJPN
    : public BaseSubspace<TwoBodySubspaceJJJPNLabels,TwoBodyStateJJJPNLabels>
    {
    
      public:

      // constructor

      TwoBodySubspaceJJJPN() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      TwoBodySubspaceJJJPN(
          const OrbitalSpacePN& orbital_space,
          TwoBodySpeciesPN two_body_species, int J, int g,
          const WeightMax& weight_max
        );
      // Set up indexing.

      // accessors
      TwoBodySpeciesPN two_body_species() const {return std::get<0>(labels_);}
      int J() const {return std::get<1>(labels_);}
      int g() const {return std::get<2>(labels_);}
      const WeightMax& weight_max() const {return weight_max_;}
      const OrbitalSubspacePN& orbital_subspace1() const {return *orbital_subspace1_ptr_;}
      const OrbitalSubspacePN& orbital_subspace2() const {return *orbital_subspace2_ptr_;}

      // diagnostic string
      std::string DebugStr() const;

      private:

      // truncation
      WeightMax weight_max_;
      
      // direct access to orbital subspaces
      const OrbitalSubspacePN* orbital_subspace1_ptr_;
      const OrbitalSubspacePN* orbital_subspace2_ptr_;
    };

  // state

  class TwoBodyStateJJJPN
    : public BaseState<TwoBodySubspaceJJJPN>
  {
    
    public:

    // pass-through constructors

    TwoBodyStateJJJPN(const SubspaceType& subspace, int index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    TwoBodyStateJJJPN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors
    TwoBodySpeciesPN two_body_species() const {return Subspace().two_body_species();}
    int J() const {return Subspace().J();}
    int g() const {return Subspace().g();}
    const OrbitalSubspacePN& orbital_subspace1() const {return Subspace().orbital_subspace1();}
    const OrbitalSubspacePN& orbital_subspace2() const {return Subspace().orbital_subspace1();}

    // state label accessors
    int index1() const {return std::get<0>(GetStateLabels());}
    int index2() const {return std::get<1>(GetStateLabels());}

  };

  // space

  class TwoBodySpaceJJJPN
    : public BaseSpace<TwoBodySubspaceJJJPN>
  {
    
    public:

    // constructor

    TwoBodySpaceJJJPN() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySpaceJJJPN(
          const OrbitalSpacePN orbital_space,
          const WeightMax& weight_max
      );
    // Enumerate subspaces.

    // accessors
      const WeightMax& weight_max() const {return weight_max_;}

    // diagnostic string
    std::string DebugStr() const;

    private:
    // truncation
    WeightMax weight_max_;

  };

  // sectors

  class TwoBodySectorsJJJPN
    : public BaseSectors<TwoBodySpaceJJJPN>
  {

    public:

    // constructor

    TwoBodySectorsJJJPN() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySectorsJJJPN(
        const TwoBodySpaceJJJPN& space,
        int J0, int g0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).
    //
    // TODO: add possibility of delta Tz

  };



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
