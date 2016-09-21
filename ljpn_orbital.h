/****************************************************************
  ljpn_orbital.h


  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  9/16/16 (mac): Created, building on code from jjjpnorb_scheme.

****************************************************************/

#ifndef LJPN_OR_H_
#define LJPN_SPACES_H_

#include "basis/indexing.h"
#include "basis/many_body.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle orbitals by (l,j,s)
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (l,j,s)
  //
  //   l (int): orbital angular momentum
  //   j (HalfInt): total angular momentum
  //   s (enum): species (kP=0 for proton, kN=1 for neutron)
  //
  // The parity label g can implicitly be deduced from l if needed.
  // We do not explicitly store it.
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered lexicographically
  // (s,l,j).
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // We will not actually enumerate states within a subspace.  We
  // provide a "dummy" int state label where a StateLabelType is
  // required by the basis template library.
  //
  // However, we keep track of the subspace dimension.  The subspace
  // dimensions must be provided to the constructor via a mapping U3S
  // -> dimension.  Then we must explicitly set the dimension_ member
  // variable (inherited from BaseSubspace).
  //
  ////////////////////////////////////////////////////////////////

  // labels

  typedef std::tuple<int,HalfInt,OrbitalSpeciesPN> OrbitalSubspaceLJPNLabels;
  // typedef std::tuple<int> OrbitalStateLJPNLabels;

  // subspace

  class OrbitalSubspaceLJPN
    : public BaseSubspace<OrbitalSubspaceLJPNLabels,int>
    {
    
      public:

      // constructor

      OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species, int l, HalfInt j, int Nmax);
      // Set up indexing and weights in traditional oscillator Nmax
      // truncation.

      // accessors

      OrbitalSpeciesPN orbital_species() const {return std::get<0>(labels_);}
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

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

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
    const OrbitalSubspacePN& orbital_subspace2() const {return Subspace().orbital_subspace2();}

    // state label accessors
    int index1() const {return std::get<0>(GetStateLabels());}
    int index2() const {return std::get<1>(GetStateLabels());}

    // state retrieval
    const OrbitalStatePN GetOrbital1() const {return OrbitalStatePN(orbital_subspace1(),index1());}
    const OrbitalStatePN GetOrbital2() const {return OrbitalStatePN(orbital_subspace2(),index2());}

    // diagnostic string
    std::string LabelStr() const;

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
          const OrbitalSpacePN& orbital_space,
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
