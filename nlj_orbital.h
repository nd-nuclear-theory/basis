/****************************************************************
  nlj_orbital.h

  Defines general single-particle orbital sets.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  7/7/16 (mac): Created (jjjpnorb_scheme), building on code from jjjt_scheme.
  7/19/16 (mac):
   - Add default constructors.
   - Use enum Rank for truncation rank.
   - Add two-body species code definitions.
   - Add GetOrbital accessors.
  7/22/16 (mac):
   - Fix reference error in TwoBodySpaceJJJPN.
   - Add debugging strings.
  9/28/16 (mac,pjf): Break out into nlj_orbital.

****************************************************************/

#ifndef NLJ_ORBITAL_H_
#define NLJ_ORBITAL_H_

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
  // subspace labels: (species)
  //
  //   species (enum): species (kP=0 for proton, kN=1 for neutron)
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
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
