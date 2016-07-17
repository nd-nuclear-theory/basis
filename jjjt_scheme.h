/****************************************************************
  jjjt_scheme.h

  Defines two-body state indexing in jjJT coupling scheme.
  Written for use in Moshinsky transformation.

  Indexing is in each case broken into a *subspace* (described by JTg
  labels and an Nmax trunction) and then indexed states within the
  subspace.

  Variant indexing schemes where the subspaces are further broken down
  by the number of oscillator quanta N are then provided, so that the
  Moshinsky transformation implementation can be more easily broken
  down by N blocks.

  See notes "Moshinsky xform for operators" (2015) for derivation
  of indexing scheme.

  Language: C++11
                                 
  Mark A. Caprio
  University of Notre Dame

  11/26/15 (mac): Created, following analogous code for LSJT
    scheme in indexing_lsjt.
  7/6/16 (mac): Overhaul to new basis module conventions.
  7/7/16 (mac): Add fixed-N subspaces in two-body basis.
  7/9/16 (mac): Add default constructors.
  7/16/16 (mac):
    - Add debug strings.
    - Move N to least significant subspace label in NLSJT basis.

****************************************************************/

#ifndef JJJT_SCHEME_H_
#define JJJT_SCHEME_H_

#include <string>

#include "am/halfint.h"

#include "basis/indexing.h"


namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJT scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (J,T,g)    P=(-)^g
  //
  //   J (int): total angular momentum
  //   T (int): isospin
  //   g (int): grade (=0,1) for the parity P
  //
  // state labels within subspace: (N1,j1,N2,j2)
  //
  //   N1 (int): oscillator quanta of particle 1
  //   j1 (HalfInt): total angular momentum of particle 1
  //   N2 (int): oscillator quanta of particle 2
  //   j2 (HalfInt): total angular momentum of particle 2
  //
  // Note that the labels l1 and l2 are not explicitly stored, as they
  // may be readily covered from the N1 and j1 or N2 and j2 labels:
  //
  //   l1 (int): orbital angular momentum of particle 1
  //   l2 (int): orbital angular momentum of particle 2
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing J (J=0,1,...,Nmax+1)
  //    -- increasing T (T=0,1)
  //    -- increasing g (g=0,1)
  // subject to:
  //    -- [implicit constraint J=Nmax+1 is excluded for T=1,
  //       but this is simply enforced by pruning to subspaces of 
  //       nonzero dimension]
  // 
  // Subspaces are pruned to those of nonzero dimension.
  //
  // Note that ordering of subspaces is lexicographic by (J,T,g).
  //
  // Truncation of the space is by the two-body Nmax.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing N (N~g)
  //   -- lexicographically increasing (N1,j1)
  //   -- lexicographically increasing (N2,j2)
  // and subject to:
  //   -- triangularity constraint on (j1,j2,J)
  //   -- parity constraint N~g
  //   -- antisymmetry constraint J+T~1 if (N1,j1)==(N2,j2)
  //
  // This basis is for *identical* particle states:
  //   -- The labels are subject to the antisymmetry constraint 
  //      (J+T~1) if the orbitals are identical.
  //   -- A canonical (lexicographic) ordering constraint is applied to the 
  //      single-particle quantum numbers.  That is, when
  //      enumerating the basis, the states
  //    
  //        |((N1,j1),(N2,j2))...>  and  |((N2,j2),(N1,j1))...>
  //
  //      would be redundant, and only the first (for (N1,j1)<=(N2,j2)) is
  //      retained.
  //
  ////////////////////////////////////////////////////////////////  

  // labels

  typedef std::tuple<int,int,int> TwoBodySubspaceJJJTLabels;
  typedef std::tuple<int,HalfInt,int,HalfInt> TwoBodyStateJJJTLabels;

  // subspace

  class TwoBodySubspaceJJJT
    : public BaseSubspace<TwoBodySubspaceJJJTLabels,TwoBodyStateJJJTLabels>
    {
    
      public:

      // constructor

      TwoBodySubspaceJJJT() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      TwoBodySubspaceJJJT(int J, int T, int g, int Nmax);
      // Set up indexing in Nmax truncation.

      // accessors

      int J() const {return std::get<0>(labels_);}
      int T() const {return std::get<1>(labels_);}
      int g() const {return std::get<2>(labels_);}
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
      int Nmax_;
    };

  // state

  class TwoBodyStateJJJT
    : public BaseState<TwoBodySubspaceJJJT>
  {
    
    public:

    // pass-through constructors

    TwoBodyStateJJJT(const SubspaceType& subspace, int index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    TwoBodyStateJJJT(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors
    int J() const {return Subspace().J();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}

    // state label accessors
    int N1() const {return std::get<0>(GetStateLabels());}
    HalfInt j1() const {return std::get<1>(GetStateLabels());}
    int N2() const {return std::get<2>(GetStateLabels());}
    HalfInt j2() const {return std::get<3>(GetStateLabels());}

    int l1() const
    {
      int N = N1();
      HalfInt j = j1();
      return (TwiceValue(j)-1)/2 + (N+(TwiceValue(j)-1)/2)%2;
    }
    int l2() const
    {
      int N = N2();
      HalfInt j = j2();
      return (TwiceValue(j)-1)/2 + (N+(TwiceValue(j)-1)/2)%2;
    }

    int N() const {return N1()+N2();}

  };

  // space

  class TwoBodySpaceJJJT
    : public BaseSpace<TwoBodySubspaceJJJT>
  {
    
    public:

    // constructor

    TwoBodySpaceJJJT() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySpaceJJJT(int Nmax);
    // Enumerate all subspaces up to a given Nmax cutoff.

    // accessors
    int Nmax() const {return Nmax_;}

    // diagnostic string
    std::string DebugStr() const;
    // Dump space contents.

    private:
    // truncation
    int Nmax_;

  };

  // sectors

  class TwoBodySectorsJJJT
    : public BaseSectors<TwoBodySpaceJJJT>
  {

    public:

    // constructor
    
    TwoBodySectorsJJJT() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySectorsJJJT(
        const TwoBodySpaceJJJT& space,
        int J0, int T0, int g0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

  };

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJT scheme -- subspaced by N
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (J,T,g,N)    P=(-)^g  (MODIFICATION for subspacing by N)
  //
  //   J (int): total angular momentum
  //   T (int): isospin
  //   g (int): grade (=0,1) for the parity P
  //
  // state labels within subspace: (N1,j1,N2,j2)
  //
  //   N1 (int): oscillator quanta of particle 1
  //   j1 (HalfInt): total angular momentum of particle 1
  //   N2 (int): oscillator quanta of particle 2
  //   j2 (HalfInt): total angular momentum of particle 2
  //
  // Note that the labels l1 and l2 are not explicitly stored, as they
  // may be readily covered from the N1 and j1 or N2 and j2 labels:
  //
  //   l1 (int): orbital angular momentum of particle 1
  //   l2 (int): orbital angular momentum of particle 2
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing J (J=0,1,...,Nmax+1)
  //    -- increasing T (T=0,1)
  //    -- increasing g (g=0,1)
  //    -- increasing N (N=0,1,...,Nmax)  (MODIFICATION for subspacing by N)
  // subject to:
  //    -- [implicit constraint J=Nmax+1 is excluded for T=1,
  //       but this is simply enforced by pruning to subspaces of 
  //       nonzero dimension]
  //    -- parity constraint N~g  (MODIFICATION for subspacing by N)
  // 
  // Subspaces are pruned to those of nonzero dimension.
  //
  // Note that ordering of subspaces is lexicographic by (J,T,g,N).  (MODIFICATION for subspacing by N)
  //
  // Truncation of the space is by the two-body Nmax.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- [increasing N  (N=N1+N2)]  (MODIFICATION for subspacing by N)
  //   -- lexicographically increasing (N1,j1)
  //   -- lexicographically increasing (N2,j2)
  // and subject to:
  //   -- triangularity constraint on (j1,j2,J)
  //   -- [parity constraint N~g]  (MODIFICATION for subspacing by N)
  //   -- antisymmetry constraint J+T~1 if (N1,j1)==(N2,j2)
  //
  // This basis is for *identical* particle states:
  //   -- The labels are subject to the antisymmetry constraint 
  //      (J+T~1) if the orbitals are identical.
  //   -- A canonical (lexicographic) ordering constraint is applied to the 
  //      single-particle quantum numbers.  That is, when
  //      enumerating the basis, the states
  //    
  //        |((N1,j1),(N2,j2))...>  and  |((N2,j2),(N1,j1))...>
  //
  //      would be redundant, and only the first (for (N1,j1)<=(N2,j2)) is
  //      retained.
  //
  ////////////////////////////////////////////////////////////////  

  // Modification for subspacing by N is by lexical replacement JJJT
  // -> NJJJT plus specific mods as flagged by MODIFICATION comments
  // in code.


  // labels

  typedef std::tuple<int,int,int,int> TwoBodySubspaceNJJJTLabels;  // (MODIFICATION for subspacing by N)
  typedef std::tuple<int,HalfInt,int,HalfInt> TwoBodyStateNJJJTLabels;

  // subspace

  class TwoBodySubspaceNJJJT
    : public BaseSubspace<TwoBodySubspaceNJJJTLabels,TwoBodyStateNJJJTLabels>
    {
    
      public:

      // constructor

      TwoBodySubspaceNJJJT() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      TwoBodySubspaceNJJJT(int J, int T, int g, int N);  // (MODIFICATION for subspacing by N)
      // Set up indexing in Nmax truncation.

      // accessors

      int J() const {return std::get<0>(labels_);}
      int T() const {return std::get<1>(labels_);}
      int g() const {return std::get<2>(labels_);}
      int N() const {return std::get<3>(labels_);}  // (MODIFICATION for subspacing by N)

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      //validation
      bool ValidLabels() const;

      // truncation
      int N_;  // (MODIFICATION for subspacing by N)
    };

  // state

  class TwoBodyStateNJJJT
    : public BaseState<TwoBodySubspaceNJJJT>
  {
    
    public:

    // pass-through constructors

    TwoBodyStateNJJJT(const SubspaceType& subspace, int index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    TwoBodyStateNJJJT(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors
    int J() const {return Subspace().J();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}

    // state label accessors
    int N1() const {return std::get<0>(GetStateLabels());}
    HalfInt j1() const {return std::get<1>(GetStateLabels());}
    int N2() const {return std::get<2>(GetStateLabels());}
    HalfInt j2() const {return std::get<3>(GetStateLabels());}

    int l1() const
    {
      int N = N1();
      HalfInt j = j1();
      return (TwiceValue(j)-1)/2 + (N+(TwiceValue(j)-1)/2)%2;
    }
    int l2() const
    {
      int N = N2();
      HalfInt j = j2();
      return (TwiceValue(j)-1)/2 + (N+(TwiceValue(j)-1)/2)%2;
    }

    int N() const {return N1()+N2();}

  };

  // space

  class TwoBodySpaceNJJJT
    : public BaseSpace<TwoBodySubspaceNJJJT>
  {
    
    public:

    // constructor

    TwoBodySpaceNJJJT() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySpaceNJJJT(int Nmax);
    // Enumerate all subspaces up to a given Nmax cutoff.

    // accessors
    int Nmax() const {return Nmax_;}

    // diagnostic string
    std::string DebugStr() const;
    // Dump space contents.

    private:
    // truncation
    int Nmax_;

  };

  // sectors

  class TwoBodySectorsNJJJT
    : public BaseSectors<TwoBodySpaceNJJJT>
  {

    public:

    // constructor

    TwoBodySectorsNJJJT() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySectorsNJJJT(
        const TwoBodySpaceNJJJT& space,
        int J0, int T0, int g0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

  };



  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
