/****************************************************************
  lsjt_scheme.h

  Defines relative and two-body state indexing in LSJT coupling scheme.
  Written for use in Moshinsky transformation.

  Indexing is in each case broken into a *subspace* (described by LSJT
  labels and an Nmax trunction) and then indexed states within
  the subspace.

  Relations: Although some of the state types are built from (Nl)
  orbitals, these are stored explicitly as (N,l) values, rather than
  making use of shell_indexing_nl.

  See notes "Moshinsky xform for operators" (2015) for derivation
  of indexing scheme.

  Language: C++11
  (for tuple and map::at)
                                 
  Mark A. Caprio, University of Notre Dame.

  11/13/15 (mac): Created, styled after shell_indexing_nlj.
  11/21/15 (mac): Rename from shell_indexing_lsjt to indexing_lsjt.
    Extracted generic class properties to templates base classes in
    indexing_base.
  11/26/15 (mac): Update iteration schemes.
  6/8/16 (mac): Update for basis package and update conventions.
    - Rename to lsjt_scheme.h.
    - Rename quantum number N to Nr.
    - Rename SpuriousRelative to RelativeCM.
    - Replace RelativeCM *state* iteration constraint lr+S+T~1 in
      constructor with *subspace* constraint Ncm+S+T+g~1 in Validate.
    - Restrict TwoBody states to canonical ordering of orbitals.
    - Replace Print() methods with Str() methods.
  6/22/16 (mac): Make explicit typedefs for label types.

****************************************************************/

#ifndef lsjt_scheme_h
#define lsjt_scheme_h

#include <string>

#include "basis/indexing.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // relative states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  // subspace labels: (L,S,J,T,g)    P=(-)^g
  // state labels within subspace: (Nr)
  //
  // Truncation of the space is by an Nmax.
  //
  // Within a subspace, states are ordered by:
  //   -- increasing Nr (Nr~g)
  // and subject to
  //   -- parity constraint Nr~g
  //
  // That is:
  //   0 <-> Nr=g
  //   1 <-> Nr=g+2
  //   ...
  //
  // Within the space, subspaces are ordered by:
  //   -- increasing L (L=0,1,...,Nmax)
  //   -- increasing S (S=0,1)
  //   -- increasing J
  //   -- [T forced by L+S+T~1]
  //   -- [g forced by g~L]
  // subject to:
  //   -- triangularity of (L,S,J)
  //   -- parity constraint L~g
  //   -- antisymmetry constraint L+S+T~1
  //
  // Assertion: sanity check on nonzero dimension for space
  //
  // Note that ordering of subspaces is therefore lexicographic by (L,S,J).

  // labels

  typedef std::tuple<int,int,int,int,int> RelativeSubspaceLSJTLabels;
  typedef std::tuple<int> RelativeStateLSJTLabels;

  // subspace

  class RelativeSubspaceLSJT
    : public BaseSubspace<RelativeSubspaceLSJTLabels,RelativeStateLSJTLabels>
    {
    
    public:

      // constructor

      RelativeSubspaceLSJT (int L, int S, int J, int T, int g, int Nmax);
      // Set up indexing in Nmax truncation.

      // accessors
 
      int L() const {return std::get<0>(labels_);}
      int S() const {return std::get<1>(labels_);}
      int J() const {return std::get<2>(labels_);}
      int T() const {return std::get<3>(labels_);}
      int g() const {return std::get<4>(labels_);}
      int Nmax() const {return Nmax_;}

    private:

      //validation
      bool ValidLabels() const;

      // truncation
      int Nmax_;

    };

  // state

  class RelativeStateLSJT
    : public BaseState<RelativeSubspaceLSJT>
  {
    
  public:

    // pass-through constructors

  RelativeStateLSJT(const SubspaceType& subspace, int index)
    // Construct state by index.
    : BaseState (subspace, index) {}

  RelativeStateLSJT(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
    : BaseState (subspace, state_labels) {}

    // pass-through accessors
    int L() const {return Subspace().L();}
    int S() const {return Subspace().S();}
    int J() const {return Subspace().J();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}

    // state label accessors
    int Nr() const {return std::get<0>(GetStateLabels());}

  };

  // space

  class RelativeSpaceLSJT
    : public BaseSpace<RelativeSubspaceLSJT>
  {

  public:
    
    // constructor
    RelativeSpaceLSJT(int Nmax);
    // Enumerates all relative LSJT subspaces of given dimension up to a
    // given Nmax cutoff.
    //
    // Arguments:
    //   Nmax (int) : one-body HO truncation on included subspaces

    // diagnostic output
    std::string Str() const;

    // accessors
    int Nmax() const {return Nmax_;}

  private:
    // truncation
    int Nmax_;

  };

  // sectors

  class RelativeSectorsLSJT
    : public BaseSectors<RelativeSpaceLSJT>
  {

  public:

    // constructor

    RelativeSectorsLSJT(const RelativeSpaceLSJT& space);
    // Enumerate all sector pairs ("all-to-all" sector enumeration).
    //
    // Sectors are enumerated in lexicographical order by (bra)(ket).
    // Sectors are included in "both directions", i.e., there is no
    // asumption of hermiticity or attempt to therefore only store
    // "half" the sectors.

    RelativeSectorsLSJT(const RelativeSpaceLSJT& space, int J0, int g0);
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

  };


  ////////////////////////////////////////////////////////////////
  // relative-cm states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  // subspace labels: (Ncm,lcm,L,S,J,T,g)
  // state labels within subspace: (Nr,lr)
  //
  // Truncation of the subspace is by the two-body Nmax.
  //
  // Subspace subject to:
  //   -- triangularity constraint on (L,S,J)
  //   -- antisymmetry constraint lr+S+T~1 on states
  //      implies Ncm+S+T+g~1 for subspace
  //
  // Within a subspace, states are ordered by:
  //   -- increasing N (N~g)
  //   -- lexicographically increasing (Nr,lr)
  // and subject to:
  //   -- triangularity constraint on (lr,lcm,L).
  //   -- parity constraint N=Nr+Ncm~g
  //   -- [antisymmetry constraint lr+S+T~1]
  //
  // A full space (and ordering of subspaces) has not been
  // implemented, as it is not needed for the Moshinsky calculations.
  //
  // However, subspaces are subject to:
  //   -- triangularity constraint on (L,S,J)

  // labels

  typedef std::tuple<int,int,int,int,int,int,int> RelativeCMSubspaceLSJTLabels;
  typedef std::tuple<int,int> RelativeCMStateLSJTLabels;

  //subspace
  
  class RelativeCMSubspaceLSJT
    : public BaseSubspace<RelativeCMSubspaceLSJTLabels,RelativeCMStateLSJTLabels>
  {
    
  public:

    // constructor

    RelativeCMSubspaceLSJT (int Ncm, int lcm, 
				  int L, int S, int J, int T, 
				  int g, int Nmax);

    // accessors

    int Ncm() const {return std::get<0>(labels_);}
    int lcm() const {return std::get<1>(labels_);}
    int L() const {return std::get<2>(labels_);}
    int S() const {return std::get<3>(labels_);}
    int J() const {return std::get<4>(labels_);}
    int T() const {return std::get<5>(labels_);}
    int g() const {return std::get<6>(labels_);}
    int Nmax() const {return Nmax_;}

  private:

    //validation
    bool ValidLabels() const;

    // truncation
    int Nmax_;

  };

  // state

  class RelativeCMStateLSJT
    : public BaseState<RelativeCMSubspaceLSJT>
  {
    
  public:

    // pass-through constructors

  RelativeCMStateLSJT(const SubspaceType& subspace, int index)
    // Construct state by index.
    : BaseState (subspace, index) {}

  RelativeCMStateLSJT(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
    : BaseState (subspace, state_labels) {}

    // pass-through accessors

    int Ncm() const {return Subspace().Ncm();}
    int lcm() const {return Subspace().lcm();}
    int L() const {return Subspace().L();}
    int S() const {return Subspace().S();}
    int J() const {return Subspace().J();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}

    // state label accessors
    int Nr() const {return std::get<0>(GetStateLabels());}
    int lr() const {return std::get<1>(GetStateLabels());}
    int N() const {return  Nr()+Ncm();}

  };

  // space -- not defined

  // sectors -- not defined



  ////////////////////////////////////////////////////////////////
  // two-body states in LSJT scheme
  ////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (L,S,J,T,g)
  //
  // state labels within subspace: (N1,l1,N2,l2)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing L (L=0,1,...,Nmax)
  //    -- increasing S (S=0,1)
  //    -- increasing J
  //    -- increasing T (T=0,1)
  //    -- increasing g (g=0,1)
  // subject to:
  //    -- triangularity of (L,S,J)
  // 
  // Subspaces are pruned to those of nonzero dimension.
  //
  // Note that ordering of subspaces is lexicographic by (L,S,J,g).
  //
  // Truncation of the space is by the two-body Nmax.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing N (N~g)
  //   -- lexicographically increasing (N1,l1)
  //   -- lexicographically increasing (N2,l2)
  // and subject to:
  //   -- triangularity constraint on (l1,l2,L)
  //   -- parity constraint N~g
  //   -- antisymmetry constraint L+S+T~1 if (N1,l1)==(N2,l2)
  //
  //
  // This basis is for *identical* particle states:
  //   -- The labels are subject to the antisymmetry constraint 
  //      (L+S+T~1) if the orbitals are identical.
  //   -- A canonical (lexicographic) ordering constraint is applied to the 
  //      single-particle quantum numbers.  That is, when
  //      enumerating the basis, the states
  //    
  //        |((N1,l1),(N2,l2))...>  and  |((N2,l2),(N1,l1))...>
  //
  //      would be redundant, and only the first (for (N1,l1)<=(N2,l2)) is
  //      retained.
  //
  // Comment: For some applications, it might be more convenient to 
  // an overcomplete basis, in which the states
  //    
  //   |(N1,l1)(N2,l2)...>  and  |(N2,l2)(N1,l1)...>
  //
  // are still counted separately in the basis.  That is, *no*
  // lexicographical ordering constraint (N1,l1)<=(N2,l2) on the two
  // single-particle states is imposed.  The basis is therefore
  // redundant.  This simplifies implementation of double contractions
  // over "particle 1" and "particle 2" indices as matrix
  // multiplication, without the calling code having to worry about
  // "swapping" single particle states within the two-body state and
  // applying the relevant phase (~L+S+T+g+1).
  //
  ////////////////////////////////////////////////////////////////  

  // labels

  typedef std::tuple<int,int,int,int,int> TwoBodySubspaceLSJTLabels;
  typedef std::tuple<int,int,int,int> TwoBodyStateLSJTLabels;

  //subspace

  class TwoBodySubspaceLSJT
    : public BaseSubspace<TwoBodySubspaceLSJTLabels,TwoBodyStateLSJTLabels>
  {
    
  public:

    // constructor

    TwoBodySubspaceLSJT (int L, int S, int J, int T, int g, int Nmax);
    // Set up indexing in Nmax truncation.

    // accessors

    int L() const {return std::get<0>(labels_);}
    int S() const {return std::get<1>(labels_);}
    int J() const {return std::get<2>(labels_);}
    int T() const {return std::get<3>(labels_);}
    int g() const {return std::get<4>(labels_);}
    int Nmax() const {return Nmax_;}

  private:

    //validation
    bool ValidLabels() const;

    // truncation
    int Nmax_;
  };

  // state

  class TwoBodyStateLSJT
    : public BaseState<TwoBodySubspaceLSJT>
  {
    
  public:

    // pass-through constructors

  TwoBodyStateLSJT(const SubspaceType& subspace, int index)
    // Construct state by index.
    : BaseState (subspace, index) {}

  TwoBodyStateLSJT(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
    : BaseState (subspace, state_labels) {}

    // pass-through accessors
    int L() const {return Subspace().L();}
    int S() const {return Subspace().S();}
    int J() const {return Subspace().J();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}

    // state label accessors
    int N1() const {return std::get<0>(GetStateLabels());}
    int l1() const {return std::get<1>(GetStateLabels());}
    int N2() const {return std::get<2>(GetStateLabels());}
    int l2() const {return std::get<3>(GetStateLabels());}

    int N() const {return  N1()+N2();}

  };

  // space

  class TwoBodySpaceLSJT
    : public BaseSpace<TwoBodySubspaceLSJT>
  {
    
  public:

    // constructor

    TwoBodySpaceLSJT(int Nmax);
    // Enumerates all relative LSJT subspaces of given dimension up to a
    // given Nmax cutoff.

    // diagnostic output
    std::string Str() const;

    // accessors
    int Nmax() const {return Nmax_;}

  private:
    // truncation
    int Nmax_;


  };

  // sectors

  class TwoBodySectorsLSJT
    : public BaseSectors<TwoBodySpaceLSJT>
  {

  public:

    // constructor

    TwoBodySectorsLSJT(const TwoBodySpaceLSJT& space);
    // Enumerate all sector pairs ("all-to-all" sector enumeration).
    //
    // Sectors are enumerated in lexicographical order by (bra)(ket).
    // Sectors are included in "both directions", i.e., there is no
    // asumption of hermiticity or attempt to therefore only store
    // "half" the sectors.

    TwoBodySectorsLSJT(const TwoBodySpaceLSJT& space, int J0, int g0);
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
