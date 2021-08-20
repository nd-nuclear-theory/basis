/************************************************************//**
  @file lsjt_scheme.h

  Defines relative and two-body state indexing in LSJT coupling scheme.
  Written for use in Moshinsky transformation.

  Indexing is in each case broken into a *subspace* (described by
  LSJTg labels and an Nmax trunction) and then indexed states within
  the subspace.

  Variant indexing schemes where the subspaces are further broken down
  by the number of oscillator quanta N are then provided, so that the
  Moshinsky transformation implementation can be more easily broken
  down by N blocks.

  See notes "Moshinsky xform for operators" (2015) for derivation
  of indexing scheme.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + 11/13/15 (mac): Created, styled after shell_indexing_nlj.
  + 11/21/15 (mac): Rename from shell_indexing_lsjt to indexing_lsjt.
    Extracted generic class properties to templates base classes in
    indexing_base.
  + 11/26/15 (mac): Update iteration schemes.
  + 6/8/16 (mac): Update for basis package and update conventions.
    - Rename to lsjt_scheme.h.
    - Rename quantum number N to Nr.
    - Rename SpuriousRelative to RelativeCM.
    - Replace RelativeCM *state* iteration constraint lr+S+T~1 in
      constructor with *subspace* constraint Ncm+S+T+g~1 in Validate.
    - Restrict TwoBody states to canonical ordering of orbitals.
    - Replace Print() methods with Str() methods.
  + 6/22/16 (mac): Make explicit typedefs for label types.
  + 6/27/16 (mac):
    - Add Jr_max cutoff on construction of relative basis.
    - Change relative-c.m. scheme from spectator (Nc,lc) to active (Nc,lc).
    - Add fixed-N subspaces in relative-c.m. basis for use with Moshinsky
      transform block structure.
    - Expand basis indexing comments.
    - Implement canonical ordering constraint on sectors.
    - Remove all-to-all sector constructors.
    - Rename Str() to DebugStr().
    - Rename labels on relative basis (e.g., J->Jr).
  + 6/30/16 (mac): Revert labels on relative basis (e.g., Jr->J).
  + 7/3/16 (mac): Add default constructors for RelativeLSJT basis.
  + 7/4/16 (mac): Add fixed-N subspaces in two-body basis for use with
    Moshinsky transform block structure.
  + 7/8/16 (mac): Add default constructors for TwoBodyLSJT basis.
  + 7/9/16 (mac):
    - Add debug strings for RelativeLSJT basis.
    - Add default constructors for remaining subspaces, spaces,
      and sectors.
  + 7/13/16 (mac): Fix relative enumeration.
  + 7/16/16 (mac):
    - Add debug strings for TwoBodyLSJT basis.
    - Move N to least significant subspace label in LSJTN bases.
  + 7/17/16 (mac):
    - Rename NLSJT to LSJTN.
    - Remove unnecessary complication of matching subspace Nmax
      to g.
    - Add one-body (square) truncation on two-body bases.
  + 7/19/16 (mac): Use enum Rank for truncation rank.
  + 7/1/17 (mac): Add n accessor to RelativeStateLSJT.
  + 05/09/19 (pjf): Use std::size_t for indices and sizes, to prevent
    integer overflow.
  + 05/27/19 (pjf): Update to initialize BaseSectors with spaces.
  + 07/03/21 (pjf): Call base class constructor for initializing labels.
  + 07/04/21 (pjf): Pass derived subspace class as template argument to
    BaseSubspace.

****************************************************************/

#ifndef BASIS_LSJT_SCHEME_H_
#define BASIS_LSJT_SCHEME_H_

#include <cstddef>
#include <string>

#include "basis.h"
#include "many_body.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // relative states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (L,S,J,T,g)    P=(-)^g
  //
  //   L (int): orbital angular momentum of *relative* motion (=lr)
  //   S (int): total spin
  //   J (int): total angular momentum of *relative* motion (=Jr)
  //            (i.e., L coupled to S)
  //   T (int): isospin
  //   g (int): grade (=0,1) for the parity P of *relative*
  //            motion (=gr)
  //
  // state labels within subspace: (N)
  //
  //   N (int): oscillator quanta of relative motion (=Nr)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
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
  // Subspaces are asserted to have nonzero dimension (as a sanity
  // check).
  //
  // Note that ordering of subspaces is lexicographic by (L,S,J).
  //
  // Truncation of the space is by the relative Nmax.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing N
  // and subject to:
  //   -- oscillator branching constraints N>=L and N~L (or,
  //      equivalently, parity constraint N~g)
  //
  // This basis is for *identical* particle states, but the
  // antisymmetry constraint is already applied at the level of
  // selecting the subspace labels L+S+T~1.
  //
  ////////////////////////////////////////////////////////////////

  // declarations
  class RelativeSubspaceLSJT;
  class RelativeStateLSJT;
  class RelativeSpaceLSJT;

  // labels

  typedef std::tuple<int,int,int,int,int> RelativeSubspaceLSJTLabels;
  typedef std::tuple<int> RelativeStateLSJTLabels;

  // subspace

  class RelativeSubspaceLSJT
    : public BaseSubspace<RelativeSubspaceLSJT,RelativeSubspaceLSJTLabels,RelativeStateLSJT,RelativeStateLSJTLabels>
    {

      public:

      // constructor

      RelativeSubspaceLSJT() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      RelativeSubspaceLSJT(int L, int S, int J, int T, int g, int Nmax);
      // Set up indexing.

      // accessors

      int L() const {return std::get<0>(labels());}
      int S() const {return std::get<1>(labels());}
      int J() const {return std::get<2>(labels());}
      int T() const {return std::get<3>(labels());}
      int g() const {return std::get<4>(labels());}
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

  class RelativeStateLSJT
    : public BaseState<RelativeSubspaceLSJT>
  {

    public:

    // pass-through constructors

    RelativeStateLSJT(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace,index) {}

    RelativeStateLSJT(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace,state_labels) {}

    // pass-through accessors
    int L() const {return subspace().L();}
    int S() const {return subspace().S();}
    int J() const {return subspace().J();}
    int T() const {return subspace().T();}
    int g() const {return subspace().g();}

    // state label accessors
    int N() const {return std::get<0>(labels());}

    // derived label
    //
    // radial number n can be recovered from N=2*n+L, but it also
    // happens to simply be state index w/in the subspace
    int n() const {return (N()-L())/2;}


  };

  // space

  class RelativeSpaceLSJT
    : public BaseSpace<RelativeSpaceLSJT, RelativeSubspaceLSJT>
  {

    public:

    // constructor

    RelativeSpaceLSJT() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    RelativeSpaceLSJT(int Nmax, int Jmax);
    // Enumerate all relative LSJT subspaces of given dimension up to
    // a given relative oscillator cutoff and relative angular
    // momentum cutoff.
    //
    // The relative angular momentum cutoff is included in recognition
    // of the common practice of truncation by highest partial wave in
    // the representation of relative interations.
    //
    // Arguments:
    //   Nmax (int) : relative oscillator truncation on included subspaces
    //   Jmax (int) : relative angular momentum truncation on included
    //     subspaces (Jmax<=Nmax+1)

    // accessors
    int Nmax() const {return Nmax_;}
    int Jmax() const {return Jmax_;}

    // diagnostic string
    std::string DebugStr() const;

    private:
    // truncation
    int Nmax_, Jmax_;

  };

  // sectors

  class RelativeSectorsLSJT
    : public BaseSectors<RelativeSpaceLSJT>
  {

    public:

    // constructor

    RelativeSectorsLSJT() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    RelativeSectorsLSJT(
        const RelativeSpaceLSJT& space,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate all sector pairs ("all-to-all" sector enumeration).
    //
    // Sectors are enumerated in lexicographical order by (bra)(ket).
    //
    // Note: This all-to-all constructor is implemented primarily for
    // the purpose of providing example code for all-to-all sector
    // enumeration.  It is not clear that there is immediate
    // application of all-to-all enumeration for this particular
    // basis.

    RelativeSectorsLSJT(
        const RelativeSpaceLSJT& space,
        int J0, int T0, int g0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

  };


  ////////////////////////////////////////////////////////////////
  // relative-cm states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (L,S,J,T,g)
  //
  //   L (int): orbital angular momentum
  //   S (int): total spin
  //   J (int): total angular momentum
  //   T (int): isospin
  //   g (int): grade (=0,1) for the parity P
  //
  // state labels within subspace: (Nr,lr,Nc,lc)
  //
  //   Nr (int): oscillator quanta of relative motion
  //   lr (int): orbital angular momentum of relative motion
  //   Nc (int): oscillator quanta of c.m. motion
  //   lc (int): orbital angular momentum of c.m. motion
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
  // Note that ordering of subspaces is lexicographic by (L,S,J,T,g).
  //
  // Truncation of the space is by the two-body Nmax.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing N  (N=Nr+Nc)
  //   -- lexicographically increasing (Nr,lr)
  //   -- lexicographically increasing (Nc,lc)
  // and subject to:
  //   -- triangularity constraint on (lr,lc,L)
  //   -- parity constraint N~g
  //   -- antisymmetry constraint lr+S+T~1 (or, equivalentsly,
  //      Nr+S+T~1)
  //
  // This basis is for *identical* particle states, as enforced by the
  // antisymmetry constraint on Nr.
  //
  ////////////////////////////////////////////////////////////////

  // declarations
  class RelativeCMSubspaceLSJT;
  class RelativeCMStateLSJT;
  class RelativeCMSpaceLSJT;

  // labels

  typedef std::tuple<int,int,int,int,int> RelativeCMSubspaceLSJTLabels;
  typedef std::tuple<int,int,int,int> RelativeCMStateLSJTLabels;

  //subspace

  class RelativeCMSubspaceLSJT
    : public BaseSubspace<RelativeCMSubspaceLSJT,RelativeCMSubspaceLSJTLabels,RelativeCMStateLSJT,RelativeCMStateLSJTLabels>
    {

      public:

      // constructor

      RelativeCMSubspaceLSJT() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      RelativeCMSubspaceLSJT(int L, int S, int J, int T, int g, int Nmax);

      // accessors

      int L() const {return std::get<0>(labels());}
      int S() const {return std::get<1>(labels());}
      int J() const {return std::get<2>(labels());}
      int T() const {return std::get<3>(labels());}
      int g() const {return std::get<4>(labels());}
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

  class RelativeCMStateLSJT
    : public BaseState<RelativeCMSubspaceLSJT>
  {

    public:

    // pass-through constructors

    RelativeCMStateLSJT(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    RelativeCMStateLSJT(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors

    int L() const {return subspace().L();}
    int S() const {return subspace().S();}
    int J() const {return subspace().J();}
    int T() const {return subspace().T();}
    int g() const {return subspace().g();}

    // state label accessors
    int Nr() const {return std::get<0>(labels());}
    int lr() const {return std::get<1>(labels());}
    int Nc() const {return std::get<2>(labels());}
    int lc() const {return std::get<3>(labels());}
    int N() const {return  Nr()+Nc();}

  };

  // space

  class RelativeCMSpaceLSJT
    : public BaseSpace<RelativeCMSpaceLSJT, RelativeCMSubspaceLSJT>
  {

    public:

    // constructor

    RelativeCMSpaceLSJT() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    RelativeCMSpaceLSJT(int Nmax);
    // Enumerate all subspaces up to a given Nmax cutoff.

    // accessors
    int Nmax() const {return Nmax_;}

    // diagnostic string
    std::string DebugStr() const;

    private:
    // truncation
    int Nmax_;

  };

  // sectors

  class RelativeCMSectorsLSJT
    : public BaseSectors<RelativeCMSpaceLSJT>
  {

    public:

    // constructor

    RelativeCMSectorsLSJT() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    RelativeCMSectorsLSJT(
        const RelativeCMSpaceLSJT& space,
        int J0, int T0, int g0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

  };


  ////////////////////////////////////////////////////////////////
  // relative-cm states in LSJT scheme -- subspaced by N
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (L,S,J,T,g,N)  (MODIFICATION for subspacing by N)
  //
  // state labels within subspace: (Nr,lr,Nc,lc)
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
  //    -- increasing N (N=0,1,...,Nmax)  (MODIFICATION for subspacing by N)
  // subject to:
  //    -- triangularity of (L,S,J)
  //    -- parity constraint N~g  (MODIFICATION for subspacing by N)
  //
  // Subspaces are pruned to those of nonzero dimension.
  //
  // Note that ordering of subspaces is lexicographic by (L,S,J,T,g,N).  (MODIFICATION for subspacing by N)
  //
  // Truncation of the space is by the two-body Nmax.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- [increasing N  (N=Nr+Nc)]  (MODIFICATION for subspacing by N)
  //   -- lexicographically increasing (Nr,lr)
  //   -- lexicographically increasing (Nc,lc)
  // and subject to:
  //   -- triangularity constraint on (lr,lc,L)
  //   -- [parity constraint N~g]  (MODIFICATION for subspacing by N)
  //   -- antisymmetry constraint lr+S+T~1 (or, equivalentsly,
  //      Nr+S+T~1)
  //
  // This basis is for *identical* particle states (see discussion
  // above for non-N version).
  //
  ////////////////////////////////////////////////////////////////

  // Modification for subspacing by N is by lexical replacement LSJT
  // -> LSJTN plus specific mods as flagged by MODIFICATION comments
  // in code.

  // declarations
  class RelativeCMSubspaceLSJTN;
  class RelativeCMStateLSJTN;
  class RelativeCMSpaceLSJTN;

  // labels

  typedef std::tuple<int,int,int,int,int,int> RelativeCMSubspaceLSJTNLabels;  // (MODIFICATION for subspacing by N)
  typedef std::tuple<int,int,int,int> RelativeCMStateLSJTNLabels;

  //subspace

  class RelativeCMSubspaceLSJTN
    : public BaseSubspace<RelativeCMSubspaceLSJTN,RelativeCMSubspaceLSJTNLabels,RelativeCMStateLSJTN,RelativeCMStateLSJTNLabels>
    {

      public:

      // constructor

      RelativeCMSubspaceLSJTN() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      RelativeCMSubspaceLSJTN(int L, int S, int J, int T, int g, int N);  // (MODIFICATION for subspacing by N)

      // accessors

      int L() const {return std::get<0>(labels());}
      int S() const {return std::get<1>(labels());}
      int J() const {return std::get<2>(labels());}
      int T() const {return std::get<3>(labels());}
      int g() const {return std::get<4>(labels());}
      int N() const {return std::get<5>(labels());}  // (MODIFICATION for subspacing by N)

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      //validation
      bool ValidLabels() const;

      private:
      // truncation
      int N_;  // (MODIFICATION for subspacing by N)


    };

  // state

  class RelativeCMStateLSJTN
    : public BaseState<RelativeCMSubspaceLSJTN>
  {

    public:

    // pass-through constructors

    RelativeCMStateLSJTN(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    RelativeCMStateLSJTN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors

    int L() const {return subspace().L();}
    int S() const {return subspace().S();}
    int J() const {return subspace().J();}
    int T() const {return subspace().T();}
    int g() const {return subspace().g();}

    // state label accessors
    int Nr() const {return std::get<0>(labels());}
    int lr() const {return std::get<1>(labels());}
    int Nc() const {return std::get<2>(labels());}
    int lc() const {return std::get<3>(labels());}
    int N() const {return  Nr()+Nc();}

  };

  // space

  class RelativeCMSpaceLSJTN
    : public BaseSpace<RelativeCMSpaceLSJTN, RelativeCMSubspaceLSJTN>
  {

    public:

    // constructor

    RelativeCMSpaceLSJTN() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    RelativeCMSpaceLSJTN(int Nmax);
    // Enumerate all subspaces up to a given Nmax cutoff.

    // accessors
    int Nmax() const {return Nmax_;}

    // diagnostic string
    std::string DebugStr() const;

    private:
    // truncation
    int Nmax_;

  };

  // sectors

  class RelativeCMSectorsLSJTN
    : public BaseSectors<RelativeCMSpaceLSJTN>
  {

    public:

    // constructor

    RelativeCMSectorsLSJTN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)


    RelativeCMSectorsLSJTN(
        const RelativeCMSpaceLSJTN& space,
        int J0, int T0, int g0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

  };


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
  // Truncation of the space is by the one-body N1max or two-body
  // N2max.
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

  // declarations
  class TwoBodySubspaceLSJT;
  class TwoBodyStateLSJT;
  class TwoBodySpaceLSJT;

  // labels

  typedef std::tuple<int,int,int,int,int> TwoBodySubspaceLSJTLabels;
  typedef std::tuple<int,int,int,int> TwoBodyStateLSJTLabels;

  //subspace

  class TwoBodySubspaceLSJT
    : public BaseSubspace<TwoBodySubspaceLSJT,TwoBodySubspaceLSJTLabels,TwoBodyStateLSJT,TwoBodyStateLSJTLabels>
    {

      public:

      // constructor

      TwoBodySubspaceLSJT() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      TwoBodySubspaceLSJT(
          int L, int S, int J, int T, int g,
          basis::Rank truncation_rank, int truncation_cutoff
        );
      // Set up indexing.

      // accessors

      int L() const {return std::get<0>(labels());}
      int S() const {return std::get<1>(labels());}
      int J() const {return std::get<2>(labels());}
      int T() const {return std::get<3>(labels());}
      int g() const {return std::get<4>(labels());}
      int N1max() const {return N1max_;}
      int N2max() const {return N2max_;}

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      //validation
      bool ValidLabels() const;

      // truncation
      int N1max_, N2max_;
    };

  // state

  class TwoBodyStateLSJT
    : public BaseState<TwoBodySubspaceLSJT>
  {

    public:

    // pass-through constructors

    TwoBodyStateLSJT(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    TwoBodyStateLSJT(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors
    int L() const {return subspace().L();}
    int S() const {return subspace().S();}
    int J() const {return subspace().J();}
    int T() const {return subspace().T();}
    int g() const {return subspace().g();}

    // state label accessors
    int N1() const {return std::get<0>(labels());}
    int l1() const {return std::get<1>(labels());}
    int N2() const {return std::get<2>(labels());}
    int l2() const {return std::get<3>(labels());}

    int N() const {return N1()+N2();}

  };

  // space

  class TwoBodySpaceLSJT
    : public BaseSpace<TwoBodySpaceLSJT, TwoBodySubspaceLSJT>
  {

    public:

    // constructor

    TwoBodySpaceLSJT() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySpaceLSJT(basis::Rank truncation_rank, int truncation_cutoff);
    // Enumerate all subspaces up to a given Nmax cutoff.

    // accessors
    int N1max() const {return N1max_;}
    int N2max() const {return N2max_;}

    // diagnostic string
    std::string DebugStr() const;

    private:
    // truncation
    int N1max_, N2max_;

  };

  // sectors

  class TwoBodySectorsLSJT
    : public BaseSectors<TwoBodySpaceLSJT>
  {

    public:

    // constructor

    TwoBodySectorsLSJT() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySectorsLSJT(
        const TwoBodySpaceLSJT& space,
        int J0, int T0, int g0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

  };

  ////////////////////////////////////////////////////////////////
  // two-body states in LSJT scheme -- subspaced by N
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (L,S,J,T,g,N)  (MODIFICATION for subspacing by N)
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
  //    -- increasing N (N=0,1,...,Nmax)  (MODIFICATION for subspacing by N)
  // subject to:
  //    -- triangularity of (L,S,J)
  //    -- parity constraint N~g  (MODIFICATION for subspacing by N)
  //
  // Subspaces are pruned to those of nonzero dimension.
  //
  // Note that ordering of subspaces is lexicographic by (L,S,J,T,g,N).  (MODIFICATION for subspacing by N)
  //
  // Truncation of the space is by the one-body N1max or two-body
  // N2max.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- [increasing N  (N=N1+N2)]  (MODIFICATION for subspacing by N)
  //   -- lexicographically increasing (N1,l1)
  //   -- lexicographically increasing (N2,l2)
  // and subject to:
  //   -- triangularity constraint on (l1,l2,L)
  //   -- [parity constraint N~g]  (MODIFICATION for subspacing by N)
  //   -- antisymmetry constraint L+S+T~1 if (N1,l1)==(N2,l2)
  //
  //
  // This basis is for *identical* particle states (see discussion
  // above for non-N version).
  //
  ////////////////////////////////////////////////////////////////

  // Modification for subspacing by N is by lexical replacement LSJT
  // -> LSJTN plus specific mods as flagged by MODIFICATION comments
  // in code.

  // declarations
  class TwoBodySubspaceLSJTN;
  class TwoBodyStateLSJTN;
  class TwoBodySpaceLSJTN;

  // labels

  typedef std::tuple<int,int,int,int,int,int> TwoBodySubspaceLSJTNLabels;  // (MODIFICATION for subspacing by N)
  typedef std::tuple<int,int,int,int> TwoBodyStateLSJTNLabels;

  //subspace

  class TwoBodySubspaceLSJTN
    : public BaseSubspace<TwoBodySubspaceLSJTN,TwoBodySubspaceLSJTNLabels,TwoBodyStateLSJTN,TwoBodyStateLSJTNLabels>
    {

      public:

      // constructor

      TwoBodySubspaceLSJTN() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      TwoBodySubspaceLSJTN(
          int L, int S, int J, int T, int g, int N,  // (MODIFICATION for subspacing by N)
          basis::Rank truncation_rank, int truncation_cutoff
        );
      // Set up indexing.

      // accessors

      int L() const {return std::get<0>(labels());}
      int S() const {return std::get<1>(labels());}
      int J() const {return std::get<2>(labels());}
      int T() const {return std::get<3>(labels());}
      int g() const {return std::get<4>(labels());}
      int N() const {return std::get<5>(labels());}  // (MODIFICATION for subspacing by N)
      int N1max() const {return N1max_;}
      int N2max() const {return N2max_;}

      // diagnostic strings
      std::string LabelStr() const;
      // Provide string representation of subspace labels.
      std::string DebugStr() const;
      // Dump subspace contents.

      private:

      //validation
      bool ValidLabels() const;

      // truncation
      int N1max_, N2max_;
      int N_;  // (MODIFICATION for subspacing by N)
    };

  // state

  class TwoBodyStateLSJTN
    : public BaseState<TwoBodySubspaceLSJTN>
  {

    public:

    // pass-through constructors

    TwoBodyStateLSJTN(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    TwoBodyStateLSJTN(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors
    int L() const {return subspace().L();}
    int S() const {return subspace().S();}
    int J() const {return subspace().J();}
    int T() const {return subspace().T();}
    int g() const {return subspace().g();}

    // state label accessors
    int N1() const {return std::get<0>(labels());}
    int l1() const {return std::get<1>(labels());}
    int N2() const {return std::get<2>(labels());}
    int l2() const {return std::get<3>(labels());}

    int N() const {return N1()+N2();}

  };

  // space

  class TwoBodySpaceLSJTN
    : public BaseSpace<TwoBodySpaceLSJTN, TwoBodySubspaceLSJTN>
  {

    public:

    // constructor

    TwoBodySpaceLSJTN() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySpaceLSJTN(basis::Rank truncation_rank, int truncation_cutoff);
    // Enumerate all subspaces up to a given Nmax cutoff.

    // accessors
    int N1max() const {return N1max_;}
    int N2max() const {return N2max_;}

    // diagnostic string
    std::string DebugStr() const;

    private:
    // truncation
    int N1max_, N2max_;

  };

  // sectors

  class TwoBodySectorsLSJTN
    : public BaseSectors<TwoBodySpaceLSJTN>
  {

    public:

    // constructor

    TwoBodySectorsLSJTN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySectorsLSJTN(
        const TwoBodySpaceLSJTN& space,
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
