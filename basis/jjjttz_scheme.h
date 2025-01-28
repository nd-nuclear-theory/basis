/************************************************************//**
  @file jjjttz_scheme.h

  Defines two-body state indexing in jjJTTz coupling scheme.

  Language: C++11

  Mark A. Caprio
  University of Notre Dame

  + Created by zz from jjjt_scheme ~12/17/23.
  + 01/28/25 (mac): Add accessors for operator properties.

****************************************************************/

#ifndef BASIS_JJJTTz_SCHEME_H_
#define BASIS_JJJTTz_SCHEME_H_

#include <cstddef>
#include <string>
#include <tuple>

#include "am/halfint.h"

#include "basis/basis.h"
#include "basis/many_body.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJTTz scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (J,T,Tz,g)    P=(-)^g
  //
  //   J (int): total angular momentum
  //   T (int): isospin
  //   g (int): grade (=0,1) for the parity P
  //   Tz (int): isospin z projection (+ for protons, - for neutrons)
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
  //    -- increasing Tz (Tz=-T,-T+1,...,T)
  //
  // subject to:
  //    -- [implicit constraint J=Nmax+1 is excluded for T=1,
  //       but this is simply enforced by pruning to subspaces of
  //       nonzero dimension]
  //
  // Subspaces are pruned to those of nonzero dimension.
  //
  // Note that ordering of subspaces is lexicographic by (J,T,Tz,g).
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

  // declarations
  class TwoBodySubspaceJJJTTz;
  class TwoBodyStateJJJTTz;
  class TwoBodySpaceJJJTTz;

  // labels

  typedef std::tuple<int,int,int,int> TwoBodySubspaceJJJTTzLabels;
  typedef std::tuple<int,HalfInt,int,HalfInt> TwoBodyStateJJJTTzLabels;

  // subspace

  class TwoBodySubspaceJJJTTz
    : public BaseSubspace<TwoBodySubspaceJJJTTz,TwoBodySubspaceJJJTTzLabels,TwoBodyStateJJJTTz,TwoBodyStateJJJTTzLabels>
    {

      public:

      // constructor

      TwoBodySubspaceJJJTTz() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      TwoBodySubspaceJJJTTz(
          int J, int T, int g, int Tz,
          basis::Rank truncation_rank, int truncation_cutoff
        );
      // Set up indexing.

      // accessors

      int J() const {return std::get<0>(labels());}
      int T() const {return std::get<1>(labels());}
      int g() const {return std::get<2>(labels());}
      int Tz() const {return std::get<3>(labels());}
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

  class TwoBodyStateJJJTTz
    : public BaseState<TwoBodySubspaceJJJTTz>
  {

    public:

    // pass-through constructors

    TwoBodyStateJJJTTz(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
      : BaseState (subspace, index) {}

    TwoBodyStateJJJTTz(const SubspaceType& subspace, const StateLabelsType& state_labels)
      // Construct state by reverse lookup on labels.
      : BaseState (subspace, state_labels) {}

    // pass-through accessors
    int J() const {return subspace().J();}
    int T() const {return subspace().T();}
    int g() const {return subspace().g();}
    int Tz() const {return subspace().Tz();}

    // state label accessors
    int N1() const {return std::get<0>(labels());}
    HalfInt j1() const {return std::get<1>(labels());}
    int N2() const {return std::get<2>(labels());}
    HalfInt j2() const {return std::get<3>(labels());}

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

  class TwoBodySpaceJJJTTz
    : public BaseSpace<TwoBodySpaceJJJTTz, TwoBodySubspaceJJJTTz>
  {

    public:

    // constructor

    TwoBodySpaceJJJTTz() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySpaceJJJTTz(basis::Rank truncation_rank, int truncation_cutoff);
    // Enumerate subspaces.

    // accessors
    int N1max() const {return N1max_;}
    int N2max() const {return N2max_;}

    // diagnostic string
    std::string DebugStr() const;
    // Dump space contents.

    private:
    // truncation
    int N1max_, N2max_;

  };

  // sectors

  class TwoBodySectorsJJJTTz
    : public BaseSectors<TwoBodySpaceJJJTTz>
  {

    public:

    // constructor

    TwoBodySectorsJJJTTz() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    TwoBodySectorsJJJTTz(
        const TwoBodySpaceJJJTTz& space,
        int J0, int g0, int Tz0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).

    // accessors
    int J0() const {return J0_;};
    int g0() const {return g0_;};
    int Tz0() const {return Tz0_;};

    private:

    // operator properties
    int J0_, g0_, Tz0_;
  };
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
