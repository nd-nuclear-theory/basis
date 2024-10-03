/****************************************************************
  jjjt_scheme.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <cassert>
#include <cstddef>
#include <iomanip>  // for debugging output
#include <sstream>
#include <utility>

#include "am/am.h"

#include "basis/jjjttz_scheme.h"
#include "basis/basis.h"
#include "basis/many_body.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJT scheme
  ////////////////////////////////////////////////////////////////

  TwoBodySubspaceJJJTTz::TwoBodySubspaceJJJTTz(
      int J, int T, int g, int Tz,
      basis::Rank truncation_rank, int truncation_cutoff
    )
    : BaseSubspace{{J,T,g,Tz}}
  {

    // set labels
    std::tie(N1max_,N2max_) = basis::TwoBodyCutoffs(truncation_rank,truncation_cutoff);

    // validate subspace labels
    assert(ValidLabels());

    // set up indexing
    // iterate over total oscillator quanta
    for (int N = g; N <= N2max_; N +=2)
      // iterate over oscillator (Nj) orbitals for particle 1
      //
      // Constraints:
      //   0 <= N1 <= N1max
      //   0 <= N2 <= N1max
      //   N1+N2=N
      // ==>  N-min(N1max,N) <= N1 <= min(N1max,N)
      {
        int N1_lower = N-std::min(N1max_,N);
        int N1_upper = std::min(N1max_,N);
        for (int N1 = N1_lower; N1 <= N1_upper; ++N1)
          for (HalfInt j1 = HalfInt(1,2); j1 <= N1+HalfInt(1,2); j1 +=1)
            {
              // iterate over oscillator (Nj) orbitals for particle 2
              // subject to given total N
              int N2 = N - N1;

              for (HalfInt j2 = HalfInt(1,2); j2 <= N2+HalfInt(1,2); j2 +=1)
                {

                  // impose canonical ordering on single-particle states
                  if (!( std::make_pair(N1,j1) <= std::make_pair(N2,j2) ))
                    continue;

                  // impose triangularity
                  if (!(am::AllowedTriangle(j1,j2,J)))
                    continue;

                  // impose antisymmetry
                  if ((N1==N2)&&(j1==j2)&&(!((J+T)%2==1)))
                    continue;

                  // keep surviving states
                  PushStateLabels(StateLabelsType(N1,j1,N2,j2));
                }
            }
      }
  }

  bool TwoBodySubspaceJJJTTz::ValidLabels() const
  {
    bool valid = true;

    // truncation
    // valid &= ((Nmax()%2)==g());

    return valid;
  }

  std::string TwoBodySubspaceJJJTTz::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << J()
       << " " << std::setw(width) << T()
       << " " << std::setw(width) << g()
       << " " << std::setw(width) << Tz()
       << " " << "]";

    return os.str();
  }


  std::string TwoBodySubspaceJJJTTz::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        TwoBodyStateJJJTTz state(*this,state_index);

        os
          << " " << "index"
          << " " << std::setw(width) << state_index
          << " " << "N1 l1 j1 N2 l2 j2"
          << " " << std::setw(width) << state.N1()
          << " " << std::setw(width) << state.l1()
          << " " << std::setw(width) << state.j1()
          << " " << std::setw(width) << state.N2()
          << " " << std::setw(width) << state.l2()
          << " " << std::setw(width) << state.j2()
          << std::endl;
      }

    return os.str();

  }


  TwoBodySpaceJJJTTz::TwoBodySpaceJJJTTz(basis::Rank truncation_rank, int truncation_cutoff)
  {

    // save truncation
    std::tie(N1max_,N2max_) = basis::TwoBodyCutoffs(truncation_rank,truncation_cutoff);

    // iterate over J
    for (int J=0; J<=N2max_+1; ++J)
        // iterate over T
        for (int T=0; T<=1; ++T)
            // iterate over g
            for (int g=0; g<=1; ++g)
              {
                  for (int Tz=-T; Tz<=T; ++Tz)
                    {
                        // downshift Nmax to match parity of subspace
                        // required to pass label validity tests
                        // int Nmax_subspace = Nmax - (Nmax-g)%2;

                        TwoBodySubspaceJJJTTz subspace(J,T,g,Tz,truncation_rank,truncation_cutoff);

                        if (subspace.size()!=0)
                          PushSubspace(subspace);
                    }
              }
  }

  std::string TwoBodySpaceJJJTTz::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        os
          << " " << "index"
          << " " << std::setw(width) << subspace_index
          << " " << "JTgTz"
          << " " << std::setw(width) << subspace.J()
          << " " << std::setw(width) << subspace.T()
          << " " << std::setw(width) << subspace.g()
          << " " << std::setw(width) << subspace.Tz()
          << " " << "N1max N2max"
          << " " << std::setw(width) << subspace.N1max()
          << " " << std::setw(width) << subspace.N2max()
          << " " << "dim"
          << " " << std::setw(width) << subspace.size()
          << " " << std::endl;
      }

    return os.str();

  }


  TwoBodySectorsJJJTTz::TwoBodySectorsJJJTTz(
      const TwoBodySpaceJJJTTz& space,
      int J0, int g0, int Tz0,
      basis::SectorDirection sector_direction
    )
    : BaseSectors(space)
  {
    for (std::size_t bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (std::size_t ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {

          // enforce canonical ordering
          if (
              (sector_direction == basis::SectorDirection::kCanonical)
              && !(bra_subspace_index<=ket_subspace_index)
            )
            continue;

          // retrieve subspaces
          const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

          // verify angular momentum, isosopin, and parity selection rules
          bool allowed = true;
          allowed &= am::AllowedTriangle(ket_subspace.J(),J0,bra_subspace.J());
          // allowed &= am::AllowedTriangle(ket_subspace.T(),T0,bra_subspace.T());
          allowed &= (ket_subspace.Tz()+Tz0==bra_subspace.Tz());
          allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);
          // allowed &= Tz0>=-T0 && Tz0<=T0;

          // push sector
          if (allowed)
            PushSector(bra_subspace_index,ket_subspace_index);
        }
  }
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace basis
