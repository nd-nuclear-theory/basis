/****************************************************************
  jjjt_scheme.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <cstddef>
#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>

#include "am/am.h"

#include "jjjt_scheme.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJT scheme
  ////////////////////////////////////////////////////////////////

  TwoBodySubspaceJJJT::TwoBodySubspaceJJJT(
      int J, int T, int g,
      basis::Rank truncation_rank, int truncation_cutoff
    )
    : BaseSubspace{{J,T,g}}
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

  bool TwoBodySubspaceJJJT::ValidLabels() const
  {
    bool valid = true;

    // truncation
    // valid &= ((Nmax()%2)==g());

    return valid;
  }

  std::string TwoBodySubspaceJJJT::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << J()
       << " " << std::setw(width) << T()
       << " " << std::setw(width) << g()
       << " " << "]";

    return os.str();
  }


  std::string TwoBodySubspaceJJJT::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        TwoBodyStateJJJT state(*this,state_index);

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


  TwoBodySpaceJJJT::TwoBodySpaceJJJT(basis::Rank truncation_rank, int truncation_cutoff)
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

                // downshift Nmax to match parity of subspace
                // required to pass label validity tests
                // int Nmax_subspace = Nmax - (Nmax-g)%2;

                TwoBodySubspaceJJJT subspace(J,T,g,truncation_rank,truncation_cutoff);

                if (subspace.size()!=0)
                  PushSubspace(subspace);
              }
  }

  std::string TwoBodySpaceJJJT::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        os
          << " " << "index"
          << " " << std::setw(width) << subspace_index
          << " " << "JTg"
          << " " << std::setw(width) << subspace.J()
          << " " << std::setw(width) << subspace.T()
          << " " << std::setw(width) << subspace.g()
          << " " << "N1max N2max"
          << " " << std::setw(width) << subspace.N1max()
          << " " << std::setw(width) << subspace.N2max()
          << " " << "dim"
          << " " << std::setw(width) << subspace.size()
          << " " << std::endl;
      }

    return os.str();

  }


  TwoBodySectorsJJJT::TwoBodySectorsJJJT(
      const TwoBodySpaceJJJT& space,
      int J0, int T0, int g0,
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
          allowed &= am::AllowedTriangle(ket_subspace.T(),T0,bra_subspace.T());
          allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);

          // push sector
          if (allowed)
            PushSector(bra_subspace_index,ket_subspace_index);
        }
  }

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJT scheme -- subspaced by N
  ////////////////////////////////////////////////////////////////

  TwoBodySubspaceJJJTN::TwoBodySubspaceJJJTN(
      int J, int T, int g, int N,
      basis::Rank truncation_rank, int truncation_cutoff
    )
    : BaseSubspace{{J,T,g,N}}, N_{N}
  {

    // set labels (MODIFICATION for subspacing by N)
    std::tie(N1max_,N2max_) = basis::TwoBodyCutoffs(truncation_rank,truncation_cutoff);

    // validate subspace labels
    assert(ValidLabels());

    // set up indexing
    // iterate over total oscillator quanta -- omit (MODIFICATION for subspacing by N)
    // for (int N = g; N <= Nmax; N +=2)
    // iterate over oscillator (Nj) orbitals for particle 1
    //
    // Constraints:
    //   0 <= N1 <= N1max
    //   0 <= N2 <= N1max
    //   N1+N2=N
    // ==>  N-min(N1max,N) <= N1 <= min(N1max,N)
    int N1_lower = N-std::min(N1max_,N);
    int N1_upper = std::min(N1max_,N);
    for (int N1 = N1_lower; N1 <= N1_upper; ++N1)
      for (HalfInt j1 = HalfInt(1,2); j1 <= N1+HalfInt(1,2); j1 += 1)
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

  bool TwoBodySubspaceJJJTN::ValidLabels() const
  {
    bool valid = true;

    // parity (MODIFICATION for subspacing by N)
    valid &= ((N()%2)==g());

    return valid;
  }

  std::string TwoBodySubspaceJJJTN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << J()
       << " " << std::setw(width) << T()
       << " " << std::setw(width) << g()
       << " " << ";"
       << " " << std::setw(width) << N()
       << " " << "]";

    return os.str();
  }

  std::string TwoBodySubspaceJJJTN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        TwoBodyStateJJJTN state(*this,state_index);

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


  TwoBodySpaceJJJTN::TwoBodySpaceJJJTN(basis::Rank truncation_rank, int truncation_cutoff)
  {
    // save truncation
    std::tie(N1max_,N2max_) = basis::TwoBodyCutoffs(truncation_rank,truncation_cutoff);

    // iterate over J
    for (int J=0; J<=N2max_+1; ++J)
        // iterate over T
        for (int T=0; T<=1; ++T)
            // iterate over g
            for (int g=0; g<=1; ++g)
              // iterate over total oscillator quanta (MODIFICATION for subspacing by N)
              for (int N = g; N <= N2max_; N +=2)
                {
                  TwoBodySubspaceJJJTN subspace(
                      J,T,g,N,truncation_rank,truncation_cutoff
                    );  // (MODIFICATION for subspacing by N)

                  if (subspace.size()!=0)
                    PushSubspace(subspace);
                }
  }

  std::string TwoBodySpaceJJJTN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        os
          << " " << "index"
          << " " << std::setw(width) << subspace_index
          << " " << "JTg"
          << " " << std::setw(width) << subspace.J()
          << " " << std::setw(width) << subspace.T()
          << " " << std::setw(width) << subspace.g()
          << " " << "N"  // (MODIFICATION for subspacing by N)
          << " " << std::setw(width) << subspace.N()  // (MODIFICATION for subspacing by N)
          << " " << "N1max N2max"
          << " " << std::setw(width) << subspace.N1max()
          << " " << std::setw(width) << subspace.N2max()
          << " " << "dim"
          << " " << std::setw(width) << subspace.size()
          << " " << std::endl;
      }

    return os.str();

  }


  TwoBodySectorsJJJTN::TwoBodySectorsJJJTN(
      const TwoBodySpaceJJJTN& space,
      int J0, int T0, int g0,
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
          allowed &= am::AllowedTriangle(ket_subspace.T(),T0,bra_subspace.T());
          allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);

          // push sector
          if (allowed)
            PushSector(bra_subspace_index,ket_subspace_index);
        }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace basis
