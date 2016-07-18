/****************************************************************
  jjjt_scheme.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>

#include "am/am.h"

#include "jjjt_scheme.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJT scheme
  ////////////////////////////////////////////////////////////////

  TwoBodySubspaceJJJT::TwoBodySubspaceJJJT(int J, int T, int g, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(J,T,g);
    Nmax_ = Nmax;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over total oscillator quanta
    for (int N = g; N <= Nmax; N +=2)
      // iterate over oscillator (Nj) orbitals for particle 1
      for (int N1 = 0; N1 <= N; ++N1)
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

    const int lw = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        TwoBodyStateJJJT state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(lw) << state_index
	  << " " << "N1 l1 j1 N2 l2 j2"
	  << " " << std::setw(lw) << state.N1()
	  << " " << std::setw(lw) << state.l1()
	  << " " << std::setw(lw) << state.j1()
	  << " " << std::setw(lw) << state.N2()
	  << " " << std::setw(lw) << state.l2()
	  << " " << std::setw(lw) << state.j2()
          << std::endl;
      }

    return os.str();

  }


  TwoBodySpaceJJJT::TwoBodySpaceJJJT(int Nmax)
  {

    // save Nmax
    Nmax_ = Nmax;

    // iterate over J
    for (int J=0; J<=Nmax+1; ++J)
	// iterate over T
	for (int T=0; T<=1; ++T)
	    // iterate over g
	    for (int g=0; g<=1; ++g)
	      {
		
		// downshift Nmax to match parity of subspace
		// required to pass label validity tests
		// int Nmax_subspace = Nmax - (Nmax-g)%2;
		    
		TwoBodySubspaceJJJT subspace(J,T,g,Nmax);

		if (subspace.size()!=0)
		  PushSubspace(subspace);
	      }
  }

  std::string TwoBodySpaceJJJT::DebugStr() const
  {

    std::ostringstream os;

    const int lw = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
	os
	  << " " << "index"
	  << " " << std::setw(lw) << subspace_index
	  << " " << "JTg"
	  << " " << std::setw(lw) << subspace.J() 
	  << " " << std::setw(lw) << subspace.T() 
	  << " " << std::setw(lw) << subspace.g()
	  << " " << "Nmax"
	  << " " << std::setw(lw) << subspace.Nmax()
	  << " " << "dim"
	  << " " << std::setw(lw) << subspace.size()
	  << " " << std::endl;
      }

    return os.str();

  }


  TwoBodySectorsJJJT::TwoBodySectorsJJJT(
      const TwoBodySpaceJJJT& space,
      int J0, int T0, int g0,
      basis::SectorDirection sector_direction
    )
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
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
            PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
        }
  }

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJT scheme -- subspaced by N
  ////////////////////////////////////////////////////////////////

  TwoBodySubspaceJJJTN::TwoBodySubspaceJJJTN(int J, int T, int g, int N)
  {

    // set values (MODIFICATION for subspacing by N)
    labels_ = SubspaceLabelsType(J,T,g,N);
    N_ = N;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over total oscillator quanta -- omit (MODIFICATION for subspacing by N)
    // for (int N = g; N <= Nmax; N +=2)
    // iterate over oscillator (Nj) orbitals for particle 1
    for (int N1 = 0; N1 <= N; ++N1)
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

    const int lw = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        TwoBodyStateJJJTN state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(lw) << state_index
	  << " " << "N1 l1 j1 N2 l2 j2"
	  << " " << std::setw(lw) << state.N1()
	  << " " << std::setw(lw) << state.l1()
	  << " " << std::setw(lw) << state.j1()
	  << " " << std::setw(lw) << state.N2()
	  << " " << std::setw(lw) << state.l2()
	  << " " << std::setw(lw) << state.j2()
          << std::endl;
      }

    return os.str();

  }


  TwoBodySpaceJJJTN::TwoBodySpaceJJJTN(int Nmax)
  {

    // save Nmax
    Nmax_ = Nmax;

    // iterate over J
    for (int J=0; J<=Nmax+1; ++J)
	// iterate over T
	for (int T=0; T<=1; ++T)
	    // iterate over g
	    for (int g=0; g<=1; ++g)
              // iterate over total oscillator quanta (MODIFICATION for subspacing by N)
              for (int N = g; N <= Nmax; N +=2)
                {	
                  TwoBodySubspaceJJJTN subspace(J,T,g,N);  // (MODIFICATION for subspacing by N)

                  if (subspace.size()!=0)
                    PushSubspace(subspace);
                }
  }

  std::string TwoBodySpaceJJJTN::DebugStr() const
  {

    std::ostringstream os;

    const int lw = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
	os
	  << " " << "index"
	  << " " << std::setw(lw) << subspace_index
	  << " " << "JTg"
	  << " " << std::setw(lw) << subspace.J() 
	  << " " << std::setw(lw) << subspace.T() 
	  << " " << std::setw(lw) << subspace.g()
	  << " " << "N"  // (MODIFICATION for subspacing by N)
	  << " " << std::setw(lw) << subspace.N()  // (MODIFICATION for subspacing by N)
	  << " " << "dim"
	  << " " << std::setw(lw) << subspace.size()
	  << " " << std::endl;
      }

    return os.str();

  }


  TwoBodySectorsJJJTN::TwoBodySectorsJJJTN(
      const TwoBodySpaceJJJTN& space,
      int J0, int T0, int g0,
      basis::SectorDirection sector_direction
    )
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
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
            PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
        }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
