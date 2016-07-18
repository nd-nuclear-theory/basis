/****************************************************************
  lsjt_scheme.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>

#include "am/am.h"

#include "lsjt_scheme.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // relative states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  RelativeSubspaceLSJT::RelativeSubspaceLSJT(int L, int S, int J, int T, int g, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(L,S,J,T,g);
    Nmax_ = Nmax;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing

    // manual version w/o lookup
    //	  dimension_ = (Nmax() - g()) / 2 + 1;

    // iterate over total oscillator quanta
    for (int N = L; N <= Nmax; N +=2)
      PushStateLabels(StateLabelsType(N));

  }

  bool RelativeSubspaceLSJT::ValidLabels() const
  {

    bool valid = true;

    // triangularity
    valid &= am::AllowedTriangle(L(),S(),J());
    // parity
    valid &= (g() == (L()%2));
    // antisymmetry
    valid &= ((L()+S()+T())%2 == 1);
    // truncation
    // valid &= ((Nmax()%2)==g());

    return valid;
  }

  std::string RelativeSubspaceLSJT::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << L() 
       << " " << std::setw(width) << S() 
       << " " << std::setw(width) << J() 
       << " " << std::setw(width) << T() 
       << " " << std::setw(width) << g()
       << " " << "]";

    return os.str();
  }

  std::string RelativeSubspaceLSJT::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        RelativeStateLSJT state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(width) << state_index
	  << " " << "N"
	  << " " << std::setw(width) << state.N() 
          << std::endl;
      }

    return os.str();

  }

  RelativeSpaceLSJT::RelativeSpaceLSJT(int Nmax, int Jmax)
    : Nmax_(Nmax), Jmax_(Jmax)
  {

    // iterate over L
    for (int L=0; L<=Nmax; ++L)
      {
        // set g
	int g = L%2;

	// iterate over S
	for (int S=0; S<=1; ++S)
	  {
            // set T
	    int T = (L+S+1)%2;

	    // iterate over J
            int J_limit = std::min(L+S,Jmax);
	    for (int J=abs(L-S); J<=J_limit; ++J)
	      {
		// downshift Nmax to match parity of subspace
		// required to pass label validity tests
		// int Nmax_subspace = Nmax - (Nmax-g)%2;

		RelativeSubspaceLSJT subspace(L,S,J,T,g,Nmax);
		assert(subspace.size()!=0);
		PushSubspace(subspace);
	      }
	  }
      }
  }

  std::string RelativeSpaceLSJT::DebugStr() const
  {
    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
	os
	  << " " << "index"
	  << " " << std::setw(width) << subspace_index
	  << " " << " (L,S,J,T,g) "
	  << " " << std::setw(width) << subspace.L() 
	  << " " << std::setw(width) << subspace.S() 
	  << " " << std::setw(width) << subspace.J() 
	  << " " << std::setw(width) << subspace.T() 
	  << " " << std::setw(width) << subspace.g()
	  << " " << "Nmax"
	  << " " << std::setw(width) << subspace.Nmax()
	  << " " << "dim"
	  << " " << std::setw(width) << subspace.size()
	  << " " << std::endl;
      }

    return os.str();
  }

  RelativeSectorsLSJT::RelativeSectorsLSJT(
      const RelativeSpaceLSJT& space,
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

          // push sector
          PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
        }
  }

  RelativeSectorsLSJT::RelativeSectorsLSJT(
      const RelativeSpaceLSJT& space,
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
  // relative-cm states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  RelativeCMSubspaceLSJT::RelativeCMSubspaceLSJT(int L, int S, int J, int T, int g, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(L,S,J,T,g);
    Nmax_ = Nmax;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over total oscillator quanta
    for (int N = g; N <= Nmax; N +=2)
      // iterate over relative oscillator (Nr,lr) orbitals
      for (int Nr = 0; Nr <= N; ++Nr)
	for (int lr = Nr%2; lr <= Nr; lr +=2) 
	  {
	    // iterate over c.m. oscillator (Nc,lc) orbitals
	    // subject to given total N
	    int Nc = N - Nr;
            
	    for (int lc = Nc%2; lc <= Nc; lc +=2) 
	      {

		// impose triangularity
		if (!(am::AllowedTriangle(lr,lc,L)))
		  continue;

		// impose antisymmetry
		if (!((lr+S+T)%2==1))
		  continue;

		// keep surviving states
		PushStateLabels(StateLabelsType(Nr,lr,Nc,lc)); 
	      }
	  }
  }

  bool RelativeCMSubspaceLSJT::ValidLabels() const
  {
    bool valid = true;
      
    // triangularity
    valid &= am::AllowedTriangle(L(),S(),J());

    // truncation
    // valid &= ((Nmax()%2)==g());

    return valid;
  }

  std::string RelativeCMSubspaceLSJT::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << L() 
       << " " << std::setw(width) << S() 
       << " " << std::setw(width) << J() 
       << " " << std::setw(width) << T() 
       << " " << std::setw(width) << g()
       << " " << "]";

    return os.str();
  }


  std::string RelativeCMSubspaceLSJT::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        RelativeCMStateLSJT state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(width) << state_index
	  << " " << "Nr lr Nc lc"
	  << " " << std::setw(width) << state.Nr() 
	  << " " << std::setw(width) << state.lr() 
	  << " " << std::setw(width) << state.Nc() 
	  << " " << std::setw(width) << state.lc()
          << std::endl;
      }

    return os.str();

  }

  RelativeCMSpaceLSJT::RelativeCMSpaceLSJT(int Nmax)
  {
    // save Nmax
    Nmax_ = Nmax;

    // iterate over L
    for (int L=0; L<=Nmax; ++L)
	// iterate over S
	for (int S=0; S<=1; ++S)
	    // iterate over J
	    // imposing triangularity (LSJ)
	    for (int J=abs(L-S); J<=L+S; ++J)
		// iterate over T
		for (int T=0; T<=1; ++T)
		    // iterate over g
		    for (int g=0; g<=1; ++g)
		      {
			
			// downshift Nmax to match parity of subspace
			// required to pass label validity tests
			// int Nmax_subspace = Nmax - (Nmax-g)%2;
		    
			RelativeCMSubspaceLSJT subspace(L,S,J,T,g,Nmax);
			if (subspace.size()!=0)
			  PushSubspace(subspace);
		      }
  }

  std::string RelativeCMSpaceLSJT::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
	os
	  << " " << "index"
	  << " " << std::setw(width) << subspace_index
	  << " " << "LSJTg"
	  << " " << std::setw(width) << subspace.L() 
	  << " " << std::setw(width) << subspace.S() 
	  << " " << std::setw(width) << subspace.J() 
	  << " " << std::setw(width) << subspace.T() 
	  << " " << std::setw(width) << subspace.g()
	  << " " << "Nmax"
	  << " " << std::setw(width) << subspace.Nmax()
	  << " " << "dim"
	  << " " << std::setw(width) << subspace.size()
	  << " " << std::endl;
      }

    return os.str();

  }

  RelativeCMSectorsLSJT::RelativeCMSectorsLSJT(
      const RelativeCMSpaceLSJT& space,
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
  // relative-cm states in LSJT scheme -- subspaced by N
  ////////////////////////////////////////////////////////////////

  RelativeCMSubspaceLSJTN::RelativeCMSubspaceLSJTN(int L, int S, int J, int T, int g, int N)
  {

    // set values (MODIFICATION for subspacing by N)
    labels_ = SubspaceLabelsType(L,S,J,T,g,N);
    N_ = N;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over total oscillator quanta -- omit (MODIFICATION for subspacing by N)
    // for (int N = g; N <= Nmax; N +=2)
    // iterate over relative oscillator (Nr,lr) orbitals
    for (int Nr = 0; Nr <= N; ++Nr)
      for (int lr = Nr%2; lr <= Nr; lr +=2) 
        {
          // iterate over c.m. oscillator (Nc,lc) orbitals
          // subject to given total N
          int Nc = N - Nr;
            
          for (int lc = Nc%2; lc <= Nc; lc +=2) 
            {

              // impose triangularity
              if (!(am::AllowedTriangle(lr,lc,L)))
                continue;

              // impose antisymmetry
              if (!((lr+S+T)%2==1))
                continue;

              // keep surviving states
              PushStateLabels(StateLabelsType(Nr,lr,Nc,lc)); 
            }
        }
  }

  bool RelativeCMSubspaceLSJTN::ValidLabels() const
  {
    bool valid = true;
      
    // triangularity
    valid &= am::AllowedTriangle(L(),S(),J());

    // parity (MODIFICATION for subspacing by N)
    valid &= ((N()%2)==g());

    return valid;
  }

  std::string RelativeCMSubspaceLSJTN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << L() 
       << " " << std::setw(width) << S() 
       << " " << std::setw(width) << J() 
       << " " << std::setw(width) << T() 
       << " " << std::setw(width) << g()
       << " " << ";"
       << " " << std::setw(width) << N()
       << " " << "]";

    return os.str();
  }


  std::string RelativeCMSubspaceLSJTN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        RelativeCMStateLSJTN state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(width) << state_index
	  << " " << "Nr lr Nc lc"
	  << " " << std::setw(width) << state.Nr() 
	  << " " << std::setw(width) << state.lr() 
	  << " " << std::setw(width) << state.Nc() 
	  << " " << std::setw(width) << state.lc()
          << std::endl;
      }

    return os.str();
  }

  RelativeCMSpaceLSJTN::RelativeCMSpaceLSJTN(int Nmax)
  {
    // save Nmax
    Nmax_ = Nmax;

   // iterate over L
    for (int L=0; L<=Nmax; ++L)
	// iterate over S
	for (int S=0; S<=1; ++S)
	    // iterate over J
	    // imposing triangularity (LSJ)
	    for (int J=abs(L-S); J<=L+S; ++J)
		// iterate over T
		for (int T=0; T<=1; ++T)
		    // iterate over g
		    for (int g=0; g<=1; ++g)
                      // iterate over total oscillator quanta (MODIFICATION for subspacing by N)
                      for (int N = g; N <= Nmax; N+=2)
                        {                      
                          RelativeCMSubspaceLSJTN subspace(L,S,J,T,g,N); // (MODIFICATION for subspacing by N)
                          if (subspace.size()!=0)
                            PushSubspace(subspace);
                        }
  }

  std::string RelativeCMSpaceLSJTN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
	os
	  << " " << "index"
	  << " " << std::setw(width) << subspace_index
	  << " " << "LSJTg"
	  << " " << std::setw(width) << subspace.L() 
	  << " " << std::setw(width) << subspace.S() 
	  << " " << std::setw(width) << subspace.J() 
	  << " " << std::setw(width) << subspace.T() 
	  << " " << std::setw(width) << subspace.g()
	  << " " << "N"  // (MODIFICATION for subspacing by N)
	  << " " << std::setw(width) << subspace.N()  // (MODIFICATION for subspacing by N)
	  << " " << "dim"
	  << " " << std::setw(width) << subspace.size()
	  << " " << std::endl;
      }

    return os.str();

  }

  RelativeCMSectorsLSJTN::RelativeCMSectorsLSJTN(
      const RelativeCMSpaceLSJTN& space,
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
  // two-body states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  TwoBodySubspaceLSJT::TwoBodySubspaceLSJT(
      int L, int S, int J, int T, int g, int truncation_cutoff, int truncation_rank
    )
  {

    // set labels
    labels_ = SubspaceLabelsType(L,S,J,T,g);

    // save truncation
    N1max_ = truncation_cutoff;
    if (truncation_rank==1)
        N2max_ = 2*truncation_cutoff;
    else if (truncation_rank==2)
        N2max_ = truncation_cutoff;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over total oscillator quanta
    for (int N = g; N <= N2max_; N +=2)
      // iterate over oscillator (Nl) orbitals for particle 1
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
          for (int l1 = N1%2; l1 <= N1; l1 +=2) 
            {
              // iterate over oscillator (Nl) orbitals for particle 2
              // subject to given total N
              int N2 = N - N1;
            
              for (int l2 = N2%2; l2 <= N2; l2 +=2) 
                {

                  // impose canonical ordering on single-particle states
                  if (!( std::make_pair(N1,l1) <= std::make_pair(N2,l2) ))
                    continue;

                  // impose triangularity
                  if (!(am::AllowedTriangle(l1,l2,L)))
                    continue;

                  // impose antisymmetry
                  if ((N1==N2)&&(l1==l2)&&(!((L+S+T)%2==1)))
                    continue;

                  // keep surviving states
                  PushStateLabels(StateLabelsType(N1,l1,N2,l2)); 
                }
            }
      }
  }

  bool TwoBodySubspaceLSJT::ValidLabels() const
  {
    bool valid = true;
      
    // triangularity
    valid &= am::AllowedTriangle(L(),S(),J());

    // truncation
    // valid &= ((Nmax()%2)==g());

    return valid;
  }

  std::string TwoBodySubspaceLSJT::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << L() 
       << " " << std::setw(width) << S() 
       << " " << std::setw(width) << J() 
       << " " << std::setw(width) << T() 
       << " " << std::setw(width) << g()
       << " " << "]";

    return os.str();
  }

  std::string TwoBodySubspaceLSJT::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        TwoBodyStateLSJT state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(width) << state_index
	  << " " << "N1 l1 N2 l2"
	  << " " << std::setw(width) << state.N1() 
	  << " " << std::setw(width) << state.l1() 
	  << " " << std::setw(width) << state.N2() 
	  << " " << std::setw(width) << state.l2()
          << std::endl;
      }

    return os.str();

  }

  TwoBodySpaceLSJT::TwoBodySpaceLSJT(int truncation_cutoff, int truncation_rank)
  {

    // save truncation
    N1max_ = truncation_cutoff;
    if (truncation_rank==1)
        N2max_ = 2*truncation_cutoff;
    else if (truncation_rank==2)
        N2max_ = truncation_cutoff;

    // iterate over L
    for (int L=0; L<=N2max_; ++L)
	// iterate over S
	for (int S=0; S<=1; ++S)
	    // iterate over J
	    // imposing triangularity (LSJ)
	    for (int J=abs(L-S); J<=L+S; ++J)
		// iterate over T
		for (int T=0; T<=1; ++T)
		    // iterate over g
		    for (int g=0; g<=1; ++g)
		      {
			
			// downshift Nmax to match parity of subspace
			// required to pass label validity tests
			// int Nmax_subspace = Nmax - (Nmax-g)%2;
		    
			TwoBodySubspaceLSJT subspace(L,S,J,T,g,truncation_cutoff,truncation_rank);
			if (subspace.size()!=0)
			  PushSubspace(subspace);
		      }
  }


  std::string TwoBodySpaceLSJT::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
	os
	  << " " << "index"
	  << " " << std::setw(width) << subspace_index
	  << " " << "LSJTg"
	  << " " << std::setw(width) << subspace.L() 
	  << " " << std::setw(width) << subspace.S() 
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

  TwoBodySectorsLSJT::TwoBodySectorsLSJT(
      const TwoBodySpaceLSJT& space,
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
  // two-body states in LSJT scheme -- subspaced by N
  ////////////////////////////////////////////////////////////////

  TwoBodySubspaceLSJTN::TwoBodySubspaceLSJTN(
      int L, int S, int J, int T, int g, int N,
      int truncation_cutoff, int truncation_rank
    )
  {

    // set labels (MODIFICATION for subspacing by N)
    labels_ = SubspaceLabelsType(L,S,J,T,g,N);
    N_ = N;

    // save truncation
    N1max_ = truncation_cutoff;
    if (truncation_rank==1)
        N2max_ = 2*truncation_cutoff;
    else if (truncation_rank==2)
        N2max_ = truncation_cutoff;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over total oscillator quanta -- omit (MODIFICATION for subspacing by N)
    // for (int N = g; N <= Nmax; N +=2)
    // iterate over oscillator (Nl) orbitals for particle 1
    //
    // Constraints:
    //   0 <= N1 <= N1max
    //   0 <= N2 <= N1max
    //   N1+N2=N
    // ==>  N-min(N1max,N) <= N1 <= min(N1max,N)
    int N1_lower = N-std::min(N1max_,N);
    int N1_upper = std::min(N1max_,N);
    for (int N1 = N1_lower; N1 <= N1_upper; ++N1)
      for (int l1 = N1%2; l1 <= N1; l1 +=2) 
        {
          // iterate over oscillator (Nl) orbitals for particle 2
          // subject to given total N
          int N2 = N - N1;
            
          for (int l2 = N2%2; l2 <= N2; l2 +=2) 
            {

              // impose canonical ordering on single-particle states
              if (!( std::make_pair(N1,l1) <= std::make_pair(N2,l2) ))
                continue;

              // impose triangularity
              if (!(am::AllowedTriangle(l1,l2,L)))
                continue;

              // impose antisymmetry
              if ((N1==N2)&&(l1==l2)&&(!((L+S+T)%2==1)))
                continue;

              // keep surviving states
              PushStateLabels(StateLabelsType(N1,l1,N2,l2)); 
            }
        }
  }

  bool TwoBodySubspaceLSJTN::ValidLabels() const
  {
    bool valid = true;
      
    // triangularity
    valid &= am::AllowedTriangle(L(),S(),J());

    // parity (MODIFICATION for subspacing by N)
    valid &= ((N()%2)==g());

    return valid;
  }

  std::string TwoBodySubspaceLSJTN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << L() 
       << " " << std::setw(width) << S() 
       << " " << std::setw(width) << J() 
       << " " << std::setw(width) << T() 
       << " " << std::setw(width) << g()
       << " " << ";"
       << " " << std::setw(width) << N()
       << " " << "]";

    return os.str();
  }

  std::string TwoBodySubspaceLSJTN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        TwoBodyStateLSJTN state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(width) << state_index
	  << " " << "N1 l1 N2 l2"
	  << " " << std::setw(width) << state.N1() 
	  << " " << std::setw(width) << state.l1() 
	  << " " << std::setw(width) << state.N2() 
	  << " " << std::setw(width) << state.l2()
          << std::endl;
      }

    return os.str();

  }


  TwoBodySpaceLSJTN::TwoBodySpaceLSJTN(int truncation_cutoff, int truncation_rank)
  {
    // save truncation
    N1max_ = truncation_cutoff;
    if (truncation_rank==1)
        N2max_ = 2*truncation_cutoff;
    else if (truncation_rank==2)
        N2max_ = truncation_cutoff;

    // iterate over L
    for (int L=0; L<=N2max_; ++L)
	// iterate over S
	for (int S=0; S<=1; ++S)
	    // iterate over J
	    // imposing triangularity (LSJ)
	    for (int J=abs(L-S); J<=L+S; ++J)
		// iterate over T
		for (int T=0; T<=1; ++T)
		    // iterate over g
		    for (int g=0; g<=1; ++g)
                      // iterate over total oscillator quanta (MODIFICATION for subspacing by N)
                      for (int N = g; N <= N2max_; N+=2)
                        {
                          TwoBodySubspaceLSJTN subspace(
                              L,S,J,T,g,N,truncation_cutoff,truncation_rank
                            );  // (MODIFICATION for subspacing by N)
                          if (subspace.size()!=0)
                            PushSubspace(subspace);
                        }
  }
  
  std::string TwoBodySpaceLSJTN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        os
          << " " << "index"
          << " " << std::setw(width) << subspace_index
          << " " << "LSJTg"
          << " " << std::setw(width) << subspace.L() 
          << " " << std::setw(width) << subspace.S() 
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

  TwoBodySectorsLSJTN::TwoBodySectorsLSJTN(
      const TwoBodySpaceLSJTN& space,
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
