/****************************************************************
  lsjt_scheme.cpp

  Mark A. Caprio, University of Notre Dame.

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
    for (int N = g; N <= Nmax; N +=2)
      PushStateLabels(StateLabelsType(N));

    // std::cout << "Done constructing " << size() << std::endl;

  }

  bool RelativeSubspaceLSJT::ValidLabels() const
  {

    bool valid = true;

    // parity
    valid &= (g() == (L()%2));
    // antisymmetry
    valid &= ((L()+S()+T())%2 == 1);
    // triangularity
    valid &= am::AllowedTriangle(L(),S(),J());
    // truncation
    valid &= ((Nmax()%2)==g());

    return valid;
  }


  RelativeSpaceLSJT::RelativeSpaceLSJT(int Nmax)
  {

    // save Nmax
    Nmax_ = Nmax;

    // iterate over L
    for (int L=0; L<=Nmax; ++L)
      {
	int g = L%2;

	// iterate over S
	for (int S=0; S<=1; ++S)
	  {
	    int T = (L+S+1)%2;

	    // iterate over J
	    for (int J=abs(L-S); J<=L+S; ++J)
	      {
		// downshift Nmax to match parity of subspace
		// required to pass label validity tests
		int Nmax_subspace = Nmax - (Nmax-g)%2;

		RelativeSubspaceLSJT subspace(L,S,J,T,g,Nmax_subspace);
		assert(subspace.size()!=0);
		PushSubspace(subspace);
	      }
	  }
      }
  }

  std::string RelativeSpaceLSJT::Str() const
  {
    std::ostringstream os;

    const int lw = 3;

    for (int s=0; s<size(); ++s)
      {
	const SubspaceType& subspace = GetSubspace(s);
	os
	  << " " << "index"
	  << " " << std::setw(lw) << s 
	  << " " << "LSJTg"
	  << " " << std::setw(lw) << subspace.L() 
	  << " " << std::setw(lw) << subspace.S() 
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

  RelativeSectorsLSJT::RelativeSectorsLSJT(const RelativeSpaceLSJT& space)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
	{
          // retrieve subspaces
          const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

          // push sector
          PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
        }
  }

  RelativeSectorsLSJT::RelativeSectorsLSJT(const RelativeSpaceLSJT& space, int J0, int g0)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
	{
          // retrieve subspaces
          const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

          // verify angular momentum and parity selection rules
          bool allowed = true;
          allowed &= am::AllowedTriangle(ket_subspace.J(),J0,bra_subspace.J());
          allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);

          // push sector
	  if (allowed)
            PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
        }
  }


  ////////////////////////////////////////////////////////////////
  // relative-cm states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  RelativeCMSubspaceLSJT::RelativeCMSubspaceLSJT (int Ncm, int lcm, 
                                                  int L, int S, int J, int T, 
                                                  int g, int Nmax)
  {
    // set values
    labels_ = SubspaceLabelsType(Ncm,lcm,L,S,J,T,g);
    Nmax_ = Nmax;

    // validate
    assert(ValidLabels()); 

    // set up indexing

    // iterate over total oscillator quanta
    for (int N = g; N <= Nmax; N +=2)
      {
	// iterate over oscillator (Nl) orbitals for relative motion
	// subject to given total N
	int Nr = N - Ncm;
	for (int lr = Nr%2; lr <= Nr; lr +=2) 
	  {
	    // impose triangularity
	    if (!(am::AllowedTriangle(lr,lcm,L)))
	      continue;

	    // // impose antisymmetry
	    // if (!((lr+S+T)%2==1))
	    //   continue;

	    // keep surviving states
	    PushStateLabels(StateLabelsType(Nr,lr)); 
	  }
      }
      

  }

  bool RelativeCMSubspaceLSJT::ValidLabels() const
  {
    bool valid = true;

    // impose antisymmetry
    valid &= ((Ncm()+S()+T()+g())%2==1);

    // triangularity
    valid &= am::AllowedTriangle(L(),S(),J());

    // truncation
    valid &= ((Nmax()%2)==g());

    return valid;
  }


  ////////////////////////////////////////////////////////////////
  // two-body states in LSJT scheme
  ////////////////////////////////////////////////////////////////

  TwoBodySubspaceLSJT::TwoBodySubspaceLSJT(int L, int S, int J, int T, int g, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(L,S,J,T,g);
    Nmax_ = Nmax;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over total oscillator quanta
    for (int N = g; N <= Nmax; N +=2)
      // iterate over oscillator (Nl) orbitals for particle 1
      for (int N1 = 0; N1 <= N; ++N1)
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

  bool TwoBodySubspaceLSJT::ValidLabels() const
  {
    bool valid = true;
      
    // triangularity
    valid &= am::AllowedTriangle(L(),S(),J());

    // truncation
    valid &= ((Nmax()%2)==g());

    return valid;
  }


  TwoBodySpaceLSJT::TwoBodySpaceLSJT(int Nmax)
  {
    // save Nmax
    Nmax_ = Nmax;

    // iterate over L
    for (int L=0; L<=Nmax; ++L)
      {
	// iterate over S
	for (int S=0; S<=1; ++S)
	  {
	    // iterate over J
	    // imposing triangularity (LSJ)
	    for (int J=abs(L-S); J<=L+S; ++J)
	      {

		// iterate over T
		for (int T=0; T<=1; ++T)
		  {

		    // iterate over g
		    for (int g=0; g<=1; ++g)

		      {
			
			// downshift Nmax to match parity of subspace
			// required to pass label validity tests
			int Nmax_subspace = Nmax - (Nmax-g)%2;
		    
			TwoBodySubspaceLSJT subspace(L,S,J,T,g,Nmax_subspace);
			// std::cout 
			//    << std::setw(3) << L 
			// 	  << std::setw(3) << S 
			// 	  << std::setw(3) << T 
			// 	  << std::setw(3) << J 
			// 	  << std::setw(3) << g 
			// 	  << std::setw(3) << Nmax_for_subspace 
			// 	  << std::setw(3) << subspace.size()
			// 	  << std::endl;
			if (subspace.size()!=0)
			  PushSubspace(subspace);
		      }
		  }
	      }
	  }
      }
  }

  std::string TwoBodySpaceLSJT::Str() const
  {

    std::ostringstream os;

    const int lw = 3;

    for (int s=0; s<size(); ++s)
      {
	const SubspaceType& subspace = GetSubspace(s);
	os
	  << " " << "index"
	  << " " << std::setw(lw) << s 
	  << " " << "LSJTg"
	  << " " << std::setw(lw) << subspace.L() 
	  << " " << std::setw(lw) << subspace.S() 
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

  TwoBodySectorsLSJT::TwoBodySectorsLSJT(const TwoBodySpaceLSJT& space)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
	{
          // retrieve subspaces
          const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

          // push sector
          PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
        }
  }

  TwoBodySectorsLSJT::TwoBodySectorsLSJT(const TwoBodySpaceLSJT& space, int J0, int g0)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
	{
          // retrieve subspaces
          const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

          // verify angular momentum and parity selection rules
          bool allowed = true;
          allowed &= am::AllowedTriangle(ket_subspace.J(),J0,bra_subspace.J());
          allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);

          // push sector
	  if (allowed)
            PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
        }
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
