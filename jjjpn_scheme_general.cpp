/****************************************************************
  jjjpn_scheme_general.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>

#include "am/am.h"

#include "jjjpn_scheme_general.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle orbitals
  ////////////////////////////////////////////////////////////////

  OrbitalSubspacePN::OrbitalSubspacePN(OrbitalSpeciesPN orbital_species, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(orbital_species);
    weight_max_ = double(Nmax);
    Nmax_ = Nmax;

    // iterate over total oscillator quanta
    for (int N = 0; N <= Nmax; ++N)
      // iterate over j within shell
      for (HalfInt j = HalfInt(1,2); j <= N+HalfInt(1,2); ++j) 
        {
          // recover derived quantum numbers (n,l) from (N,j)
          int l = (TwiceValue(j)-1)/2 + (N+(TwiceValue(j)-1)/2)%2;
          int n = (N-l)/2;
            
          // save state
          PushStateLabels(StateLabelsType(n,l,j)); 
            
          // save oscillator quantum number as weight
          weights_.push_back(double(N));
        }
  }

  std::string OrbitalSubspacePN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        OrbitalStatePN state(*this,state_index);

        os
	  << " " << "index"
	  << " " << std::setw(width) << state_index
	  << " " << "nlj"
	  << " " << std::setw(width) << state.n()
	  << " " << std::setw(width) << state.l()
	  << " " << std::setw(width+2) << state.j().Str()
	  << " " << "weight"
          << " " << state.weight()
          << std::endl;
      }

    return os.str();

  }


  OrbitalSpacePN::OrbitalSpacePN(int Nmax)
  {

    // save truncation
    weight_max_ = double(Nmax);
    Nmax_ = Nmax;

    // iterate over species
    for (OrbitalSpeciesPN orbital_species : {OrbitalSpeciesPN::kP,OrbitalSpeciesPN::kN})
      {
        OrbitalSubspacePN subspace(orbital_species,Nmax);
        PushSubspace(subspace);
      }
  }

  std::string OrbitalSpacePN::DebugStr() const
  {

    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);
        
	os
	  << " " << "index"
	  << " " << std::setw(width) << subspace_index
	  << " " << "species"
	  << " " << std::setw(width) << int(subspace.orbital_species()) 
	  << " " << "dim"
	  << " " << std::setw(width) << subspace.size()
	  << " " << std::endl;
      }

    return os.str();

  }


  std::string OrbitalSpacePN::OrbitalDefinitionStr() const
  {

    std::ostringstream os;

    // header comments
    os << "# MFDn SPorbital file" << std::endl;
    os << "#   version" << std::endl;
    os << "#   norb_p norb_n" << std::endl;
    os << "#   index n l 2*j species weight" << std::endl;

    // header line 1: version
    int version = 15055;
    os << version << std::endl;

    // header line 2: dimensions
    os << GetSubspace(0).size() << " " << GetSubspace(1).size() << std::endl;
        
    // data
    const int width = 3;
    const int precision = 8;
    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        // retrieve subspace
	const SubspaceType& subspace = GetSubspace(subspace_index);

        // iterate over states
        for (int state_index=0; state_index<subspace.size(); ++state_index)
          {
            OrbitalStatePN state(subspace,state_index);
            
            os
              << " " << std::setw(width) << state.index()+1  // 1-based
              << " " << std::setw(width) << state.n()
              << " " << std::setw(width) << state.l()
              << " " << std::setw(width) << TwiceValue(state.j())
              << " " << std::setw(width) << int(state.orbital_species())+1 // 1-based
              << " " << std::fixed << std::setw(width+1+precision)
              << std::setprecision(precision) << state.weight()
              << std::endl;
          }
      }
    
    return os.str();

  }

  //  ////////////////////////////////////////////////////////////////
  //  // two-body states in jjJT scheme
  //  ////////////////////////////////////////////////////////////////
  //
  //  TwoBodySubspaceJJJT::TwoBodySubspaceJJJT(int J, int T, int g, int Nmax)
  //  {
  //
  //    // set values
  //    labels_ = SubspaceLabelsType(J,T,g);
  //    Nmax_ = Nmax;
  //
  //    // validate subspace labels
  //    assert(ValidLabels()); 
  //
  //    // set up indexing
  //    // iterate over total oscillator quanta
  //    for (int N = g; N <= Nmax; N +=2)
  //      // iterate over oscillator (Nj) orbitals for particle 1
  //      for (int N1 = 0; N1 <= N; ++N1)
  //	for (HalfInt j1 = HalfInt(1,2); j1 <= N1+HalfInt(1,2); j1 +=1) 
  //	  {
  //	    // iterate over oscillator (Nj) orbitals for particle 2
  //	    // subject to given total N
  //	    int N2 = N - N1;
  //
  //	    for (HalfInt j2 = HalfInt(1,2); j2 <= N2+HalfInt(1,2); j2 +=1) 
  //	      {
  //
  //                // impose canonical ordering on single-particle states
  //                if (!( std::make_pair(N1,j1) <= std::make_pair(N2,j2) ))
  //                  continue;
  //
  //		// impose triangularity
  //		if (!(am::AllowedTriangle(j1,j2,J)))
  //		  continue;
  //
  //		// impose antisymmetry
  //		if ((N1==N2)&&(j1==j2)&&(!((J+T)%2==1)))
  //		  continue;
  //
  //		// keep surviving states
  //		PushStateLabels(StateLabelsType(N1,j1,N2,j2)); 
  //	      }
  //	  }
  //  }
  //
  //  bool TwoBodySubspaceJJJT::ValidLabels() const
  //  {
  //    bool valid = true;
  //      
  //    // truncation
  //    valid &= ((Nmax()%2)==g());
  //
  //    return valid;
  //  }
  //
  //  std::string TwoBodySubspaceJJJT::DebugStr() const
  //  {
  //
  //    std::ostringstream os;
  //
  //    const int width = 3;
  //
  //    for (int state_index=0; state_index<size(); ++state_index)
  //      {
  //        TwoBodyStateJJJT state(*this,state_index);
  //
  //        os
  //	  << " " << "index"
  //	  << " " << std::setw(width) << state_index
  //	  << " " << "N1 l1 j1 N2 l2 j2"
  //	  << " " << std::setw(width) << state.N1()
  //	  << " " << std::setw(width) << state.l1()
  //	  << " " << std::setw(width) << state.j1()
  //	  << " " << std::setw(width) << state.N2()
  //	  << " " << std::setw(width) << state.l2()
  //	  << " " << std::setw(width) << state.j2()
  //          << std::endl;
  //      }
  //
  //    return os.str();
  //
  //  }
  //
  //
  //  TwoBodySpaceJJJT::TwoBodySpaceJJJT(int Nmax)
  //  {
  //
  //    // save Nmax
  //    Nmax_ = Nmax;
  //
  //    // iterate over J
  //    for (int J=0; J<=Nmax+1; ++J)
  //      {
  //	// iterate over T
  //	for (int T=0; T<=1; ++T)
  //	  {
  //	    // iterate over g
  //	    for (int g=0; g<=1; ++g)
  //	      {
  //		
  //		// downshift Nmax to match parity of subspace
  //		// required to pass label validity tests
  //		int Nmax_subspace = Nmax - (Nmax-g)%2;
  //		    
  //		TwoBodySubspaceJJJT subspace(J,T,g,Nmax_subspace);
  //
  //		if (subspace.size()!=0)
  //		  PushSubspace(subspace);
  //	      }
  //	  }
  //      }
  //  }
  //
  //  std::string TwoBodySpaceJJJT::DebugStr() const
  //  {
  //
  //    std::ostringstream os;
  //
  //    const int width = 3;
  //
  //    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
  //      {
  //	const SubspaceType& subspace = GetSubspace(subspace_index);
  //	os
  //	  << " " << "index"
  //	  << " " << std::setw(width) << subspace_index
  //	  << " " << "JTg"
  //	  << " " << std::setw(width) << subspace.J() 
  //	  << " " << std::setw(width) << subspace.T() 
  //	  << " " << std::setw(width) << subspace.g()
  //	  << " " << "Nmax"
  //	  << " " << std::setw(width) << subspace.Nmax()
  //	  << " " << "dim"
  //	  << " " << std::setw(width) << subspace.size()
  //	  << " " << std::endl;
  //      }
  //
  //    return os.str();
  //
  //  }
  //
  //
  //  TwoBodySectorsJJJT::TwoBodySectorsJJJT(
  //      const TwoBodySpaceJJJT& space,
  //      int J0, int T0, int g0,
  //      basis::SectorDirection sector_direction
  //    )
  //  {
  //    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
  //      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
  //	{
  //
  //          // enforce canonical ordering
  //          if (
  //              (sector_direction == basis::SectorDirection::kCanonical)
  //              && !(bra_subspace_index<=ket_subspace_index)
  //            )
  //            continue;
  //
  //          // retrieve subspaces
  //          const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
  //          const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);
  //
  //          // verify angular momentum, isosopin, and parity selection rules
  //          bool allowed = true;
  //          allowed &= am::AllowedTriangle(ket_subspace.J(),J0,bra_subspace.J());
  //          allowed &= am::AllowedTriangle(ket_subspace.T(),T0,bra_subspace.T());
  //          allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);
  //
  //          // push sector
  //	  if (allowed)
  //            PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
  //        }
  //  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
