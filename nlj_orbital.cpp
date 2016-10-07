/****************************************************************
  nlj_orbital.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>

#include "am/am.h"

#include "nlj_orbital.h"

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

  std::string OrbitalSubspacePN::LabelStr() const
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << int(orbital_species())
       << " " << "]";

    return os.str();
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

  ////////////////////////////////////////////////////////////////
  // single-particle orbitals - lj subspaces
  ////////////////////////////////////////////////////////////////

  OrbitalSubspaceLJPN::OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species,
                                           int l, HalfInt j, int Nmax) {
    // set values
    labels_ = SubspaceLabelsType(orbital_species,l,j);
    weight_max_ = double(Nmax);
    Nmax_ = Nmax;

    // iterate over radial quantum number
    for (int n = 0; (2*n+l) <= Nmax; ++n) {
      // save state
      PushStateLabels(StateLabelsType(n));

      // save oscillator quantum number as weight
      weights_.push_back(double(2*n+l));
    }
  }

  OrbitalSubspaceLJPN::OrbitalSubspaceLJPN(OrbitalSpeciesPN orbital_species,
                                           int l, HalfInt j,
                                           std::vector<OrbitalPNInfo> const &states) {
    labels_ = SubspaceLabelsType(orbital_species,l,j);
    weight_max_ = 0.0;
    Nmax_ = -1;
    for (const OrbitalPNInfo& state : states) {
      if (state.orbital_species == orbital_species
          && state.l == l && state.j == j) {
        PushStateLabels(StateLabelsType(state.n));
        weights_.push_back(state.weight);
        if (state.weight>weight_max_) weight_max_ = state.weight;
      }
    }
  }


  std::string OrbitalSubspaceLJPN::LabelStr() const {
    std::ostringstream os;

    const int width = 3;

    os << "["
       << " " << std::setw(width) << int(orbital_species())
       << " " << std::setw(width) << l()
       << " " << std::setw(width+2) << j().Str()
       << " " << "]";

    return os.str();
  }

  std::string OrbitalSubspaceLJPN::DebugStr() const {
    std::ostringstream os;

    const int width = 3;

    for (int state_index=0; state_index<size(); ++state_index) {
      OrbitalStateLJPN state(*this,state_index);

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

  OrbitalSpaceLJPN::OrbitalSpaceLJPN(int Nmax) {
    // save truncation
    weight_max_ = double(Nmax);
    Nmax_ = Nmax;

    // iterate over species
    for (OrbitalSpeciesPN orbital_species :
          {OrbitalSpeciesPN::kP,OrbitalSpeciesPN::kN}) {
      for (int l=0; l<=Nmax; ++l) {
        for (HalfInt j = l-HalfInt(1,2); j<=(l+HalfInt(1,2)); ++j) {
          if (j<0) continue;
          OrbitalSubspaceLJPN subspace(orbital_species,l,j,Nmax);
          PushSubspace(subspace);
        }
      }
    }
  }

  OrbitalSpaceLJPN::OrbitalSpaceLJPN(std::vector<OrbitalPNInfo> const &states) {
    for (const OrbitalPNInfo& state : states) {
      OrbitalSubspaceLJPNLabels labels(state.orbital_species,state.l,state.j);
      if (!ContainsSubspace(labels)) {
        OrbitalSubspaceLJPN subspace(state.orbital_species,state.l,state.j,states);
        PushSubspace(subspace);
      }
    }
  }

  std::string OrbitalSpaceLJPN::DebugStr() const {

    std::ostringstream os;

    const int width = 3;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index) {
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

  std::string OrbitalSpaceLJPN::OrbitalDefinitionStr() const {

    std::ostringstream header;
    std::ostringstream body;
    std::ostringstream os;

    // header comments
    header << "# MFDn SPorbital file" << std::endl;
    header << "#   version" << std::endl;
    header << "#   norb_p norb_n" << std::endl;
    header << "#   index n l 2*j species weight" << std::endl;

    // header line 1: version
    int version = 15055;
    header << version << std::endl;

    // data
    const int width = 3;
    const int precision = 8;
    // orbital indices
    int p_index = 0;
    int n_index = 0;
    int output_index = 0;
    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
    {
      // retrieve subspace
      const SubspaceType& subspace = GetSubspace(subspace_index);

      // iterate over states
      for (int state_index=0; state_index<subspace.size(); ++state_index) {
        OrbitalStateLJPN state(subspace,state_index);
        if (state.orbital_species() == OrbitalSpeciesPN::kP) {
          output_index = ++p_index;
        } else if (state.orbital_species() == OrbitalSpeciesPN::kN) {
          output_index = ++n_index;
        }

        body << " " << std::setw(width) << output_index
             << " " << std::setw(width) << state.n()
             << " " << std::setw(width) << state.l()
             << " " << std::setw(width) << TwiceValue(state.j())
             << " " << std::setw(width) << int(state.orbital_species())+1 // 1-based
             << " " << std::fixed << std::setw(width+1+precision)
             << std::setprecision(precision) << state.weight()
             << std::endl;
      }
    }

    // header line 2: dimensions
    header << p_index << " " << n_index << std::endl;

    // assemble file
    os << header.str() << body.str();

    return os.str();

  }

  OrbitalSectorsLJPN::OrbitalSectorsLJPN(
      const OrbitalSpaceLJPN& space,
      basis::SectorDirection sector_direction)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index) {
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index) {
        if ((sector_direction == basis::SectorDirection::kCanonical)
            && (bra_subspace_index>ket_subspace_index)) {
          continue;
        }

        // retrieve subspaces
        const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
        const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

        // push sector
        PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
      }
    }
  }

  OrbitalSectorsLJPN::OrbitalSectorsLJPN(
      const OrbitalSpaceLJPN& space,
      int l0, int j0, int g0,
      basis::SectorDirection sector_direction)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index) {
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index) {
        if ((sector_direction == basis::SectorDirection::kCanonical)
            && (bra_subspace_index>ket_subspace_index)) {
          continue;
        }

        // retrieve subspaces
        const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
        const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

        bool allowed = true;
        allowed &= am::AllowedTriangle(ket_subspace.l(),l0,bra_subspace.l());
        allowed &= am::AllowedTriangle(ket_subspace.j(),j0,bra_subspace.j());
        allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);

        // push sector
        if (allowed) {
          PushSector(SectorType(bra_subspace_index,ket_subspace_index,
                                bra_subspace,ket_subspace));
        }
      }
    }
  }

  std::string OrbitalSectorsLJPN::DebugStr() const {
    std::ostringstream os;
    int width = 3;

    for (int sector_index=0; sector_index<size(); ++sector_index) {
      // const basis::OrbitalSectorLJPN& sector = GetSector(sector_index);
      // os << sector_index+1 << sector.DebugStr();
      const basis::BaseSector<basis::OrbitalSubspaceLJPN>& sector =
        GetSector(sector_index);
      os << std::setw(width) << sector_index+1
         << " bra " << std::setw(width) << sector.bra_subspace_index() + 1
            << " (" << sector.bra_subspace().l()
            << ", " << sector.bra_subspace().j().Str()
            << ", " << int(sector.bra_subspace().orbital_species())+1 << ")"
         << " ket " << std::setw(width) << sector.ket_subspace_index() + 1
            << " (" << sector.ket_subspace().l()
            << ", " << sector.ket_subspace().j().Str()
            << ", " << int(sector.ket_subspace().orbital_species())+1 << ")"
         << std::endl;
    }
    return os.str();
  }

  // std::string OrbitalSectorLJPN::DebugStr() const {
  //   std::ostringstream os;
  //
  //   const int width = 3;
  //
  //   os << " bra " << std::setw(width) << bra_subspace_index() + 1
  //        << " (" << bra_subspace().l()
  //        << ", " << bra_subspace().j().Str()
  //        << ", " << int(bra_subspace().orbital_species())+1 << ")"
  //      << " ket " << std::setw(width) << ket_subspace_index() + 1
  //        << " (" << ket_subspace().l()
  //        << ", " << ket_subspace().j().Str()
  //        << ", " << int(ket_subspace().orbital_species())+1 << ")"
  //     << std::endl;
  //
  //   return os.str();
  // }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
