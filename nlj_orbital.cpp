/****************************************************************
  nlj_orbital.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <iomanip>  // for debugging output
#include <iostream>
#include <sstream>
#include <set>

#include "am/am.h"
#include "mcutils/parsing.h"

#include "nlj_orbital.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle definition file parsing
  ////////////////////////////////////////////////////////////////

  std::vector<OrbitalPNInfo> ParseOrbitalPNStream(std::istream& is, bool standalone)
  {

    // set up line counter for use in error messages
    std::string line;
    int line_count = 0;


    int num_orbitals_p, num_orbitals_n;
    if (standalone)
      {
        // line 1: version -- but first gobble any comment lines
        while (std::getline(is,line), line[0]=='#') {++line_count;};
        {
          ++line_count;
          std::istringstream line_stream(line);
          int version;
          line_stream >> version;
          ParsingCheck(line_stream,line_count,line);
          assert(version==15055);
        }

        // line 2: number of p,n orbitals
        {
          ++line_count;
          std::getline(is,line);
          std::istringstream line_stream(line);
          line_stream >> num_orbitals_p >> num_orbitals_n;
          ParsingCheck(line_stream,line_count,line);
        }
      }

    // lines 3+: orbital definitions
    std::vector<OrbitalPNInfo> states;
    int num_orbitals_p_extracted=0, num_orbitals_n_extracted=0;
    while ( getline(is,line))
      {
        // count line
        ++line_count;

        // set up for parsing
        std::istringstream line_stream(line);
        if (line.size() == 0)
          continue;

        int index, n, l, twice_j, orbital_species_raw;
        double weight;
        line_stream >> index
                    >> n >> l >> twice_j
                    >> orbital_species_raw >> weight;
        ParsingCheck(line_stream,line_count,line);

        HalfInt j(twice_j,2);
        OrbitalSpeciesPN orbital_species = static_cast<OrbitalSpeciesPN>(orbital_species_raw-1);

        // count orbitals by type
        num_orbitals_p_extracted += int(orbital_species == OrbitalSpeciesPN::kP);
        num_orbitals_n_extracted += int(orbital_species == OrbitalSpeciesPN::kN);

        OrbitalPNInfo state(orbital_species,n,l,j,weight);
        states.push_back(state);
      }

    if (standalone)
      {
        assert(num_orbitals_p==num_orbitals_p_extracted);
        assert(num_orbitals_n==num_orbitals_n_extracted);
      }

    return states;
  }

  std::string OrbitalDefinitionStr(const std::vector<OrbitalPNInfo>& orbitals, bool standalone)
  {

    std::ostringstream header;
    std::ostringstream body;
    std::ostringstream os;

    // construct body
    const int width = 3;
    const int precision = 8;
    // orbital indices
    int p_index = 0;
    int n_index = 0;
    int output_index = 0;
    for (const OrbitalPNInfo& state: orbitals)
      {
        // iterate over states
        if (state.orbital_species == OrbitalSpeciesPN::kP) {
          output_index = ++p_index;
        } else if (state.orbital_species == OrbitalSpeciesPN::kN) {
          output_index = ++n_index;
        }

        body << " " << std::setw(width) << output_index
             << " " << std::setw(width) << state.n
             << " " << std::setw(width) << state.l
             << " " << std::setw(width) << TwiceValue(state.j)
             << " " << std::setw(width) << int(state.orbital_species)+1 // 1-based
             << " " << std::fixed << std::setw(width+1+precision)
             << std::setprecision(precision) << state.weight
             << std::endl;
      }

    // construct header
    //
    // We defer constructing the header until after constructing the
    // body, since we need statistics obtained while constructing the
    // body.
    if (standalone)
      {
        // header comments
        header << "# MFDn SPorbital file" << std::endl;
        header << "#   version" << std::endl;
        header << "#   norb_p norb_n" << std::endl;
        header << "#   index n l 2*j species weight" << std::endl;

        // header line 1: version
        int version = 15055;
        header << version << std::endl;

        // header line 2: dimensions
        header << p_index << " " << n_index << std::endl;
      }

    // assemble file
    os << header.str() << body.str();

    return os.str();

  }

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

  OrbitalSubspacePN::OrbitalSubspacePN(OrbitalSpeciesPN orbital_species,
                                       const std::vector<OrbitalPNInfo>& states)
  {
    labels_ = SubspaceLabelsType(orbital_species);
    weight_max_ = 0.0;
    Nmax_ = -1; // TODO call IsOscillatorLike_()
    for (const OrbitalPNInfo& state : states) {
      if (state.orbital_species == orbital_species) {
        PushStateLabels(StateLabelsType(state.n,state.l,state.j));
        weights_.push_back(state.weight);
        weight_max_ = std::max(weight_max_,state.weight);
      }
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

  std::vector<OrbitalPNInfo> OrbitalSubspacePN::OrbitalInfo() const
  {
    std::vector<OrbitalPNInfo> orbitals;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        OrbitalStatePN state(*this,state_index);
        orbitals.push_back(state.OrbitalInfo());
      }

    return orbitals;

  }

  OrbitalPNInfo OrbitalStatePN::OrbitalInfo() const
  {
    OrbitalPNInfo orbital;

    orbital.orbital_species = orbital_species();
    orbital.n = n();
    orbital.l = l();
    orbital.j = j();
    orbital.weight = weight();

    return orbital;
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

  OrbitalSpacePN::OrbitalSpacePN(const std::vector<OrbitalPNInfo>& states)
  {
    weight_max_ = 0.0;
    Nmax_ = -1; // TODO call IsOscillatorLike_()
    // collect orbital_species subspace labels sorted in canonical order
    std::set<OrbitalSubspacePNLabels> subspace_labels_set;
    for (int state_index=0; state_index<states.size(); ++state_index)
      {
        OrbitalPNInfo state = states[state_index];
        OrbitalSubspacePNLabels labels(state.orbital_species);
        subspace_labels_set.insert(labels);
      }

    // construct subspaces
    for (const OrbitalSubspacePNLabels& labels : subspace_labels_set)
      {
        OrbitalSpeciesPN orbital_species;
        std::tie(orbital_species) = labels;
        OrbitalSubspacePN subspace(orbital_species,states);
        PushSubspace(subspace);
        weight_max_ = std::max(weight_max_,subspace.weight_max());
      }

  }

  bool OrbitalSpacePN::IsOscillatorLike()
  {

    bool is_oscillator_like = true;

    // first check if subspaces are populated, else weight_max is nonsense
    is_oscillator_like &= ((GetSubspace(0).size()!=0)&&(GetSubspace(1).size()!=0));
    if (!is_oscillator_like)
      return false;

    // see if can extract viable Nmax
    is_oscillator_like &= ((GetSubspace(0).size()!=0)&&(GetSubspace(1).size()!=0));
    double weight_max_p = GetSubspace(0).weight_max();
    double weight_max_n = GetSubspace(1).weight_max();
    int Nmax = int(weight_max_p);
    is_oscillator_like &= (Nmax == weight_max_p); // check integer weight
    is_oscillator_like &= (Nmax >=0 ); // check positive weight
    is_oscillator_like &= (weight_max_p == weight_max_n); // check proton and neutron consistent
    if (!is_oscillator_like)
      return false;

    // TODO check all orbital labels and weights against Nmax-constructed space

    // orbitals are oscillator like!
    Nmax_ = Nmax;
    return true;
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

  std::vector<OrbitalPNInfo> OrbitalSpacePN::OrbitalInfo() const
  {
    std::vector<OrbitalPNInfo> orbitals;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        std::vector<OrbitalPNInfo> subspace_orbitals;

        // get orbitals for subspace and append to vector
        subspace_orbitals = subspace.OrbitalInfo();
        orbitals.insert(orbitals.end(),
                        subspace_orbitals.begin(),
                        subspace_orbitals.end()
                       );

      }

    return orbitals;

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
                                           const std::vector<OrbitalPNInfo>& states) {
    labels_ = SubspaceLabelsType(orbital_species,l,j);
    weight_max_ = 0.0;
    Nmax_ = -1;
    for (auto&& state : states) {
      if (state.orbital_species == orbital_species
          && state.l == l && state.j == j) {
        PushStateLabels(StateLabelsType(state.n));
        weights_.push_back(state.weight);
        weight_max_ = std::max(weight_max_,state.weight);
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

  std::vector<OrbitalPNInfo> OrbitalSubspaceLJPN::OrbitalInfo() const
  {
    std::vector<OrbitalPNInfo> orbitals;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        OrbitalStateLJPN state(*this,state_index);
        orbitals.push_back(state.OrbitalInfo());
      }

    return orbitals;

  }

  OrbitalPNInfo OrbitalStateLJPN::OrbitalInfo() const
  {
    OrbitalPNInfo orbital;

    orbital.orbital_species = orbital_species();
    orbital.n = n();
    orbital.l = l();
    orbital.j = j();
    orbital.weight = weight();

    return orbital;
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

  OrbitalSpaceLJPN::OrbitalSpaceLJPN(const std::vector<OrbitalPNInfo>& states)
  {
    weight_max_ = 0.0;
    Nmax_ = -1; // TODO call IsOscillatorLike_()

    // collect (l,j) subspace labels sorted in canonical order
    std::set<OrbitalSubspaceLJPNLabels> subspace_labels_set;
    for (int state_index=0; state_index<states.size(); ++state_index)
      {
        OrbitalPNInfo state = states[state_index];
        OrbitalSubspaceLJPNLabels labels(state.orbital_species,state.l,state.j);
        subspace_labels_set.insert(labels);
      }

    // construct subspaces
    for (const OrbitalSubspaceLJPNLabels& labels : subspace_labels_set)
      {
        OrbitalSpeciesPN orbital_species;
        int l;
        HalfInt j;
        std::tie(orbital_species,l,j) = labels;
        OrbitalSubspaceLJPN subspace(orbital_species,l,j,states);
        PushSubspace(subspace);
        weight_max_ = std::max(weight_max_,subspace.weight_max());
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

  std::vector<OrbitalPNInfo> OrbitalSpaceLJPN::OrbitalInfo() const
  {
    std::vector<OrbitalPNInfo> orbitals;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
        std::vector<OrbitalPNInfo> subspace_orbitals;

        // get orbitals for subspace and append to vector
        subspace_orbitals = subspace.OrbitalInfo();
        orbitals.insert(orbitals.end(),
                        subspace_orbitals.begin(),
                        subspace_orbitals.end()
                       );

      }

    return orbitals;

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

  std::string OrbitalSectorsLJPN::DebugStr() const
  {
    std::ostringstream os;
    int width = 3;

    for (int sector_index=0; sector_index<size(); ++sector_index) {
      // const basis::OrbitalSectorLJPN& sector = GetSector(sector_index);
      // os << sector_index+1 << sector.DebugStr();
      const basis::BaseSector<basis::OrbitalSubspaceLJPN>& sector =
        GetSector(sector_index);
      os << std::setw(width) << sector_index
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
