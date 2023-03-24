/****************************************************************
  oscillator_orbital.cpp

  Mark A. Caprio, Patrick J. Fasano, Jakub Herko, Shwetha L. Vittal
  University of Notre Dame

****************************************************************/


// #include <cstddef>
// #include <algorithm>
#include <iomanip>
// #include <set>
// #include <istream>
#include <sstream>

#include "am/am.h"

#include "oscillator_orbital.h"
#include "basis.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // states (single-particle orbitals)
  ////////////////////////////////////////////////////////////////

  OscillatorOrbitalSubspace::OscillatorOrbitalSubspace(int l, int Nmax)
    : BaseSubspace{{l}}, Nmax_{Nmax}
  {
    
    // validate subspace labels
    assert(ValidLabels());

    // iterate over total oscillator quanta
    for (int n = 0; n <= (Nmax-l)/2; ++n)
      // DEBUGGING: Must use arguments Nmax and l (or data members Nmax_ and l_)
      // here, instead of (hidden) accessors Nmax() and l().
      PushStateLabels(StateLabelsType(n));
  }

  bool OscillatorOrbitalSubspace::ValidLabels() const
  {

    bool valid = true;

    // nonnegativity
    valid &= l()>=0;

    // nonemptiness in given Nmax truncation
    valid &= l()<=Nmax();
    
    return valid;
  }

  std::string OscillatorOrbitalSubspace::LabelStr() const
  // TODO (mac): reimplement in terms of fmt library
  {
    std::ostringstream os;

    const int width = 0;  // for now, no fixed width

    os << "["
       << " " << std::setw(width) << int(l())
       << " " << "]";

    return os.str();
  }

  std::string OscillatorOrbitalSubspace::DebugStr() const
  // TODO (mac): reimplement in terms of fmt library
  {

    std::ostringstream os;

    const int width = 3;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        StateType state(*this,state_index);

        os
          << " " << "index"
          << " " << std::setw(width) << state_index
          << " " << "n"
          << " " << std::setw(width) << state.n()
          << " " << "=> " << "N"
          << " " << std::setw(width) << state.N()
          << std::endl;
      }

    return os.str();

  }

  std::string OscillatorOrbitalState::LabelStr() const
  // TODO (mac): reimplement in terms of fmt library
  {
    std::ostringstream os;
  
    const int width = 0;  // for now, no fixed width
  
    os << "["
       << " " << "l"
       << " " << std::setw(width) << l()
       << " " << "index"
       << " " << std::setw(width) << index()
       << " :"
       << " " << "n"
       << " " << std::setw(width) << n()
       << " " << "]";
  
    return os.str();
  }
  
  OscillatorOrbitalSpace::OscillatorOrbitalSpace(int Nmax)
    : Nmax_(Nmax)
  {
  
    // iterate over l
    for (int l=0; l<=Nmax_; ++l)
      {
        SubspaceType subspace(l,Nmax);
        PushSubspace(subspace);
      }
  }
  std::string OscillatorOrbitalSpace::DebugStr() const
  {
  
    std::ostringstream os;
  
    const int width = 3;
  
    os << " Nmax " << Nmax()
       << std::endl;
  
    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);
  
        os
          << " " << "index"
          << " " << std::setw(width) << subspace_index
          << " " << "dim"
          << " " << std::setw(width) << subspace.dimension()
          << " " << "labels"
          << " " << subspace.LabelStr()
          << " " << std::endl;
      }
  
    return os.str();
  
  }


  OscillatorOrbitalSectors::OscillatorOrbitalSectors(
      const OscillatorOrbitalSpace& space,
      int L0, int g0,
      basis::SectorDirection sector_direction
    )
    : BaseSectors(space), L0_(L0), g0_(g0)
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
          allowed &= am::AllowedTriangle(ket_subspace.l(),L0,bra_subspace.l());
          allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);

          // push sector
          if (allowed)
            PushSector(bra_subspace_index,ket_subspace_index);
        }
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
