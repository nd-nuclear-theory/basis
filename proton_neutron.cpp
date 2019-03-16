/****************************************************************
  proton_neutron.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "proton_neutron.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // single-particle orbitals
  ////////////////////////////////////////////////////////////////

  // notational definitions for orbital species
  const std::array<HalfInt,2> kOrbitalSpeciesPNCodeTz({HalfInt(+1,2),HalfInt(-1,2)});
  const std::array<int,2> kOrbitalSpeciesPNCodeDecimal({1,2});
  const std::array<const char*,2> kOrbitalSpeciesPNCodeChar({"p","n"});

  const std::unordered_map<HalfInt, OrbitalSpeciesPN> kTzCodeOrbitalSpeciesPN({
      {HalfInt(+1,2),OrbitalSpeciesPN::kP},
      {HalfInt(-1,2),OrbitalSpeciesPN::kN}
    });
  const std::unordered_map<int, OrbitalSpeciesPN> kIndexCodeOrbitalSpeciesPN({
      {0,OrbitalSpeciesPN::kP},
      {1,OrbitalSpeciesPN::kN}
    });
  const std::unordered_map<int, OrbitalSpeciesPN> kDecimalCodeOrbitalSpeciesPN({
      {1,OrbitalSpeciesPN::kP},
      {2,OrbitalSpeciesPN::kN}
    });
  const std::unordered_map<std::string, OrbitalSpeciesPN> kCharCodeOrbitalSpeciesPN({
    {"p",OrbitalSpeciesPN::kP},
    {"n",OrbitalSpeciesPN::kN}
    });

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJpn scheme with general orbitals
  ////////////////////////////////////////////////////////////////

  // notational definitions for two-body state species -- static array definitions
  const std::array<int,3> kTwoBodySpeciesPNCodeTz({+1,-1,0});
  const std::array<int,3> kTwoBodySpeciesPNCodeDecimal({11,22,12});
  const std::array<const char*,3> kTwoBodySpeciesPNCodeChar({"pp","nn","pn"});

  ////////////////////////////////////////////////////////////////
  // proton-neutron operator types
  ////////////////////////////////////////////////////////////////

  // notational definitions for proton-neutron operator types
  const std::array<const char*,3> kOperatorTypePNCodeChar({"p","n","total"});


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
