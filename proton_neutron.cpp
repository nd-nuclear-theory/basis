/****************************************************************
  proton_neutron.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "proton_neutron.h"

namespace basis {

  ////////////////////////////////////////////////////////////////
  // two-body states in jjJpn scheme with general orbitals
  ////////////////////////////////////////////////////////////////

  // notational definitions for two-body state species -- static array definitions
  const std::array<int,3> kTwoBodySpeciesPNCodeTz({+1,-1,0});
  const std::array<int,3> kTwoBodySpeciesPNCodeDecimal({11,22,12});
  const std::array<const char*,3> kTwoBodySpeciesPNCodeChar({"pp","nn","pn"});

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
