$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# units
module_units_h += basis degenerate hypersector jt_operator map operator product sector space state subspace type_traits
module_units_cpp-h += jjjpn_operator jjjpn_scheme jjjt_operator jjjt_scheme jjjttz_operator jjjttz_scheme lsjt_operator lsjt_scheme many_body m_scheme nlj_operator nlj_orbital oscillator_orbital proton_neutron

# programs
module_programs_cpp +=

# test programs
## module_programs_cpp_test += 

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
