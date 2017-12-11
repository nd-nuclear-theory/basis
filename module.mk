$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# units
module_units_h += basis degenerate hypersector operator
module_units_cpp-h += many_body proton_neutron
module_units_cpp-h += lsjt_scheme lsjt_operator jjjt_scheme jjjt_operator
module_units_cpp-h += jjjpn_scheme jjjpn_operator nlj_orbital nlj_operator
module_units_cpp-h += m_scheme

# programs
module_programs_cpp += 

# test programs
module_programs_cpp_test += degenerate_test hypersector_test
module_programs_cpp_test += lsjt_scheme_test lsjt_operator_test jjjt_scheme_test
module_programs_cpp_test += jjjpn_scheme_test jjjpn_operator_test
module_programs_cpp_test += nlj_orbital_test m_scheme_test

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
