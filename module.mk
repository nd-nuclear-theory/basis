$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := indexing
module_units_cpp-h := many_body lsjt_scheme lsjt_operator jjjt_scheme jjjt_operator
module_units_cpp-h += jjjpnorb_scheme jjjpnorb_operator
# module_units_f := 
module_programs_cpp := lsjt_scheme_test lsjt_operator_test jjjt_scheme_test
module_programs_cpp += jjjpnorb_scheme_test jjjpnorb_operator_test
# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
