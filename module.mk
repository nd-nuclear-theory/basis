$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := indexing
module_units_cpp-h := lsjt_scheme lsjt_operator
# module_units_f := 
module_programs_cpp := lsjt_scheme_test lsjt_operator_test
# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
