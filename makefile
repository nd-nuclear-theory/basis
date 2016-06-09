################################
# project configuration
################################

libname = basis

# modules -- header-only
modules_h = indexing

# modules -- header-plus-object 
##modules_ho = jt_scheme
modules_ho = lsjt_scheme lsjt_operator

# programs
##programs = jt_scheme_test
programs = lsjt_scheme_test lsjt_operator_test
##programs += write_lsjt_relative
CC := $(CXX)

CXXFLAGS = -std=c++11

# set flag for linking to FORTRAN
# need_fortran = 

################################
# common definitions
################################

COMMON_MAKE_DIR ?= .
include $(COMMON_MAKE_DIR)/common.mk

################################
# options and dependencies
################################

# program linking
CC := $(CXX)

# external libraries
LDLIBS +=  -lhalfint -lgsl

