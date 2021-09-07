cmake_minimum_required(VERSION 3.20)

# Determine if basis is built as a subproject (using add_subdirectory) or if it
# is the master project.
if(NOT DEFINED BASIS_MASTER_PROJECT)
  set(BASIS_MASTER_PROJECT OFF)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(BASIS_MASTER_PROJECT ON)
    message(STATUS "CMake version: ${CMAKE_VERSION}")
  endif()
endif()

project(basis CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)

option(BASIS_ENABLE_HASH "use Boost hash for unordered_map" ON)
option(BASIS_ENABLE_DEPRECATED "enable deprecated basis APIs" OFF)

# add headers and sources
set(basis_headers
    basis.h
    degenerate.h
    hypersector.h
    operator.h
    jt_operator.h
    many_body.h
    proton_neutron.h
    lsjt_scheme.h
    lsjt_operator.h
    jjjt_scheme.h
    jjjt_operator.h
    jjjpn_scheme.h
    jjjpn_operator.h
    nlj_orbital.h
    nlj_operator.h
    m_scheme.h
)
set(basis_sources
    jjjpn_operator.cpp
    jjjpn_scheme.cpp
    jjjt_operator.cpp
    jjjt_scheme.cpp
    lsjt_operator.cpp
    lsjt_scheme.cpp
    m_scheme.cpp
    many_body.cpp
    nlj_operator.cpp
    nlj_orbital.cpp
    proton_neutron.cpp
)

# define basis library and add source files
add_library(basis ${basis_sources})
add_library(basis::basis ALIAS basis)

# ##############################################################################
# find external projects/dependencies
# ##############################################################################

if(NOT TARGET Eigen3::Eigen)
  find_package(Eigen3 REQUIRED NO_MODULE)
endif()

if(NOT TARGET mcutils::mcutils)
  find_package(mcutils REQUIRED)
endif()

if(NOT TARGET am::am)
  find_package(am REQUIRED)
endif()

# ##############################################################################
# configure optional features
# ##############################################################################

if(BASIS_ENABLE_HASH)
  if(NOT TARGET Boost::headers)
    find_package(Boost REQUIRED)
  endif()
  target_compile_definitions(basis PUBLIC BASIS_HASH)
  target_link_libraries(basis PUBLIC Boost::headers)
endif()
if(BASIS_ENABLE_DEPRECATED)
  target_compile_definitions(basis BASIS_ALLOW_DEPRECATED)
endif()

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(basis PUBLIC m am::am Eigen3::Eigen mcutils::mcutils)

# ##############################################################################
# define include directory
# ##############################################################################
cmake_path(GET basis_SOURCE_DIR PARENT_PATH basis_PARENT_DIR)
target_include_directories(
  basis INTERFACE $<BUILD_INTERFACE:${basis_PARENT_DIR}>
                  $<INSTALL_INTERFACE:include>
)

# ##############################################################################
# define installation rules
# ##############################################################################
install(
  TARGETS basis
  DESTINATION lib
  EXPORT basisTargets
)
install(FILES ${basis_headers} DESTINATION include/basis)

install(
  EXPORT basisTargets
  NAMESPACE basis::
  FILE basisTargets.cmake
  DESTINATION lib/cmake/basis
)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/basisConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/basis"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/basisConfigVersion.cmake"
  VERSION "${basis_VERSION_MAJOR}.${basis_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)

# install the configuration file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/basisConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/basisConfigVersion.cmake
        DESTINATION lib/cmake/basis
)

export(
  EXPORT basisTargets
  NAMESPACE basis::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/basisTargets.cmake"
)
