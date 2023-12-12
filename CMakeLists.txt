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

if(BASIS_MASTER_PROJECT AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)

option(BASIS_ENABLE_HASH "use Boost hash for unordered_map" ON)
option(BASIS_ENABLE_DEPRECATED "enable deprecated basis APIs" OFF)

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

if(BASIS_ENABLE_HASH)
  if(NOT TARGET Boost::headers)
    find_package(Boost REQUIRED COMPONENTS headers)
  endif()
endif()

# ##############################################################################
# define basis headers and sources
# ##############################################################################

# define units
set(basis_UNITS_H basis degenerate hypersector operator jt_operator)
set(basis_UNITS_H_CPP
    many_body
    proton_neutron
    lsjt_scheme
    lsjt_operator
    jjjt_scheme
    jjjt_operator
    jjjpn_scheme
    jjjpn_operator
    nlj_orbital
    nlj_operator
    m_scheme
    oscillator_orbital
)

# construct lists of headers and sources
set(basis_HEADERS ${basis_UNITS_H} ${basis_UNITS_H_CPP})
list(TRANSFORM basis_HEADERS APPEND .h)
set(basis_SOURCES ${basis_UNITS_H_CPP})
list(TRANSFORM basis_SOURCES APPEND .cpp)

# define basis library and add source/header files
add_library(basis ${basis_SOURCES})
add_library(basis::basis ALIAS basis)
foreach(header IN LISTS basis_HEADERS)
  target_sources(
    basis INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${header}>
                    $<INSTALL_INTERFACE:include/basis/${header}>
  )
endforeach()

# ##############################################################################
# configure optional features
# ##############################################################################

if(BASIS_ENABLE_HASH)
  target_compile_definitions(basis PUBLIC BASIS_HASH)
endif()
if(BASIS_ENABLE_DEPRECATED)
  target_compile_definitions(basis PUBLIC BASIS_ALLOW_DEPRECATED)
endif()

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(basis PUBLIC m am::am Eigen3::Eigen mcutils::mcutils)
if(BASIS_ENABLE_HASH)
  target_link_libraries(basis PUBLIC Boost::headers)
endif()

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
install(FILES ${basis_HEADERS} DESTINATION include/basis)

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

# ##############################################################################
# define tests
# ##############################################################################

set(basis_UNITS_TEST
    degenerate_test
    hypersector_test
    lsjt_scheme_test
    lsjt_operator_test
    jjjt_scheme_test
    jjjpn_scheme_test
    jjjpn_operator_test
    m_scheme_test
    nlj_orbital_test
    oscillator_orbital_test
)

add_custom_target(basis_tests)
foreach(test_name IN LISTS basis_UNITS_TEST)
  add_executable(${test_name} EXCLUDE_FROM_ALL ${test_name}.cpp)
  target_link_libraries(${test_name} basis::basis)
  add_dependencies(basis_tests ${test_name})
endforeach()

if(BASIS_MASTER_PROJECT)
  add_custom_target(tests)
  add_dependencies(tests basis_tests)
endif()
