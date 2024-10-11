cmake_minimum_required(VERSION 3.23)

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

option(BASIS_ENABLE_HASH "use unordered_map for lookup tables" ON)
option(BASIS_ENABLE_BOOST_HASH "use Boost hash for unordered_map" ON)
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
  if(BASIS_ENABLE_BOOST_HASH AND NOT TARGET Boost::headers)
    find_package(Boost REQUIRED COMPONENTS headers)
  endif()
elseif(BASIS_ENABLE_BOOST_HASH)
  message(WARNING "BASIS_ENABLE_BOOST_HASH without BASIS_ENABLE_HASH does nothing.")
endif()

# ##############################################################################
# define headers and sources
# ##############################################################################

# define units
set(${PROJECT_NAME}_UNITS_H
    basis type_traits map subspace state space sector
    degenerate hypersector operator product
    jt_operator
)
set(${PROJECT_NAME}_UNITS_H_CPP
    many_body
    proton_neutron
    lsjt_scheme
    lsjt_operator
    jjjt_scheme
    jjjt_operator
    jjjttz_scheme
    jjjttz_operator
    jjjpn_scheme
    jjjpn_operator
    nlj_orbital
    nlj_operator
    m_scheme
    oscillator_orbital
)

# construct lists of headers and sources
set(${PROJECT_NAME}_HEADERS ${${PROJECT_NAME}_UNITS_H} ${${PROJECT_NAME}_UNITS_H_CPP})
list(TRANSFORM ${PROJECT_NAME}_HEADERS PREPEND "${PROJECT_NAME}/")
list(TRANSFORM ${PROJECT_NAME}_HEADERS APPEND .h)
set(${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_UNITS_H_CPP})
list(TRANSFORM ${PROJECT_NAME}_SOURCES PREPEND "${PROJECT_NAME}/")
list(TRANSFORM ${PROJECT_NAME}_SOURCES APPEND .cpp)

# define basis library and add source/header files
add_library(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
target_sources(
  ${PROJECT_NAME} PUBLIC
  FILE_SET HEADERS
  TYPE HEADERS
  FILES ${${PROJECT_NAME}_HEADERS}
)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

# ##############################################################################
# configure optional features
# ##############################################################################

if(BASIS_ENABLE_HASH)
  target_compile_definitions(${PROJECT_NAME} PUBLIC BASIS_HASH)
  if(BASIS_ENABLE_BOOST_HASH)
    target_compile_definitions(${PROJECT_NAME} PUBLIC BASIS_BOOST_HASH)
  else()
    target_compile_definitions(${PROJECT_NAME} PUBLIC BASIS_STD_HASH)
  endif()
endif()
if(BASIS_ENABLE_DEPRECATED)
  target_compile_definitions(${PROJECT_NAME} PUBLIC BASIS_ALLOW_DEPRECATED)
endif()

# ##############################################################################
# link dependencies
# ##############################################################################

target_link_libraries(${PROJECT_NAME} PUBLIC m am::am Eigen3::Eigen mcutils::mcutils)
if(BASIS_ENABLE_HASH)
  target_link_libraries(${PROJECT_NAME} PUBLIC Boost::headers)
endif()

# ##############################################################################
# define installation rules
# ##############################################################################
include(GNUInstallDirs)

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  FILE_SET HEADERS
)

install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}Targets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)

# install the configuration file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

export(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake"
)

# ##############################################################################
# define tests
# ##############################################################################

set(${PROJECT_NAME}_UNITS_TEST
    degenerate_test
    product_test
    hypersector_test
    lsjt_scheme_test
    lsjt_operator_test
    jjjt_scheme_test
    jjjttz_scheme_test
    jjjttz_operator_test
    jjjpn_scheme_test
    jjjpn_operator_test
    m_scheme_test
    nlj_orbital_test
    oscillator_orbital_test
)

add_custom_target(${PROJECT_NAME}_tests)
foreach(test_name IN LISTS ${PROJECT_NAME}_UNITS_TEST)
  add_executable(${test_name} EXCLUDE_FROM_ALL test/${test_name}.cpp)
  target_link_libraries(${test_name} ${PROJECT_NAME}::${PROJECT_NAME})
  add_dependencies(${PROJECT_NAME}_tests ${test_name})
endforeach()

if(BASIS_MASTER_PROJECT)
  add_custom_target(tests)
  add_dependencies(tests ${PROJECT_NAME}_tests)
endif()
