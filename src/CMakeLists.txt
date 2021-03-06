cmake_minimum_required ( VERSION 3.0 )
cmake_policy ( VERSION 3.0 )

project ( maxwell CXX )
enable_testing ()




# RPATH settings
set ( CMAKE_SKIP_BUILD_RPATH  FALSE )
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
list ( FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir )
if ( "${isSystemDir}" STREQUAL "-1" )
  set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )
endif ( "${isSystemDir}" STREQUAL "-1" )




# TBB
include_directories( ${CMAKE_SOURCE_DIR}/../../TBB/include )
link_directories( ${CMAKE_SOURCE_DIR}/../../TBB/build/linux_intel64_gcc_cc8.3.0_libc2.28_kernel4.19.0_release )




# CXX LIBRARY
add_library ( maxwell SHARED SimManager.cxx UnsteadyMaxwellKernel.cxx Mesh3D.cxx Cell.cxx DataKeeper.cxx DataProcessor.cxx Vec3D.cxx ScalarField.cxx VectorField.cxx Matrix.cxx MatrixFEM.cxx ElectrostaticSimManager.cxx ElectrostaticKernel.cxx ElectrostaticDataKeeper.cxx )
target_link_libraries(maxwell tbb)
install (
  TARGETS maxwell
  LIBRARY DESTINATION lib
  )




# TESTS
add_executable( test_matrix test_matrix.cxx )
target_link_libraries( test_matrix maxwell )

add_executable( test_vec3D test_vec3D.cxx )
target_link_libraries( test_vec3D maxwell )




# REFERENCE DOCUMENTATION
find_program ( DOXYGEN_EXECUTABLE
  NAMES doxygen
  DOC /usr/ 
  HINTS $ENV{DOXYGEN_ROOT}
  PATH_SUFFIXES bin
)

FIND_PROGRAM (DOXYGEN_DOT_EXECUTABLE
  NAMES dot
  DOC /usr/bin
  PATHS $ENV{DOT_HOME}
)

add_custom_target ( reference_doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/Doxyfile
)




# PYTHON LIBRARY
find_package ( PythonInterp 3.0 REQUIRED )
find_package ( pybind11 REQUIRED )

#include_directories ( /chemin/vers/pybind11/pybind11.h/si/installation/hors/path/std )
pybind11_add_module ( _maxwell MODULE Interface.cxx )
target_link_libraries ( _maxwell PUBLIC maxwell )
install (
  TARGETS _maxwell
  LIBRARY DESTINATION lib/python${PYTHON_VERSION}
  )


# Local Variables:
# mode: cmake
# coding: utf-8
# End:
