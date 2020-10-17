# Compilation du projet sous CMake
#
# 1. Créer un répertoire de construction (=builddir) hors du répertoire de sources (=srcdir):
#    mkdir build
#    cd build
#
# 2. Construire le projet (on référence le répertoire des sources en absolu ou en relatif):
#    cmake <srcdir>
#
#   Si on envisage d'installer le projet définitivement (parce que terminé) après la construction,
#   on précise tout de suite le futur répertoire d'installation (=installdir):
#   cmake <srcdir> -DCMAKE_INSTALL_PREFIX=<installdir>
#
#   Pour activer le mode optimisé ou debug, ajouter la variable CMAKE_BUILD_TYPE avec les valeurs
#   Release ou Debug ou RelWithDebInfo. Ex:
#   cmake <srcdir> -DCMAKE_INSTALL_PREFIX=<installdir> -DCMAKE_INSTALL_PREFIX=Debug
#
#   Ici:
#    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/install
#
# NOTA: cette opération n'est à faire qu'une seule fois (bootstrapping). Ensuite, chaque fichier est
# suivi par CMake et reconstruit automatiquement en cas de modification par la commande 'make' de l'étape 3.
#
# 3. Compiler le projet:
#    make
#
# 4. Installer le projet:
#    make install
#
# 5. Tester le projet:
#    make test
#
#   Ici:
#    make test PYTHONPATH=$PWD/install/lib/python
#
# NOTA: en temps normal il devrait y avoir une batterie de tests à lancer avant l'installation, puis une autre
#       à lancer après l'installation (obligatoire pour les test Python). Ici j'ai tout réuni par simplicité.
#
# That's all folks !
#
#
# Dans ton cas:
#  mkdir build
#  cd build
#  cmake .. -CMAKE_INSTALL_PREFIX=$PWD/install
#  make install
#  export PYTHONPATH=$PWD/install/lib/python
#  python3 ../scriptVtx.py
#


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


add_library ( maxwell SHARED SimManager.cxx UnsteadyMaxwellKernel.cxx SparseMatrix.cxx Mesh3D.cxx BoundaryCondition.cxx Cell.cxx Node.cxx DataKeeper.cxx )
install (
  TARGETS maxwell
  LIBRARY DESTINATION lib
  )


# test : nom=SimManager_Test commande=./SimManagerTest
# add_executable ( SimManagerTest SimManagerTest.cxx )
# target_link_libraries ( SimManagerTest vortex )
# add_test (
#   NAME SimManager_Test
#   COMMAND SimManagerTest )


#find_program ( DOXYGEN_EXECUTABLE
#  NAMES doxygen
#  DOC /usr/ 
#  HINTS $ENV{DOXYGEN_ROOT}
#  PATH_SUFFIXES bin
#)

#FIND_PROGRAM (DOXYGEN_DOT_EXECUTABLE
#  NAMES dot
#  DOC /usr/bin
#  PATHS $ENV{DOT_HOME}
#)

#add_custom_target ( reference_doc
#    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/Doxyfile
#)


# Détection des prérequis
find_package ( PythonInterp 3.0 REQUIRED )
find_package ( pybind11 REQUIRED )


# Construction du module Python (import _vortex) avec Pybind11
# NOTA: décommenter au besoin si installation de Pybind11 hors système
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
