cmake_minimum_required(VERSION 3.10)

project(mpm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_definitions(-DBOOST_BIND_GLOBAL_PLACEHOLDERS)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# SET(CMAKE_BUILD_TYPE "Debug")

# General compile settings

IF (NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

# Set optimization flags based on build type
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")

# GNU specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -w")
endif()

# Intel specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# Clang specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
endif()

# Check if Ninja
set(USED_CMAKE_GENERATOR "${CMAKE_GENERATOR}" CACHE STRING "Expose CMAKE_GENERATOR" FORCE)
if (USED_CMAKE_GENERATOR MATCHES "Ninja")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
endif()

# CMake seems to have no way to enable/disable testing per subproject,
# so we provide an option similar to BUILD_TESTING, but just for MPM.
option(MPM_BUILD_TESTING "enable testing for mpm" ON)

# CMake option for building mpm as a shared library for CI
option(MPM_BUILD_LIB "enable libmpm" ON)

# Halo exchange
option(HALO_EXCHANGE "Enable halo exchange" OFF)

# CMake Modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Boost Archive
find_package(Boost REQUIRED COMPONENTS filesystem system)
include_directories(${BOOST_INCLUDE_DIRS})
link_libraries(${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# MPI
#find_package(MPI)
#if (MPI_FOUND)
#  if(HALO_EXCHANGE)
#    add_definitions(-DUSE_HALO_EXCHANGE)
#  endif()
#  add_definitions("-DUSE_MPI")
#  include_directories(${MPI_CXX_INCLUDE_DIRS})
#  link_libraries(${MPI_CXX_LIBRARIES})
#endif()

# HDF5
ENABLE_LANGUAGE(C)
find_package(HDF5 COMPONENTS CXX HL)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  link_libraries(${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})
  add_definitions(${HDF5_DEFINITIONS})
endif()

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

# TBB
find_package(TBB REQUIRED)
add_definitions(${TBB_DEFINITIONS})
include_directories(${TBB_INCLUDE_DIRS})
link_libraries(${TBB_LIBRARIES})
if(NOT EXISTS "${TBB_INCLUDE_DIRS}/tbb/global_control.h")
  message(FATAL_ERROR
    "MPM requires TBB global_control support so -p can bound solver "
    "parallelism (TBB 2017 preview / 2019 U4 or newer)")
endif()

# pthreads
find_package (Threads)
link_libraries (${CMAKE_THREAD_LIBS_INIT})

# MKL
find_package(MKL)
if (MKL_FOUND)
  if (NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    link_libraries(iomp5)
    include_directories(${MKL_INCLUDE_DIR})
    link_libraries(${MKL_LIBRARIES})
    add_definitions("-DUSE_MKL")
  endif()
endif()


# KaHIP
if (MPI_FOUND)
  if (NO_KAHIP)
    find_package(KaHIP)
  else()
    find_package(KaHIP REQUIRED)
  endif()
if (KAHIP_FOUND)
  add_definitions("-DUSE_GRAPH_PARTITIONING")
  include_directories(${KAHIP_INCLUDE_DIRS})
  link_libraries(${KAHIP_LIBRARIES})
endif()
endif()

# VTK
find_package(VTK REQUIRED)
if (VTK_FOUND)
  include_directories(${VTK_INCLUDE_DIRS})
  link_libraries(${VTK_LIBRARIES})
  add_definitions("-DUSE_VTK")
  set(mpm_vtk ${mpm_SOURCE_DIR}/src/io/vtk_writer.cc)
  message(STATUS "FOUND VTK, DIRS= ${VTK_INCLUDE_DIRS}")
endif()

# Partio
find_package(Partio)
if (PARTIO_FOUND)
  add_definitions("-DUSE_PARTIO")
  include_directories(${PARTIO_INCLUDE_DIRS})
  link_libraries(${PARTIO_LIBRARIES})
endif()

# Include directories
include_directories(BEFORE
  ${mpm_SOURCE_DIR}/include/
  ${mpm_SOURCE_DIR}/include/elements/
  ${mpm_SOURCE_DIR}/include/elements/2d
  ${mpm_SOURCE_DIR}/include/elements/3d
  ${mpm_SOURCE_DIR}/include/functions/
  ${mpm_SOURCE_DIR}/include/generators/
  ${mpm_SOURCE_DIR}/include/io/
  ${mpm_SOURCE_DIR}/include/loads_bcs/
  ${mpm_SOURCE_DIR}/include/material/
  ${mpm_SOURCE_DIR}/include/matrix/
  ${mpm_SOURCE_DIR}/include/particles/
  ${mpm_SOURCE_DIR}/include/mesh/
  ${mpm_SOURCE_DIR}/include/solvers/
  ${mpm_SOURCE_DIR}/external/
)

# mpm executable
SET(mpm_src
  ${mpm_SOURCE_DIR}/src/affine_transform.cc
  ${mpm_SOURCE_DIR}/src/cell.cc
  ${mpm_SOURCE_DIR}/src/element.cc
  ${mpm_SOURCE_DIR}/src/functions/functions.cc
  ${mpm_SOURCE_DIR}/src/functions/linear_function.cc
  ${mpm_SOURCE_DIR}/src/functions/sin_function.cc
  ${mpm_SOURCE_DIR}/src/geometry.cc
  ${mpm_SOURCE_DIR}/src/hdf5_particle.cc
  ${mpm_SOURCE_DIR}/src/io/io.cc
  ${mpm_SOURCE_DIR}/src/io/io_mesh.cc
  ${mpm_SOURCE_DIR}/src/io/logger.cc
  ${mpm_SOURCE_DIR}/src/io/partio_writer.cc
  ${mpm_SOURCE_DIR}/src/material.cc
  ${mpm_SOURCE_DIR}/src/matrix.cc
  ${mpm_SOURCE_DIR}/src/mpm.cc
  ${mpm_SOURCE_DIR}/src/node.cc
  ${mpm_SOURCE_DIR}/src/particle.cc
  ${mpm_SOURCE_DIR}/src/quadrature.cc
)

# add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc ${mpm_src} ${mpm_vtk})

# if(MPM_BUILD_LIB)
# set(CMAKE_SHARED_LIBRARY_PREFIX "")
# add_library(lmpm SHARED ${mpm_src} ${mpm_vtk})
# add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc)
# target_link_libraries(mpm lmpm ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
# target_include_directories(lmpm PRIVATE ${PYTHON_LIBRARIES})
# else()
# add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc ${mpm_src} ${mpm_vtk})
# endif()

# if(MPM_BUILD_LIB)
#   find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter Development)
#   find_package(Boost REQUIRED COMPONENTS python)

# message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
# message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
# message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
# message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
# message(STATUS "BOOST_INCLUDE_DIRS = ${BOOST_INCLUDE_DIRS}")
# message(STATUS "EIGEN3_INCLUDE_DIR = $ENV{EIGEN3_INCLUDE_DIR}")

# set(CMAKE_SHARED_LIBRARY_PREFIX "")
# add_library(lmpm SHARED ${mpm_SOURCE_DIR}/src/main_mpmdem.cc ${mpm_SOURCE_DIR}/src/mpm_wrapper.cc ${mpm_src} ${mpm_vtk})
# # add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc ${mpm_src} ${mpm_vtk} ${PYTHON_INCLUDE_DIRS})
# # add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc)
# target_link_libraries(lmpm ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} )
# target_include_directories( lmpm PRIVATE ${PYTHON_INCLUDE_DIRS} 
#   ${mpm_SOURCE_DIR}/include/
#   ${mpm_SOURCE_DIR}/include/elements/
#   ${mpm_SOURCE_DIR}/include/elements/2d
#   ${mpm_SOURCE_DIR}/include/elements/3d
#   ${mpm_SOURCE_DIR}/include/functions/
#   ${mpm_SOURCE_DIR}/include/generators/
#   ${mpm_SOURCE_DIR}/include/io/
#   ${mpm_SOURCE_DIR}/include/loads_bcs/
#   ${mpm_SOURCE_DIR}/include/material/
#   ${mpm_SOURCE_DIR}/include/matrix/
#   ${mpm_SOURCE_DIR}/include/particles/
#   ${mpm_SOURCE_DIR}/include/mesh/
#   ${mpm_SOURCE_DIR}/include/solvers/
#   ${mpm_SOURCE_DIR}/external/
#   ${mpm_SOURCE_DIR}/tests/include/)

# else()
# # add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc ${mpm_src} ${mpm_vtk})
# endif()

# 关闭 Python / boost_python 版本，只编译纯 C++ 的 mpm 可执行文件

add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc ${mpm_src} ${mpm_vtk})

target_include_directories(mpm PRIVATE 
  ${mpm_SOURCE_DIR}/include/
  ${mpm_SOURCE_DIR}/include/elements/
  ${mpm_SOURCE_DIR}/include/elements/2d
  ${mpm_SOURCE_DIR}/include/elements/3d
  ${mpm_SOURCE_DIR}/include/functions/
  ${mpm_SOURCE_DIR}/include/generators/
  ${mpm_SOURCE_DIR}/include/io/
  ${mpm_SOURCE_DIR}/include/loads_bcs/
  ${mpm_SOURCE_DIR}/include/material/
  ${mpm_SOURCE_DIR}/include/matrix/
  ${mpm_SOURCE_DIR}/include/particles/
  ${mpm_SOURCE_DIR}/include/mesh/
  ${mpm_SOURCE_DIR}/include/solvers/
  ${mpm_SOURCE_DIR}/external/
)

if(MPM_BUILD_TESTING)
  enable_testing()
  # SANISAND regression tests are a separate executable and are not linked
  # into the production mpm executable.
  add_executable(sanisand_material_test
    ${mpm_SOURCE_DIR}/tests/materials/sanisand_test.cc
    ${mpm_SOURCE_DIR}/src/material.cc
    ${mpm_SOURCE_DIR}/src/io/logger.cc)
  add_test(NAME sanisand_material_test COMMAND sanisand_material_test)

  add_executable(sanisand_low_pressure_driver
    ${mpm_SOURCE_DIR}/tests/materials/sanisand_low_pressure_driver.cc
    ${mpm_SOURCE_DIR}/src/material.cc
    ${mpm_SOURCE_DIR}/src/io/logger.cc)

  # Audited single-point mechanism comparison. Its CSV stdout is consumed by
  # the phase-lag manuscript workflow; the executable returns nonzero if the
  # shared path, constitutive states, or pressure/tangent checks fail.
  add_executable(cyclic_constitutive_comparison_driver
    ${mpm_SOURCE_DIR}/tests/materials/cyclic_constitutive_comparison_driver.cc
    ${mpm_SOURCE_DIR}/src/material.cc
    ${mpm_SOURCE_DIR}/src/io/logger.cc)
  add_test(NAME cyclic_constitutive_comparison_driver
    COMMAND cyclic_constitutive_comparison_driver)

  add_executable(mohr_coulomb_handoff_test
    ${mpm_SOURCE_DIR}/tests/materials/mohr_coulomb_handoff_test.cc
    ${mpm_SOURCE_DIR}/src/material.cc
    ${mpm_SOURCE_DIR}/src/io/logger.cc)
  add_test(NAME mohr_coulomb_handoff_test COMMAND mohr_coulomb_handoff_test)

  if(HDF5_FOUND)
    add_executable(hdf5_particle_test
      ${mpm_SOURCE_DIR}/tests/io/hdf5_particle_test.cc
      ${mpm_SOURCE_DIR}/src/hdf5_particle.cc)
    add_test(NAME hdf5_particle_test COMMAND hdf5_particle_test)

    add_executable(hdf5_checkpoint_audit
      ${mpm_SOURCE_DIR}/src/hdf5_checkpoint_audit.cc
      ${mpm_SOURCE_DIR}/src/hdf5_particle.cc)
    add_executable(hdf5_checkpoint_audit_fixture
      ${mpm_SOURCE_DIR}/tests/io/hdf5_checkpoint_audit_fixture.cc
      ${mpm_SOURCE_DIR}/src/hdf5_particle.cc)
    foreach(audit_mode IN ITEMS valid duplicate nonfinite bad_porosity false_status
                                bad_nstate nonfinite_svar legacy legacy_allowed)
      if(audit_mode STREQUAL valid OR audit_mode STREQUAL legacy_allowed)
        set(audit_expect_success ON)
      else()
        set(audit_expect_success OFF)
      endif()
      add_test(
        NAME hdf5_checkpoint_audit_${audit_mode}
        COMMAND ${CMAKE_COMMAND}
          -DFIXTURE=$<TARGET_FILE:hdf5_checkpoint_audit_fixture>
          -DAUDITOR=$<TARGET_FILE:hdf5_checkpoint_audit>
          -DMODE=${audit_mode}
          -DEXPECT_SUCCESS=${audit_expect_success}
          -P ${mpm_SOURCE_DIR}/tests/io/hdf5_checkpoint_audit_test.cmake)
    endforeach()
  endif()

  add_executable(prescribed_pressure_validation_test
    ${mpm_SOURCE_DIR}/tests/solvers/prescribed_pressure_validation_test.cc)
  add_test(NAME prescribed_pressure_validation_test
    COMMAND prescribed_pressure_validation_test)

  # Three-DOF rigid-pipeline section and contact regression tests.
  add_executable(rigid_pipeline2d_test
    ${mpm_SOURCE_DIR}/tests/pipeline/rigid_pipeline2d_test.cc)
  add_test(NAME rigid_pipeline2d_test COMMAND rigid_pipeline2d_test)

  add_executable(free_surface_node_test
    ${mpm_SOURCE_DIR}/tests/mesh/free_surface_node_test.cc
    ${mpm_SOURCE_DIR}/src/affine_transform.cc)
  add_test(NAME free_surface_node_test COMMAND free_surface_node_test)

  add_executable(seabed_pressure_boundary_test
    ${mpm_SOURCE_DIR}/tests/particles/seabed_pressure_boundary_test.cc
    ${mpm_SOURCE_DIR}/src/affine_transform.cc)
  add_test(NAME seabed_pressure_boundary_test
    COMMAND seabed_pressure_boundary_test)

  add_executable(parallelism_test
    ${mpm_SOURCE_DIR}/tests/runtime/parallelism_test.cc
    ${mpm_SOURCE_DIR}/src/io/io.cc
    ${mpm_SOURCE_DIR}/src/io/logger.cc)
  target_compile_definitions(parallelism_test PRIVATE
    MPM_TEST_FIXTURE_DIR="${mpm_SOURCE_DIR}/tests/runtime/fixtures/")
  add_test(NAME parallelism_test COMMAND parallelism_test)

endif()
