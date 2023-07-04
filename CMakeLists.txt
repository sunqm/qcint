cmake_minimum_required (VERSION 3.5)
project (qcint C)
set(qcint_VERSION_MAJOR "5")
set(qcint_VERSION_MINOR "4")
set(qcint_VERSION_PATCH "0")
set(qcint_VERSION_TWEAK "0")
set(qcint_VERSION "${qcint_VERSION_MAJOR}.${qcint_VERSION_MINOR}.${qcint_VERSION_PATCH}")
set(qcint_SOVERSION "${qcint_VERSION_MAJOR}")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 -DNDEBUG")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()
set(CMAKE_VERBOSE_MAKEFILE OFF)
option(BUILD_MARCH_NATIVE "gcc flag -march=native" on)
if(BUILD_MARCH_NATIVE)
  message("Build with -march=native")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -O2")
else()
  # This appears to be the lowest version that compiles
  message("Build with -march=sse3")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3 -O2")
endif()
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -unroll-aggressive -ipo ")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-math-errno")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -ftree-vectorize -ffast-math")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-math-errno")
endif()

cmake_policy(SET CMP0042 NEW)
include(GNUInstallDirs)
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(QUADMATH)
if(QUADMATH_FOUND)
  set(HAVE_QUADMATH_H 1 CACHE INTERNAL "Have QUADMATH")
endif()

if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.12)
  cmake_policy(SET CMP0075 NEW)
endif()

include(CheckFunctionExists)
set(CMAKE_REQUIRED_INCLUDES math.h)
set(CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(expl HAVE_EXPL)
check_function_exists(sqrtl HAVE_SQRTL)
#set(CMAKE_REQUIRED_INCLUDES quadmath.h)
#set(CMAKE_REQUIRED_LIBRARIES quadmath)
#check_function_exists(fabsq HAVE_QUADMATH_H)

include_directories(
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/src)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/src/config.h")

configure_file(
  "${PROJECT_SOURCE_DIR}/include/cint.h.in"
  "${PROJECT_BINARY_DIR}/include/cint.h")

set(cintSrc
  src/breit.c src/c2f.c src/cart2sph.c src/cint1e.c src/cint2c2e.c
  src/cint2e.c src/cint3c1e.c src/cint3c2e.c src/cint_bas.c src/fblas.c
  src/g1e.c src/g2e.c src/g2e_simd1.c src/g3c1e.c src/g3c2e.c
  src/gout2e.c src/misc.c src/optimizer.c
  src/fmt.c src/rys_wheeler.c src/eigh.c src/rys_roots.c src/polyfits.c
  src/cint1e_a.c src/cint3c1e_a.c
  src/cint1e_grids.c src/g1e_grids.c
  src/autocode/breit1.c src/autocode/dkb.c src/autocode/gaunt1.c
  src/autocode/grad1.c src/autocode/grad2.c src/autocode/hess.c
  src/autocode/int3c1e.c src/autocode/int3c2e.c src/autocode/intor1.c
  src/autocode/intor2.c src/autocode/intor3.c src/autocode/intor4.c
  src/autocode/int1e_grids1.c)

#option(WITH_RANGE_COULOMB "Range separated Coulomb" on)
if(WITH_RANGE_COULOMB)
# defined in config.h
#  add_definitions(-DWITH_RANGE_COULOMB)
  message("Enabled WITH_RANGE_COULOMB")
endif(WITH_RANGE_COULOMB)

if(WITH_COULOMB_ERF)
  set(cintSrc ${cintSrc} src/g2e_coulerf.c src/cint2e_coulerf.c)
  add_definitions(-DWITH_COULOMB_ERF)
  message("Enabled WITH_COULOMB_ERF")
endif(WITH_COULOMB_ERF)

if(WITH_F12)
  set(cintSrc ${cintSrc} src/g2e_f12.c src/stg_roots.c src/cint2e_f12.c)
  add_definitions(-DWITH_F12)
  message("Enabled WITH_F12")
endif(WITH_F12)

if(WITH_GTG)
  set(cintSrc ${cintSrc} src/g2e_gtg.c src/cint2e_gtg.c src/cint3c2e_gtg.c
    src/cint2c2e_gtg.c)
  add_definitions(-DWITH_GTG)
  message("Enabled WITH_GTG")
  message("Enabled WITH_GTG. Note there are bugs in gtg type integrals")
endif(WITH_GTG)

if(PYPZPX)
  add_definitions(-DPYPZPX)
  message("P orbitals convention (py, pz, px)")
else()
  message("P orbitals convention (px, py, pz)")
endif(PYPZPX)

if(MIN_EXPCUTOFF)
  message("Minimal cutoff in exponent for integral screening ${MIN_EXPCUTOFF}")
  add_definitions(-DMIN_EXPCUTOFF=${MIN_EXPCUTOFF})
endif()

if(KEEP_GOING)
  message("Do not trigger hard exit for numerical issues in Rys quadrature")
  add_definitions(-DKEEP_GOING)
endif()

option(WITH_FORTRAN "Fortran interface" on)
if(WITH_FORTRAN)
  add_definitions(-DWITH_FORTRAN)
  message("Enable Fortran interface")
else()
  message("Exclude Fortran interface")
endif(WITH_FORTRAN)

option(WITH_CINT2_INTERFACE "Old libcint (version 2) interface" on)
if(WITH_CINT2_INTERFACE)
  add_definitions(-DWITH_CINT2_INTERFACE)
  message("Enable old cint (version 2) interface")
else()
  message("Exclude old cint (version 2) interface")
endif(WITH_CINT2_INTERFACE)

option(BUILD_SHARED_LIBS "build shared libraries" 1)
option(ENABLE_STATIC "Enforce static library build" 0)
if(ENABLE_STATIC)
  set(BUILD_SHARED_LIBS 0)
endif()

add_library(cint ${cintSrc})
set_target_properties(cint PROPERTIES
  VERSION ${qcint_VERSION}
  SOVERSION ${qcint_SOVERSION}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
if(QUADMATH_FOUND)
  target_link_libraries(cint quadmath)
endif()
target_link_libraries(cint "-lm")

set(CintHeaders
  ${PROJECT_SOURCE_DIR}/include/cint_funcs.h
  ${PROJECT_BINARY_DIR}/include/cint.h)

install(TARGETS cint DESTINATION "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" COMPONENT "lib")
install(FILES ${CintHeaders} DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR} COMPONENT "dev")
