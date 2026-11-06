cmake_minimum_required (VERSION 3.16)
project (qcint C)
set(qcint_VERSION_MAJOR "6")
set(qcint_VERSION_MINOR "2")
set(qcint_VERSION_PATCH "0")
set(qcint_VERSION_TWEAK "0")
set(qcint_VERSION "${qcint_VERSION_MAJOR}.${qcint_VERSION_MINOR}.${qcint_VERSION_PATCH}")
set(cint_SOVERSION "${qcint_VERSION_MAJOR}")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 -DNDEBUG")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
  message("Build type RELWITHDEBINFO")
else()
  message("Build type ${CMAKE_BUILD_TYPE}")
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

include(GNUInstallDirs)
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(QUADMATH)
if(QUADMATH_FOUND)
  set(HAVE_QUADMATH_H 1 CACHE INTERNAL "Have QUADMATH")
endif()

include(CheckFunctionExists)
set(CMAKE_REQUIRED_INCLUDES math.h)
set(CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(expl HAVE_EXPL)
check_function_exists(sqrtl HAVE_SQRTL)
#set(CMAKE_REQUIRED_INCLUDES quadmath.h)
#set(CMAKE_REQUIRED_LIBRARIES quadmath)
#check_function_exists(fabsq HAVE_QUADMATH_H)

include(CheckNativeVectorization)
check_native_vectorization()
message("sse3 works: ${SSE3_WORKS}")
message("avx works: ${AVX_WORKS}")
message("avx512f works: ${AVX512F_WORKS}")
if(AVX512F_WORKS)
  set(QCINT_VECTOR_LEVEL "avx512f")
  set(QCINT_SIMDD 8)
elseif(AVX_WORKS)
  set(QCINT_VECTOR_LEVEL "avx")
  set(QCINT_SIMDD 4)
elseif(SSE3_WORKS)
  set(QCINT_VECTOR_LEVEL "sse3")
  set(QCINT_SIMDD 2)
else()
  message(FATAL_ERROR "SSE3 is the minimum supported instruction set")
endif()

message("Will compile with vectorization width: ${QCINT_VECTOR_LEVEL}, SIMDD=${QCINT_SIMDD}")

configure_file(
  "${PROJECT_SOURCE_DIR}/src/cint_config.h.in"
  "${PROJECT_BINARY_DIR}/src/cint_config.h")

# I8 and CACHE_SIZE_I8 are consumed by #cmakedefine in include/cint.h.in
option(I8 "Use int64_t instead of int for FINT, the integer type of the API (atm, bas, shls, natm, nbas, ...)" off)
if(I8)
  message("Use 64-bit integers (int64_t) in the API")
endif()

# CACHE_SIZE_T falls back to FINT, so it is int64_t whenever I8 is set
include(CMakeDependentOption)
cmake_dependent_option(CACHE_SIZE_I8 "Use int64_t for CACHE_SIZE_T, the return type of the integral functions, without changing FINT" off "NOT I8" on)
if(CACHE_SIZE_I8)
  message("Use 64-bit integers (int64_t) for the cache size returned by integral functions")
endif()

configure_file(
  "${PROJECT_SOURCE_DIR}/include/cint.h.in"
  "${PROJECT_BINARY_DIR}/include/cint.h")

set(cintSrc
  src/breit.c src/c2f.c src/cart2sph.c src/cint1e.c src/cint2c2e.c
  src/cint2e.c src/cint3c1e.c src/cint3c2e.c src/cint_bas.c src/fblas.c
  src/g1e.c src/g2e.c src/g2e_simd1.c src/g3c1e.c src/g3c2e.c
  src/gout2e.c src/misc.c src/optimizer.c
  src/fmt.c src/rys_wheeler.c src/eigh.c src/rys_roots.c src/find_roots.c
  src/polyfits.c
  src/cint1e_a.c src/cint3c1e_a.c
  src/cint1e_grids.c src/g1e_grids.c
  src/autocode/breit1.c src/autocode/dkb.c src/autocode/gaunt1.c
  src/autocode/grad1.c src/autocode/grad2.c src/autocode/hess.c
  src/autocode/int3c1e.c src/autocode/int3c2e.c src/autocode/intor1.c
  src/autocode/intor2.c src/autocode/intor3.c src/autocode/intor4.c
  src/autocode/deriv3.c src/autocode/int1e_grids1.c
  src/autocode/deriv4.c src/autocode/lresc.c)

#option(WITH_RANGE_COULOMB "Range separated Coulomb" on)
if(WITH_RANGE_COULOMB)
# defined in cint_config.h
#  add_definitions(-DWITH_RANGE_COULOMB)
#  message("Enabled WITH_RANGE_COULOMB")
endif(WITH_RANGE_COULOMB)

if(WITH_POLYNOMIAL_FIT)
  set(cintSrc ${cintSrc} src/sr_rys_polyfits.c)
  add_definitions(-DWITH_POLYNOMIAL_FIT)
  message("Enabled WITH_POLYNOMIAL_FIT")
endif(WITH_POLYNOMIAL_FIT)

if(WITH_COULOMB_ERF)
  message("WITH_COULOMB_ERF is deprecated since v6.0")
endif(WITH_COULOMB_ERF)

if(WITH_F12)
  set(cintSrc ${cintSrc} src/g2e_f12.c src/stg_roots.c src/cint2e_f12.c)
  add_definitions(-DWITH_F12)
  message("Enabled WITH_F12")
endif(WITH_F12)

if(WITH_4C1E)
  set(cintSrc ${cintSrc} src/cint4c1e.c src/g4c1e.c)
  message("Enabled WITH_4C1E")
endif(WITH_4C1E)

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

target_include_directories(cint
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  PRIVATE
  ${PROJECT_BINARY_DIR}/src
  ${PROJECT_SOURCE_DIR}/src
)

set_target_properties(cint PROPERTIES
  C_STANDARD 99
  VERSION ${qcint_VERSION}
  SOVERSION ${cint_SOVERSION}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
# PRIVATE: a shared libcint resolves these itself, so consumers of cint::cint
# shouldn't link them. (A static libcint still passes them on, as LINK_ONLY.)
if(QUADMATH_FOUND)
  target_link_libraries(cint PRIVATE quadmath)
endif()
target_link_libraries(cint PRIVATE m)

# for consistency between FetchContent / add_subdirectory
# and find_package, so one can always link to cint::cint.
add_library(cint::cint ALIAS cint)

set(CintHeaders
  ${PROJECT_SOURCE_DIR}/include/cint_funcs.h
  ${PROJECT_BINARY_DIR}/include/cint.h)

# install headers
install(FILES ${CintHeaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT "dev")

# install the library and associate the installed targets
# with export set "cint-targets"
install(TARGETS cint
        EXPORT cint-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# qcint is a drop-in replacement for libcint, so it is exported under the same
# package name (cint) and namespace (cint::) for find_package and pkg-config.

# generate cintTargets.cmake and prepend the namespace "cint::"
# to exported targets.
export(EXPORT cint-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/cintTargets.cmake
        NAMESPACE cint::
)

# install cintTargets.cmake and associate this file with
# the export set "cint-targets"
install(EXPORT cint-targets
        FILE cintTargets.cmake
        NAMESPACE cint::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cint
)

include(CMakePackageConfigHelpers)

# configure cintConfig.cmake
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/cintConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cintConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cint)

# generate cintConfigVersion.cmake
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cintConfigVersion.cmake"
    VERSION ${qcint_VERSION}
    COMPATIBILITY SameMajorVersion
)

# and install both of these two files...
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/cintConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/cintConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cint
)

# configure the pkg-config file...
configure_file(${PROJECT_SOURCE_DIR}/cmake/cint.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/cint.pc @ONLY)

# and install it.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/cint.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
