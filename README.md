qcint (quick libcint) 
======================

An optimized libcint branch for X86 platform

version 6.2.0
2026-10-16


What is qcint
-------------
(GPL licensed) qcint is a branch of [libcint](https://github.com/sunqm/libcint.git)
library.  It provides exactly the same APIs as libcint.  However, the code
is optimized against SSE3, AVX, AVX2 and AVX-512 instructions.  On x86_64 platform, qcint can
be 5 ~ 50% faster than libcint.  Please refer to libcint for
more details of the features and installation instructions.

Like libcint, qcint can be built with 64-bit integers in the API
(`-DI8=1`), which makes `FINT` (the type of the `atm`, `bas`, `shls` and
`dims` arrays, `natm`, `nbas` and the return values) `int64_t` instead of
`int`, e.g. to match Fortran codes compiled with 8-byte default integers.
Callers must then pass 64-bit integer arrays; the library is not ABI
compatible with a default build. `-DCACHE_SIZE_I8=1` instead makes only
`CACHE_SIZE_T`, the return type of the integral functions (which return the
required cache size when `out` is NULL), `int64_t`; it is implied by
`-DI8=1`.


Bug report
----------
Qiming Sun <osirpt.sun@gmail.com>

