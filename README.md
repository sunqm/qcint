qcint (quick libcint) 
======================

An optimized libcint branch for X86 platform

version 3.0.16

2019-03-23


What is qcint
-------------
(GPL licensed) qcint is a branch of [libcint](https://github.com/sunqm/libcint.git)
library.  It provides exactly the same APIs as libcint.  However, the code
is optimized against SSE3, AVX, AVX2 and AVX-512 instructions.  On x86_64 platform, qcint can
be 5 ~ 50% faster than libcint.  Please refer to libcint for
more details of the features and installation instructions.


Bug report
----------
Qiming Sun <osirpt.sun@gmail.com>

