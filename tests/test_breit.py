'''
Compare the Breit integrals of qcint with those of libcint (through PySCF).

int2e_breit_r1p2_spinor read part of its output buffer uninitialized
(the stride between SIMD lanes in its gout function was nf*9 instead of
nf), so its results depended on what was in memory. Run with
MALLOC_PERTURB_ set, e.g. MALLOC_PERTURB_=85, so that such reads show up as
large errors rather than happening to look plausible.
'''

import os
import sys
import ctypes
import numpy
from pyscf import gto

_cint = ctypes.CDLL(os.path.abspath(os.path.join(__file__, '../../build/libcint.so')))

mol = gto.M(atom='O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692',
            basis='cc-pvtz', verbose=0)

def run(intor, thr=1e-12):
    fn = getattr(_cint, intor)
    atm = mol._atm.ctypes.data_as(ctypes.c_void_p)
    bas = mol._bas.ctypes.data_as(ctypes.c_void_p)
    env = mol._env.ctypes.data_as(ctypes.c_void_p)
    natm = ctypes.c_int(mol.natm)
    nbas = ctypes.c_int(mol.nbas)
    failed = False
    maxerr = 0
    for i in range(mol.nbas):
        for j in range(i+1):
            for k in range(mol.nbas):
                for l in range(k+1):
                    shls = [i, j, k, l]
                    ref = mol.intor_by_shell(intor, shls)
                    buf = numpy.empty_like(ref, order='F')
                    buf[:] = numpy.nan
                    fn(buf.ctypes.data_as(ctypes.c_void_p), None,
                       (ctypes.c_int*4)(*shls), atm, natm, bas, nbas, env,
                       None, None)
                    err = abs(buf - ref).max()
                    maxerr = max(maxerr, err)
                    if not err <= thr:
                        if not failed:
                            print(intor, 'shls', shls, 'max|err|', err)
                        failed = True
    print(intor, 'failed' if failed else 'pass', 'max|err| %.3g' % maxerr)
    return failed

failed = False
for intor in ('int2e_breit_r1p2_spinor', 'int2e_breit_r2p2_spinor',
              'int2e_breit_ssp1ssp2_spinor'):
    failed |= run(intor)
if failed:
    sys.exit(1)
