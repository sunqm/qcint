import ctypes
import numpy
from rys_roots import *

def fp(a):
    a = numpy.asarray(a)
    return a.ravel().dot(numpy.cos(numpy.arange(a.size)))

SIMDD = 8
def cint_call(fname, nroots, x, low=None):
    cint = ctypes.CDLL('../build/libcint.so')
    buf = numpy.zeros(nroots*SIMDD*2+128)
    start = SIMDD - buf.ctypes.data % (SIMDD * 8) // buf.data.itemsize
    end = start + SIMDD*nroots
    r = numpy.ndarray((nroots,SIMDD), buffer=buf.data[start:end])
    start, end = end, end + SIMDD*nroots
    w = numpy.ndarray((nroots,SIMDD), buffer=buf.data[start:end])
    start, end = end, end + SIMDD
    x0 = numpy.ndarray(SIMDD, buffer=buf.data[start:end])
    x0[0] = x
    fun = getattr(cint, fname)
    if low is None:
        fun(ctypes.c_int(nroots),
            x0.ctypes.data_as(ctypes.c_void_p),
            r.ctypes.data_as(ctypes.c_void_p),
            w.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(1))
    else:
        start, end = end, end + SIMDD
        low0 = numpy.ndarray(SIMDD, buffer=buf.data[start:end])
        low0[0] = low
        fun(ctypes.c_int(nroots),
            x0.ctypes.data_as(ctypes.c_void_p),
            low0.ctypes.data_as(ctypes.c_void_p),
            r.ctypes.data_as(ctypes.c_void_p),
            w.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(1))
    return r[:,0], w[:,0]

def check_turnover_point(fun1, fun2):
    es = 2**numpy.arange(-2, 5, .5)
    turnover_points = []
    for i in range(1, 14):
        for x in es:
            f1 = fun1(x, i)
            f2 = fun2(x, i)
            diff = (f1[i] - f2[i])/f1[i]
            if diff < 2.2e-16:
                turnover_points.append(x)
                break
    return turnover_points

def test_turnover_point():
    import functools
    import mpmath
    mpmath.mp.dps = 15
    print(check_turnover_point(fmt1, fmt2))
    tps = []
    for low in 2**numpy.arange(-4, 5, .5):
        tp = check_turnover_point(functools.partial(fmt1_erfc, low=low),
                                  functools.partial(fmt2_erfc, low=low))
        print(low, tp)
        tps.append(tp)
    print(numpy.array(tps).max(axis=0))
    print(fmt1_erfc(0.9, 2, .6))
    print(fmt1_erfc(1.1, 2, .2))
    print(fmt2_erfc(1.1, 2, .2))
    print(fmt1_erfc(8.1, 2, .2))
    print(fmt2_erfc(8.1, 2, .2))
    print(fmt1(1.1, 2))
    print(fmt2(1.1, 2))
    print(fmt1(8.1, 2))
    print(fmt2(8.1, 2))

    print(abs(fmt(1, 2) - [0.746824132812427, 0.18947234582049235, 0.10026879814501737]).max())
    print(abs(fmt(7, 2) - [0.3349010581765593, 0.023856369729357483, 0.005046944801608424]).max())

def test_rys_roots_mpmath():
    numpy.random.seed(4)
    a = numpy.random.rand(4,4)
    a = a.T.dot(a)
    cs = R_dsmit(a, 4)
    print(fp(cs) - -1.1566333099933561)
    r, w = rys_roots_weights(6, .4)
    print(fp(r) - 2.898038115730877)
    print(fp(w) - 0.1193288112731978)

    print(fp(rys_roots_weights(1, .3)) - 0.93475892362307533)
    print(fp(rys_roots_weights(2, .8)) - 0.94983728587499205)
    print(fp(rys_roots_weights(3, .5)) - -2.7165192846520481)
    print(fp(rys_roots_weights(4, .5)) - -11.442155925744938)

    ff = fmt_erfc(.4, 6*2, 0.2)
    r, w = rys_roots_weights(6, .4, ff)
    print(fp(r) - 3.0557596858378053)
    print(fp(w) - 0.00070334763896715)

def test_rys_roots_weights():
    cint = ctypes.CDLL('../build/libcint.so')
    def check(nroots, x):
        r_ref, w_ref = rys_roots_weights(nroots, x)
        r, w = cint_call('CINTrys_roots', nroots, x)
        return np.array([abs(r-r_ref).max(), abs(w-w_ref).max()]).astype(float)

    failed = False
    es = 2**numpy.arange(-3, 7, .25)
    for i in range(1, 12):
        for x in es:
            diffs = check(i, x)
            if any(s > 1e-8 for s in diffs):
                print(i, x, diffs)
                failed |= any(s > 1e-5 for s in diffs)
    es = [50, 40, 36, 24, 20, 16, 12, 8, 4, 2, 0.5, 6e-8]
    for i in range(1, 5):
        for x in es:
            diffs = check(i, x)
            if any(s > 1e-8 for s in diffs):
                print(i, x, diffs)
                failed |= any(s > 1e-5 for s in diffs)
    if failed:
        print('test_rys_roots_weights .. failed')
    else:
        print('test_rys_roots_weights .. pass')


def test_rys_roots_weights_erfc():
    cint = ctypes.CDLL('../build/libcint.so')
    def check(nroots, x, lower):
        r_ref, w_ref = rys_roots_weights(nroots, x, lower)
        r, w = cint_call('CINTerfc_rys_roots', nroots, x, lower)
        return np.array([abs(r-r_ref).max(), abs(w-w_ref).max()]).astype(float)

    es = 2**numpy.arange(-3, 6, .25)
    failed = False
    for i in range(1, 12):
        for x in es:
            for low in [.1, .2, .3, .4, .5, .6, .7, .8, .9]:
                diffs = check(i, x, low)
                if any(s > 1e-7 for s in diffs):
                    print(i, x, low, diffs)
                    failed |= any(s > 1e-4 for s in diffs)
    if failed:
        print('test_rys_roots_weights_erfc .. failed')
    else:
        print('test_rys_roots_weights_erfc .. pass')

def test_stg_roots():
    cint = ctypes.CDLL('../build/libcint.so')
    def stg(nroots, t, u):
        r = numpy.zeros(nroots)
        w = numpy.zeros(nroots)
        cint.CINTstg_roots(ctypes.c_int(nroots),
                           ctypes.c_double(t), ctypes.c_double(u),
                           r.ctypes.data_as(ctypes.c_void_p),
                           w.ctypes.data_as(ctypes.c_void_p))
        return r, w

    print(fp(stg(1, 2.2, 1.5)) - 0.653262713484748 )
    print(fp(stg(2, 0.2, 8.5)) - 1.2362174210548105)
    print(fp(stg(4, 1.0, 0.5)) - -0.6907781084439245)


if __name__ == '__main__':
    # test_rys_roots_mpmath()
    # test_polyfit()
    # test_rys_roots_vs_polyfit()
    # test_rys_roots_weights()
    test_rys_roots_weights_erfc()
