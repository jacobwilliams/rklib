"""
Generate the dense output (interpolant) routines for the methods whose interpolants are
published as coefficient files.

This script writes the `<name>_dense_coefficients` and `<name>_dense_extra_stages`
routines for these methods into `src/rklib_variable_dense.f90`, replacing the block between the
`BEGIN GENERATED DENSE OUTPUT` and `END GENERATED DENSE OUTPUT` marker lines. The rest of that
file (e.g., the `rkdp54` and `rkdp85` routines) is written by hand and is not changed.

Usage (from the root of the repository):

    python scripts/generate_dense.py            # download (if needed), check, and update the source file
    python scripts/generate_dense.py --check    # also verify all the order conditions of the interpolants
    python scripts/generate_dense.py --stdout   # print the generated code instead of updating the source file

Only the Python standard library is needed.

Data sources
------------
1. Jim Verner's coefficient files ("Jim Verner's Refuge for Runge--Kutta Pairs",
   https://www.sfu.ca/~jverner/), which accompany:

       J.H. Verner, "Numerically optimal Runge--Kutta pairs with interpolants",
       Numerical Algorithms, 53, (2010) pp. 383-396. https://doi.org/10.1007/s11075-009-9290-3

   for `rkv65e`, `rkv65r`, `rkv76e`, `rkv76r`, `rkv87e`, `rkv87r`, `rkv98e`, `rkv98r`.
   Verner's site is no longer online, so the files are downloaded from the Internet Archive
   (specific snapshots). For the 7(6) pairs, the files are the July 2024 versions, in which
   Verner corrected the order 7 interpolants (the originally posted ones "satisfied some but
   not all conditions of order 7"); the step coefficients of these pairs are unchanged.

2. Peter Stone's Maple worksheets ("Maple worksheets on interpolation for Runge-Kutta schemes",
   http://www.peterstone.name/Maplepgs/RKinterp.html), in which he constructs interpolants for
   `rkdp87` (Prince-Dormand 8(7)), `rkev87` (Enright-Verner 8(7)), and `dverk78` (Maple's
   dverk78). The worksheets contain several intermediate results of the construction; the
   coefficient sets that are used (the final ones) are specified by their position in the
   worksheet. The coefficients are exact rationals. The worksheets are only parsed as data:
   only the numeric values of assignments like `a[3,2] = 1/16` are read (the expressions are
   checked to contain only numbers and arithmetic operators before they are evaluated).

   For `rkev87` and `dverk78`, the 13th stage of Stone's tableau is `f` at the end of the step,
   rather than the 13th stage of the method (which is only used for the error estimate, and is
   not needed by the interpolant), so the stages of the interpolant after the first 12 are
   stored after the 13 stages of the method.

   For `rkdp87`, rklib uses the rational coefficients published by Prince and Dormand, which
   are approximations (to about 18 digits) of the coefficients in Stone's worksheet.

All the files are downloaded into a cache directory (`build/dense_coefficients` by default,
which is not under version control), and each one is checked against a SHA-256 hash, so that
the generated code does not silently change if a different file is served.

What is generated
-----------------
For each method:

* `<name>_dense_extra_stages`: the additional function evaluations needed by the interpolant.
  If one of them is `f` at the end of the step (`c=1`, with `a` equal to the weights `b`),
  it is computed at the final state of the step and saved in the function evaluation cache
  (`set_fsal_cache`), so that the next step does not compute it again.

* `<name>_dense_coefficients`: the interpolant coefficients. The interpolants are published in
  the polynomial form, b_i(theta) = sum_k b[i,k] theta^k. For these coefficient sets, evaluating
  that form in double precision has a large roundoff error (up to ~1e-9 for `rkv87e`, because
  of cancellation between coefficients as large as ~1e6), so they are converted to the form used
  in DOPRI5/DOP853 (`dense_form_dopri`, see `dense_coefficients` in `src/rklib_module.F90`),
  which has roundoff errors 100-1000 times smaller. The conversion is exact (high precision
  arithmetic on the published values).

All the coefficients are written with 40 significant digits (enough for quadruple precision).

Checks
------
Before writing anything, the script verifies (in high precision arithmetic) that:

* the step coefficients in the coefficient file (the nodes, the `a` matrix, and the weights
  given by the interpolant at theta=1) are the same as the ones in rklib's step function
  (in `src/rklib_variable_steps.f90`), to within a tolerance for each method,
* each node `c[i]` equals the sum of row `i` of the `a` matrix, and
* the conversion to the DOPRI form is exact.

With `--check`, it also verifies that each interpolant satisfies all the (continuous) order
conditions up to its order (all the rooted trees, e.g., 486 trees for order 9).

The accuracy and order of the generated routines are tested in `test/rk_test_dense.f90`.
"""

import argparse
import hashlib
import os
import re
import sys
import time
import urllib.request
from decimal import Decimal as D, getcontext
from functools import lru_cache
from math import comb

getcontext().prec = 120

VERNER_URL = 'https://web.archive.org/web/{timestamp}id_/https://www.sfu.ca/~jverner/{file}'
VERNER_LINK = 'https://web.archive.org/web/https://www.sfu.ca/~jverner/{file}'  # for the comments
STONE_URL = 'http://www.peterstone.name/Maplepgs/Maple/nmthds/RKinterp/{file}'
STONE_PAGE = 'http://www.peterstone.name/Maplepgs/RKinterp.html'

# file name: (download URL, link for the comments, SHA-256)
def _verner(file, timestamp, sha256):
    return file, (VERNER_URL.format(timestamp=timestamp, file=file), VERNER_LINK.format(file=file), sha256)
def _stone(file, sha256):
    return file, (STONE_URL.format(file=file), STONE_URL.format(file=file), sha256)
FILES = dict([
    _verner('RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT', '2023',
            'a2ba2da870d50c43e351dab4729672774ca1f46b9002e680ed19716722b39a3a'),
    _verner('RKV65.IIIXb.Robust.00010102836.081204.CoeffsOnlyFLOAT', '2023',
            '75f849f81fd130329cda3bcab37b02d9d97307009ff798adfdbd6335e7e2e72d'),
    _verner('RKV76.IIa.Efficient.00001675585.240712.CoeffsOnlyFLOAT', '20250314130621',
            'b603d2cd21d58a61cb275c10a7cfb7e3e6e86231520f3162d834b3ab1b2ed4f8'),
    _verner('RKV76.IIa.Robust.000027015646.240702.CoeffsOnlyFLOAT', '20251208072439',
            'd9ee55c31c1ac88f9cf7ec1a842a2520716c81b79de9483cd1d8ec475942caa3'),
    _verner('RKV87.IIa.Efficient.000000282866.081208.CoeffsOnlyFLOAT', '2023',
            'ec0898bff0f5e1fcedfd1afc6776d307a28cebf1a7c7578193e56da54a9ae52d'),
    _verner('RKV87.IIa.Robust.00000754677.081208.CoeffsOnlyFLOAT', '2023',
            'fb91a5c967b58b71d894622d38adf64623205b0e95ce12cac70cddc1dfb7c1c1'),
    _verner('RKV98.IIa.Efficient.000000349.081209.FLOAT6040OnWeb', '2022',
            'b14441346716e86f916a1d58bdb3142abf0877a114a55663c8a50b0c3710638b'),
    _verner('RKV98.IIa.Efficient.000000349.081210.CoeffsOnlyRADandFLOATS', '2022',
            '1e9fafa7915000528ed951896113229b9248fc303cf5d0fdcdf9e67edd25e46d'),
    _verner('RKV98.IIa.Robust.000000351.081209.CoeffsOnlyFLOAT6040', '2022',
            '2dca4e3d0365392f17e90743527fc8ed035c960bd1fa8f9861d718c510c0a580'),
    _verner('RKV98.IIa.Robust.000000351.081209.CoeffsOnlyRADandFLOATS40', '2022',
            'fbf1d7c14521ae68e58ca951f54e3030069851420a58691237dac74484c7e1f4'),
    _stone('RKinterp8a.mws', '09394af3dd097ac5d7991313e5b1be106b55d43a445dd13d22781c754d978b20'),
    _stone('RKinterp8b.mws', 'd3d55b5d36a9d188edfb157737d7aff7566c5098bdbfd78e040608936cd72d88'),
    _stone('RKinterp8c.mws', '18213dcde5edef227e088a606cd04b3f849f9befbbf593dd8537f4deb09165ee'),
])

VERNER_REF = ('!  * J.H. Verner, "Numerically optimal Runge--Kutta pairs with interpolants",\n'
              '!    Numerical Algorithms, 53, (2010) pp. 383-396. [doi](https://doi.org/10.1007/s11075-009-9290-3)\n')
STONE_REF = (f'!  * P. Stone, [Maple worksheets on interpolation for Runge-Kutta schemes]({STONE_PAGE})\n')

# The methods:
#  * name, description
#  * source: 'verner' (coefficient files [the first one has the step coefficients], interpolant name), or
#            'stone' (worksheet, position of the stage coefficients set, position of the interpolant set,
#                     number of stages that are the same as the method's)
#  * tolerance for comparing the step coefficients with rklib's
METHODS = [
    dict(name='rkv65e', desc='Verner\'s "most efficient" 6(5)', source='verner',
         files=['RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT'], interp='bi6', tol='1e-35'),
    dict(name='rkv65r', desc='Verner\'s "most robust" 6(5)', source='verner',
         files=['RKV65.IIIXb.Robust.00010102836.081204.CoeffsOnlyFLOAT'], interp='bi6', tol='1e-35'),
    dict(name='rkv76e', desc='Verner\'s "most efficient" 7(6)', source='verner',
         files=['RKV76.IIa.Efficient.00001675585.240712.CoeffsOnlyFLOAT'], interp='bi7', tol='1e-35'),
    dict(name='rkv76r', desc='Verner\'s "most robust" 7(6)', source='verner',
         files=['RKV76.IIa.Robust.000027015646.240702.CoeffsOnlyFLOAT'], interp='bi7', tol='1e-35'),
    dict(name='rkv87e', desc='Verner\'s "most efficient" 8(7)', source='verner',
         files=['RKV87.IIa.Efficient.000000282866.081208.CoeffsOnlyFLOAT'], interp='bi8', tol='1e-35'),
    dict(name='rkv87r', desc='Verner\'s "most robust" 8(7)', source='verner',
         files=['RKV87.IIa.Robust.00000754677.081208.CoeffsOnlyFLOAT'], interp='bi8', tol='1e-35'),
    dict(name='rkv98e', desc='Verner\'s "most efficient" 9(8)', source='verner',
         files=['RKV98.IIa.Efficient.000000349.081209.FLOAT6040OnWeb',
                'RKV98.IIa.Efficient.000000349.081210.CoeffsOnlyRADandFLOATS'], interp='bi9', tol='1e-35'),
    dict(name='rkv98r', desc='Verner\'s "most robust" 9(8)', source='verner',
         files=['RKV98.IIa.Robust.000000351.081209.CoeffsOnlyFLOAT6040',
                'RKV98.IIa.Robust.000000351.081209.CoeffsOnlyRADandFLOATS40'], interp='bi9', tol='1e-35'),
    dict(name='rkdp87', desc='the Prince-Dormand 8(7) (RK8(7)13M)', source='stone',
         files=['RKinterp8c.mws'], stages_set=124, interp_set=125, n_same=13,
         tol='1e-16'),  # (rklib uses the published rational approximations of the coefficients)
    dict(name='rkev87', desc='the Enright-Verner 8(7)', source='stone',
         files=['RKinterp8b.mws'], stages_set=363, interp_set=329, n_same=12, tol='1e-35'),
    dict(name='dverk78', desc='the Verner 8(7) (Maple\'s `dverk78`)', source='stone',
         files=['RKinterp8a.mws'], stages_set=663, interp_set=665, n_same=12, tol='1e-35'),
]

SOURCE_FILE = 'src/rklib_variable_dense.f90'
STEPS_FILE = 'src/rklib_variable_steps.f90'
BEGIN_MARKER = '! BEGIN GENERATED DENSE OUTPUT'
END_MARKER = '! END GENERATED DENSE OUTPUT'
SEP = '!' + '*' * 89 + '\n'

################################################################################################
# Downloading the coefficient files

def download(cache, file, allow_changed):
    """Download a coefficient file into the cache directory (if it's not already there) and check its hash."""
    url, link, sha256 = FILES[file]
    path = os.path.join(cache, file)
    if not os.path.isfile(path):
        os.makedirs(cache, exist_ok=True)
        for attempt in range(6):  # the archive is sometimes busy or rate limited
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    data = r.read()
                break
            except OSError as e:
                print(f'  download failed ({e}), retrying...', file=sys.stderr)
                time.sleep(20 * (attempt + 1))
        else:
            sys.exit(f'Error: could not download {url}')
        if b'<html' in data[:1000].lower():
            sys.exit(f'Error: {url} returned an html page (not a coefficient file)')
        with open(path, 'wb') as f:
            f.write(data)
        print(f'  downloaded {file}', file=sys.stderr)
    with open(path, 'rb') as f:
        h = hashlib.sha256(f.read()).hexdigest()
    if h != sha256:
        msg = f'{path} has SHA-256 {h}, expected {sha256}'
        if not allow_changed:
            sys.exit(f'Error: {msg}\n  (delete it to download it again, or use --allow-changed)')
        print(f'  Warning: {msg}', file=sys.stderr)
    return path

################################################################################################
# Reading the coefficients

# Verner's files: assignments like `a[3,2] = .7669e-1` or `c[18]:=  .742` (the files use both forms)
VERNER_COEFF = re.compile(r'^\s*([A-Za-z]+[0-9]*)\s*\[\s*(\d+)\s*(?:,\s*(\d+)\s*)?\]\s*:?=\s*([-+]?[0-9.]+(?:[eE][-+]?\d+)?)\s*$')

def read_verner(paths):
    """Read the coefficients: {(name, i, j or None): published string}. The first file takes precedence."""
    d = {}
    for p in paths:
        with open(p, errors='replace') as f:
            for line in f:
                m = VERNER_COEFF.match(line.replace('\r', ''))
                if m:
                    n, i, j, v = m.groups()
                    d.setdefault((n, int(i), int(j) if j else None), v)
    return d

def stone_inputs(path):
    """The Maple input lines of a worksheet (`{MPLTEXT ...}` blocks)."""
    with open(path, errors='replace') as f:
        t = f.read()
    out = []
    for c in re.findall(r'\{MPLTEXT 1 0 \d+ "((?:[^"\\]|\\.)*)"', t):
        # the text is wrapped at a fixed width (the line breaks can be anywhere, even in numbers
        # and names), so remove them; the escaped "\n" are the actual line breaks of the input:
        c = c.replace('\r', '').replace('\n', '')
        c = c.replace('\\n', ' ').replace('\\+', '').replace('\\"', '"').replace('\\{', '{').replace('\\}', '}')
        c = c.replace('\\\\', '\\')
        out.append(re.sub(r'[ \t]+', ' ', c))
    return out

def numeric(expr):
    """The value of a numeric Maple expression (numbers, `+ - * / ^` and parentheses), or None."""
    e = expr.replace(' ', '')
    if not re.fullmatch(r'[-+*/^().0-9eE]+', e):
        return None
    e = e.replace('^', '**')
    e = re.sub(r'(?<![\w.])(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?', lambda m: f'D("{m.group(0)}")', e)
    try:
        return eval(e, {'D': D, '__builtins__': {}})  # (only numbers and operators, checked above)
    except Exception:
        return None

def stone_set(path, position):
    """The coefficients of the assignment `name := {x[i,j] = value, ...}` at a position (input number)
    in a worksheet: {(var, i, j or None): Decimal}."""
    c = stone_inputs(path)[position]
    m = re.match(r'\s*(\w+)\s*:=\s*\{(.*)\}\s*[:;]?\s*$', c)
    if not m:
        sys.exit(f'Error: input {position} of {path} is not a set of coefficients')
    d = {}
    for em in re.finditer(r'(\w+)\[\s*(\d+)\s*(?:,\s*(\d+)\s*)?\]\s*=\s*([^,{}]+?)(?=,\s*\w+\[|$)', m.group(2)):
        var, i, j, val = em.groups()
        v = numeric(val)
        if v is None:
            sys.exit(f'Error: could not evaluate {var}[{i},{j}] in input {position} of {path}')
        d[(var, int(i), int(j) if j else None)] = v
    return d

def rklib_tableau(name):
    """The step coefficients of a method from rklib's source (`aI`: nodes, `bIJ`: matrix, `cI`: weights)."""
    with open(STEPS_FILE) as f:
        src = f.read()
    blk = re.search(rf"    module procedure {name}\n(.*?)    end procedure {name}\n", src, re.S).group(1)
    env = {}
    for nm, expr in re.findall(r"parameter\s*::\s*(\w+)\s*=\s*([^\n!]+)", blk):
        e = expr.strip().replace('_wp', '')
        e = re.sub(r'(?<![\w.])(\d+\.?\d*|\.\d+)([eEdD][-+]?\d+)?',
                   lambda m: 'D("%s")' % m.group(0).replace('d', 'e').replace('D', 'e'), e)
        env[nm] = eval(e, {'D': D, '__builtins__': {}}, env)  # (rklib's own source)
    s = len(re.findall(r'call me%(?:f|check_fsal_cache|set_fsal_cache)\(', blk))  # number of stages
    A = {}
    for k, v in env.items():
        m = re.fullmatch(r'b(\d+)', k)
        if m:
            dg = m.group(1)
            ij = [(int(dg[:n]), int(dg[n:])) for n in range(1, len(dg)) if dg[n] != '0']
            ij = [(i, j) for i, j in ij if 2 <= i <= s and 1 <= j < i]
            assert len(ij) == 1, f'{name}: ambiguous coefficient name {k}'
            A[ij[0]] = v
    b = {int(k[1:]): v for k, v in env.items() if re.fullmatch(r'c\d+', k)}
    return s, A, b

################################################################################################
# Conversion and checks

def to_dopri(B):
    """Convert polynomial-form coefficients (row i: coefficients of theta^1..theta^deg in b_i)
    to the DOPRI form: x(t+theta*h) = x + theta*(r_1 + (1-theta)*(r_2 + theta*(r_3 + ...))),
    where r_j = h*sum_i R[i][j] f_i. The r_j term is multiplied by
    P_j = theta^ceil(j/2) (1-theta)^floor(j/2), whose theta^j coefficient is (-1)^floor(j/2),
    so the highest power is eliminated first."""
    deg = len(B[0])
    R = [[D(0)] * deg for _ in B]
    for i, row in enumerate(B):
        bt = list(row)
        for j in range(deg, 0, -1):
            p, q = (j + 1) // 2, j // 2
            pj = [D(0)] * (deg + 1)
            for m in range(q + 1):
                pj[p + m] = D(comb(q, m) * (-1) ** m)
            R[i][j - 1] = bt[j - 1] / pj[j]
            for k in range(1, j + 1):
                bt[k - 1] -= R[i][j - 1] * pj[k]
        assert all(abs(x) < D('1e-90') for x in bt), 'conversion to the DOPRI form is not exact'
    # (the values that are zero to within the roundoff of the high precision arithmetic are set to zero)
    return [[x if abs(x) > D('1e-90') else D(0) for x in row] for row in R]

_trees = {1: [()]}
def trees(n):
    """All rooted trees with n nodes (a tree is a sorted tuple of (order, subtree) pairs)."""
    if n not in _trees:
        res = set()
        def parts(rem, maxt, cur):
            if rem == 0:
                res.add(tuple(sorted(cur)))
                return
            for m in range(1, rem + 1):
                for t in trees(m):
                    if (m, t) <= maxt:
                        parts(rem - m, (m, t), cur + [(m, t)])
        parts(n - 1, (n, ()), [])
        _trees[n] = sorted(res)
    return _trees[n]

@lru_cache(None)
def gamma(t):
    g = 1 + sum(m for m, _ in t)
    for _, ch in t:
        g *= gamma(ch)
    return g

def check_order(A, B, p):
    """Max residual of the continuous order conditions sum_i B[i][k-1] Phi_i(t) = delta_{k,|t|}/gamma(t),
    for all trees t with up to p nodes."""
    n = len(A)
    phi_cache = {}
    def phi(t):
        if t not in phi_cache:
            v = [D(1)] * n
            for _, ch in t:
                pc = phi(ch)
                s = [sum((A[i][j] * pc[j] for j in range(i)), D(0)) for i in range(n)]
                v = [x * y for x, y in zip(v, s)]
            phi_cache[t] = v
        return phi_cache[t]
    worst = D(0)
    for order in range(1, p + 1):
        for t in trees(order):
            ph = phi(t)
            for k in range(1, len(B[0]) + 1):
                lhs = sum((B[i][k - 1] * ph[i] for i in range(n)), D(0))
                rhs = D(1) / gamma(t) if k == order else D(0)
                worst = max(worst, abs(lhs - rhs))
    return worst

def close(terms, value, tol):
    """sum(terms) == value, to within a tolerance relative to the magnitude of the terms."""
    return abs(sum(terms, D(0)) - value) <= tol * max(D(1), sum((abs(x) for x in terms), D(0)))

################################################################################################
# The data for each method

class Method:
    """The data for a method, from its source:

    * `s`: number of stages of rklib's step,
    * `A`, `B`: the tableau and the interpolant (polynomial form) of the source, with the source's
      numbering of the stages (`ns` stages; used for the order check),
    * `reg`: the register (column of `me%funcs`) of each stage of the source,
    * `nr`: the number of registers,
    * `extra`: the extra stages: (source stage, literal of c, [(source stage j, literal of a)]),
    * `b`: the weights of the step (for the end-of-step check)."""

def load(m, cache, allow_changed):
    paths = [download(cache, f, allow_changed) for f in m['files']]
    M = Method()
    M.name, M.desc = m['name'], m['desc']
    tol = D(m['tol'])
    s_rk, A_rk, b_rk = rklib_tableau(M.name)
    M.s = s_rk
    if m['source'] == 'verner':
        d = read_verner(paths)
        bi = m['interp']
        s = max(i for (n, i, _) in d if n == 'b')
        assert s == s_rk, f"{M.name}: rklib's step has {s_rk} stages, the coefficient file has {s}"
        ns = max(i for (n, i, _) in d if n == bi)
        deg = max(j for (n, _, j) in d if n == bi)
        M.p = int(bi[2:])
        M.A = [[D(d.get(('a', i, j), '0')) for j in range(1, ns + 1)] for i in range(1, ns + 1)]
        M.B = [[D(d.get((bi, i, k), '0')) for k in range(1, deg + 1)] for i in range(1, ns + 1)]
        M.c = [D(0)] + [D(d[('c', i, None)]) for i in range(2, ns + 1)]
        n_same = s
        lit_c = lambda i: lit(d[('c', i, None)])
        lit_a = lambda i, j: lit(d.get(('a', i, j), '0'))
        M.what = f'the interpolant of order {M.p} (`{bi}`) from the coefficient file'
        M.refs = VERNER_REF
    else:
        da = stone_set(paths[0], m['stages_set'])
        dd = stone_set(paths[0], m['interp_set'])
        ns = max(i for (_, i, _) in dd)
        deg = max(j for (_, _, j) in dd)
        M.p = deg  # (the order is verified with --check)
        M.A = [[da.get(('a', i, j), D(0)) for j in range(1, ns + 1)] for i in range(1, ns + 1)]
        M.B = [[dd.get(('d', i, k), D(0)) for k in range(1, deg + 1)] for i in range(1, ns + 1)]
        M.c = [D(0)] + [da.get(('c', i, None), sum(M.A[i-1], D(0))) for i in range(2, ns + 1)]
        n_same = m['n_same']
        lit_c = lambda i: lit40(M.c[i-1])
        lit_a = lambda i, j: lit40(M.A[i-1][j-1])
        M.what = f'an interpolant of order {M.p} constructed by P. Stone'
        M.refs = STONE_REF
    M.ns, M.deg = ns, deg
    M.links = [FILES[f][1] for f in m['files']]
    # stages 1..n_same of the source are the stages of the method; the other ones are stored after
    # the stages of the method:
    M.reg = {i: (i if i <= n_same else i + (M.s - n_same)) for i in range(1, ns + 1)}
    M.nr = max(M.reg.values())

    # the step coefficients must be the same as rklib's:
    for i in range(2, n_same + 1):
        for j in range(1, i):
            diff = abs(M.A[i-1][j-1] - A_rk.get((i, j), D(0)))
            assert diff <= tol * max(D(1), abs(A_rk.get((i, j), D(0)))), \
                f"{M.name}: a[{i},{j}] is not the same as rklib's ({float(diff):.2e})"
    # the weights (the interpolant at theta=1):
    M.b = ([D(d.get(('b', i, None), '0')) for i in range(1, ns + 1)] if m['source'] == 'verner'
           else [sum(row, D(0)) for row in M.B])
    for i in range(1, ns + 1):
        bw = b_rk.get(i, D(0)) if i <= n_same else D(0)
        assert close(M.B[i-1], bw, tol), f"{M.name}: the interpolant does not reproduce rklib's b[{i}] at theta=1"
    for i in range(n_same + 1, M.s + 1):
        # a stage of the method that is not in the source must not be used by the method's solution
        assert b_rk.get(i, D(0)) == 0, f"{M.name}: rklib's stage {i} is not in the coefficient file"
    # the nodes:
    for i in range(2, ns + 1):
        assert close(M.A[i-1], M.c[i-1], max(tol, D('1e-38'))), f'{M.name}: c[{i}] is not the row sum of a'

    # the extra stages (one of them may be f at the end of the step):
    M.extra = []
    M.end_stage = None
    for i in range(n_same + 1, ns + 1):
        if M.end_stage is None and M.c[i-1] == 1 and all(
                abs(M.A[i-1][j-1] - (M.b[j-1] if j <= n_same else D(0))) <= D('1e-38') * max(D(1), abs(M.b[j-1]))
                for j in range(1, i)):
            M.end_stage = i
        M.extra.append((i, lit_c(i), [(j, lit_a(i, j)) for j in range(1, i) if M.A[i-1][j-1] != 0]))
    return M

################################################################################################
# Fortran output

def lit(v):
    """Fortran literal from a published value (digits unchanged)."""
    s = v.strip().lower()
    if not re.search(r'[.e]', s):
        s += '.0'
    return s + '_wp'

def lit40(x):
    """Fortran literal with 40 significant digits."""
    if x == 0:
        return '0.0_wp'
    m, e = f'{x:.39e}'.split('e')
    return f'{m}e{int(e)}_wp'

def fortran(M):
    out = []
    w = out.append
    nextra = len(M.extra)
    links = ''.join(f'!    * [{os.path.basename(l)}]({l})\n' for l in M.links)
    R = to_dopri(M.B)

    # coefficients (in the order of the registers):
    w(f"""{SEP}!>
!  Dense output for {M.desc} method [[{M.name}]].
!  This is {M.what}, which requires {nextra} additional
!  function evaluations (see [[{M.name}_dense_extra_stages]]).
!
!@note The interpolant is published in the polynomial form. Here it has been
!      converted (exactly, using high precision arithmetic)
!      to the form used in `DOPRI5` and `DOP853` (see [[dense_coefficients]]), which has
!      much smaller roundoff error for these coefficients. This code is generated by
!      `scripts/generate_dense.py`.
!
!### References
{M.refs}!  * Coefficients:
{links}
    module procedure {M.name}_dense_coefficients

    form = dense_form_dopri
    allocate(b({M.nr},{M.deg}))
    b = 0.0_wp
""")
    for i in sorted(M.reg, key=lambda i: M.reg[i]):
        row = R[i - 1]
        if not any(x != 0 for x in row):
            continue
        r = M.reg[i]
        vals = (', &\n' + ' ' * (14 + len(str(r)))).join(lit40(x) for x in row)
        w(f'    b({r},:) = [ {vals} ]\n')
    w(f"""
    end procedure {M.name}_dense_coefficients
{SEP}
""")

    # extra stages:
    regs = sorted(M.reg[i] for i, _, _ in M.extra)
    note = ''
    if M.end_stage:
        note = (f"Stage {M.reg[M.end_stage]} is `f` at the end of the step,\n!  which is cached for the first stage "
                f"of the next step (see [[set_fsal_cache]]).\n")
    w(f"""{SEP}!>
!  The additional function evaluations for the dense output of [[{M.name}]].
!  These are stored in `me%funcs(:,{regs[0]}:{regs[-1]})`.{(' ' + note) if note else chr(10)}
    module procedure {M.name}_dense_extra_stages

""")
    for i, c, terms in M.extra:
        if i == M.end_stage:
            continue
        r = M.reg[i]
        w(f"    real(wp),parameter :: c{r} = {c}\n")
        for j, a in terms:
            w(f'    real(wp),parameter :: a{r:02d}{M.reg[j]:02d} = {a}\n')
        w('\n')
    w('    associate (t => me%dense_t, h => me%dense_h, x => me%dense_x, k => me%funcs)\n\n')
    for n_, (i, c, terms) in enumerate(M.extra):
        r = M.reg[i]
        last = n_ == len(M.extra) - 1
        if i == M.end_stage:
            w('        ! f at the end of the step (this is also the first stage of the next step,\n')
            w('        ! so it is saved in the function evaluation cache):\n')
            w(f'        call me%set_fsal_cache(t+h, me%dense_xf, k(:,{r}))\n')
        else:
            tt = [f'a{r:02d}{M.reg[j]:02d}*k(:,{M.reg[j]})' for j, a in terms]
            body = (' + &\n' + ' ' * 24).join(' + '.join(tt[m:m + 4]) for m in range(0, len(tt), 4))
            w(f'        call me%f(t+c{r}*h, x+h*({body}), k(:,{r}))\n')
        if not last:
            w('        if (me%stopped) return\n')
    w(f"""
    end associate

    end procedure {M.name}_dense_extra_stages
{SEP}
""")
    return ''.join(out)

def generate(cache, check, allow_changed):
    out = []
    for m in METHODS:
        M = load(m, cache, allow_changed)
        if check:
            res = check_order(M.A, M.B, M.p)
            assert res < D('1e-30'), f'{M.name}: the interpolant does not satisfy the order {M.p} conditions ({res:.2e})'
            print(f'  {M.name}: order {M.p} conditions satisfied (max residual {float(res):.1e})', file=sys.stderr)
        out.append(fortran(M))
        print(f'  {M.name}: {M.s} stages, {M.nr - M.s} extra stages (registers {M.s+1}..{M.nr}), '
              f'interpolant of degree {M.deg}', file=sys.stderr)
    return ''.join(out)

################################################################################################

def main():
    parser = argparse.ArgumentParser(description='Generate the dense output routines in '
                                                 f'{SOURCE_FILE} (see the docstring of this script).')
    parser.add_argument('--cache', default=os.path.join('build', 'dense_coefficients'),
                        help='directory for the downloaded coefficient files (default: %(default)s)')
    parser.add_argument('--check', action='store_true',
                        help='also verify all the order conditions of the interpolants')
    parser.add_argument('--stdout', action='store_true',
                        help='print the generated code instead of updating the source file')
    parser.add_argument('--allow-changed', action='store_true',
                        help='only warn (rather than stop) if a coefficient file does not match its SHA-256 hash')
    args = parser.parse_args()

    code = generate(args.cache, args.check, args.allow_changed)
    if args.stdout:
        sys.stdout.write(code)
        return
    with open(SOURCE_FILE) as f:
        src = f.read()
    i, j = src.find(BEGIN_MARKER), src.find(END_MARKER)
    if i < 0 or j < i:
        sys.exit(f'Error: the {BEGIN_MARKER} / {END_MARKER} lines were not found in {SOURCE_FILE}')
    i = src.index('\n', i) + 1  # keep the marker lines
    src = src[:i] + code + src[j:]
    with open(SOURCE_FILE, 'w') as f:
        f.write(src)
    print(f'updated {SOURCE_FILE}', file=sys.stderr)

if __name__ == '__main__':
    main()
