"""
Generate the dense output (interpolant) routines for the Verner methods.

This script writes the `<name>_dense_coefficients` and `<name>_dense_extra_stages`
routines for the eight Verner methods (`rkv65e`, `rkv65r`, `rkv76e`, `rkv76r`,
`rkv87e`, `rkv87r`, `rkv98e`, `rkv98r`) into `src/rklib_variable_dense.f90`,
replacing the block between the `BEGIN GENERATED VERNER DENSE OUTPUT` and
`END GENERATED VERNER DENSE OUTPUT` marker lines. The rest of that file
(e.g., the `rkdp54` and `rkdp85` routines) is written by hand and is not changed.

Usage (from the root of the repository):

    python scripts/generate_verner_dense.py            # download (if needed), check, and update the source file
    python scripts/generate_verner_dense.py --check    # also verify all the order conditions of the interpolants
    python scripts/generate_verner_dense.py --stdout   # print the generated code instead of updating the file

Only the Python standard library is needed.

Data source
-----------
The coefficients are from Jim Verner's coefficient files ("Jim Verner's Refuge for
Runge--Kutta Pairs", https://www.sfu.ca/~jverner/), which accompany:

    J.H. Verner, "Numerically optimal Runge--Kutta pairs with interpolants",
    Numerical Algorithms, 53, (2010) pp. 383-396. https://doi.org/10.1007/s11075-009-9290-3

Verner's site is no longer online, so the files are downloaded from the Internet Archive
(specific snapshots, listed in `METHODS`) into a cache directory (`build/verner_coefficients`
by default, which is not under version control). Each file is checked against a SHA-256 hash,
so that the generated code does not silently change if the archive serves a different file.
The files are only parsed as data (lines like `a[3,2] = .7669e-1`); nothing in them is executed.

For the 7(6) pairs, the files are the July 2024 versions, in which Verner corrected the order 7
interpolants (the originally posted ones "satisfied some but not all conditions of order 7").
The step coefficients of these pairs are unchanged from the original versions.

What is generated
-----------------
For each method, using the highest order interpolant in Verner's file (`bi6` ... `bi9`):

* `<name>_dense_extra_stages`: the additional function evaluations needed by the interpolant
  (stages `s+1` to `ns`, where `s` is the number of stages of the step), using Verner's
  published 40-digit coefficients as they are.

* `<name>_dense_coefficients`: the interpolant coefficients. Verner publishes these in the
  polynomial form, b_i(theta) = sum_k bi[i,k] theta^k. For these coefficient sets, evaluating
  that form in double precision has a large roundoff error (up to ~1e-9 for `rkv87e`, because of
  cancellation between coefficients as large as ~1e6), so they are converted to the form used in
  DOPRI5/DOP853 (`dense_form_dopri`, see `dense_coefficients` in `src/rklib_module.F90`), which has
  roundoff errors 100-1000 times smaller. The conversion is exact (60-digit arithmetic on the
  published 40-digit values), and the results are rounded once, to 40 digits.

Checks
------
Before writing anything, the script verifies (in 60-digit arithmetic) that:

* each node `c[i]` equals the sum of row `i` of the `a` matrix,
* the interpolant reproduces the step solution at theta=1 (sum_k bi[i,k] = b[i]), and
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

getcontext().prec = 60

ARCHIVE = 'https://web.archive.org/web/{timestamp}id_/https://www.sfu.ca/~jverner/{file}'
ARCHIVE_LINK = 'https://web.archive.org/web/https://www.sfu.ca/~jverner/{file}'  # for the comments

# (file name, archive snapshot timestamp, SHA-256)
FILES = {
    'RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT':
        ('2023', 'a2ba2da870d50c43e351dab4729672774ca1f46b9002e680ed19716722b39a3a'),
    'RKV65.IIIXb.Robust.00010102836.081204.CoeffsOnlyFLOAT':
        ('2023', '75f849f81fd130329cda3bcab37b02d9d97307009ff798adfdbd6335e7e2e72d'),
    'RKV76.IIa.Efficient.00001675585.240712.CoeffsOnlyFLOAT':
        ('20250314130621', 'b603d2cd21d58a61cb275c10a7cfb7e3e6e86231520f3162d834b3ab1b2ed4f8'),
    'RKV76.IIa.Robust.000027015646.240702.CoeffsOnlyFLOAT':
        ('20251208072439', 'd9ee55c31c1ac88f9cf7ec1a842a2520716c81b79de9483cd1d8ec475942caa3'),
    'RKV87.IIa.Efficient.000000282866.081208.CoeffsOnlyFLOAT':
        ('2023', 'ec0898bff0f5e1fcedfd1afc6776d307a28cebf1a7c7578193e56da54a9ae52d'),
    'RKV87.IIa.Robust.00000754677.081208.CoeffsOnlyFLOAT':
        ('2023', 'fb91a5c967b58b71d894622d38adf64623205b0e95ce12cac70cddc1dfb7c1c1'),
    'RKV98.IIa.Efficient.000000349.081209.FLOAT6040OnWeb':
        ('2022', 'b14441346716e86f916a1d58bdb3142abf0877a114a55663c8a50b0c3710638b'),
    'RKV98.IIa.Efficient.000000349.081210.CoeffsOnlyRADandFLOATS':
        ('2022', '1e9fafa7915000528ed951896113229b9248fc303cf5d0fdcdf9e67edd25e46d'),
    'RKV98.IIa.Robust.000000351.081209.CoeffsOnlyFLOAT6040':
        ('2022', '2dca4e3d0365392f17e90743527fc8ed035c960bd1fa8f9861d718c510c0a580'),
    'RKV98.IIa.Robust.000000351.081209.CoeffsOnlyRADandFLOATS40':
        ('2022', 'fbf1d7c14521ae68e58ca951f54e3030069851420a58691237dac74484c7e1f4'),
}

# (method, description, coefficient files [the first one has the step coefficients], interpolant)
METHODS = [
    ('rkv65e', '"most efficient" 6(5)', ['RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT'], 'bi6'),
    ('rkv65r', '"most robust" 6(5)',    ['RKV65.IIIXb.Robust.00010102836.081204.CoeffsOnlyFLOAT'], 'bi6'),
    ('rkv76e', '"most efficient" 7(6)', ['RKV76.IIa.Efficient.00001675585.240712.CoeffsOnlyFLOAT'], 'bi7'),
    ('rkv76r', '"most robust" 7(6)',    ['RKV76.IIa.Robust.000027015646.240702.CoeffsOnlyFLOAT'], 'bi7'),
    ('rkv87e', '"most efficient" 8(7)', ['RKV87.IIa.Efficient.000000282866.081208.CoeffsOnlyFLOAT'], 'bi8'),
    ('rkv87r', '"most robust" 8(7)',    ['RKV87.IIa.Robust.00000754677.081208.CoeffsOnlyFLOAT'], 'bi8'),
    ('rkv98e', '"most efficient" 9(8)', ['RKV98.IIa.Efficient.000000349.081209.FLOAT6040OnWeb',
                                         'RKV98.IIa.Efficient.000000349.081210.CoeffsOnlyRADandFLOATS'], 'bi9'),
    ('rkv98r', '"most robust" 9(8)',    ['RKV98.IIa.Robust.000000351.081209.CoeffsOnlyFLOAT6040',
                                         'RKV98.IIa.Robust.000000351.081209.CoeffsOnlyRADandFLOATS40'], 'bi9'),
]

SOURCE_FILE = 'src/rklib_variable_dense.f90'
BEGIN_MARKER = '! BEGIN GENERATED VERNER DENSE OUTPUT'
END_MARKER = '! END GENERATED VERNER DENSE OUTPUT'
SEP = '!' + '*' * 89 + '\n'
TOL = D('1e-38')  # relative tolerance for the consistency checks (the published values have 40 significant digits)

################################################################################################
# Downloading and reading the coefficient files

def download(cache, file, allow_changed):
    """Download a coefficient file into the cache directory (if it's not already there) and check its hash."""
    timestamp, sha256 = FILES[file]
    path = os.path.join(cache, file)
    if not os.path.isfile(path):
        os.makedirs(cache, exist_ok=True)
        url = ARCHIVE.format(timestamp=timestamp, file=file)
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

# assignments like `a[3,2] = .7669e-1` or `c[18]:=  .742` (the files use both forms)
COEFF = re.compile(r'^\s*([A-Za-z]+[0-9]*)\s*\[\s*(\d+)\s*(?:,\s*(\d+)\s*)?\]\s*:?=\s*([-+]?[0-9.]+(?:[eE][-+]?\d+)?)\s*$')

def read(paths):
    """Read the coefficients: {(name, i, j or None): published string}. The first file takes precedence."""
    d = {}
    for p in paths:
        with open(p, errors='replace') as f:
            for line in f:
                m = COEFF.match(line.replace('\r', ''))
                if m:
                    n, i, j, v = m.groups()
                    d.setdefault((n, int(i), int(j) if j else None), v)
    return d

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
        assert all(abs(x) < D('1e-45') for x in bt), 'conversion to the DOPRI form is not exact'
    return R

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

def generate(cache, check, allow_changed):
    out = []
    w = out.append
    for name, desc, files, bi in METHODS:
        paths = [download(cache, f, allow_changed) for f in files]
        d = read(paths)
        s = max(i for (n, i, _) in d if n == 'b')    # number of stages of the step
        ns = max(i for (n, i, _) in d if n == bi)    # number of stages used by the interpolant
        deg = max(j for (n, _, j) in d if n == bi)   # degree of the interpolant
        p = int(bi[2:])                              # order of the interpolant
        A = [[D(d.get(('a', i, j), '0')) for j in range(1, ns + 1)] for i in range(1, ns + 1)]
        B = [[D(d.get((bi, i, k), '0')) for k in range(1, deg + 1)] for i in range(1, ns + 1)]

        # consistency checks:
        # (the tolerance is relative to the magnitude of the terms, since the values are rounded to 40 digits)
        def close(terms, value):
            return abs(sum(terms, D(0)) - value) <= TOL * max(D(1), sum((abs(x) for x in terms), D(0)))
        for i in range(2, ns + 1):
            assert close(A[i - 1], D(d[('c', i, None)])), f'{name}: c[{i}] is not the row sum of a'
        for i in range(1, ns + 1):
            bi_ = D(d.get(('b', i, None), '0')) if i <= s else D(0)
            assert close(B[i - 1], bi_), f'{name}: the interpolant does not reproduce b[{i}] at theta=1'
        if check:
            res = check_order(A, B, p)
            assert res < D('1e-30'), f'{name}: the interpolant does not satisfy the order {p} conditions ({res:.2e})'
            print(f'  {name}: order {p} conditions satisfied (max residual {float(res):.1e})', file=sys.stderr)

        R = to_dopri(B)
        refs = ''.join(f'!    * [{f}]({ARCHIVE_LINK.format(file=f)})\n' for f in files)

        # coefficients:
        w(f"""{SEP}!>
!  Dense output for Verner's {desc} method [[{name}]].
!  This is the interpolant of order {p} (`{bi}`), which requires {ns - s} additional
!  function evaluations (see [[{name}_dense_extra_stages]]).
!
!@note Verner publishes these coefficients in the polynomial form. Here they have been
!      converted (exactly, using the published 40-digit values and 60-digit arithmetic)
!      to the form used in `DOPRI5` and `DOP853` (see [[dense_coefficients]]), which has
!      much smaller roundoff error for these coefficients. This code is generated by
!      `scripts/generate_verner_dense.py`.
!
!### References
!  * J.H. Verner, "Numerically optimal Runge--Kutta pairs with interpolants",
!    Numerical Algorithms, 53, (2010) pp. 383-396. [doi](https://doi.org/10.1007/s11075-009-9290-3)
!  * Coefficients (40 digits):
{refs}
    module procedure {name}_dense_coefficients

    form = dense_form_dopri
    allocate(b({ns},{deg}))
    b = 0.0_wp
""")
        for i in range(1, ns + 1):
            row = R[i - 1]
            if not any(x != 0 for x in row):
                continue
            vals = (', &\n' + ' ' * (14 + len(str(i)))).join(lit40(x) for x in row)
            w(f'    b({i},:) = [ {vals} ]\n')
        w(f"""
    end procedure {name}_dense_coefficients
{SEP}
""")

        # extra stages:
        w(f"""{SEP}!>
!  The additional function evaluations for the dense output of [[{name}]].
!  These are stored in `me%funcs(:,{s + 1}:{ns})`.

    module procedure {name}_dense_extra_stages

""")
        for i in range(s + 1, ns + 1):
            w(f"    real(wp),parameter :: c{i} = {lit(d[('c', i, None)])}\n")
            for j in range(1, i):
                v = d.get(('a', i, j), '0')
                if D(v) != 0:
                    w(f'    real(wp),parameter :: a{i:02d}{j:02d} = {lit(v)}\n')
            w('\n')
        w('    associate (t => me%dense_t, h => me%dense_h, x => me%dense_x, k => me%funcs)\n\n')
        for i in range(s + 1, ns + 1):
            terms = [f'a{i:02d}{j:02d}*k(:,{j})' for j in range(1, i) if D(d.get(('a', i, j), '0')) != 0]
            body = (' + &\n' + ' ' * 24).join(' + '.join(terms[m:m + 4]) for m in range(0, len(terms), 4))
            w(f'        call me%f(t+c{i}*h, x+h*({body}), k(:,{i}))\n')
            if i < ns:
                w('        if (me%stopped) return\n')
        w(f"""
    end associate

    end procedure {name}_dense_extra_stages
{SEP}
""")
        print(f'  {name}: {s} stages + {ns - s} extra, interpolant {bi} (degree {deg})', file=sys.stderr)
    return ''.join(out)

################################################################################################

def main():
    parser = argparse.ArgumentParser(description='Generate the dense output routines for the Verner methods '
                                                 f'in {SOURCE_FILE} (see the docstring of this script).')
    parser.add_argument('--cache', default=os.path.join('build', 'verner_coefficients'),
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
