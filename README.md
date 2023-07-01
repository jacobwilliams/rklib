![rklib](media/rklib.png)
============

**rklib**: A modern Fortran library of fixed and variable-step Runge-Kutta solvers.

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/rklib.svg)](https://github.com/jacobwilliams/rklib/releases/latest)
[![CI Status](https://github.com/jacobwilliams/rklib/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/rklib/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/rklib/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/rklib)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/rklib)](https://github.com/jacobwilliams/rklib/commits/master)

### Description

**This is a work in progress!**

The focus of this library is single-step, explicit Runge-Kutta solvers for 1st order differential equations.

#### Novel features:

 * The library includes a wide range of both fixed and variable-step Runge-Kutta methods, from very low to very high order.
 * It is object-oriented and written in modern Fortran.
 * It allows for defining a variable-step size integrator with a custom-tuned step size selection method. See `stepsize_class` in the code.
 * The `real` kind is selectable via a compiler directive (`REAL32`, `REAL64`, or `REAL128`).
 * Integration to an event is also supported.

#### Available methods:

Method name | Type | Order | Number of Stages | Reference
--- | --- | --- | --- | ---
`euler`    | Fixed-step    | 1  | 1        | [Euler (1768)](https://archive.org/details/institutionescal020326mbp)
`midpoint` | Fixed-step    | 2  | 2        | ?
`heun`     | Fixed-step    | 2  | 2        | ?
`rk3`      | Fixed-step    | 3  | 3        | ?
`rk4`      | Fixed-step    | 4  | 4        | [Kutta (1901)](https://archive.org/stream/zeitschriftfrma12runggoog#page/n449/mode/2up)
`rks4`     | Fixed-step    | 4  | 4        | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rks5`     | Fixed-step    | 5  | 5        | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkb6`     | Fixed-step    | 6  | 7        | [Butcher (1963)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/40DFE501CAB781C9AAE1439B6B8F481A/S1446788700023387a.pdf/div-class-title-on-runge-kutta-processes-of-high-order-div.pdf)
`rk7`      | Fixed-step    | 7  | 9        | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rk8-10`   | Fixed-step    | 8  | 10       | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkcv8`    | Fixed-step    | 8  | 11       | [Cooper & Verner (1972)](https://epubs.siam.org/doi/abs/10.1137/0709037)
`rk8-12`   | Fixed-step    | 8  | 12       | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkz10`    | Fixed-step    | 10 | 16       | [Zhang (2019)](https://arxiv.org/abs/1911.00318)
`rkbs32`   | Variable-step | 3  | 4 (FSAL) | [Bogacki & Shampine (1989)](https://www.sciencedirect.com/science/article/pii/0893965989900797)
`rkf45`    | Variable-step | 4  | 6        | [Fehlberg (1969)](https://ntrs.nasa.gov/api/citations/19690021375/downloads/19690021375.pdf)
`rkck54`   | Variable-step | 5  | 6        | [Cash & Karp (1990)](http://www.elegio.it/mc2/rk/doc/p201-cash-karp.pdf)
`rkdp54`   | Variable-step | 5  | 7 (FSAL) | [Dormand & Prince (1980)](https://www.sciencedirect.com/science/article/pii/0771050X80900133?via%3Dihub)
`rkt54`    | Variable-step | 5  | 7 (FSAL) | [Tsitouras (2011)](https://www.sciencedirect.com/science/article/pii/S0898122111004706/pdf)
`rkc65`    | Variable-step | 6  | 9        | [Calvo (1990)](https://www.sciencedirect.com/science/article/pii/089812219090064Q)
`rktp64`   | Variable-step | 6  | 7        | [Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)
`rkv65e`   | Variable-step | 6  | 9        | [Verner (1994)](https://www.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT)
`rktp75`   | Variable-step | 7  | 9        | [Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)
`rkv76e`   | Variable-step | 7  | 10       | [Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)
`rkf78`    | Variable-step | 7  | 13       | [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)
`rkv78`    | Variable-step | 7  | 13       | [Verner (1978)](https://www.jstor.org/stable/2156853)
`rktp86`   | Variable-step | 8  | 12       | [Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)
`rkdp87`   | Variable-step | 8  | 13       | [Prince & Dormand (1981)](https://www.sciencedirect.com/science/article/pii/0771050X81900103)
`rkv87e`   | Variable-step | 8  | 13       | [Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)
`rkf89`    | Variable-step | 8  | 17       | [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)
`rkv89`    | Variable-step | 8  | 16       | [Verner (1978)](https://www.jstor.org/stable/2156853)
`rkt98a`   | Variable-step | 9  | 16       | [Tsitouras (2001)](https://www.sciencedirect.com/science/article/abs/pii/S0168927401000253)
`rkv98e`   | Variable-step | 9  | 16       | [Verner (1978)](https://www.jstor.org/stable/2156853)
`rkf108`   | Variable-step | 10 | 17       | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk108.txt)
`rkc108`   | Variable-step | 10 | 21       | [Curtis (1975)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10a(8)_2.pdf)
`rkf1210`  | Variable-step | 12 | 25       | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1210.txt)
`rko129`   | Variable-step | 12 | 29       | [Ono (2006)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK12/RKcoeff12h(9)_1.pdf)
`rkf1412`  | Variable-step | 14 | 35       | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1412.txt)

### Compiling

A `fmp.toml` file is provided for compiling `rklib` with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

```
fpm build --profile release
```

### Documentation

The latest API documentation for the `master` branch can be found [here](https://jacobwilliams.github.io/rklib/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### Notes

This code was split off from the [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit) in September 2022.

#### TODO

 Some things that might be added at some point:

  * Dense output (for the ones that have interpolants)
  * Add [Verner's](https://www.sfu.ca/~jverner/) "more robust" methods.
  * Add max number of steps optional input.
  * Add a fixed-step option to `stepsize_class`, so that the variable step ones can be run as fixed-step integrators.
  * Comprehensive accuracy testing of all the methods.

### License

The `rklib` source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/rklib/blob/master/LICENSE.md) (BSD-3).

### References

  * E. B. Shanks, "[Higher Order Approximations of Runge-Kutta Type](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)", NASA Technical Note, NASA TN D-2920, Sept. 1965.
  * E. B. Shanks, "[Solutions of Differential Equations by Evaluations of Functions](https://www.ams.org/journals/mcom/1966-20-093/S0025-5718-1966-0187406-1/S0025-5718-1966-0187406-1.pdf)" Math. Comp. 20 (1966).
  * E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control", [NASA TR R-2870](https://ntrs.nasa.gov/citations/19680027281), 1968.
  * E. Fehlberg, "[Low-order classical Runge-Kutta formulas with stepsize control and their application to some heat transfer problems](https://ntrs.nasa.gov/api/citations/19690021375/downloads/19690021375.pdf)", NASA Technical Report R-315, July 1, 1969.
  * J. H. Verner, "Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error", SIAM Journal on Numerical Analysis, 15(4), 772-790, 1978.
  * T. Feagin, "[High-Order Explicit Runge-Kutta Methods](https://sce.uhcl.edu/rungekutta/)"
  * J. C. Butcher, "[A history of Runge-Kutta methods](https://www.sciencedirect.com/science/article/abs/pii/0168927495001085)", Applied Numerical Mathematics, Volume 20, Issue 3, March 1996, Pages 247-260
  * J. C. Butcher, "[On Runge-Kutta Processes of High Order](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/40DFE501CAB781C9AAE1439B6B8F481A/S1446788700023387a.pdf)", Oct. 28, 1963.
  * G. E. Müllges & F. Uhlig, "Numerical Algorithms with Fortran", Springer, 1996.
  * K. Fox, "[Numerical Integration of the Equations of Motion of Celestial Mechanics](https://adsabs.harvard.edu/full/1984CeMec..33..127F)", Celestial Mechanics 33, p 127-142, 1984.
  * [Mathematics Source Library](http://www.mymathlib.com/diffeq/)
  * [Maple worksheets on the derivation of Runge-Kutta schemes](http://www.peterstone.name/Maplepgs/RKcoeff.html)
  * [Index of numerical integrators](http://ketch.github.io/numipedia/index.html)
  * J. Williams, [Fehlberg's Runge-Kutta Methods](https://degenerateconic.com/fehlbergs-runge-kutta-methods.html), Feb. 10, 2018.

### See also

* [FOODIE](https://github.com/Fortran-FOSS-Programmers/FOODIE)
* [FLINT](https://github.com/princemahajan/FLINT)
* [DDEABM](https://github.com/jacobwilliams/ddeabm)
* [DOP853](https://github.com/jacobwilliams/dop853)
* [DVODE](https://github.com/jacobwilliams/dvode)
* [libode](https://github.com/markmbaum/libode)