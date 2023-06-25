![rklib](media/rklib.png)
============

**rklib**: A modern Fortran library of fixed and variable-step Runge-Kutta solvers.

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/rklib.svg)](https://github.com/jacobwilliams/rklib/releases/latest)
[![CI Status](https://github.com/jacobwilliams/rklib/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/rklib/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/rklib/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/rklib)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/rklib)](https://github.com/jacobwilliams/rklib/commits/master)

### Description

**This is a work in progress!**

#### Novel features:

 * The library includes a wide range of both fixed and variable-step Runge-Kutta methods, from very low to very high order.
 * It is object-oriented and written in modern Fortran.
 * It allows for defining a variable-step size integrator with a custom-tuned step size selection method. See `stepsize_class` in the code.
 * The `real` kind is selectable via a compiler directive (`REAL32`, `REAL64`, or `REAL128`).

#### Available methods:

Method name | Type | Order | Number of Points | Reference
--- | --- | --- | --- | ---
`euler`    | Fixed-step    | 1  | 1  |  [Euler (1768)](https://archive.org/details/institutionescal020326mbp)
`midpoint` | Fixed-step    | 2  | 2  | ?
`heun`     | Fixed-step    | 2  | 2  | ?
`rk3`      | Fixed-step    | 3  | 3  | ?
`rk4`      | Fixed-step    | 4  | 4  | [Kutta (1901)](https://archive.org/stream/zeitschriftfrma12runggoog#page/n449/mode/2up)
`rks4`     | Fixed-step    | 4  | 4  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rks5`     | Fixed-step    | 5  | 5  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rk7`      | Fixed-step    | 7  | 9  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rk8-10`   | Fixed-step    | 8  | 10 | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rk8-12`   | Fixed-step    | 8  | 12 | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rktp64`   | Variable-step | 6  | 7  | [Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)
`rktp75`   | Variable-step | 7  | 9  | [Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)
`rkf78`    | Variable-step | 7  | 13 | [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)
`rkv78`    | Variable-step | 7  | 13 | [Verner (1978)](https://www.jstor.org/stable/2156853)
`rktp86`   | Variable-step | 8  | 12 | [Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)
`rkdp87`   | Variable-step | 8  | 13 | [Prince & Dormand (1981)](https://www.sciencedirect.com/science/article/pii/0771050X81900103)
`rkf89`    | Variable-step | 8  | 17 | [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)
`rkv89`    | Variable-step | 8  | 16 | [Verner (1978)](https://www.jstor.org/stable/2156853)
`rkf108`   | Variable-step | 10 | 17 | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk108.txt)
`rkf1210`  | Variable-step | 12 | 25 | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1210.txt)
`rkf1412`  | Variable-step | 14 | 35 | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1412.txt)

### Compiling

A `fmp.toml` file is provided for compiling `rklib` with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

```
fpm build --profile release
```

### Documentation

The latest API documentation for the `master` branch can be found [here](https://jacobwilliams.github.io/rklib/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### Notes

This code was split off from the [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit) in September 2022.

### License

The `rklib` source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/rklib/blob/master/LICENSE.md) (BSD-3).

### References

  * E. B. Shanks, "[Higher Order Approximations of Runge-Kutta Type](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)", NASA Technical Note, NASA TN D-2920, Sept. 1965.
  * E. B. Shanks, "[Solutions of Differential Equations by Evaluations of Functions](https://www.ams.org/journals/mcom/1966-20-093/S0025-5718-1966-0187406-1/S0025-5718-1966-0187406-1.pdf)" Math. Comp. 20 (1966).
  * E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control", [NASA TR R-2870](https://ntrs.nasa.gov/citations/19680027281), 1968.
  * J. H. Verner, "Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error", SIAM Journal on Numerical Analysis, 15(4), 772-790, 1978.
  * T. Feagin, "[High-Order Explicit Runge-Kutta Methods](https://sce.uhcl.edu/rungekutta/)"
  * J. C. Butcher, "[A history of Runge-Kutta methods](https://www.sciencedirect.com/science/article/abs/pii/0168927495001085)", Applied Numerical Mathematics, Volume 20, Issue 3, March 1996, Pages 247-260
  * J. C. Butcher, "[On Runge-Kutta Processes of High Order](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/40DFE501CAB781C9AAE1439B6B8F481A/S1446788700023387a.pdf)", Oct. 28, 1963.
  * G. E. Müllges & F. Uhlig, "Numerical Algorithms with Fortran", Springer, 1996.
  * K. Fox, "[Numerical Integration of the Equations of Motion of Celestial Mechanics](https://adsabs.harvard.edu/full/1984CeMec..33..127F)", Celestial Mechanics 33, p 127-142, 1984.
  * [Mathematics Source Library](http://www.mymathlib.com/diffeq/)

### See also

* [FOODIE](https://github.com/Fortran-FOSS-Programmers/FOODIE)
* [FLINT](https://github.com/princemahajan/FLINT)
* [DDEABM](https://github.com/jacobwilliams/ddeabm)
* [DOP853](https://github.com/jacobwilliams/dop853)
* [DVODE](https://github.com/jacobwilliams/dvode)
