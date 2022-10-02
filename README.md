**runge-kutta-fortran**: Fixed and variable-step Runge-Kutta solvers in Modern Fortran

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/runge-kutta-fortran.svg)](https://github.com/jacobwilliams/runge-kutta-fortran/releases/latest)
[![CI Status](https://github.com/jacobwilliams/runge-kutta-fortran/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/runge-kutta-fortran/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/runge-kutta-fortran/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/runge-kutta-fortran)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/runge-kutta-fortran)](https://github.com/jacobwilliams/runge-kutta-fortran/commits/master)

### Description


**This is a work in progress!**

Available methods:

Method name | Type | Order | Number of Points | Reference
--- | --- | --- | --- | ---
`rk4` | Fixed-step | 4 | 4 | [Kutta (1901)](https://archive.org/stream/zeitschriftfrma12runggoog#page/n449/mode/2up)
`rk7` | Fixed-step | 7 | 9 | [Shanks (1966)](https://www.ams.org/journals/mcom/1966-20-093/S0025-5718-1966-0187406-1/S0025-5718-1966-0187406-1.pdf)
`rk8-10` | Fixed-step | 8 | 10 | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkf78` | Variable-step | 7 | 13 | [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)
`rkf89`   | Variable-step | 8  | 17 | [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)
`rkv89`   | Variable-step | 8  | 16 | [Verner (1978)](https://www.jstor.org/stable/2156853)
`rkf108`  | Variable-step | 10 | 17 | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk108.txt)
`rkf1210` | Variable-step | 12 | 25 | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1210.txt)
`rkf1412` | Variable-step | 14 | 35 | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1412.txt)


### Compiling

A `fmp.toml` file is provided for compiling `runge-kutta-fortran` with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

```
fpm build --profile release
```

### Documentation

The latest API documentation for the `master` branch can be found [here](https://jacobwilliams.github.io/runge-kutta-fortran/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### Notes

This code was split off from the [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit) in September 2022.

### License

The `runge-kutta-fortran` source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/runge-kutta-fortran/blob/master/LICENSE.md) (BSD-3).

### References

  * E. B. Shanks, "[Higher Order Approximations of Runge-Kutta Type](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)", NASA Technical Note, NASA TN D-2920, Sept. 1965.
  * E. B. Shanks, "[Solutions of Differential Equations by Evaluations of Functions](https://www.ams.org/journals/mcom/1966-20-093/S0025-5718-1966-0187406-1/S0025-5718-1966-0187406-1.pdf)" Math. Comp. 20 (1966).
  * E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas with Stepsize Control", [NASA TR R-2870](https://ntrs.nasa.gov/citations/19680027281), 1968.
  * J. H. Verner, "Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error", SIAM Journal on Numerical Analysis, 15(4), 772-790, 1978.
  * T. Feagin, "[High-Order Explicit Runge-Kutta Methods](https://sce.uhcl.edu/rungekutta/)"
  * J. C. Butcher, "[A history of Runge-Kutta methods](https://www.sciencedirect.com/science/article/abs/pii/0168927495001085)", Applied Numerical Mathematics, Volume 20, Issue 3, March 1996, Pages 247-260


### See also

* [FOODIE](https://github.com/Fortran-FOSS-Programmers/FOODIE)
* [FLINT](https://github.com/princemahajan/FLINT)
* [DDEABM](https://github.com/jacobwilliams/ddeabm)
* [DOP853](https://github.com/jacobwilliams/dop853)