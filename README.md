
![rklib](media/rklib.png)
============

**rklib**: A modern Fortran library of fixed and variable-step Runge-Kutta solvers.

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![GitHub release](https://img.shields.io/github/release/jacobwilliams/rklib.svg)](https://github.com/jacobwilliams/rklib/releases/latest)
[![CI Status](https://github.com/jacobwilliams/rklib/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/rklib/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/rklib/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/rklib)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/rklib)](https://github.com/jacobwilliams/rklib/commits/master)

### Description

**This is a work in progress!**

The focus of this library is single-step, explicit Runge-Kutta solvers for 1st order differential equations.

### Novel features:

  * The library includes a wide range of both fixed and variable-step Runge-Kutta methods, from very low to very high order.
  * It is object-oriented and written in modern Fortran.
  * It allows for defining a variable-step size integrator with a custom-tuned step size selection method. See `stepsize_class` in the code.
  * The `real` kind is selectable via a compiler directive (`REAL32`, `REAL64`, or `REAL128`).
  * Integration to an event is also supported. The root-finding method is also selectable (via the [roots-fortran](https://github.com/jacobwilliams/roots-fortran) library).

### Available Runge-Kutta methods:

  * Number of fixed-step methods:    27
  * Number of variable-step methods: 48
  * Total number of methods:         75

### Fixed-step methods:

Name       | Description| Properties | Order | Stages   | Registers | CFL  | Reference
---        | ---        | ---        | ---   | ---      | ---       | ---  | ---
`euler` | Euler |  | 1 | 1 | 1 | 1.0 | [Euler (1768)](https://archive.org/details/institutionescal020326mbp)
`midpoint` | Midpoint |  | 2 | 2 | 2 |  | ?
`heun` | Heun |  | 2 | 2 | 2 |  | ?
`rkssp22` | 2-stage, 2nd order TVD Runge-Kutta Shu-Osher | SSP | 2 | 2 | 1 | 1.0 | [Shu & Oscher (1988)](https://ntrs.nasa.gov/api/citations/19880014833/downloads/19880014833.pdf)
`rk3` | 3th order Runge-Kutta |  | 3 | 3 | 3 |  | ?
`rkssp33` | 3-stage, 3rd order TVD Runge-Kutta Shu-Osher | SSP | 3 | 3 | 1 | 1.0 | [Shu & Oscher (1988)](https://ntrs.nasa.gov/api/citations/19880014833/downloads/19880014833.pdf)
`rkssp53` | 5-stage, 3rd order SSP Runge-Kutta Spiteri-Ruuth | SSP | 3 | 5 | 2 | 2.65 | [Ruuth (2006)](https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf)
`rk4` | Classic 4th order Runge-Kutta |  | 4 | 4 | 4 |  | [Kutta (1901)](https://archive.org/stream/zeitschriftfrma12runggoog#page/n449/mode/2up)
`rks4` | 4th order Runge-Kutta Shanks |  | 4 | 4 | 4 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkr4` | 4th order Runge-Kutta Ralston |  | 4 | 4 | 4 |  | [Ralston (1962)](https://doi.org/10.1090%2FS0025-5718-1962-0150954-0)
`rkls44` | 4-stage, 4th order low storage non-TVD Runge-Kutta Jiang-Shu | LS | 4 | 4 | 2 |  | [Jiang and Shu (1988)](https://ntrs.nasa.gov/api/citations/19960007052/downloads/19960007052.pdf)
`rkls54` | 5-stage, 4th order low storage Runge-Kutta Carpenter-Kennedy | LS | 4 | 5 | 2 | 0.32 | [Carpenter & Kennedy (1994)](https://ntrs.nasa.gov/api/citations/19940028444/downloads/19940028444.pdf)
`rkssp54` | 5-stage, 4th order SSP Runge-Kutta Spiteri-Ruuth | SSP | 4 | 5 | 4 | 1.51 | [Ruuth (2006)](https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf)
`rks5` | 5th order Runge-Kutta Shanks |  | 5 | 5 | 5 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rk5` | 5th order Runge-Kutta |  | 5 | 6 | 6 |  | ?
`rkc5` | 5th order Runge-Kutta Cassity |  | 5 | 6 | 6 |  | [Cassity (1966)](https://epubs.siam.org/doi/10.1137/0703052)
`rkl5` | 5th order Runge-Kutta Lawson |  | 5 | 6 | 6 |  | [Lawson (1966)](https://epubs.siam.org/doi/abs/10.1137/0703051)
`rklk5a` | 5th order Runge-Kutta Luther-Konen 1 |  | 5 | 6 | 6 |  | [Luther & Konen (1965)](https://epubs.siam.org/doi/abs/10.1137/1007112)
`rklk5b` | 5th order Runge-Kutta Luther-Konen 2 |  | 5 | 6 | 6 |  | [Luther & Konen (1965)](https://epubs.siam.org/doi/abs/10.1137/1007112)
`rkb6` | 6th order Runge-Kutta Butcher |  | 6 | 7 | 7 |  | [Butcher (1963)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/40DFE501CAB781C9AAE1439B6B8F481A/S1446788700023387a.pdf/div-class-title-on-runge-kutta-processes-of-high-order-div.pdf)
`rk7` | 7th order Runge-Kutta Shanks |  | 7 | 9 | 9 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rk8_10` | 10-stage, 8th order Runge-Kutta Shanks |  | 8 | 10 | 10 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkcv8` | 11-stage, 8th order Runge-Kutta Cooper-Verner |  | 8 | 11 | 11 |  | [Cooper & Verner (1972)](https://epubs.siam.org/doi/abs/10.1137/0709037)
`rk8_12` | 12-stage, 8th order Runge-Kutta Shanks |  | 8 | 12 | 12 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkz10` | 10th order Runge-Kutta Zhang |  | 10 | 16 | 16 |  | [Zhang (2019)](https://arxiv.org/abs/1911.00318)
`rko10` | 10th order Runge-Kutta Ono |  | 10 | 17 | 17 |  | [Ono (2003)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10f_1.pdf)
`rkh10` | 10th order Runge-Kutta Hairer |  | 10 | 17 | 17 |  | [Hairer (1978)](https://www.researchgate.net/publication/31221486_A_Runge-Kutta_Method_of_Order_10)



### Variable-step methods:

Name       | Description| Properties | Order | Stages   | Registers | CFL  | Reference
---        | ---        | ---        | ---   | ---      | ---       | ---  | ---
`rkbs32` | Bogacki & Shampine 3(2) | FSAL | 3 | 4 | 4 |  | [Bogacki & Shampine (1989)](https://www.sciencedirect.com/science/article/pii/0893965989900797)
`rkssp43` | 4-stage, 3rd order SSP | SSP, LS | 3 | 4 | 2 | 2.0 | [Kraaijevanger (1991)](https://doi.org/10.1007/BF01933264), [Conde et al. (2018)](https://doi.org/10.48550/arXiv.1806.08693)
`rkf45` | Fehlberg 4(5) |  | 4 | 6 | 6 |  | [Fehlberg (1969)](https://ntrs.nasa.gov/api/citations/19690021375/downloads/19690021375.pdf)
`rkck54` | Cash & Karp 5(4) |  | 5 | 6 | 6 |  | [Cash & Karp (1990)](http://www.elegio.it/mc2/rk/doc/p201-cash-karp.pdf)
`rkdp54` | Dormand-Prince 5(4) | FSAL | 5 | 7 | 7 |  | [Dormand & Prince (1980)](https://www.sciencedirect.com/science/article/pii/0771050X80900133?via%3Dihub)
`rkt54` | Tsitouras 5(4) | FSAL | 5 | 7 | 7 |  | [Tsitouras (2011)](https://www.sciencedirect.com/science/article/pii/S0898122111004706/pdf)
`rks54` | Stepanov 5(4) | FSAL | 5 | 7 | 7 |  | [Stepanov (2022)](https://arxiv.org/pdf/2108.12590.pdf)
`rkpp54` | Papakostas-PapaGeorgiou 5(4) | FSAL | 5 | 7 | 7 |  | [Papakostas & Papageorgiou (1996)](https://www.jstor.org/stable/2153797)
`rkpp54b` | Papakostas-PapaGeorgiou 5(4) b | FSAL | 5 | 7 | 7 |  | [Papakostas & Papageorgiou (1996)](https://www.jstor.org/stable/2153797)
`rkbs54` | Bogacki & Shampine 5(4) |  | 5 | 8 | 8 |  | [Bogacki & Shampine (1996)](https://www.sciencedirect.com/science/article/pii/0898122196001411)
`rkss54` | Sharp & Smart 5(4) |  | 5 | 7 | 7 |  | [Sharp & Smart (1993)](https://epubs.siam.org/doi/10.1137/0914021)
`rkdp65` | Dormand-Prince 6(5) |  | 6 | 8 | 8 |  | [Dormand & Prince (1981)](https://www.sciencedirect.com/science/article/pii/0771050X81900103)
`rkc65` | Calvo 6(5) |  | 6 | 9 | 9 |  | [Calvo (1990)](https://www.sciencedirect.com/science/article/pii/089812219090064Q)
`rktp64` | Tsitouras & Papakostas NEW6(4) |  | 6 | 7 | 7 |  | [Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)
`rkv65e` | Verner efficient (9,6(5)) | FSAL | 6 | 9 | 9 |  | [Verner (1994)](https://www.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT)
`rkv65r` | Verner robust (9,6(5)) | FSAL | 6 | 9 | 9 |  | [Verner (1994)](https://www.sfu.ca/~jverner/RKV65.IIIXb.Robust.00010102836.081204.RATOnWeb)
`rkv65` | Verner 6(5) |  | 6 | 8 | 8 |  | [Verner (2006)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK6/RKcoeff6e_3.pdf)
`dverk65` | Verner 6(5) "DVERK" |  | 6 | 8 | 8 |  | Verner (?)
`rktf65` | Tsitouras & Famelis 6(5) | FSAL | 6 | 9 | 9 |  | [Tsitouras & Famelis (2006)](http://users.uoa.gr/~tsitourasc/ModifiedRK-ICNAAM2006.pdf)
`rktp75` | Tsitouras & Papakostas NEW7(5) |  | 7 | 9 | 9 |  | [Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)
`rktmy7` | 7th order Tanaka-Muramatsu-Yamashita |  | 7 | 10 | 10 |  | [Tanaka, Muramatsu & Yamashita (1992)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK7/RKcoeff7d_4.pdf)
`rktmy7s` | 7th order Stable Tanaka-Muramatsu-Yamashita |  | 7 | 10 | 10 |  | [Tanaka, Muramatsu & Yamashita (1992)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK7/RKcoeff7d_3.pdf)
`rkv76e` | Verner efficient (10:7(6)) |  | 7 | 10 | 10 |  | [Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)
`rkv76r` | Verner robust (10:7(6)) |  | 7 | 10 | 10 |  | [Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)
`rkss76` | Sharp & Smart 7(6) |  | 7 | 11 | 11 |  | [Sharp & Smart (1993)](https://epubs.siam.org/doi/10.1137/0914021)
`rkf78` | Fehlberg 7(8) |  | 7 | 13 | 13 |  | [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)
`rkv78` | Verner 7(8) |  | 7 | 13 | 13 |  | [Verner (1978)](https://www.jstor.org/stable/2156853)
`dverk78` | Verner "Maple" 7(8) |  | 7 | 13 | 13 |  | [Verner (?)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK8/RKcoeff8c_2.pdf)
`rkdp85` | Dormand-Prince 8(5) |  | 8 | 12 | 12 |  | [Hairer (1993)](https://github.com/jacobwilliams/dop853)
`rktp86` | Tsitouras & Papakostas NEW8(6) |  | 8 | 12 | 12 |  | [Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)
`rkdp87` | Dormand & Prince RK8(7)13M |  | 8 | 13 | 13 |  | [Prince & Dormand (1981)](https://www.sciencedirect.com/science/article/pii/0771050X81900103)
`rkv87e` | Verner efficient (8)7 |  | 8 | 13 | 13 |  | [Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)
`rkv87r` | Verner robust (8)7 |  | 8 | 13 | 13 |  | [Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)
`rkev87` | Enright-Verner (8)7 |  | 8 | 13 | 13 |  | [Enright (1993)](https://epubs.siam.org/doi/10.1137/0730074)
`rkk87` | Kovalnogov-Fedorov-Karpukhina-Simos-Tsitouras 8(7) |  | 8 | 13 | 13 |  | [Kovalnogov, Fedorov, Karpukhina, Simos, Tsitouras (2022)](https://www.researchgate.net/publication/363396601_Runge-Kutta_Embedded_Methods_of_Orders_87_for_Use_in_Quadruple_Precision_Computations)
`rkf89` | Fehlberg 8(9) |  | 8 | 17 | 17 |  | [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)
`rkv89` | Verner 8(9) |  | 8 | 16 | 16 |  | [Verner (1978)](https://www.jstor.org/stable/2156853)
`rkt98a` | Tsitouras 9(8) A |  | 9 | 16 | 16 |  | [Tsitouras (2001)](https://www.sciencedirect.com/science/article/abs/pii/S0168927401000253)
`rkv98e` | Verner efficient (16:9(8)) |  | 9 | 16 | 16 |  | [Verner (1978)](https://www.jstor.org/stable/2156853)
`rkv98r` | Verner robust (16:9(8)) |  | 9 | 16 | 16 |  | [Verner (1978)](https://www.jstor.org/stable/2156853)
`rks98` | Sharp 9(8) |  | 9 | 16 | 16 |  | [Sharp (2000)](https://www.hindawi.com/journals/ads/2000/853972/)
`rkf108` | Feagin 8(10) |  | 10 | 17 | 17 |  | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk108.txt)
`rkc108` | Curtis 10(8) |  | 10 | 21 | 21 |  | [Curtis (1975)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10a(8)_2.pdf)
`rkb109` | Baker 10(9) |  | 10 | 21 | 21 |  | [Baker (?)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10c_1.pdf)
`rks1110a` | Stone 11(10) |  | 11 | 26 | 26 |  | [Stone (2015)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK11/RKcoeff11_a.pdf)
`rkf1210` | Feagin 12(10) |  | 12 | 25 | 25 |  | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1210.txt)
`rko129` | Ono 12(9) |  | 12 | 29 | 29 |  | [Ono (2006)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK12/RKcoeff12h(9)_1.pdf)
`rkf1412` | Feagin 14(12) |  | 14 | 35 | 35 |  | [Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1412.txt)



#### Properties key:
  * LS = Low storage
  * SSP = Strong stability preserving
  * FSAL = First same as last
  * CFL = Courant-Friedrichs-Lewy

### Example use case

Basic use of the library is shown here (this uses the `rktp86` method):

```fortran
program rklib_example

  use rklib_module, wp => rk_module_rk
  use iso_fortran_env, only: output_unit

  implicit none

  integer,parameter :: n = 2 !! dimension of the system
  real(wp),parameter :: tol = 1.0e-12_wp !! integration tolerance
  real(wp),parameter :: t0 = 0.0_wp !! initial t value
  real(wp),parameter :: dt = 1.0_wp !! initial step size
  real(wp),parameter :: tf = 100.0_wp !! endpoint of integration
  real(wp),dimension(n),parameter :: x0 = [0.0_wp,0.1_wp] !! initial x value

  real(wp),dimension(n) :: xf !! final x value
  type(rktp86_class) :: prop
  character(len=:),allocatable :: message

  call prop%initialize(n=n,f=fvpol,rtol=[tol],atol=[tol])
  call prop%integrate(t0,x0,dt,tf,xf)
  call prop%status(message=message)

  write (output_unit,'(A)') message
  write (output_unit,'(A,F7.2/,A,2E18.10)') &
              'tf =',tf ,'xf =',xf(1),xf(2)

contains

  subroutine fvpol(me,t,x,f)
    !! Right-hand side of van der Pol equation

    class(rk_class),intent(inout)     :: me
    real(wp),intent(in)               :: t
    real(wp),dimension(:),intent(in)  :: x
    real(wp),dimension(:),intent(out) :: f

    f(1) = x(2)
    f(2) = 0.2_wp*(1.0_wp-x(1)**2)*x(2) - x(1)

  end subroutine fvpol

end program rklib_example
```

The result is:

```
Success
tf = 100.00
xf = -0.1360372426E+01  0.1325538438E+01
```

### Example performance comparison

Running the unit tests will generate some performance plots. The following is for the variable-step methods compiled with quadruple precision (e.g, `fpm test rk_test_variable_step --compiler ifort --flag "-DREAL128"`): [rk_test_variable_step_R16.pdf](media/rk_test_variable_step_R16.pdf)

### Compiling

A [Fortran Package Manager](https://github.com/fortran-lang/fpm) manifest file is included, so that the library and test cases can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

To use `rklib` within your FPM project, add the following to your `fpm.toml` file:
```toml
[dependencies]
rklib = { git="https://github.com/jacobwilliams/rklib.git" }
```

By default, the library is built with double precision (`real64`) real values. Explicitly specifying the real kind can be done using the following processor flags:

Preprocessor flag | Kind  | Number of bytes
----------------- | ----- | ---------------
`REAL32`  | `real(kind=real32)`  | 4
`REAL64`  | `real(kind=real64)`  | 8
`REAL128` | `real(kind=real128)` | 16

For example, to build a single precision version of the library, use:

```
fpm build --profile release --flag "-DREAL32"
```

To generate the documentation using [FORD](https://github.com/Fortran-FOSS-Programmers/ford), run:

```
ford ford.md
```

### 3rd Party Dependencies

  * The library requires [roots-fortran](https://github.com/jacobwilliams/roots-fortran).
  * The unit tests require [pyplot-fortran](https://github.com/jacobwilliams/pyplot-fortran), to generate the performance plots.
  * The `coefficients` app (not required to use the library, but used to generate some of the code) requires the [mpfun2020-var1](https://github.com/jacobwilliams/mpfun2020-var1.git) arbitrary precision library.

All of these will be automatically downloaded by FPM.

### Documentation

The latest API documentation for the `master` branch can be found [here](https://jacobwilliams.github.io/rklib/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### Notes

The original version of this code was split off from the [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit) in September 2022.

### For developers

To add a new method to this library:

  * Update the tables (either the fixed or variable one in `scripts/generate_files.py`)
  * Run `python scripts/generate_files.py` to update all the include files. This script will generate all the boilerplate code for all the methods. It will also this `README` file.
  * Add a step function (either in `rklib_fixed_steps.f90` or `rklib_variable_steps.f90`). Note that you can generate a template of an RK step function using the `scripts/generate_rk_code.py` script. The two command line arguments are the number of function evaluations required and the method name (e.g., `'rk4'`). Edit the template accordingly (note at the FSAL ones have a slightly different format).
  * Update the unit tests.

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
  * C.-W. Shu, S. Osher, "[Efficient implementation of essentially non-oscillatory shock-capturing schemes](https://doi.org/10.1016/0021-9991(88)90177-5)", Journal of Computational Physics, 77(2), 439-471, 1988.
  * S. Ruuth, "[Global optimization of explicit strong-stability-preserving Runge-Kutta methods.](https://doi.org/10.1090/S0025-5718-05-01772-2)" Mathematics of Computation 75.253 (2006): 183-207.
  * Jiang, Guang-Shan, and Chi-Wang Shu. "[Efficient implementation of weighted ENO schemes.](https://doi.org/10.1006/jcph.1996.0130)" Journal of computational physics 126.1 (1996): 202-228.

### See also

  * [FOODIE](https://github.com/Fortran-FOSS-Programmers/FOODIE) Fortran Object-Oriented Differential-equations Integration Environment
  * [FLINT](https://github.com/princemahajan/FLINT) Fortran Library for numerical INTegration of differential equations
  * [DDEABM](https://github.com/jacobwilliams/ddeabm) Modern Fortran implementation of the DDEABM Adams-Bashforth algorithm
  * [DOP853](https://github.com/jacobwilliams/dop853) Modern Fortran Edition of Hairer's DOP853 ODE Solver. An explicit Runge-Kutta method of order 8(5,3) for problems y'=f(x,y); with dense output of order 7
  * [DVODE](https://github.com/jacobwilliams/dvode) Modern Fortran Edition of the DVODE ODE Solver
  * [ODEPACK](https://github.com/jacobwilliams/odepack) Work in Progress to refactor and modernize the ODEPACK Library
  * [libode](https://github.com/markmbaum/libode) Easy-to-compile, high-order ODE solvers as C++ classes


