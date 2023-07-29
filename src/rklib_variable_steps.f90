!*****************************************************************************************
!>
!  Variable-step RK formulas.

    submodule(rklib_module) rklib_variable_steps

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Bogacki-Shampine 3(2) method
!
!### Reference
!  * Bogacki, P and Shampine, L. F. (1989),
!    "A 3(2) pair of Runge-Kutta formulas",
!    Applied Mathematics Letters, 2 (4): 321-325
!
!@note This is a first-same-as-last (FSAL) step.

    module procedure rkbs32

    real(wp),parameter :: a2 =  1.0_wp / 2.0_wp
    real(wp),parameter :: a3 =  3.0_wp / 4.0_wp

    real(wp),parameter :: b21 = 1.0_wp / 2.0_wp
    real(wp),parameter :: b32 = 3.0_wp / 4.0_wp
    ! real(wp),parameter :: b41 = 2.0_wp / 9.0_wp ! FSAL
    ! real(wp),parameter :: b42 = 1.0_wp / 3.0_wp
    ! real(wp),parameter :: b43 = 4.0_wp / 9.0_wp

    real(wp),parameter :: c1  = 2.0_wp / 9.0_wp  ! 3th order
    real(wp),parameter :: c2  = 1.0_wp / 3.0_wp
    real(wp),parameter :: c3  = 4.0_wp / 9.0_wp

    real(wp),parameter :: d1  = 7.0_wp / 24.0_wp  ! 2nd order
    real(wp),parameter :: d2  = 1.0_wp / 4.0_wp
    real(wp),parameter :: d3  = 1.0_wp / 3.0_wp
    real(wp),parameter :: d4  = 1.0_wp / 8.0_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e2  = c2  - d2
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  =     - d4

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4))

        ! check the cached function eval of the last step:
        call me%check_fsal_cache(t,x,f1)

        call me%f(t+a2*h, x+h*(b21*f1),f2)
        call me%f(t+a3*h, x+h*(b32*f2),f3)

        ! last point is cached for the next step:
        xf = x + h*(c1*f1 + c2*f2 + c3*f3)
        call me%set_fsal_cache(t+h,xf,f4)

        xerr = h*(e1*f1 + e2*f2 + e3*f3 + e4*f4)

    end associate

    end procedure rkbs32
!*****************************************************************************************

!*****************************************************************************************
!>
!  4-stage, 3rd order SSP Runge-Kutta method.
!
! @note
!    On the current performance tests, this method shows 3rd order convergence, but its
! relative error is one order of magnitude higher than that obtained with other 3rd order
! methods. An implementation error cannot be excluded.
!
!### References
!
!   * Main method: Kraaijevanger, Johannes Franciscus Bernardus Maria. "Contractivity of
!   runge-kutta methods." BIT Numerical Mathematics 31.3 (1991): 482-528.
!   [DOI: 10.1007/BF01933264](https://doi.org/10.1007/BF01933264).
!
!   * Embedded method: Conde, Sidafa, Imre Fekete, and John N. Shadid. "Embedded error
!   estimation and adaptive step-size control for optimal explicit strong stability preserving
!   Runge--Kutta methods." arXiv preprint arXiv:1806.08693 (2018).
!   [arXiv: 1806.08693](https://arXiv.org/abs/1806.08693)
!
!   * Implementation:  Ranocha, Hendrik, et al. "Optimized Runge-Kutta methods with automatic
!   step size control for compressible computational fluid dynamics." Communications on Applied
!   Mathematics and Computation 4.4 (2022): 1191-1228.
!   [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
!
    module procedure rkssp43

    associate (xtilde => me%funcs(:,1), &
               fs => me%funcs(:,2), &
               half_h => h / 2.0_wp)

        call me%f(t, x, fs)
        xf = x + half_h*fs
        call me%f(t + half_h, xf, fs)
        xf = xf + half_h*fs
        call me%f(t + h, xf, fs)
        xf = xf + half_h*fs
        xtilde = (x + 2.0_wp*xf) / 3.0_wp
        xf = (2.0_wp*x + xf) / 3.0_wp
        call me%f(t + half_h, xf, fs)

        xf = xf + half_h*fs

        xerr = (xtilde - xf) / 2.0_wp

    end associate

    end procedure rkssp43
!*****************************************************************************************

!*****************************************************************************************
!>
!  Fehlberg's 4(5) method.
!
!  This is Table III, RK4(5), Formula 2 in the reference.
!
!### References
!  * E. Fehlberg, "Low-order classical Runge-Kutta formulas with stepsize control and
!    their application to some heat transfer problems", NASA Technical Report R-315,
!    July 1, 1969.

    module procedure rkf45

    real(wp),parameter :: a2 = 1.0_wp  / 4.0_wp
    real(wp),parameter :: a3 = 3.0_wp  / 8.0_wp
    real(wp),parameter :: a4 = 12.0_wp / 13.0_wp
    real(wp),parameter :: a6 = 1.0_wp  / 2.0_wp

    real(wp),parameter :: b21 =  1.0_wp / 4.0_wp
    real(wp),parameter :: b31 =  3.0_wp / 32.0_wp
    real(wp),parameter :: b32 =  9.0_wp / 32.0_wp
    real(wp),parameter :: b41 =  1932.0_wp  / 2197.0_wp
    real(wp),parameter :: b42 =  -7200.0_wp / 2197.0_wp
    real(wp),parameter :: b43 =  7296.0_wp  / 2197.0_wp
    real(wp),parameter :: b51 = 439.0_wp  / 216.0_wp
    real(wp),parameter :: b52 = -8.0_wp
    real(wp),parameter :: b53 = 3680.0_wp / 513.0_wp
    real(wp),parameter :: b54 = -845.0_wp / 4104.0_wp
    real(wp),parameter :: b61 = -8.0_wp     / 27.0_wp
    real(wp),parameter :: b62 = 2.0_wp
    real(wp),parameter :: b63 = -3544.0_wp  / 2565.0_wp
    real(wp),parameter :: b64 =  1859.0_wp  / 4104.0_wp
    real(wp),parameter :: b65 = -11.0_wp    / 40.0_wp

    real(wp),parameter :: c1 = 16.0_wp    / 135.0_wp ! 5th order
    real(wp),parameter :: c3 = 6656.0_wp  / 12825.0_wp
    real(wp),parameter :: c4 = 28561.0_wp / 56430.0_wp
    real(wp),parameter :: c5 = -9.0_wp    / 50.0_wp
    real(wp),parameter :: c6 = 2.0_wp     / 55.0_wp

    real(wp),parameter :: d1 = 25.0_wp / 216.0_wp    ! 4th order
    real(wp),parameter :: d3 = 1408.0_wp / 2565.0_wp
    real(wp),parameter :: d4 = 2197.0_wp / 4104.0_wp
    real(wp),parameter :: d5 = -1.0_wp  / 5.0_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6))

        call me%f(t,       x,f1)
        call me%f(t+a2*h,  x+h*(b21*f1),f2)
        call me%f(t+a3*h,  x+h*(b31*f1 + b32*f2),f3)
        call me%f(t+a4*h,  x+h*(b41*f1 + b42*f2 + b43*f3),f4)
        call me%f(t+h,     x+h*(b51*f1 + b52*f2 + b53*f3 + b54*f4),f5)
        call me%f(t+a6*h,  x+h*(b61*f1 + b62*f2 + b63*f3 + b64*f4 + b65*f5),f6)

        xf = x + h*(d1*f1 + d3*f3 + d4*f4 + d5*f5)

        xerr =   h*(e1*f1 + e3*f3 + e4*f4 + e5*f5 + e6*f6)

    end associate

    end procedure rkf45
!*****************************************************************************************

!*****************************************************************************************
!>
!  Runge Kutta Cash-Karp.
!
!### Reference
!  * J. R. Cash, A. H. Karp. "A variable order Runge-Kutta method
!    for initial value problems with rapidly varying right-hand sides",
!    ACM Transactions on Mathematical Software 16: 201-222, 1990

    module procedure rkck54

    real(wp),parameter :: a2 =  1.0_wp / 5.0_wp
    real(wp),parameter :: a3 =  3.0_wp / 10.0_wp
    real(wp),parameter :: a4 =  3.0_wp / 5.0_wp
    real(wp),parameter :: a6 =  7.0_wp / 8.0_wp

    real(wp),parameter :: b21 = 1.0_wp / 5.0_wp
    real(wp),parameter :: b31 = 3.0_wp / 40.0_wp
    real(wp),parameter :: b32 = 9.0_wp / 40.0_wp
    real(wp),parameter :: b41 = 3.0_wp / 10.0_wp
    real(wp),parameter :: b42 = -9.0_wp / 10.0_wp
    real(wp),parameter :: b43 = 6.0_wp / 5.0_wp
    real(wp),parameter :: b51 = -11.0_wp / 54.0_wp
    real(wp),parameter :: b52 = 5.0_wp / 2.0_wp
    real(wp),parameter :: b53 = -70.0_wp / 27.0_wp
    real(wp),parameter :: b54 = 35.0_wp / 27.0_wp
    real(wp),parameter :: b61 = 1631.0_wp / 55296.0_wp
    real(wp),parameter :: b62 = 175.0_wp / 512.0_wp
    real(wp),parameter :: b63 = 575.0_wp / 13824.0_wp
    real(wp),parameter :: b64 = 44275.0_wp / 110592.0_wp
    real(wp),parameter :: b65 = 253.0_wp / 4096.0_wp

    real(wp),parameter :: c1  = 37.0_wp / 378.0_wp ! 5th order
    real(wp),parameter :: c3  = 250.0_wp / 621.0_wp
    real(wp),parameter :: c4  = 125.0_wp / 594.0_wp
    real(wp),parameter :: c6  = 512.0_wp / 1771.0_wp

    real(wp),parameter :: d1  = 2825.0_wp / 27648.0_wp ! 4th order
    real(wp),parameter :: d3  = 18575.0_wp / 48384.0_wp
    real(wp),parameter :: d4  = 13525.0_wp / 55296.0_wp
    real(wp),parameter :: d5  = 277.0_wp / 14336.0_wp
    real(wp),parameter :: d6  = 1.0_wp / 4.0_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  =     - d5
    real(wp),parameter :: e6  = c6  - d6

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6))

        call me%f(t,       x,f1)
        call me%f(t+a2*h,  x+h*(b21*f1),f2)
        call me%f(t+a3*h,  x+h*(b31*f1 + b32*f2),f3)
        call me%f(t+a4*h,  x+h*(b41*f1 + b42*f2 + b43*f3),f4)
        call me%f(t+h,     x+h*(b51*f1 + b52*f2 + b53*f3 + b54*f4),f5)
        call me%f(t+a6*h,  x+h*(b61*f1 + b62*f2 + b63*f3 + b64*f4 + b65*f5),f6)

        xf = x + h*(c1*f1 + c3*f3 + c4*f4 + c6*f6  )

        xerr =   h*(e1*f1 + e3*f3 + e4*f4 + e5*f5 + e6*f6 )

    end associate

    end procedure rkck54
!*****************************************************************************************

!*****************************************************************************************
!>
!   Dormand-Prince 5(4) method (RKDP or DOPRI)
!
!### Reference
!  * Dormand, J. R.; Prince, P. J. (1980), "A family of embedded Runge-Kutta formulae",
!    Journal of Computational and Applied Mathematics, 6 (1): 19-26
!  * https://en.wikipedia.org/wiki/Dormand-Prince_method
!
!@note This is a first-same-as-last (FSAL) step.

    module procedure rkdp54

    real(wp),parameter :: a2 =  1.0_wp / 5.0_wp
    real(wp),parameter :: a3 =  3.0_wp / 10.0_wp
    real(wp),parameter :: a4 =  4.0_wp / 5.0_wp
    real(wp),parameter :: a5 =  8.0_wp / 9.0_wp

    real(wp),parameter :: b21 =  1.0_wp / 5.0_wp
    real(wp),parameter :: b31 =  3.0_wp / 40.0_wp
    real(wp),parameter :: b32 =  9.0_wp / 40.0_wp
    real(wp),parameter :: b41 =  44.0_wp / 45.0_wp
    real(wp),parameter :: b42 =  -56.0_wp / 15.0_wp
    real(wp),parameter :: b43 =  32.0_wp / 9.0_wp
    real(wp),parameter :: b51 =  19372.0_wp / 6561.0_wp
    real(wp),parameter :: b52 =  -25360.0_wp / 2187.0_wp
    real(wp),parameter :: b53 =  64448.0_wp / 6561.0_wp
    real(wp),parameter :: b54 =  -212.0_wp / 729.0_wp
    real(wp),parameter :: b61 =  9017.0_wp / 3168.0_wp
    real(wp),parameter :: b62 =  -355.0_wp / 33.0_wp
    real(wp),parameter :: b63 =  46732.0_wp / 5247.0_wp
    real(wp),parameter :: b64 =  49.0_wp / 176.0_wp
    real(wp),parameter :: b65 =  -5103.0_wp / 18656.0_wp
    real(wp),parameter :: b71 =  35.0_wp / 384.0_wp
    real(wp),parameter :: b73 =  500.0_wp / 1113.0_wp
    real(wp),parameter :: b74 =  125.0_wp / 192.0_wp
    real(wp),parameter :: b75 =  -2187.0_wp / 6784.0_wp
    real(wp),parameter :: b76 =  11.0_wp / 84.0_wp

    real(wp),parameter :: c1  = b71 ! FSAL
    real(wp),parameter :: c3  = b73
    real(wp),parameter :: c4  = b74
    real(wp),parameter :: c5  = b75
    real(wp),parameter :: c6  = b76

    real(wp),parameter :: d1  = 5179.0_wp / 57600.0_wp
    real(wp),parameter :: d3  = 7571.0_wp / 16695.0_wp
    real(wp),parameter :: d4  = 393.0_wp / 640.0_wp
    real(wp),parameter :: d5  = -92097.0_wp / 339200.0_wp
    real(wp),parameter :: d6  = 187.0_wp / 2100.0_wp
    real(wp),parameter :: d7  = 1.0_wp / 40.0_wp

    real(wp),parameter :: e1  = c1 - d1
    real(wp),parameter :: e3  = c3 - d3
    real(wp),parameter :: e4  = c4 - d4
    real(wp),parameter :: e5  = c5 - d5
    real(wp),parameter :: e6  = c6 - d6
    real(wp),parameter :: e7  =    - d7

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               tf => t + h)

        ! check the cached function eval of the last step:
        call me%check_fsal_cache(t,x,f1)

        call me%f(t+a2*h,  x+h*(b21*f1),f2)
        call me%f(t+a3*h,  x+h*(b31*f1 + b32*f2),f3)
        call me%f(t+a4*h,  x+h*(b41*f1 + b42*f2 + b43*f3),f4)
        call me%f(t+a5*h,  x+h*(b51*f1 + b52*f2 + b53*f3 + b54*f4),f5)
        call me%f(tf,      x+h*(b61*f1 + b62*f2 + b63*f3 + b64*f4 + b65*f5),f6)

        ! last point is cached for the next step:
        xf = x + h*(c1*f1 + c3*f3 + c4*f4 + c5*f5 + c6*f6)
        call me%set_fsal_cache(tf,xf,f7)

        xerr = h*(e1*f1 + e3*f3 + e4*f4 + e5*f5 + e6*f6 + e7*f7)

    end associate

    end procedure rkdp54
!*****************************************************************************************

!*****************************************************************************************
!>
!  Tsitouras 5(4) method
!
!### Reference
!  * Ch. Tsitouras, "Runge-Kutta pairs of order 5(4) satisfying only the
!    first column simplifying assumption", Computers & Mathematics with
!    Applications, Vol. 62, No. 2, pp. 770 - 775, 2011.
!  * [Higher-precision coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK5/RKcoeff5n_1.pdf)

    module procedure rkt54

    real(wp),parameter :: a2  = 0.161_wp
    real(wp),parameter :: a3  = 0.327_wp
    real(wp),parameter :: a4  = 0.9_wp
    real(wp),parameter :: a5  = 0.9800255409045096857298102862870245954942137979563024768854764293221195950761080302604_wp
    real(wp),parameter :: a6  = 1.0_wp
    real(wp),parameter :: b32 = 0.3354806554923569885444268742502307746751211773934303915373692342452941929761641411569_wp
    real(wp),parameter :: b42 = -6.359448489975074843148159912383825625952700647415626703305928850207288721235210244366_wp
    real(wp),parameter :: b43 = 4.362295432869581411017727318190886861027813359713760212991062156752264926097707165077_wp
    real(wp),parameter :: b52 = -11.74888356406282787774717033978577296188744178259862899288666928009020615663593781589_wp
    real(wp),parameter :: b53 = 7.495539342889836208304604784564358155658679161518186721010132816213648793440552049753_wp
    real(wp),parameter :: b54 = -0.9249506636175524925650207933207191611349983406029535244034750452930469056411389539635e-1_wp
    real(wp),parameter :: b62 = -12.92096931784710929170611868178335939541780751955743459166312250439928519268343184452_wp
    real(wp),parameter :: b63 = 8.159367898576158643180400794539253485181918321135053305748355423955009222648673734986_wp
    real(wp),parameter :: b64 = -0.7158497328140099722453054252582973869127213147363544882721139659546372402303777878835e-1_wp
    real(wp),parameter :: b65 = -0.2826905039406838290900305721271224146717633626879770007617876201276764571291579142206e-1_wp
    real(wp),parameter :: b21 = a2
    real(wp),parameter :: b31 = a3 - (b32)
    real(wp),parameter :: b41 = a4 - (b42+b43)
    real(wp),parameter :: b51 = a5 - (b52+b53+b54)
    real(wp),parameter :: b61 = a6 - (b62+b63+b64+b65)
    real(wp),parameter :: c1  = 0.9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1_wp
    real(wp),parameter :: c2  = 0.1e-1_wp
    real(wp),parameter :: c3  = 0.4798896504144995747752495322905965199130404621990332488332634944254542060153074523509_wp
    real(wp),parameter :: c4  = 1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331_wp
    real(wp),parameter :: c5  = -3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677_wp
    real(wp),parameter :: c6  = 2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841_wp
    real(wp),parameter :: d1  = 0.9468075576583945807478876255758922856117527357724631226139574065785592789071067303271e-1_wp
    real(wp),parameter :: d2  = 0.9183565540343253096776363936645313759813746240984095238905939532922955247253608687270e-2_wp
    real(wp),parameter :: d3  = 0.4877705284247615707855642599631228241516691959761363774365216240304071651579571959813_wp
    real(wp),parameter :: d4  = 1.234297566930478985655109673884237654035539930748192848315425833500484878378061439761_wp
    real(wp),parameter :: d5  = -2.707712349983525454881109975059321670689605166938197378763992255714444407154902012702_wp
    real(wp),parameter :: d6  = 1.866628418170587035753719399566211498666255505244122593996591602841258328965767580089_wp
    real(wp),parameter :: d7  = 1.0_wp / 66.0_wp
    real(wp),parameter :: e1  = c1 - d1
    real(wp),parameter :: e2  = c2 - d2
    real(wp),parameter :: e3  = c3 - d3
    real(wp),parameter :: e4  = c4 - d4
    real(wp),parameter :: e5  = c5 - d5
    real(wp),parameter :: e6  = c6 - d6
    real(wp),parameter :: e7  =    - d7

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               tf => t + h)

        ! check the cached function eval of the last step:
        call me%check_fsal_cache(t,x,f1)

        call me%f(t+a2*h, x+h*(b21*f1),f2)
        call me%f(t+a3*h, x+h*(b31*f1 + b32*f2),f3)
        call me%f(t+a4*h, x+h*(b41*f1 + b42*f2 + b43*f3),f4)
        call me%f(t+a5*h, x+h*(b51*f1 + b52*f2 + b53*f3 + b54*f4),f5)
        call me%f(tf,     x+h*(b61*f1 + b62*f2 + b63*f3 + b64*f4 + b65*f5),f6)

        ! last point is cached for the next step:
        xf = x + h*(c1*f1 + c2*f2 + c3*f3 + c4*f4 + c5*f5 + c6*f6)
        call me%set_fsal_cache(tf,xf,f7)

        xerr = h*(e1*f1 + e2*f2 + e3*f3 + e4*f4 + e5*f5 + e6*f6 + e7*f7)

    end associate

    end procedure rkt54
!*****************************************************************************************

!*****************************************************************************************
!>
!  Stepanov 5(4) method.
!
!### References
!  * Misha Stepanov, Calcolo volume 59, Article number: 41 (2022),
!    [Embedded (4,5) pairs of explicit 7-stage Runge-Kutta methods with FSAL property](https://link.springer.com/article/10.1007/s10092-022-00486-1).
!    [arxiv](https://arxiv.org/pdf/2108.12590.pdf) (see Table 4) Note that there is also
!    a 4th order continuously differential interpolant given.

    module procedure rks54

    real(wp),parameter :: a2 = 1.0_wp / 5.0_wp
    real(wp),parameter :: a3 = 21.0_wp / 65.0_wp
    real(wp),parameter :: a4 = 9.0_wp / 10.0_wp
    real(wp),parameter :: a5 = 39.0_wp / 40.0_wp

    real(wp),parameter :: b21 = 1.0_wp / 5.0_wp
    real(wp),parameter :: b31 = 21.0_wp / 338.0_wp
    real(wp),parameter :: b32 = 441.0_wp / 1690.0_wp
    real(wp),parameter :: b41 = 639.0_wp / 392.0_wp
    real(wp),parameter :: b42 = -729.0_wp / 140.0_wp
    real(wp),parameter :: b43 = 1755.0_wp / 392.0_wp
    real(wp),parameter :: b51 = 4878991.0_wp / 1693440.0_wp
    real(wp),parameter :: b52 = -16601.0_wp / 1792.0_wp
    real(wp),parameter :: b53 = 210067.0_wp / 28224.0_wp
    real(wp),parameter :: b54 = -1469.0_wp / 17280.0_wp
    real(wp),parameter :: b61 = 13759919.0_wp / 4230954.0_wp
    real(wp),parameter :: b62 = -2995.0_wp / 287.0_wp
    real(wp),parameter :: b63 = 507312091.0_wp / 61294590.0_wp
    real(wp),parameter :: b64 = -22.0_wp / 405.0_wp
    real(wp),parameter :: b65 = -7040.0_wp / 180687.0_wp

    real(wp),parameter :: c1 = 1441.0_wp / 14742.0_wp
    real(wp),parameter :: c3 = 114244.0_wp / 234927.0_wp
    real(wp),parameter :: c4 = 118.0_wp / 81.0_wp
    real(wp),parameter :: c5 = -12800.0_wp / 4407.0_wp
    real(wp),parameter :: c6 = 41.0_wp / 22.0_wp

    real(wp),parameter :: e1 = -1.0_wp / 273.0_wp
    real(wp),parameter :: e3 = 2197.0_wp / 174020.0_wp
    real(wp),parameter :: e4 = -4.0_wp / 15.0_wp
    real(wp),parameter :: e5 = 1280.0_wp / 1469.0_wp
    real(wp),parameter :: e6 = -33743.0_wp / 52712.0_wp
    real(wp),parameter :: e7 = 127.0_wp / 4792.0_wp

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7))

        ! check the cached function eval of the last step:
        call me%check_fsal_cache(t,x,f1)

        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+h,   x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)

        ! last point is cached for the next step:
        xf = x+h*(c1*f1+c3*f3+c4*f4+c5*f5+c6*f6)
        call me%set_fsal_cache(t+h,xf,f7)

        xerr = h*(e1*f1+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7)

    end associate

    end procedure rks54
!*****************************************************************************************

!*****************************************************************************************
!>
!  Papakostas-PapaGeorgiou 5(4) method (more stable, higher error)
!
!### Reference
!  * S.N. Papakostas and G. PapaGeorgiou, "A Family of Fifth-order Runge-Kutta Pairs",
!    Mathematics of Computation, Vol. 65, No. 215 (Jul., 1996), pp. 1165-1181
!
!@note This is a first-same-as-last (FSAL) step.

    module procedure rkpp54

    real(wp),parameter :: a2 = 2.03174603174603174603174603174603174603174603174603174603174603174603174603174603175e-1_wp
    real(wp),parameter :: a3 = 3.01837270341207349081364829396325459317585301837270341207349081364829396325459317585e-1_wp
    real(wp),parameter :: a4 = 8.14973262032085561497326203208556149732620320855614973262032085561497326203208556150e-1_wp
    real(wp),parameter :: a5 = 8.92857142857142857142857142857142857142857142857142857142857142857142857142857142857e-1_wp

    real(wp),parameter :: b21 = 2.03174603174603174603174603174603174603174603174603174603174603174603174603174603175e-1_wp
    real(wp),parameter :: b31 = 7.76317438051542769752206171079008824684316035298737264141194949056564779796226259119e-2_wp
    real(wp),parameter :: b32 = 2.24205526536053072106144212288424576849153698307396614793229586459172918345836691673e-1_wp
    real(wp),parameter :: b41 = 9.97054225910649851955403847106996156261144080556828071531226332913609633807524002826e-1_wp
    real(wp),parameter :: b42 = -3.9229583098359766207115693643000730264502583412156151006542354260709699970151262085e0_wp
    real(wp),parameter :: b43 = 3.74087734595741233025349172040163301992173458151440200238504117871885768941081076180e0_wp
    real(wp),parameter :: b51 = 2.41107195213820140711520163970784555943075679078314370129890823084996830515582884769e0_wp
    real(wp),parameter :: b52 = -9.6949301103878662103154820500485383234086722205716227648448127788897312407254087039e0_wp
    real(wp),parameter :: b53 = 8.37096577602619810632388905477462929171641908259741309319801335095185921382851168188e0_wp
    real(wp),parameter :: b54 = -1.9425047491939044598075150157679367059564650995179117250925166005495342111607468280e-1_wp
    real(wp),parameter :: b61 = 2.71680576250618168562042732618844645610969561203458029063121928183564071601176708510e0_wp
    real(wp),parameter :: b62 = -1.0546836035293105493140447781471570084309037601955971025884876060757459965008137190e1_wp
    real(wp),parameter :: b63 = 8.82438454669452608734254952418112316791646628404483378649096846259826474086840654822e0_wp
    real(wp),parameter :: b64 = 3.30231934937991900184188201559235426225444020547394478857972406538849522552348656011e-1_wp
    real(wp),parameter :: b65 = -3.2458620884559418000671727045723496594256831467083753009528409021529501442438509930e-1_wp
    ! real(wp),parameter :: b71 = 9.15102514359618091216858762219939898816995701624253490052873825554414393852942295257e-2_wp ! FSAL
    ! real(wp),parameter :: b72 = 0
    ! real(wp),parameter :: b73 = 4.53088687963682563226735078441923499983686103419218207034508391590507417156754767246e-1_wp
    ! real(wp),parameter :: b74 = 7.70212354683768742558003248195356972937967314339813992394948819785937766391608532599e-1_wp
    ! real(wp),parameter :: b75 = -4.6993115056469338687729550442725993708636301781314244468307254289871602785250570955e-1_wp
    ! real(wp),parameter :: b76 = 1.55119856481280271970871301567985474283010029891684896248327948966829404918848180181e-1_wp

    real(wp),parameter :: c1 = 9.15102514359618091216858762219939898816995701624253490052873825554414393852942295257e-2_wp
    real(wp),parameter :: c3 = 4.53088687963682563226735078441923499983686103419218207034508391590507417156754767246e-1_wp
    real(wp),parameter :: c4 = 7.70212354683768742558003248195356972937967314339813992394948819785937766391608532599e-1_wp
    real(wp),parameter :: c5 = -4.6993115056469338687729550442725993708636301781314244468307254289871602785250570955e-1_wp
    real(wp),parameter :: c6 = 1.55119856481280271970871301567985474283010029891684896248327948966829404918848180181e-1_wp

    real(wp),parameter :: d1 = 9.00540754825040114941655607289341515504490822369589853929344566160995864634775952099e-2_wp
    real(wp),parameter :: d3 = 4.58093041355081442291135932317120322249435302966675363526235643007610910215270397680e-1_wp
    real(wp),parameter :: d4 = 7.17141736822002520962978129215068133609337518199443217706780208321061003994969882396e-1_wp
    real(wp),parameter :: d5 = -3.9730095904319871050513227050199810426175540655366523944169265633781699285936171784e-1_wp
    real(wp),parameter :: d6 = 1.08202581574086926233043124431351687328723979341063863291932824583521682661834318740e-1_wp
    real(wp),parameter :: d7 = 2.38095238095238095238095238095238095238095238095238095238095238095238095238095238095e-2_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  =     - d7

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7))

        ! check the cached function eval of the last step:
        call me%check_fsal_cache(t,x,f1)

        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+h,   x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)

        ! last point is cached for the next step:
        xf = x+h*(c1*f1+c3*f3+c4*f4+c5*f5+c6*f6)
        call me%set_fsal_cache(t+h,xf,f7)

        xerr = h*(e1*f1+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7)

    end associate

    end procedure rkpp54
!*****************************************************************************************

!*****************************************************************************************
!>
!  Papakostas-PapaGeorgiou 5(4) method (less stable, lower error).
!
!### Reference
!  * S.N. Papakostas and G. PapaGeorgiou, "A Family of Fifth-order Runge-Kutta Pairs",
!    Mathematics of Computation, Vol. 65, No. 215 (Jul., 1996), pp. 1165-1181
!
!@note This is a first-same-as-last (FSAL) step.

    module procedure rkpp54b

    real(wp),parameter :: a2 = 2.201257861635220125786163522012578616352201257861635220125786163522012578616352201257e-1_wp
    real(wp),parameter :: a3 = 3.206106870229007633587786259541984732824427480916030534351145038167938931297709923664e-1_wp
    real(wp),parameter :: a4 = 9.160839160839160839160839160839160839160839160839160839160839160839160839160839160839e-1_wp
    real(wp),parameter :: a5 = 9.545454545454545454545454545454545454545454545454545454545454545454545454545454545454e-1_wp

    real(wp),parameter :: b21 = 2.201257861635220125786163522012578616352201257861635220125786163522012578616352201257e-1_wp
    real(wp),parameter :: b31 = 8.712778975584173416467571819824019579278596818367228017015325447235009614824310937591e-2_wp
    real(wp),parameter :: b32 = 2.334828972670590291941029077559582774896567799079307732649612493444437969815278829905e-1_wp
    real(wp),parameter :: b41 = 1.229131120244552267075198689039393789189448076016053493121879281174531378872492304514e0_wp
    real(wp),parameter :: b42 = -5.174619722997219626464507457640199024688660926548912688647080129723360618339765067163e0_wp
    real(wp),parameter :: b43 = 4.861572518836583443305392684684721319415296766616775279441284764632745323383356678732e0_wp
    real(wp),parameter :: b51 = 1.635985585148686070813241441537203183216255340359063210526474302100593166648405715598e0_wp
    real(wp),parameter :: b52 = -6.911673476406501672554353133581688231701838702998878887272831441157207088140814883585e0_wp
    real(wp),parameter :: b53 = 6.264598153851822604079697377077162805912462110097874626859336947971237445955455470077e0_wp
    real(wp),parameter :: b54 = -3.436480804855245688404023048722321197233329291260440465843435436916897900850084754553e-2_wp
    real(wp),parameter :: b61 = 2.080199455181079060589945391266072650640389894799843423576531143602673095720452243744e0_wp
    real(wp),parameter :: b62 = -8.741042793278798222268721650787838365080091632936417117311519288907488660315135346032e0_wp
    real(wp),parameter :: b63 = 7.705947570635899833548328481096646034173724896740120966754373004394344077008935321798e0_wp
    real(wp),parameter :: b64 = 8.908088473605147947645976479787347297467362259887092033560387203210318720083482332807e-2_wp
    real(wp),parameter :: b65 = -1.341851172742321513460119863727537927086967812024181933549887311216316996150870428391e-1_wp
    ! real(wp),parameter :: b71 = 9.660368812495311950055679608569466802836486587849728525846301200717776508398100546410e-2_wp ! FSAL
    ! real(wp),parameter :: b72 = 0
    ! real(wp),parameter :: b73 = 4.847050840514515173851124515232483753229089551347247274647212546244695209679216658037e-1_wp
    ! real(wp),parameter :: b74 = 2.985721504066570633782817413202582623066368759430782759246016403146502612390979295903e0_wp
    ! real(wp),parameter :: b75 = -3.882058366130615720106688908002536902372698760219285670845605164160172370352994326721e0_wp
    ! real(wp),parameter :: b76 = 1.315028089887640449438202247191011235955056179775280898876404494382022471910112359550e0_wp

    real(wp),parameter :: c1 = 9.660368812495311950055679608569466802836486587849728525846301200717776508398100546410e-2_wp
    real(wp),parameter :: c3 = 4.847050840514515173851124515232483753229089551347247274647212546244695209679216658037e-1_wp
    real(wp),parameter :: c4 = 2.985721504066570633782817413202582623066368759430782759246016403146502612390979295903e0_wp
    real(wp),parameter :: c5 = -3.882058366130615720106688908002536902372698760219285670845605164160172370352994326721e0_wp
    real(wp),parameter :: c6 = 1.315028089887640449438202247191011235955056179775280898876404494382022471910112359550e0_wp

    real(wp),parameter :: d1 = 9.532926220038256951865552620720509864047974971891244940167035171716827175894410248606e-2_wp
    real(wp),parameter :: d3 = 4.890503958633595206242721536021852509608016393506331933941960604717211861520491115111e-1_wp
    real(wp),parameter :: d4 = 2.782787369422591128285131424607412688536477162495912788909265685055380447336451073683e0_wp
    real(wp),parameter :: d5 = -3.544320847545574336397313295146726421831558091536665262636276708453748484122524958198e0_wp
    real(wp),parameter :: d6 = 1.168225248630669689397825619301351955122370968542635402359716039780907150303652099088e0_wp
    real(wp),parameter :: d7 = 8.928571428571428571428571428571428571428571428571428571428571428571428571428571428571e-3_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  =     - d7

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7))

        ! check the cached function eval of the last step:
        call me%check_fsal_cache(t,x,f1)

        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+h,   x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)

        ! last point is cached for the next step:
        xf = x+h*(c1*f1+c3*f3+c4*f4+c5*f5+c6*f6)
        call me%set_fsal_cache(t+h,xf,f7)

        xerr = h*(e1*f1+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7)

    end associate

    end procedure rkpp54b
!*****************************************************************************************

!*****************************************************************************************
!>
!  Bogacki & Shampine 5(4)
!
!### Reference
!  * P. Bogacki & L.F. Shampine, "An efficient Runge-Kutta (4,5) pair"
!    Computers & Mathematics with Applications,
!    Volume 32, Issue 6, September 1996, Pages 15-28
!  * [Coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK5/RKcoeff5p_1.pdf)
!  * This pair is in [rksuite](https://github.com/jacobwilliams/rksuite) (`METHD = 2`)

    module procedure rkbs54

    real(wp),parameter :: a2 = 1.0_wp / 6.0_wp
    real(wp),parameter :: a3 = 2.0_wp / 9.0_wp
    real(wp),parameter :: a4 = 3.0_wp / 7.0_wp
    real(wp),parameter :: a5 = 2.0_wp / 3.0_wp
    real(wp),parameter :: a6 = 3.0_wp / 4.0_wp

    real(wp),parameter :: b21 = 1.0_wp / 6.0_wp
    real(wp),parameter :: b31 = 2.0_wp / 27.0_wp
    real(wp),parameter :: b32 = 4.0_wp / 27.0_wp
    real(wp),parameter :: b41 = 183.0_wp / 1372.0_wp
    real(wp),parameter :: b42 = -162.0_wp / 343.0_wp
    real(wp),parameter :: b43 = 1053.0_wp / 1372.0_wp
    real(wp),parameter :: b51 = 68.0_wp / 297.0_wp
    real(wp),parameter :: b52 = -4.0_wp / 11.0_wp
    real(wp),parameter :: b53 = 42.0_wp / 143.0_wp
    real(wp),parameter :: b54 = 1960.0_wp / 3861.0_wp
    real(wp),parameter :: b61 = 597.0_wp / 22528.0_wp
    real(wp),parameter :: b62 = 81.0_wp / 352.0_wp
    real(wp),parameter :: b63 = 63099.0_wp / 585728.0_wp
    real(wp),parameter :: b64 = 58653.0_wp / 366080.0_wp
    real(wp),parameter :: b65 = 4617.0_wp / 20480.0_wp
    real(wp),parameter :: b71 = 174197.0_wp / 959244.0_wp
    real(wp),parameter :: b72 = -30942.0_wp / 79937.0_wp
    real(wp),parameter :: b73 = 8152137.0_wp / 19744439.0_wp
    real(wp),parameter :: b74 = 666106.0_wp / 1039181.0_wp
    real(wp),parameter :: b75 = -29421.0_wp / 29068.0_wp
    real(wp),parameter :: b76 = 482048.0_wp / 414219.0_wp
    real(wp),parameter :: b81 = 587.0_wp / 8064.0_wp
    real(wp),parameter :: b83 = 4440339.0_wp / 15491840.0_wp
    real(wp),parameter :: b84 = 24353.0_wp / 124800.0_wp
    real(wp),parameter :: b85 = 387.0_wp / 44800.0_wp
    real(wp),parameter :: b86 = 2152.0_wp / 5985.0_wp
    real(wp),parameter :: b87 = 7267.0_wp / 94080.0_wp

    real(wp),parameter :: c1 = 587.0_wp / 8064.0_wp
    real(wp),parameter :: c3 = 4440339.0_wp / 15491840.0_wp
    real(wp),parameter :: c4 = 24353.0_wp / 124800.0_wp
    real(wp),parameter :: c5 = 387.0_wp / 44800.0_wp
    real(wp),parameter :: c6 = 2152.0_wp / 5985.0_wp
    real(wp),parameter :: c7 = 7267.0_wp / 94080.0_wp

    ! real(wp),parameter :: d1 = 6059 / 80640
    ! real(wp),parameter :: d2 = 0
    ! real(wp),parameter :: d3 = 8559189 / 30983680
    ! real(wp),parameter :: d4 = 26411 / 124800
    ! real(wp),parameter :: d5 = -927 / 89600
    ! real(wp),parameter :: d6 = 443 / 1197
    ! real(wp),parameter :: d7 = 7267 / 94080

    real(wp),parameter :: d1 = 2479.0_wp / 34992.0_wp
    real(wp),parameter :: d3 = 123.0_wp / 416.0_wp
    real(wp),parameter :: d4 = 612941.0_wp / 3411720.0_wp
    real(wp),parameter :: d5 = 43.0_wp / 1440.0_wp
    real(wp),parameter :: d6 = 2272.0_wp / 6561.0_wp
    real(wp),parameter :: d7 = 79937.0_wp / 1113912.0_wp
    real(wp),parameter :: d8 = 3293.0_wp / 556956.0_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  =     - d8

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+h,   x+h*(b71*f1+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+h,   x+h*(b81*f1+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)

        xf = x+h*(c1*f1+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7)

        xerr = h*(e1*f1+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8)

    end associate

    end procedure rkbs54
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sharp & Smart 5(4) method.
!
!### Reference
!  * P.W. Sharp and E. Smart, "Explicit Runge-Kutta Pairs with One More Derivative Evaluation than the Minimum",
!    Siam Journal of Scientific Computing, Vol. 14, No. 2, pages. 338-348, March 1993.

    module procedure rkss54

    real(wp),parameter :: a2 = 16.0_wp / 105.0_wp
    real(wp),parameter :: a3 = 8.0_wp / 35.0_wp
    real(wp),parameter :: a4 = 9.0_wp / 20.0_wp
    real(wp),parameter :: a5 = 2.0_wp / 3.0_wp
    real(wp),parameter :: a6 = 7.0_wp / 9.0_wp

    real(wp),parameter :: b21 = 16.0_wp / 105.0_wp
    real(wp),parameter :: b31 = 2.0_wp / 35.0_wp
    real(wp),parameter :: b32 = 6.0_wp / 35.0_wp
    real(wp),parameter :: b41 = 8793.0_wp / 40960.0_wp
    real(wp),parameter :: b42 = -5103.0_wp / 8192.0_wp
    real(wp),parameter :: b43 = 17577.0_wp / 20480.0_wp
    real(wp),parameter :: b51 = 347.0_wp / 1458.0_wp
    real(wp),parameter :: b52 = -7.0_wp / 20.0_wp
    real(wp),parameter :: b53 = 3395.0_wp / 10044.0_wp
    real(wp),parameter :: b54 = 49792.0_wp / 112995.0_wp
    real(wp),parameter :: b61 = -1223224109959.0_wp / 9199771214400.0_wp
    real(wp),parameter :: b62 = 1234787701.0_wp / 2523942720.0_wp
    real(wp),parameter :: b63 = 568994101921.0_wp / 3168810084960.0_wp
    real(wp),parameter :: b64 = -105209683888.0_wp / 891227836395.0_wp
    real(wp),parameter :: b65 = 9.0_wp / 25.0_wp
    real(wp),parameter :: b71 = 2462504862877.0_wp / 8306031988800.0_wp
    real(wp),parameter :: b72 = -123991.0_wp / 287040.0_wp
    real(wp),parameter :: b73 = 106522578491.0_wp / 408709510560.0_wp
    real(wp),parameter :: b74 = 590616498832.0_wp / 804646848915.0_wp
    real(wp),parameter :: b75 = -319138726.0_wp / 534081275.0_wp
    real(wp),parameter :: b76 = 52758.0_wp / 71449.0_wp

    real(wp),parameter :: c1 = 1093.0_wp / 15120.0_wp
    real(wp),parameter :: c3 = 60025.0_wp / 190992.0_wp
    real(wp),parameter :: c4 = 3200.0_wp / 20709.0_wp
    real(wp),parameter :: c5 = 1611.0_wp / 11960.0_wp
    real(wp),parameter :: c6 = 712233.0_wp / 2857960.0_wp
    real(wp),parameter :: c7 = 3.0_wp / 40.0_wp

    real(wp),parameter :: d1 = 84018211.0_wp / 991368000.0_wp
    real(wp),parameter :: d3 = 92098979.0_wp / 357791680.0_wp
    real(wp),parameter :: d4 = 17606944.0_wp / 67891005.0_wp
    real(wp),parameter :: d5 = 3142101.0_wp / 235253200.0_wp
    real(wp),parameter :: d6 = 22004596809.0_wp / 70270091500.0_wp
    real(wp),parameter :: d7 = 9.0_wp / 125.0_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+h,   x+h*(b71*f1+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)

        xf = x+h*(c1*f1+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7)

        xerr = h*(e1*f1+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7)

    end associate

    end procedure rkss54
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dormand-Prince 6(5) method.
!  This is `RK6(5)8M` from the reference.
!
!### Reference
!  * P.J. Prince, J.R. Dormand, "High order embedded Runge-Kutta formulae",
!    Journal of Computational and Applied Mathematics,
!    Volume 7, Issue 1, March 1981, Pages 67-75.

    module procedure rkdp65

    real(wp),parameter :: a2 = 1.0_wp / 10.0_wp
    real(wp),parameter :: a3 = 2.0_wp / 9.0_wp
    real(wp),parameter :: a4 = 3.0_wp / 7.0_wp
    real(wp),parameter :: a5 = 3.0_wp / 5.0_wp
    real(wp),parameter :: a6 = 4.0_wp / 5.0_wp

    real(wp),parameter :: b21 = 1.0_wp / 10.0_wp
    real(wp),parameter :: b31 = -2.0_wp / 81.0_wp
    real(wp),parameter :: b32 = 20.0_wp / 81.0_wp
    real(wp),parameter :: b41 = 615.0_wp / 1372.0_wp
    real(wp),parameter :: b42 = -270.0_wp / 343.0_wp
    real(wp),parameter :: b43 = 1053.0_wp / 1372.0_wp
    real(wp),parameter :: b51 = 3243.0_wp / 5500.0_wp
    real(wp),parameter :: b52 = -54.0_wp / 55.0_wp
    real(wp),parameter :: b53 = 50949.0_wp / 71500.0_wp
    real(wp),parameter :: b54 = 4998.0_wp / 17875.0_wp
    real(wp),parameter :: b61 = -26492.0_wp / 37125.0_wp
    real(wp),parameter :: b62 = 72.0_wp / 55.0_wp
    real(wp),parameter :: b63 = 2808.0_wp / 23375.0_wp
    real(wp),parameter :: b64 = -24206.0_wp / 37125.0_wp
    real(wp),parameter :: b65 = 338.0_wp / 459.0_wp
    real(wp),parameter :: b71 = 5561.0_wp / 2376.0_wp
    real(wp),parameter :: b72 = -35.0_wp / 11.0_wp
    real(wp),parameter :: b73 = -24117.0_wp / 31603.0_wp
    real(wp),parameter :: b74 = 899983.0_wp / 200772.0_wp
    real(wp),parameter :: b75 = -5225.0_wp / 1836.0_wp
    real(wp),parameter :: b76 = 3925.0_wp / 4056.0_wp
    real(wp),parameter :: b81 = 465467.0_wp / 266112.0_wp
    real(wp),parameter :: b82 = -2945.0_wp / 1232.0_wp
    real(wp),parameter :: b83 = -5610201.0_wp / 14158144.0_wp
    real(wp),parameter :: b84 = 10513573.0_wp / 3212352.0_wp
    real(wp),parameter :: b85 = -424325.0_wp / 205632.0_wp
    real(wp),parameter :: b86 = 376225.0_wp / 454272.0_wp

    real(wp),parameter :: c1 = 61.0_wp / 864.0_wp
    real(wp),parameter :: c3 = 98415.0_wp / 321776.0_wp
    real(wp),parameter :: c4 = 16807.0_wp / 146016.0_wp
    real(wp),parameter :: c5 = 1375.0_wp / 7344.0_wp
    real(wp),parameter :: c6 = 1375.0_wp / 5408.0_wp
    real(wp),parameter :: c7 = -37.0_wp / 1120.0_wp
    real(wp),parameter :: c8 = 1.0_wp / 10.0_wp

    real(wp),parameter :: d1 = 821.0_wp / 10800.0_wp
    real(wp),parameter :: d3 = 19683.0_wp / 71825.0_wp
    real(wp),parameter :: d4 = 175273.0_wp / 912600.0_wp
    real(wp),parameter :: d5 = 395.0_wp / 3672.0_wp
    real(wp),parameter :: d6 = 785.0_wp / 2704.0_wp
    real(wp),parameter :: d7 = 3.0_wp / 50.0_wp

    real(wp),parameter :: e1  = c1 - d1
    real(wp),parameter :: e3  = c3 - d3
    real(wp),parameter :: e4  = c4 - d4
    real(wp),parameter :: e5  = c5 - d5
    real(wp),parameter :: e6  = c6 - d6
    real(wp),parameter :: e7  = c7 - d7
    real(wp),parameter :: e8  = c8

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+h,   x+h*(b71*f1+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+h,   x+h*(b81*f1+b82*f2+b83*f3+b84*f4+b85*f5+b86*f6),f8)

        xf = x+h*(c1*f1+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8)

        xerr = h*(e1*f1+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8)

    end associate

    end procedure rkdp65
!*****************************************************************************************

!*****************************************************************************************
!>
!  Calvo 6(5) method.
!
!### Reference
!  * M. Calvo, J.I. Montijano, L. Randez,
!    "A new embedded pair of Runge-Kutta formulas of orders 5 and 6",
!    Computers & Mathematics with Applications, Volume 20, Issue 1, 1990, Pages 15-24

    module procedure rkc65

    real(wp),parameter :: a2  = 2.0_wp/15.0_wp
    real(wp),parameter :: a3  = 1.0_wp/5.0_wp
    real(wp),parameter :: a4  = 3.0_wp/10.0_wp
    real(wp),parameter :: a5  = 14.0_wp/25.0_wp
    real(wp),parameter :: a6  = 19.0_wp/25.0_wp
    real(wp),parameter :: a7  = 35226607.0_wp/35688279.0_wp

    real(wp),parameter :: b21 = 2.0_wp/15.0_wp
    real(wp),parameter :: b31 = 1.0_wp/20.0_wp
    real(wp),parameter :: b32 = 3.0_wp/20.0_wp
    real(wp),parameter :: b41 = 3.0_wp/40.0_wp
    real(wp),parameter :: b43 = 9.0_wp/40.0_wp
    real(wp),parameter :: b51 = 86727015.0_wp/196851553.0_wp
    real(wp),parameter :: b52 = -60129073.0_wp/52624712.0_wp
    real(wp),parameter :: b53 = 957436434.0_wp/1378352377.0_wp
    real(wp),parameter :: b54 = 83886832.0_wp/147842441.0_wp
    real(wp),parameter :: b61 = -86860849.0_wp/45628967.0_wp
    real(wp),parameter :: b62 = 111022885.0_wp/25716487.0_wp
    real(wp),parameter :: b63 = 108046682.0_wp/101167669.0_wp
    real(wp),parameter :: b64 = -141756746.0_wp/36005461.0_wp
    real(wp),parameter :: b65 = 73139862.0_wp/60170633.0_wp
    real(wp),parameter :: b71 = 77759591.0_wp/16096467.0_wp
    real(wp),parameter :: b72 = -49252809.0_wp/6452555.0_wp
    real(wp),parameter :: b73 = -381680111.0_wp/51572984.0_wp
    real(wp),parameter :: b74 = 879269579.0_wp/66788831.0_wp
    real(wp),parameter :: b75 = -90453121.0_wp/33722162.0_wp
    real(wp),parameter :: b76 = 111179552.0_wp/157155827.0_wp
    real(wp),parameter :: b81 = 237564263.0_wp/39280295.0_wp
    real(wp),parameter :: b82 = -100523239.0_wp/10677940.0_wp
    real(wp),parameter :: b83 = -265574846.0_wp/27330247.0_wp
    real(wp),parameter :: b84 = 317978411.0_wp/18988713.0_wp
    real(wp),parameter :: b85 = -124494385.0_wp/35453627.0_wp
    real(wp),parameter :: b86 = 86822444.0_wp/100138635.0_wp
    real(wp),parameter :: b87 = -12873523.0_wp/724232625.0_wp
    real(wp),parameter :: b91 = 17572349.0_wp/289262523.0_wp
    real(wp),parameter :: b93 = 57513011.0_wp/201864250.0_wp
    real(wp),parameter :: b94 = 15587306.0_wp/354501571.0_wp
    real(wp),parameter :: b95 = 71783021.0_wp/234982865.0_wp
    real(wp),parameter :: b96 = 29672000.0_wp/180480167.0_wp
    real(wp),parameter :: b97 = 65567621.0_wp/127060952.0_wp
    real(wp),parameter :: b98 = -79074570.0_wp/210557597.0_wp

    real(wp),parameter :: c1 = 17572349.0_wp/289262523.0_wp
    real(wp),parameter :: c3 = 57513011.0_wp/201864250.0_wp
    real(wp),parameter :: c4 = 15587306.0_wp/354501571.0_wp
    real(wp),parameter :: c5 = 71783021.0_wp/234982865.0_wp
    real(wp),parameter :: c6 = 29672000.0_wp/180480167.0_wp
    real(wp),parameter :: c7 = 65567621.0_wp/127060952.0_wp
    real(wp),parameter :: c8 = -79074570.0_wp/210557597.0_wp

    real(wp),parameter :: d1 = 15231665.0_wp/510830334.0_wp
    real(wp),parameter :: d3 = 59452991.0_wp/116050448.0_wp
    real(wp),parameter :: d4 = -28398517.0_wp/122437738.0_wp
    real(wp),parameter :: d5 = 56673824.0_wp/137010559.0_wp
    real(wp),parameter :: d6 = 68003849.0_wp/426673583.0_wp
    real(wp),parameter :: d7 = 7097631.0_wp/37564021.0_wp
    real(wp),parameter :: d8 = -71226429.0_wp/583093742.0_wp
    real(wp),parameter :: d9 = 1.0_wp/20.0_wp

    real(wp),parameter :: e1  = c1 - d1
    real(wp),parameter :: e3  = c3 - d3
    real(wp),parameter :: e4  = c4 - d4
    real(wp),parameter :: e5  = c5 - d5
    real(wp),parameter :: e6  = c6 - d6
    real(wp),parameter :: e7  = c7 - d7
    real(wp),parameter :: e8  = c8 - d8
    real(wp),parameter :: e9  =    - d9

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9))

        call me%f(t,       x,f1)
        call me%f(t+a2*h,  x+h*(b21*f1),f2)
        call me%f(t+a3*h,  x+h*(b31*f1 + b32*f2),f3)
        call me%f(t+a4*h,  x+h*(b41*f1          + b43*f3),f4)
        call me%f(t+a5*h,  x+h*(b51*f1 + b52*f2 + b53*f3 + b54*f4),f5)
        call me%f(t+a6*h,  x+h*(b61*f1 + b62*f2 + b63*f3 + b64*f4 + b65*f5),f6)
        call me%f(t+a7*h,  x+h*(b71*f1 + b72*f2 + b73*f3 + b74*f4 + b75*f5 + b76*f6),f7)
        call me%f(t+h,     x+h*(b81*f1 + b82*f2 + b83*f3 + b84*f4 + b85*f5 + b86*f6 + b87*f7),f8)
        call me%f(t+h,     x+h*(b91*f1          + b93*f3 + b94*f4 + b95*f5 + b96*f6 + b97*f7 + b98*f8),f9)

        xf = x + h*(c1*f1 + c3*f3 + c4*f4 + c5*f5 + c6*f6 + c7*f7 + c8*f8)

        xerr =   h*(e1*f1 + e3*f3 + e4*f4 + e5*f5 + e6*f6 + e7*f7 + e8*f8 + e9*f9)

    end associate

    end procedure rkc65
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's "most efficient" Runge-Kutta (9,6(5)) pair.
!
!### Reference
!  * J.H. Verner, "Strategies for deriving new explicit Runge-Kutta pairs",
!    Annals of Num. Math 1 1994, 225-244.
!  * [Jim Verner's Refuge for Runge-Kutta Pairs](https://www.sfu.ca/~jverner/)
!  * [Rational coefficients](https://www.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.RATOnWeb)
!  * [Floating point coefficients](https://www.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT)
!
!@note This is a first-same-as-last (FSAL) step.

    module procedure rkv65e

    real(wp),parameter :: a2 =  3.0_wp    / 50.0_wp
    real(wp),parameter :: a3 =  1439.0_wp / 15000.0_wp
    real(wp),parameter :: a4 =  1439.0_wp / 10000.0_wp
    real(wp),parameter :: a5 =  4973.0_wp / 10000.0_wp
    real(wp),parameter :: a6 =  389.0_wp  / 400.0_wp
    real(wp),parameter :: a7 =  1999.0_wp / 2000.0_wp

    real(wp),parameter :: b21 =  .6e-1_wp
    real(wp),parameter :: b31 =  .1923996296296296296296296296296296296296e-1_wp
    real(wp),parameter :: b32 =  .7669337037037037037037037037037037037037e-1_wp
    real(wp),parameter :: b41 =  .35975e-1_wp
    real(wp),parameter :: b43 =  .107925_wp
    real(wp),parameter :: b51 =  1.318683415233148260919747276431735612861_wp
    real(wp),parameter :: b53 = -5.042058063628562225427761634715637693344_wp
    real(wp),parameter :: b54 =  4.220674648395413964508014358283902080483_wp
    real(wp),parameter :: b61 = -41.87259166432751461803757780644346812905_wp
    real(wp),parameter :: b63 =  159.4325621631374917700365669070346830453_wp
    real(wp),parameter :: b64 = -122.1192135650100309202516203389242140663_wp
    real(wp),parameter :: b65 =  5.531743066200053768252631238332999150076_wp
    real(wp),parameter :: b71 = -54.43015693531650433250642051294142461271_wp
    real(wp),parameter :: b73 =  207.0672513650184644273657173866509835987_wp
    real(wp),parameter :: b74 = -158.6108137845899991828742424365058599469_wp
    real(wp),parameter :: b75 =  6.991816585950242321992597280791793907096_wp
    real(wp),parameter :: b76 = -.1859723106220323397765171799549294623692e-1_wp
    real(wp),parameter :: b81 = -54.66374178728197680241215648050386959351_wp
    real(wp),parameter :: b83 =  207.9528062553893734515824816699834244238_wp
    real(wp),parameter :: b84 = -159.2889574744995071508959805871426654216_wp
    real(wp),parameter :: b85 =  7.018743740796944434698170760964252490817_wp
    real(wp),parameter :: b86 = -.1833878590504572306472782005141738268361e-1_wp
    real(wp),parameter :: b87 = -.5119484997882099077875432497245168395840e-3_wp
    ! real(wp),parameter :: b91 =  .3438957868357036009278820124728322386520e-1_wp ! FSAL
    ! real(wp),parameter :: b94 =  .2582624555633503404659558098586120858767_wp
    ! real(wp),parameter :: b95 =  .4209371189673537150642551514069801967032_wp
    ! real(wp),parameter :: b96 =  4.405396469669310170148836816197095664891_wp
    ! real(wp),parameter :: b97 = -176.4831190242986576151740942499002125029_wp
    ! real(wp),parameter :: b98 =  172.3641334014150730294022582711902413315_wp

    real(wp),parameter :: c1  =  .3438957868357036009278820124728322386520e-1_wp
    real(wp),parameter :: c4  =  .2582624555633503404659558098586120858767_wp
    real(wp),parameter :: c5  =  .4209371189673537150642551514069801967032_wp
    real(wp),parameter :: c6  =  4.405396469669310170148836816197095664891_wp
    real(wp),parameter :: c7  = -176.4831190242986576151740942499002125029_wp
    real(wp),parameter :: c8  =  172.3641334014150730294022582711902413315_wp

    real(wp),parameter :: d1 =  .4909967648382489730906854927971225836479e-1_wp
    real(wp),parameter :: d4 =  .2251112229516524153401395320539875329485_wp
    real(wp),parameter :: d5 =  .4694682253029562039431948525047387412553_wp
    real(wp),parameter :: d6 =  .8065792249988867707634161808995217981443_wp
    real(wp),parameter :: d8 = -.6071194891777959797672951465256217122488_wp
    real(wp),parameter :: d9 =  .5686113944047569241147603178766138153594e-1_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  =     - d9

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9))

        ! check the cached function eval of the last step:
        call me%check_fsal_cache(t,x,f1)

        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+h,   x+h*(b81*f1+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)

        ! last point is cached for the next step:
        xf = x+h*(c1*f1+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8)
        call me%set_fsal_cache(t+h,xf,f9)

        xerr =   h*(e1*f1 + e4*f4 + e5*f5 + e6*f6 + e7*f7 + e8*f8 + e9*f9)

    end associate

    end procedure rkv65e
!*****************************************************************************************

!*****************************************************************************************
!>
!  Tsitouras & Famelis Runge-Kutta 6(5) method.
!
!### Reference
!  * Ch. Tsitouras and I. Th. Famelis,
!    [Phase-Fitted modified Runge-Kutta pairs of orders 6(5)](https://www.researchgate.net/publication/251740152_Phase-Fitted_modified_Runge-Kutta_pairs_of_orders_65),
!    ICNAAM 2006, Crete, 2006.
!  * [More precise rational coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK6/RKcoeff6g_5.pdf)
!
!@note The floating point coefficients here were generated from the rational ones from the 2nd reference.
!
!@note This is a first-same-as-last (FSAL) step.

    module procedure rktf65

    real(wp),parameter :: a2 = 9.28961748633879781420765027322404371584699453551912568306010928961748633879781420765e-2_wp
    real(wp),parameter :: a3 = 1.44578313253012048192771084337349397590361445783132530120481927710843373493975903614e-1_wp
    real(wp),parameter :: a4 = 2.16867469879518072289156626506024096385542168674698795180722891566265060240963855422e-1_wp
    real(wp),parameter :: a5 = 5.68000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a6 = 7.11864406779661016949152542372881355932203389830508474576271186440677966101694915254e-1_wp
    real(wp),parameter :: a7 = 9.95000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp

    real(wp),parameter :: b21 = 9.28961748633879781420765027322404371584699453551912568306010928961748633879781420765e-2_wp
    real(wp),parameter :: b31 = 3.20715889781663863106572284887245651635599805316233039884555941697335052470690700435e-2_wp
    real(wp),parameter :: b32 = 1.12506724274845661882113855848624832426801465251509226132026333541109868246906833571e-1_wp
    real(wp),parameter :: b41 = 5.42168674698795180722891566265060240963855421686746987951807228915662650602409638554e-2_wp
    real(wp),parameter :: b43 = 1.62650602409638554216867469879518072289156626506024096385542168674698795180722891566e-1_wp
    real(wp),parameter :: b51 = 6.56598126617283950617283950617283950617283950617283950617283950617283950617283950617e-1_wp
    real(wp),parameter :: b53 = -2.4972770465185185185185185185185185185185185185185185185185185185185185185185185185e0_wp
    real(wp),parameter :: b54 = 2.40867891990123456790123456790123456790123456790123456790123456790123456790123456790e0_wp
    real(wp),parameter :: b61 = -1.7212338830272756350112289300338046256920059840520750991645177506321959252208613366e0_wp
    real(wp),parameter :: b63 = 7.22310751192021126636920155870658124799459070542113115974697356535677039490201029989e0_wp
    real(wp),parameter :: b64 = -5.4959191727393518000813340726320306436524484390690989489777480799761402436788850115e0_wp
    real(wp),parameter :: b65 = 7.05909950626077185672513986332135377282067107530551362971563451692243740099430963481e-1_wp
    real(wp),parameter :: b71 = 4.12867189716161527469018071566157538750756968750429124243842122794843113666860197270e0_wp
    real(wp),parameter :: b73 = -1.6914025304289438220281346004343739703895247887189472122417356767305063273262951589e1_wp
    real(wp),parameter :: b74 = 1.43289813217405868844323023574174855877382144987680899683823277400972814702596244738e1_wp
    real(wp),parameter :: b75 = -1.5533550903735309941352947451105003919180320410408681830389632181656533578398372284e0_wp
    real(wp),parameter :: b76 = 1.00472717576076705529415767637517912056749574195795909463557101742500402417456237099e0_wp
    real(wp),parameter :: b81 = 4.46860784284429599022330889710676757001193443328341367781153268526799285074027050787e0_wp
    real(wp),parameter :: b83 = -1.8345441869429650066499014099391358577583747335752433112108627235199546393109090809e1_wp
    real(wp),parameter :: b84 = 1.55237707233543199854527260736987323920859616527670939779256214402670438227794675418e1_wp
    real(wp),parameter :: b85 = -1.7228800213316940986071454499455146268823602061289992205388429896934101951354882736e0_wp
    real(wp),parameter :: b86 = 1.08151571740214517180934640319081789154027660669157008569533046194596896725782768317e0_wp
    real(wp),parameter :: b87 = -5.5723928394169823792218246594446491720651508606454087850143625880490525329866504633e-3_wp
    ! real(wp),parameter :: b91   = 6.42309093721083246173486105540744934898802830008152394733186013531402447815498986849e-2_wp ! FSAL
    ! real(wp),parameter :: b94   = 3.32861824699421109595374929094856793240150741775741341700388991831307262888094726546e-1_wp
    ! real(wp),parameter :: b95   = 2.67859229165778080498089104036336501876704241846808785499560813816788199727204331902e-1_wp
    ! real(wp),parameter :: b96   = 1.79863899670938711343662526258687438581919457377762276997401953795124106054682426881e-1_wp
    ! real(wp),parameter :: b97   = 1.51075784805762160255845148219838318965956604578826959670696217369891977986727762107e0_wp
    ! real(wp),parameter :: b98   = -1.3555737109658678286129266521423384168482207697893972403776325344952795933188090051e0_wp

    real(wp),parameter :: c1 = 6.42309093721083246173486105540744934898802830008152394733186013531402447815498986849e-2_wp
    real(wp),parameter :: c4 = 3.32861824699421109595374929094856793240150741775741341700388991831307262888094726546e-1_wp
    real(wp),parameter :: c5 = 2.67859229165778080498089104036336501876704241846808785499560813816788199727204331902e-1_wp
    real(wp),parameter :: c6 = 1.79863899670938711343662526258687438581919457377762276997401953795124106054682426881e-1_wp
    real(wp),parameter :: c7 = 1.51075784805762160255845148219838318965956604578826959670696217369891977986727762107e0_wp
    real(wp),parameter :: c8 = -1.3555737109658678286129266521423384168482207697893972403776325344952795933188090051e0_wp

    real(wp),parameter :: d1 = 6.22980954171238844288765383876964441563622264193747564518688416336360055072506207144e-2_wp
    real(wp),parameter :: d4 = 3.40203519635783669970356076069438053031579043631460538976870898652474922063170157204e-1_wp
    real(wp),parameter :: d5 = 2.35997541364108988248916956975923713134336326558729708330181566604979639253608698633e-1_wp
    real(wp),parameter :: d6 = 2.20640276396042514942275961685396046511990069015617783025687162561106239011510198458e-1_wp
    real(wp),parameter :: d7 = 1.15044133952868915480448599961883762022583301353216360785466101429534020477553812807e0_wp
    real(wp),parameter :: d8 = -1.0029141056750815457282448660706252103934340124906797279726028170808703439444111364e0_wp
    real(wp),parameter :: d9 = -6.6666666666666666666666666666666666666666666666666666666666666666666666666666666667e-3_wp

    real(wp),parameter :: e1 = c1 - d1
    real(wp),parameter :: e4 = c4 - d4
    real(wp),parameter :: e5 = c5 - d5
    real(wp),parameter :: e6 = c6 - d6
    real(wp),parameter :: e7 = c7 - d7
    real(wp),parameter :: e8 = c8 - d8
    real(wp),parameter :: e9 =    - d9

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9))

        ! check the cached function eval of the last step:
        call me%check_fsal_cache(t,x,f1)

        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+h,   x+h*(b81*f1+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)

        ! last point is cached for the next step:
        xf = x+h*(c1*f1+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8)
        call me%set_fsal_cache(t+h,xf,f9)

        xerr = h*(e1*f1+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8+e9*f9)

    end associate

    end procedure rktf65
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's "most robust" Runge-Kutta (9,6(5)) pair.
!
!### Reference
!  * J.H. Verner, "Strategies for deriving new explicit Runge-Kutta pairs",
!    Annals of Num. Math 1 1994, 225-244.
!  * [Jim Verner's Refuge for Runge-Kutta Pairs](https://www.sfu.ca/~jverner/)
!  * [Rational coefficients](https://www.sfu.ca/~jverner/RKV65.IIIXb.Robust.00010102836.081204.CoeffsOnlyRAT)
!  * [Floating point coefficients](https://www.sfu.ca/~jverner/RKV65.IIIXb.Robust.00010102836.081204.CoeffsOnlyFLOAT)
!
!@note This is a first-same-as-last (FSAL) step.

    module procedure rkv65r

    real(wp),parameter :: a2 = .18_wp
    real(wp),parameter :: a3 = .1666666666666666666666666666666666666667_wp
    real(wp),parameter :: a4 = .25_wp
    real(wp),parameter :: a5 = .53_wp
    real(wp),parameter :: a6 = .6_wp
    real(wp),parameter :: a7 = .8_wp

    real(wp),parameter :: b21 =  .18_wp
    real(wp),parameter :: b31 =  .8950617283950617283950617283950617283951e-1_wp
    real(wp),parameter :: b32 =  .7716049382716049382716049382716049382716e-1_wp
    real(wp),parameter :: b41 =  .625e-1_wp
    real(wp),parameter :: b43 =  .1875_wp
    real(wp),parameter :: b51 =  .316516_wp
    real(wp),parameter :: b53 = -1.044948_wp
    real(wp),parameter :: b54 =  1.258432_wp
    real(wp),parameter :: b61 =  .2723261273648562625722506556667430550251_wp
    real(wp),parameter :: b63 = -.8251336032388663967611336032388663967611_wp
    real(wp),parameter :: b64 =  1.048091767881241565452091767881241565452_wp
    real(wp),parameter :: b65 =  .1047157079927685687367911796908817762840_wp
    real(wp),parameter :: b71 = -.1669941859971651431432960727896179733320_wp
    real(wp),parameter :: b73 =  .6317085020242914979757085020242914979757_wp
    real(wp),parameter :: b74 =  .1746104455277387608214675883848816179643_wp
    real(wp),parameter :: b75 = -1.066535645908606612252519473401868067778_wp
    real(wp),parameter :: b76 =  1.227210884353741496598639455782312925170_wp
    real(wp),parameter :: b81 =  .3642375168690958164642375168690958164642_wp
    real(wp),parameter :: b83 = -.2040485829959514170040485829959514170040_wp
    real(wp),parameter :: b84 = -.3488373781606864313631230924464007170774_wp
    real(wp),parameter :: b85 =  3.261932303285686744333360874714258172905_wp
    real(wp),parameter :: b86 = -2.755102040816326530612244897959183673469_wp
    real(wp),parameter :: b87 =  .6818181818181818181818181818181818181818_wp
    ! real(wp),parameter :: b91 =  .7638888888888888888888888888888888888889e-1_wp  ! FSAL
    ! real(wp),parameter :: b92 = 0
    ! real(wp),parameter :: b93 = 0
    ! real(wp),parameter :: b94 =  .3694083694083694083694083694083694083694_wp
    ! real(wp),parameter :: b95 = 0
    ! real(wp),parameter :: b96 =  .2480158730158730158730158730158730158730_wp
    ! real(wp),parameter :: b97 =  .2367424242424242424242424242424242424242_wp
    ! real(wp),parameter :: b98 =  .6944444444444444444444444444444444444444e-1_wp

    real(wp),parameter :: c1 = .7638888888888888888888888888888888888889e-1_wp
    real(wp),parameter :: c4 = .3694083694083694083694083694083694083694_wp
    real(wp),parameter :: c6 = .2480158730158730158730158730158730158730_wp
    real(wp),parameter :: c7 = .2367424242424242424242424242424242424242_wp
    real(wp),parameter :: c8 = .6944444444444444444444444444444444444444e-1_wp

    real(wp),parameter :: d1 =  .5870020964360587002096436058700209643606e-1_wp
    real(wp),parameter :: d4 =  .4807256235827664399092970521541950113379_wp
    real(wp),parameter :: d5 = -.8534124207691908557883209486122831308356_wp
    real(wp),parameter :: d6 =  1.204648526077097505668934240362811791383_wp
    real(wp),parameter :: d8 = -.5924237307216030620285939434875605088371e-1_wp
    real(wp),parameter :: d9 =  .1685804345378813463919846898570302825622_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  =     - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  =     - d9

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9))

        ! check the cached function eval of the last step:
        call me%check_fsal_cache(t,x,f1)

        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+h,   x+h*(b81*f1+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)

        ! last point is cached for the next step:
        xf = x+h*(c1*f1+c4*f4+c6*f6+c7*f7+c8*f8)
        call me%set_fsal_cache(t+h,xf,f9)

        xerr = h*(e1*f1+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8+e9*f9)

    end associate

    end procedure rkv65r
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner 6(5) method.
!
!### References
!  * A Contrast of a New RK56 pair with DP56, by Jim Verner,
!    Department of Mathematics. Simon Fraser University, Burnaby, Canada, 2006.
!  * [Coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK6/RKcoeff6e_3.pdf)

    module procedure rkv65

    real(wp),parameter :: a2 = 1.0_wp / 7.0_wp
    real(wp),parameter :: a3 = 2.0_wp / 9.0_wp
    real(wp),parameter :: a4 = 3.0_wp / 7.0_wp
    real(wp),parameter :: a5 = 2.0_wp / 3.0_wp
    real(wp),parameter :: a6 = 3.0_wp / 4.0_wp

    real(wp),parameter :: b21 = 1.0_wp / 7.0_wp
    real(wp),parameter :: b31 = 4.0_wp / 81.0_wp
    real(wp),parameter :: b32 = 14.0_wp / 81.0_wp
    real(wp),parameter :: b41 = 291.0_wp / 1372.0_wp
    real(wp),parameter :: b42 = -27.0_wp / 49.0_wp
    real(wp),parameter :: b43 = 1053.0_wp / 1372.0_wp
    real(wp),parameter :: b51 = 86.0_wp / 297.0_wp
    real(wp),parameter :: b52 = -14.0_wp / 33.0_wp
    real(wp),parameter :: b53 = 42.0_wp / 143.0_wp
    real(wp),parameter :: b54 = 1960.0_wp / 3861.0_wp
    real(wp),parameter :: b61 = -267.0_wp / 22528.0_wp
    real(wp),parameter :: b62 = 189.0_wp / 704.0_wp
    real(wp),parameter :: b63 = 63099.0_wp / 585728.0_wp
    real(wp),parameter :: b64 = 58653.0_wp / 366080.0_wp
    real(wp),parameter :: b65 = 4617.0_wp / 20480.0_wp
    real(wp),parameter :: b71 = 10949.0_wp / 6912.0_wp
    real(wp),parameter :: b72 = -69.0_wp / 32.0_wp
    real(wp),parameter :: b73 = -90891.0_wp / 68096.0_wp
    real(wp),parameter :: b74 = 112931.0_wp / 25920.0_wp
    real(wp),parameter :: b75 = -69861.0_wp / 17920.0_wp
    real(wp),parameter :: b76 = 26378.0_wp / 10773.0_wp
    real(wp),parameter :: b81 = 1501.0_wp / 19008.0_wp
    real(wp),parameter :: b82 = -21.0_wp / 88.0_wp
    real(wp),parameter :: b83 = 219519.0_wp / 347776.0_wp
    real(wp),parameter :: b84 = 163807.0_wp / 926640.0_wp
    real(wp),parameter :: b85 = -417.0_wp / 640.0_wp
    real(wp),parameter :: b86 = 1544.0_wp / 1539.0_wp

    real(wp),parameter :: c1 = 79.0_wp / 1080.0_wp
    real(wp),parameter :: c3 = 19683.0_wp / 69160.0_wp
    real(wp),parameter :: c4 = 16807.0_wp / 84240.0_wp
    real(wp),parameter :: c6 = 2816.0_wp / 7695.0_wp
    real(wp),parameter :: c7 = 1.0_wp / 100.0_wp
    real(wp),parameter :: c8 = 187.0_wp / 2800.0_wp

    real(wp),parameter :: d1 = 763.0_wp / 10800.0_wp
    real(wp),parameter :: d3 = 59049.0_wp / 197600.0_wp
    real(wp),parameter :: d4 = 88837.0_wp / 526500.0_wp
    real(wp),parameter :: d5 = 243.0_wp / 4000.0_wp
    real(wp),parameter :: d6 = 12352.0_wp / 38475.0_wp
    real(wp),parameter :: d8 = 2.0_wp / 25.0_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  =     - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7
    real(wp),parameter :: e8  = c8  - d8

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+h,   x+h*(b71*f1+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+h,   x+h*(b81*f1+b82*f2+b83*f3+b84*f4+b85*f5+b86*f6),f8)

        xf = x+h*(c1*f1+c3*f3+c4*f4+c6*f6+c7*f7+c8*f8)

        xerr = h*(e1*f1+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8)

    end associate

    end procedure rkv65
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner 6(5) 'DVERK' method.

    module procedure dverk65

    real(wp),parameter :: a2 = 1.0_wp / 6.0_wp
    real(wp),parameter :: a3 = 4.0_wp / 15.0_wp
    real(wp),parameter :: a4 = 2.0_wp / 3.0_wp
    real(wp),parameter :: a5 = 5.0_wp / 6.0_wp
    real(wp),parameter :: a7 = 1.0_wp / 15.0_wp

    real(wp),parameter :: b21 = 1.0_wp      / 6.0_wp
    real(wp),parameter :: b31 = 4.0_wp      / 75.0_wp
    real(wp),parameter :: b32 = 16.0_wp     / 75.0_wp
    real(wp),parameter :: b41 = 5.0_wp      / 6.0_wp
    real(wp),parameter :: b42 = -8.0_wp     / 3.0_wp
    real(wp),parameter :: b43 = 5.0_wp      / 2.0_wp
    real(wp),parameter :: b51 = -165.0_wp   / 64.0_wp
    real(wp),parameter :: b52 = 55.0_wp     / 6.0_wp
    real(wp),parameter :: b53 = -425.0_wp   / 64.0_wp
    real(wp),parameter :: b54 = 85.0_wp     / 96.0_wp
    real(wp),parameter :: b61 = 12.0_wp     / 5.0_wp
    real(wp),parameter :: b62 = -8.0_wp
    real(wp),parameter :: b63 = 4015.0_wp   / 612.0_wp
    real(wp),parameter :: b64 = -11.0_wp    / 36.0_wp
    real(wp),parameter :: b65 = 88.0_wp     / 255.0_wp
    real(wp),parameter :: b71 = -8263.0_wp  / 15000.0_wp
    real(wp),parameter :: b72 = 124.0_wp    / 75.0_wp
    real(wp),parameter :: b73 = -643.0_wp   / 680.0_wp
    real(wp),parameter :: b74 = -81.0_wp    / 250.0_wp
    real(wp),parameter :: b75 = 2484.0_wp   / 10625.0_wp
    real(wp),parameter :: b81 = 3501.0_wp   / 1720.0_wp
    real(wp),parameter :: b82 = -300.0_wp   / 43.0_wp
    real(wp),parameter :: b83 = 297275.0_wp / 52632.0_wp
    real(wp),parameter :: b84 = -319.0_wp   / 2322.0_wp
    real(wp),parameter :: b85 = 24068.0_wp  / 84065.0_wp
    real(wp),parameter :: b87 = 3850.0_wp   / 26703.0_wp

    real(wp),parameter :: c1 = 3.0_wp   / 40.0_wp
    real(wp),parameter :: c3 = 875.0_wp / 2244.0_wp
    real(wp),parameter :: c4 = 23.0_wp  / 72.0_wp
    real(wp),parameter :: c5 = 264.0_wp / 1955.0_wp
    real(wp),parameter :: c7 = 125.0_wp / 11592.0_wp
    real(wp),parameter :: c8 = 43.0_wp  / 616.0_wp

    real(wp),parameter :: d1 = 13.0_wp   / 160.0_wp
    real(wp),parameter :: d3 = 2375.0_wp / 5984.0_wp
    real(wp),parameter :: d4 = 5.0_wp    / 16.0_wp
    real(wp),parameter :: d5 = 12.0_wp   / 85.0_wp
    real(wp),parameter :: d6 = 3.0_wp    / 44.0_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  =     - d6
    real(wp),parameter :: e7  = c7
    real(wp),parameter :: e8  = c8

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+h,   x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b72*f2+b73*f3+b74*f4+b75*f5),f7)
        call me%f(t+h,   x+h*(b81*f1+b82*f2+b83*f3+b84*f4+b85*f5+b87*f7),f8)

        xf = x+h*(c1*f1+c3*f3+c4*f4+c5*f5+c7*f7+c8*f8)

        xerr = h*(e1*f1+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8)

    end associate

    end procedure dverk65
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's "most efficient" Runge-Kutta (10:7(6)) pair.
!
!### Reference
!  * J.H. Verner, "Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error",
!    SIAM Journal on Numerical Analysis, Vol. 15, Issue 4 (1978)
!  * [Jim Verner's Refuge for Runge-Kutta Pairs](https://www.sfu.ca/~jverner/)
!  * [Coefficients](https://www.sfu.ca/~jverner/RKV76.IIa.Efficient.00001675585.081206.OnWeb)

    module procedure rkv76e

    real(wp),parameter :: a2 = .5e-2_wp
    real(wp),parameter :: a3 = .1088888888888888888888888888888888888889_wp
    real(wp),parameter :: a4 = .1633333333333333333333333333333333333333_wp
    real(wp),parameter :: a5 = .4555_wp
    real(wp),parameter :: a6 = .6095094489978381317087004421486024949638_wp
    real(wp),parameter :: a7 = .884_wp
    real(wp),parameter :: a8 = .925_wp

    real(wp),parameter :: b21 =  .5e-2_wp
    real(wp),parameter :: b31 = -1.076790123456790123456790123456790123457_wp
    real(wp),parameter :: b32 =  1.185679012345679012345679012345679012346_wp
    real(wp),parameter :: b41 =  .4083333333333333333333333333333333333333e-1_wp
    real(wp),parameter :: b43 =  .1225_wp
    real(wp),parameter :: b51 =  .6389139236255726780508121615993336109954_wp
    real(wp),parameter :: b53 = -2.455672638223656809662640566430653894211_wp
    real(wp),parameter :: b54 =  2.272258714598084131611828404831320283215_wp
    real(wp),parameter :: b61 = -2.661577375018757131119259297861818119279_wp
    real(wp),parameter :: b63 =  10.80451388645613769565396655365532838482_wp
    real(wp),parameter :: b64 = -8.353914657396199411968048547819291691541_wp
    real(wp),parameter :: b65 =  .8204875949566569791420417341743839209619_wp
    real(wp),parameter :: b71 =  6.067741434696770992718360183877276714679_wp
    real(wp),parameter :: b73 = -24.71127363591108579734203485290746001803_wp
    real(wp),parameter :: b74 =  20.42751793078889394045773111748346612697_wp
    real(wp),parameter :: b75 = -1.906157978816647150624096784352757010879_wp
    real(wp),parameter :: b76 =  1.006172249242068014790040335899474187268_wp
    real(wp),parameter :: b81 =  12.05467007625320299509109452892778311648_wp
    real(wp),parameter :: b83 = -49.75478495046898932807257615331444758322_wp
    real(wp),parameter :: b84 =  41.14288863860467663259698416710157354209_wp
    real(wp),parameter :: b85 = -4.461760149974004185641911603484815375051_wp
    real(wp),parameter :: b86 =  2.042334822239174959821717077708608543738_wp
    real(wp),parameter :: b87 = -0.9834843665406107379530801693870224403537e-1_wp
    real(wp),parameter :: b91 =  10.13814652288180787641845141981689030769_wp
    real(wp),parameter :: b93 = -42.64113603171750214622846006736635730625_wp
    real(wp),parameter :: b94 =  35.76384003992257007135021178023160054034_wp
    real(wp),parameter :: b95 = -4.348022840392907653340370296908245943710_wp
    real(wp),parameter :: b96 =  2.009862268377035895441943593011827554771_wp
    real(wp),parameter :: b97 =  .3487490460338272405953822853053145879140_wp
    real(wp),parameter :: b98 = -.2714390051048312842371587140910297407572_wp
    real(wp),parameter :: b101 = -45.03007203429867712435322405073769635151_wp
    real(wp),parameter :: b103 =  187.3272437654588840752418206154201997384_wp
    real(wp),parameter :: b104 = -154.0288236935018690596728621034510402582_wp
    real(wp),parameter :: b105 =  18.56465306347536233859492332958439136765_wp
    real(wp),parameter :: b106 = -7.141809679295078854925420496823551192821_wp
    real(wp),parameter :: b107 =  1.308808578161378625114762706007696696508_wp

    real(wp),parameter :: c1 =  .4715561848627222170431765108838175679569e-1_wp
    real(wp),parameter :: c4 =  .2575056429843415189596436101037687580986_wp
    real(wp),parameter :: c5 =  .2621665397741262047713863095764527711129_wp
    real(wp),parameter :: c6 =  .1521609265673855740323133199165117535523_wp
    real(wp),parameter :: c7 =  .4939969170032484246907175893227876844296_wp
    real(wp),parameter :: c8 = -.2943031171403250441557244744092703429139_wp
    real(wp),parameter :: c9 =  .8131747232495109999734599440136761892478e-1_wp

    real(wp),parameter :: d1 =  .4460860660634117628731817597479197781432e-1_wp
    real(wp),parameter :: d4 =  .2671640378571372680509102260943837899738_wp
    real(wp),parameter :: d5 =  .2201018300177293019979715776650753096323_wp
    real(wp),parameter :: d6 =  .2188431703143156830983120833512893824578_wp
    real(wp),parameter :: d7 =  .2289871705411202883378173889763552365362_wp
    real(wp),parameter :: d10 =  .2029518466335628222767054793810430358554e-1_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8
    real(wp),parameter :: e9  = c9
    real(wp),parameter :: e10 =    - d10

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10))

        call me%f(t,       x,f1)
        call me%f(t+a2*h,  x+h*(b21*f1),f2)
        call me%f(t+a3*h,  x+h*(b31*f1 + b32*f2),f3)
        call me%f(t+a4*h,  x+h*(b41*f1          + b43*f3),f4)
        call me%f(t+a5*h,  x+h*(b51*f1          + b53*f3 + b54*f4),f5)
        call me%f(t+a6*h,  x+h*(b61*f1          + b63*f3 + b64*f4 + b65*f5),f6)
        call me%f(t+a7*h,  x+h*(b71*f1          + b73*f3 + b74*f4 + b75*f5 + b76*f6),f7)
        call me%f(t+a8*h,  x+h*(b81*f1          + b83*f3 + b84*f4 + b85*f5 + b86*f6 + b87*f7),f8)
        call me%f(t+h,     x+h*(b91*f1          + b93*f3 + b94*f4 + b95*f5 + b96*f6 + b97*f7 + b98*f8),f9)
        call me%f(t+h,     x+h*(b101*f1         + b103*f3+ b104*f4+ b105*f5+ b106*f6+ b107*f7 ),f10)

        xf = x + h*(c1*f1 + c4*f4 + c5*f5 + c6*f6 + c7*f7 + c8*f8 + c9*f9)

        xerr =   h*(e1*f1 + e4*f4 + e5*f5 + e6*f6 + e7*f7 + e8*f8 + e9*f9 + e10*f10)

    end associate

    end procedure rkv76e
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's "most robust" Runge-Kutta (10:7(6)) pair.
!
!### Reference
!  * J.H. Verner, "Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error",
!    SIAM Journal on Numerical Analysis, Vol. 15, Issue 4 (1978)
!  * [Jim Verner's Refuge for Runge-Kutta Pairs](https://www.sfu.ca/~jverner/)
!  * [Coefficients](https://www.sfu.ca/~jverner/RKV76.IIa.Robust.000027015646.081206.CoeffsOnlyFLOAT)

    module procedure rkv76r

    real(wp),parameter :: a2  = .5e-2_wp
    real(wp),parameter :: a3  = .1088888888888888888888888888888888888889_wp
    real(wp),parameter :: a4  = .1633333333333333333333333333333333333333_wp
    real(wp),parameter :: a5  = .455_wp
    real(wp),parameter :: a6  = .6059617471462913245758145021744683294809_wp
    real(wp),parameter :: a7  = .835_wp
    real(wp),parameter :: a8  = .915_wp

    real(wp),parameter :: b21  =  .5e-2_wp
    real(wp),parameter :: b31  = -1.076790123456790123456790123456790123457_wp
    real(wp),parameter :: b32  =  1.185679012345679012345679012345679012346_wp
    real(wp),parameter :: b41  =  .4083333333333333333333333333333333333333e-1_wp
    real(wp),parameter :: b43  =  .1225_wp
    real(wp),parameter :: b51  =  .6360714285714285714285714285714285714286_wp
    real(wp),parameter :: b53  = -2.444464285714285714285714285714285714286_wp
    real(wp),parameter :: b54  =  2.263392857142857142857142857142857142857_wp
    real(wp),parameter :: b61  = -2.535121107934924522925638355466021548721_wp
    real(wp),parameter :: b63  =  10.29937465444926792043851446075602491361_wp
    real(wp),parameter :: b64  = -7.951303288599057994949321745826687653648_wp
    real(wp),parameter :: b65  =  .7930114892310059220122601427111526182380_wp
    real(wp),parameter :: b71  =  1.001876581252463296196919658309499980821_wp
    real(wp),parameter :: b73  = -4.166571282442379833131393800547097145319_wp
    real(wp),parameter :: b74  =  3.834343292912864241255266521825137866520_wp
    real(wp),parameter :: b75  = -.5023333356071084754746433022861176561240_wp
    real(wp),parameter :: b76  =  .6676847438841607711538509226985769541026_wp
    real(wp),parameter :: b81  =  27.25501835463076713033396381917500571735_wp
    real(wp),parameter :: b83  = -42.00461727841063835531864544390929536961_wp
    real(wp),parameter :: b84  = -10.53571312661948991792108160054652610372_wp
    real(wp),parameter :: b85  =  80.49553671141193714798365215892682663420_wp
    real(wp),parameter :: b86  = -67.34388227179051346854907596321297564093_wp
    real(wp),parameter :: b87  =  13.04865761077793746347118702956696476271_wp
    real(wp),parameter :: b91  = -3.039737805711496514694365865875576322688_wp
    real(wp),parameter :: b93  =  10.13816141032980111185794619070970015044_wp
    real(wp),parameter :: b94  = -6.429305674864721572146282562955529806444_wp
    real(wp),parameter :: b95  = -1.586437148340827658711531285379861057947_wp
    real(wp),parameter :: b96  =  1.892178184196842441086430890913135336502_wp
    real(wp),parameter :: b97  =  .1969933540760886906129236016333644283801e-1_wp
    real(wp),parameter :: b98  =  .5441698982793323546510272424795257297790e-2_wp
    real(wp),parameter :: b101 = -1.444951891677773513735100317935571236052_wp
    real(wp),parameter :: b103 =  8.031891385995591922411703322301956043504_wp
    real(wp),parameter :: b104 = -7.583174166340134682079888302367158860498_wp
    real(wp),parameter :: b105 =  3.581616935319007421124768544245287869686_wp
    real(wp),parameter :: b106 = -2.436972263219952941118380906569375238373_wp
    real(wp),parameter :: b107 =  .8515899999232617933968976603248614217339_wp

    real(wp),parameter :: c1  = .4742583783370675608356917271757453469893e-1_wp
    real(wp),parameter :: c4  = .2562236165937056265996172745827462344816_wp
    real(wp),parameter :: c5  = .2695137683307420661947381725807595288676_wp
    real(wp),parameter :: c6  = .1268662240909278284598913836473917324788_wp
    real(wp),parameter :: c7  = .2488722594206007162204644942764749276729_wp
    real(wp),parameter :: c8  = .3074483740820063133530438847909918476864e-2_wp
    real(wp),parameter :: c9  = .4802380998949694330818906334714312332321e-1_wp

    real(wp),parameter :: d1  = .4748524769929963103753127380572796155227e-1_wp
    real(wp),parameter :: d4  = .2559941258869063329715491824590539387050_wp
    real(wp),parameter :: d5  = .2705847808106768872253089109926813573239_wp
    real(wp),parameter :: d6  = .1250561868442599291363882232374691792045_wp
    real(wp),parameter :: d7  = .2520446872374386050718404382019744256218_wp
    real(wp),parameter :: d10 = .4883497152141861455738197130309313759259e-1_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8
    real(wp),parameter :: e9  = c9
    real(wp),parameter :: e10 =     - d10

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+h,   x+h*(b91*f1+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+h,   x+h*(b101*f1+b103*f3+b104*f4+b105*f5+b106*f6+b107*f7),f10)

        xf = x+h*(c1*f1+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8+c9*f9)

        xerr = h*(e1*f1+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10)

    end associate

    end procedure rkv76r
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sharp & Smart 7(6)) method.
!
!### Reference
!  * P.W. Sharp and E. Smart, "Explicit Runge-Kutta Pairs with One More Derivative Evaluation than the Minimum",
!    Siam Journal of Scientific Computing, Vol. 14, No. 2, pages. 338-348, March 1993.

    module procedure rkss76

    real(wp),parameter :: a2    = 2.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: a3    = 2.16000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a4    = 4.10000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a5    = 5.70000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a6    = 8.60000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a7    = 1.89734415567554714447412739459850829747319955773081476233291234325271903458210730096e-1_wp
    real(wp),parameter :: a8    = 7.20000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a9    = 8.33333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1_wp

    real(wp),parameter :: b21   = 2.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: b31   = -9.5040000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: b32   = 1.16640000000000000000000000000000000000000000000000000000000000000000000000000000000e0_wp
    real(wp),parameter :: b41   = 2.08796296296296296296296296296296296296296296296296296296296296296296296296296296296e-2_wp
    real(wp),parameter :: b43   = 3.89120370370370370370370370370370370370370370370370370370370370370370370370370370370e-1_wp
    real(wp),parameter :: b51   = 1.18750000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: b53   = 1.16301546391752577319587628865979381443298969072164948453608247422680412371134020619e-1_wp
    real(wp),parameter :: b54   = 3.34948453608247422680412371134020618556701030927835051546391752577319587628865979381e-1_wp
    real(wp),parameter :: b61   = -1.5500000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: b63   = 1.26419952627798163356436989147106218572232900382472285320441881685984429300871726173e0_wp
    real(wp),parameter :: b64   = -1.4923539518900343642611683848797250859106529209621993127147766323024054982817869416e0_wp
    real(wp),parameter :: b65   = 1.24315442561205273069679849340866290018832391713747645951035781544256120527306967985e0_wp
    real(wp),parameter :: b71   = 5.32404818460847244880465705884559948010013316415508489398809027494008261005225044158e-2_wp
    real(wp),parameter :: b73   = 2.63190905005944241395058267679896015089060637096436681624641838124193103498139523636e-1_wp
    real(wp),parameter :: b74   = -2.2992382584611251312628854152112504080421217996536744225053579460154724416916797118e-1_wp
    real(wp),parameter :: b75   = 1.15020646834503298021503942059822317395934020502096307719511499416916110119759292822e-1_wp
    real(wp),parameter :: b76   = -1.1793792272865036330907499347198456734463853501634919800207211363690892091042619602e-2_wp
    real(wp),parameter :: b81   = 1.01833250525779410925056563242748813162631092358712093910139906641055533787858600752e-1_wp
    real(wp),parameter :: b83   = 3.67650402264021539690171893675768442020160413827130555257229267948121503927894530089e0_wp
    real(wp),parameter :: b84   = -1.3855064886371957262483408369918407532890732638371784299323099278110242610274153464e0_wp
    real(wp),parameter :: b85   = 1.10585536543242802574001937437432027003273412482516272090263172771172002666187126910e0_wp
    real(wp),parameter :: b86   = -1.0442483660597752261508675907091877183025101004713070424083201191960255004658089509e-1_wp
    real(wp),parameter :: b87   = -2.6742613133552495847033672783119939782776450815708712332119223741033637886546789293e0_wp
    real(wp),parameter :: b91   = -7.0447278725441551253173322525085182085211898876535016082680181826378260058107213334e-2_wp
    real(wp),parameter :: b93   = 2.57849022691031881978588351318784339446713208794819222430391001640703602812013523694e0_wp
    real(wp),parameter :: b94   = -1.7350648742705534137364675443051006373149953263240864527861457323256158060409724050e0_wp
    real(wp),parameter :: b95   = 1.42076925842656275743055431879937793726979032667289381452783078677285932226590451543e0_wp
    real(wp),parameter :: b96   = -2.9994002356151422690421464715731127083992859674624280646321971569589440552603013442e-2_wp
    real(wp),parameter :: b97   = -1.2738969184577295296576255258337856595600446902141816152083573928383427523639249936e0_wp
    real(wp),parameter :: b98   = -5.6523078193672326545416641274185392359344306198325340774902191286635758037098793739e-2_wp
    real(wp),parameter :: b101  = 1.12188929566497595110309434330843411185625863434326301999225595935353619019786372693e-1_wp
    real(wp),parameter :: b103  = -9.4478179712463597825861610011719624575041188851623316508747159935931144002457572402e0_wp
    real(wp),parameter :: b104  = 4.72928685969506458483541566908828462557623316160618665528745761537386487691574498692e0_wp
    real(wp),parameter :: b105  = -2.1535862973784519387089628138790863673524951190739517405566618338182156484450238408e0_wp
    real(wp),parameter :: b106  = -1.5003176034455324417891524003264793431625218068048301751161484381435733494665524146e0_wp
    real(wp),parameter :: b107  = 7.32228493516381472239660644916433114018057892752031963313864026830512965942191992469e0_wp
    real(wp),parameter :: b108  = -1.0171881348813826762245691938848731969668479025856205360601503764094927919654990399e0_wp
    real(wp),parameter :: b109  = 2.95514928252634993696651385667894218804354576106590151218235316235004803476538125126e0_wp
    real(wp),parameter :: b111  = 1.26777165953424622263709262280041924195783046360695933368685779852909039382377368327e-2_wp
    real(wp),parameter :: b113  = -1.0069960582593041001624297555860950454903273549350063548770708526802347773053466643e1_wp
    real(wp),parameter :: b114  = 4.62534242829213472690316528870584469183038900200033492988915271866480429986075212694e0_wp
    real(wp),parameter :: b115  = -2.1610283337552771303871477895705714199299596209426922872912146988024763221546690181e0_wp
    real(wp),parameter :: b116  = -1.6781218510451494897112741151843512565528970972689550849978046265492823347869148500e0_wp
    real(wp),parameter :: b117  = 8.10720558296219407284479196100506119564876298707858160850972380908658055161496706221e0_wp
    real(wp),parameter :: b118  = -9.9034046697481677826772383527914286511512256496144470729835169541686418859152250889e-1_wp
    real(wp),parameter :: b119  = 3.15422550651861313801611511995610591660252253880816949662233444183429486317261609431e0_wp

    real(wp),parameter :: c1    = 5.49699932022185175241442181883597507655769828815946565807124447154669050002194891238e-2_wp
    real(wp),parameter :: c5    = 5.27017899545901255716693676497458451526562246965950291388111745594043275344112358966e-1_wp
    real(wp),parameter :: c6    = -5.7640010460253178146313005308307447324854309728047707823605247362054750175374841292e-1_wp
    real(wp),parameter :: c7    = 2.99473766075857495911574410718252874676717226006797736321538429237884690787951256718e-1_wp
    real(wp),parameter :: c8    = -4.3511367243178106005606352407152054011086352330218025896985334826437712480556570901e-1_wp
    real(wp),parameter :: c9    = 1.06972120095994204935876677034564950812746194295206963084848042699967266056102464470e0_wp
    real(wp),parameter :: c10   = 6.03309172503935230080145014048744282630882217762450220670627753378570948660063724257e-2_wp

    real(wp),parameter :: d1    = 5.51948596604973952360558600324258025823235723557842084969195331387651217650789597054e-2_wp
    real(wp),parameter :: d5    = 5.38771472959210564343392307526371633391618262286599376534904900758161647057057668278e-1_wp
    real(wp),parameter :: d6    = -7.2002415738222743602899270061101901552733320251540477745765896611188358469552535128e-1_wp
    real(wp),parameter :: d7    = 2.98535385184885222079230265390696358266728601821898916299256942162002807672289159564e-1_wp
    real(wp),parameter :: d8    = -4.8443496217457505118337730236228616963449840154387472665426235404619128753787634184e-1_wp
    real(wp),parameter :: d9    = 1.24611378035303235082118128195796776952198421286248671471499632269996834100646561750e0_wp
    real(wp),parameter :: d11   = 6.58436213991769547325102880658436213991769547325102880658436213991769547325102880658e-2_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10  = c10
    real(wp),parameter :: e11  =    - d11

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+h,   x+h*(b101*f1+b103*f3+b104*f4+b105*f5+b106*f6+b107*f7+b108*f8+b109*f9),f10)
        call me%f(t+h,   x+h*(b111*f1+b113*f3+b114*f4+b115*f5+b116*f6+b117*f7+b118*f8+b119*f9),f11)

        xf = x+h*(c1*f1+c5*f5+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10)

        xerr = h*(e1*f1+e5*f5+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11)

    end associate

    end procedure rkss76
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's "most efficient" Runge-Kutta (13:8(7)) pair.
!
!### Reference
!  * J.H. Verner, "Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error",
!    SIAM Journal on Numerical Analysis, Vol. 15, Issue 4 (1978)
!  * [Jim Verner's Refuge for Runge-Kutta Pairs](https://www.sfu.ca/~jverner/)
!  * [Coefficients](https://www.sfu.ca/~jverner/RKV87.IIa.Efficient.000000282866.081208.FLOAT40OnWeb)

    module procedure rkv87e

    real(wp),parameter :: a2    = .5e-1_wp
    real(wp),parameter :: a3    = .1065625_wp
    real(wp),parameter :: a4    = .15984375_wp
    real(wp),parameter :: a5    = .39_wp
    real(wp),parameter :: a6    = .465_wp
    real(wp),parameter :: a7    = .155_wp
    real(wp),parameter :: a8    = .943_wp
    real(wp),parameter :: a9    = .9018020417358569582597079406783721499560_wp
    real(wp),parameter :: a10   = .909_wp
    real(wp),parameter :: a11   = .94_wp

    real(wp),parameter :: b21   =  .5e-1_wp
    real(wp),parameter :: b31   = -.69931640625e-2_wp
    real(wp),parameter :: b32   =  .1135556640625_wp
    real(wp),parameter :: b41   =  .399609375e-1_wp
    real(wp),parameter :: b43   =  .1198828125_wp
    real(wp),parameter :: b51   =  .3613975628004575124052940721184028345129_wp
    real(wp),parameter :: b53   = -1.341524066700492771819987788202715834917_wp
    real(wp),parameter :: b54   =  1.370126503900035259414693716084313000404_wp
    real(wp),parameter :: b61   =  .4904720279720279720279720279720279720280e-1_wp
    real(wp),parameter :: b64   =  .2350972042214404739862988335493427143122_wp
    real(wp),parameter :: b65   =  .1808555929813567288109039636534544884850_wp
    real(wp),parameter :: b71   =  .6169289044289044289044289044289044289044e-1_wp
    real(wp),parameter :: b74   =  .1123656831464027662262557035130015442303_wp
    real(wp),parameter :: b75   = -.3885046071451366767049048108111244567456e-1_wp
    real(wp),parameter :: b76   =  .1979188712522045855379188712522045855379e-1_wp
    real(wp),parameter :: b81   = -1.767630240222326875735597119572145586714_wp
    real(wp),parameter :: b84   = -62.5_wp
    real(wp),parameter :: b85   = -6.061889377376669100821361459659331999758_wp
    real(wp),parameter :: b86   =  5.650823198222763138561298030600840174201_wp
    real(wp),parameter :: b87   =  65.62169641937623283799566054863063741227_wp
    real(wp),parameter :: b91   = -1.180945066554970799825116282628297957882_wp
    real(wp),parameter :: b94   = -41.50473441114320841606641502701994225874_wp
    real(wp),parameter :: b95   = -4.434438319103725011225169229846100211776_wp
    real(wp),parameter :: b96   =  4.260408188586133024812193710744693240761_wp
    real(wp),parameter :: b97   =  43.75364022446171584987676829438379303004_wp
    real(wp),parameter :: b98   =  .7871425489912310687446475044226307550860e-2_wp
    real(wp),parameter :: b101  = -1.281405999441488405459510291182054246266_wp
    real(wp),parameter :: b104  = -45.04713996013986630220754257136007322267_wp
    real(wp),parameter :: b105  = -4.731362069449576477311464265491282810943_wp
    real(wp),parameter :: b106  =  4.514967016593807841185851584597240996214_wp
    real(wp),parameter :: b107  =  47.44909557172985134869022392235929015114_wp
    real(wp),parameter :: b108  =  .1059228297111661135687393955516542875228e-1_wp
    real(wp),parameter :: b109  = -.5746842263844616254432318478286296232021e-2_wp
    real(wp),parameter :: b111  = -1.724470134262485191756709817484481861731_wp
    real(wp),parameter :: b114  = -60.92349008483054016518434619253765246063_wp
    real(wp),parameter :: b115  = -5.951518376222392455202832767061854868290_wp
    real(wp),parameter :: b116  =  5.556523730698456235979791650843592496839_wp
    real(wp),parameter :: b117  =  63.98301198033305336837536378635995939281_wp
    real(wp),parameter :: b118  =  .1464202825041496159275921391759452676003e-1_wp
    real(wp),parameter :: b119  =  .6460408772358203603621865144977650714892e-1_wp
    real(wp),parameter :: b1110 = -.7930323169008878984024452548693373291447e-1_wp
    real(wp),parameter :: b121  = -3.301622667747079016353994789790983625569_wp
    real(wp),parameter :: b124  = -118.0112723597525085666923303957898868510_wp
    real(wp),parameter :: b125  = -10.14142238845611248642783916034510897595_wp
    real(wp),parameter :: b126  =  9.139311332232057923544012273556827000619_wp
    real(wp),parameter :: b127  =  123.3759428284042683684847180986501894364_wp
    real(wp),parameter :: b128  =  4.623244378874580474839807625067630924792_wp
    real(wp),parameter :: b129  = -3.383277738068201923652550971536811240814_wp
    real(wp),parameter :: b1210 =  4.527592100324618189451265339351129035325_wp
    real(wp),parameter :: b1211 = -5.828495485811622963193088019162985703755_wp
    real(wp),parameter :: b131  = -3.039515033766309030040102851821200251056_wp
    real(wp),parameter :: b134  = -109.2608680894176254686444192322164623352_wp
    real(wp),parameter :: b135  = -9.290642497400293449717665542656897549158_wp
    real(wp),parameter :: b136  =  8.430504981764911142134299253836167803454_wp
    real(wp),parameter :: b137  =  114.2010010378331313557424041095523427476_wp
    real(wp),parameter :: b138  = -.9637271342145479358162375658987901652762_wp
    real(wp),parameter :: b139  = -5.034884088802189791198680336183332323118_wp
    real(wp),parameter :: b1310 =  5.958130824002923177540402165388172072794_wp

    real(wp),parameter :: c1    =  .4427989419007951074716746668098518862111e-1_wp
    real(wp),parameter :: c6    =  .3541049391724448744815552028733568354121_wp
    real(wp),parameter :: c7    =  .2479692154956437828667629415370663023884_wp
    real(wp),parameter :: c8    = -15.69420203883808405099207034271191213468_wp
    real(wp),parameter :: c9    =  25.08406496555856261343930031237186278518_wp
    real(wp),parameter :: c10   = -31.73836778626027646833156112007297739997_wp
    real(wp),parameter :: c11   =  22.93828327398878395231483560344797018313_wp
    real(wp),parameter :: c12   = -.2361324633071542145259900641263517600737_wp

    real(wp),parameter :: d1    =  .4431261522908979212486436510209029764893e-1_wp
    real(wp),parameter :: d6    =  .3546095642343226447863179350895055038855_wp
    real(wp),parameter :: d7    =  .2478480431366653069619986721504458660016_wp
    real(wp),parameter :: d8    =  4.448134732475784492725128317159648871312_wp
    real(wp),parameter :: d9    =  19.84688636611873369930932399297687935291_wp
    real(wp),parameter :: d10   = -23.58162337746561841969517960870394965085_wp
    real(wp),parameter :: d13   = -.3601679437289775162124536737746202409110_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10 = c10 - d10
    real(wp),parameter :: e11 = c11
    real(wp),parameter :: e12 = c12
    real(wp),parameter :: e13 =     - d13

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13))

        call me%f(t,       x,f1)
        call me%f(t+a2*h,  x+h*(b21*f1),f2)
        call me%f(t+a3*h,  x+h*(b31*f1 + b32*f2),f3)
        call me%f(t+a4*h,  x+h*(b41*f1 + b43*f3),f4)
        call me%f(t+a5*h,  x+h*(b51*f1 + b53*f3 + b54*f4),f5)
        call me%f(t+a6*h,  x+h*(b61*f1 + b64*f4 + b65*f5),f6)
        call me%f(t+a7*h,  x+h*(b71*f1 + b74*f4 + b75*f5 + b76*f6),f7)
        call me%f(t+a8*h,  x+h*(b81*f1 + b84*f4 + b85*f5 + b86*f6 + b87*f7),f8)
        call me%f(t+a9*h,  x+h*(b91*f1 + b94*f4 + b95*f5 + b96*f6 + b97*f7 + b98*f8),f9)
        call me%f(t+a10*h, x+h*(b101*f1+ b104*f4+ b105*f5+ b106*f6+ b107*f7+ b108*f8+ &
                                b109*f9),f10)
        call me%f(t+a11*h, x+h*(b111*f1+ b114*f4+ b115*f5+ b116*f6+ b117*f7+ b118*f8+ &
                                b119*f9+ b1110*f10),f11)
        call me%f(t+h,     x+h*(b121*f1+ b124*f4+ b125*f5+ b126*f6+ b127*f7+ b128*f8+ &
                                b129*f9+ b1210*f10+ b1211*f11),f12)
        call me%f(t+h,     x+h*(b131*f1+ b134*f4+ b135*f5+ b136*f6+ b137*f7+ b138*f8+ &
                                b139*f9+ b1310*f10),f13)

        xf = x + h*(c1*f1 + c6*f6 + c7*f7 + c8*f8 + c9*f9 + c10*f10 + c11*f11 + c12*f12)

        xerr =   h*(e1*f1 + e6*f6 + e7*f7 + e8*f8 + e9*f9 + e10*f10 + e11*f11 + e12*f12 + e13*f13)

    end associate

    end procedure rkv87e
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's "most robust" Runge-Kutta (13:8(7)) pair.
!
!### Reference
!  * J.H. Verner, "Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error",
!    SIAM Journal on Numerical Analysis, Vol. 15, Issue 4 (1978)
!  * [Jim Verner's Refuge for Runge-Kutta Pairs](https://www.sfu.ca/~jverner/)
!  * [Coefficients](https://www.sfu.ca/~jverner/RKV87.IIa.Robust.00000754677.081208.CoeffsOnlyFLOAT)

    module procedure rkv87r

    real(wp),parameter :: a2  = .25_wp
    real(wp),parameter :: a3  = .1128884514435695538057742782152230971129_wp
    real(wp),parameter :: a4  = .1693326771653543307086614173228346456693_wp
    real(wp),parameter :: a5  = .424_wp
    real(wp),parameter :: a6  = .509_wp
    real(wp),parameter :: a7  = .867_wp
    real(wp),parameter :: a8  = .15_wp
    real(wp),parameter :: a9  = .7090680365138684008060140010282474786750_wp
    real(wp),parameter :: a10 = .32_wp
    real(wp),parameter :: a11 = .45_wp

    real(wp),parameter :: b21   =   .25_wp
    real(wp),parameter :: b31   =   .8740084650491523205268632759487741197705e-1_wp
    real(wp),parameter :: b32   =   .2548760493865432175308795062034568513581e-1_wp
    real(wp),parameter :: b41   =   .4233316929133858267716535433070866141732e-1_wp
    real(wp),parameter :: b43   =   .1269995078740157480314960629921259842520_wp
    real(wp),parameter :: b51   =   .4260950588874226149488144523757227409094_wp
    real(wp),parameter :: b53   =  -1.598795284659152326542773323065718111709_wp
    real(wp),parameter :: b54   =   1.596700225771729711593958870689995370799_wp
    real(wp),parameter :: b61   =   .5071933729671392951509061813851363923933e-1_wp
    real(wp),parameter :: b64   =   .2543337726460040758275471440887777803137_wp
    real(wp),parameter :: b65   =   .2039468900572819946573622377727085804470_wp
    real(wp),parameter :: b71   =  -.2900037471752311097038837928542589612409_wp
    real(wp),parameter :: b74   =   1.344187391026078988943868110941433700318_wp
    real(wp),parameter :: b75   =  -2.864777943361442730961110382703656282947_wp
    real(wp),parameter :: b76   =   2.677594299510594851721126064616481543870_wp
    real(wp),parameter :: b81   =   .9853501133799354646974040298072701428476e-1_wp
    real(wp),parameter :: b85   =   .2219268063075138484202403649819738790358_wp
    real(wp),parameter :: b86   =  -.1814062291180699431269033828807395245747_wp
    real(wp),parameter :: b87   =   .1094441147256254823692261491803863125415e-1_wp
    real(wp),parameter :: b91   =   .3871105254573114467944461816516637340565_wp
    real(wp),parameter :: b94   =  -1.442445497485527757125674555307792776717_wp
    real(wp),parameter :: b95   =   2.905398189069950931769134644923384844174_wp
    real(wp),parameter :: b96   =  -1.853771069630105929084333267581197802518_wp
    real(wp),parameter :: b97   =   .1400364809872815426949732510977124147922_wp
    real(wp),parameter :: b98   =   .5727394081149581657574677462444770648875_wp
    real(wp),parameter :: b101  =  -.1612440344443930810063001619791348059544_wp
    real(wp),parameter :: b104  =  -.1733960295735898408357840447396256789490_wp
    real(wp),parameter :: b105  =  -1.301289281406514740601681274517249252974_wp
    real(wp),parameter :: b106  =   1.137950375173861730855879213143100347212_wp
    real(wp),parameter :: b107  =  -.3174764966396688010692352113804302469898e-1_wp
    real(wp),parameter :: b108  =   .9335129382493366643981106448605688485659_wp
    real(wp),parameter :: b109  =  -.8378631833473385270330085562961643320150e-1_wp
    real(wp),parameter :: b111  =  -.1919944488158953328151080465148357607314e-1_wp
    real(wp),parameter :: b114  =   .2733085726526428490794232625401612427562_wp
    real(wp),parameter :: b115  =  -.6753497320694437291969161121094238085624_wp
    real(wp),parameter :: b116  =   .3415184981384601607173848997472838271198_wp
    real(wp),parameter :: b117  =  -.6795006480337577247892051619852462939191e-1_wp
    real(wp),parameter :: b118  =   .9659175224762387888426558649121637650975e-1_wp
    real(wp),parameter :: b119  =   .1325308251118210118072103846654538995123_wp
    real(wp),parameter :: b1110 =   .3685495936038611344690632995153166681295_wp
    real(wp),parameter :: b121  =   .6091877403645289867688841211158881778458_wp
    real(wp),parameter :: b124  =  -2.272569085898001676899980093141308839972_wp
    real(wp),parameter :: b125  =   4.757898342694029006815525588191478549755_wp
    real(wp),parameter :: b126  =  -5.516106706692758482429468966784424824484_wp
    real(wp),parameter :: b127  =   .2900596369680119270909581856594617437818_wp
    real(wp),parameter :: b128  =   .5691423963359036822910985845480184914563_wp
    real(wp),parameter :: b129  =   .7926795760332167027133991620589332757995_wp
    real(wp),parameter :: b1210 =   .1547372045328882289412619077184989823205_wp
    real(wp),parameter :: b1211 =   1.614970895662181624708321510633454443497_wp
    real(wp),parameter :: b131  =   .8873576220853471966321169405198102270488_wp
    real(wp),parameter :: b134  =  -2.975459782108536755851363280470930158198_wp
    real(wp),parameter :: b135  =   5.600717009488163059799039254835009892383_wp
    real(wp),parameter :: b136  =  -5.915607450536674468001493018994165735184_wp
    real(wp),parameter :: b137  =   .2202968915613492701687914254080763833124_wp
    real(wp),parameter :: b138  =   .1015509782446221666614327134090299699755_wp
    real(wp),parameter :: b139  =   1.151434564738605590978039775212585055356_wp
    real(wp),parameter :: b1310 =   1.929710166527123939613436190080584365307_wp

    real(wp),parameter :: c1  =  .4472956466669571420301584042904938246647e-1_wp
    real(wp),parameter :: c6  =  .1569103352770819981336869801072664540918_wp
    real(wp),parameter :: c7  =  .1846097340815163774070245187352627789204_wp
    real(wp),parameter :: c8  =  .2251638060208699104247941940035072197092_wp
    real(wp),parameter :: c9  =  .1479461565197023468700517988544914175374_wp
    real(wp),parameter :: c10 =  .7605554244495582526979836191033649101273e-1_wp
    real(wp),parameter :: c11 =  .1227729023501861961082434631592143738854_wp
    real(wp),parameter :: c12 =  .4181195863899163158338484280087188237679e-1_wp

    real(wp),parameter :: d1  =  .4584711140049592587866473012201028209588e-1_wp
    real(wp),parameter :: d6  =  .2623189140415238743744335658484580339239_wp
    real(wp),parameter :: d7  =  .1916937233785261190448573863568842900803_wp
    real(wp),parameter :: d8  =  .2170917232790261833097840742290644856820_wp
    real(wp),parameter :: d9  =  .1273818962483370679680316945065673786790_wp
    real(wp),parameter :: d10 =  .1151053038536532625824051575004319214889_wp
    real(wp),parameter :: d13 =  .4056132779843756684182339143658360805005e-1_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10 = c10 - d10
    real(wp),parameter :: e11 = c11
    real(wp),parameter :: e12 = c12
    real(wp),parameter :: e13 =     - d13

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13))

        call me%f(t+h,x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b104*f4+b105*f5+b106*f6+b107*f7+&
                            b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b114*f4+b115*f5+b116*f6+b117*f7+&
                            b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+h,    x+h*(b121*f1+b124*f4+b125*f5+b126*f6+b127*f7+&
                            b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+h,    x+h*(b131*f1+b134*f4+b135*f5+b136*f6+b137*f7+&
                            b138*f8+b139*f9+b1310*f10),f13)

        xf = x+h*(c1*f1+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12)

        xerr = h*(e1*f1+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13)

    end associate

    end procedure rkv87r
!*****************************************************************************************

!*****************************************************************************************
!>
!  Kovalnogov, Fedorov, Karpukhina, Simos, Tsitouras 8(7) method
!
!### References
!  * Kovalnogov, Fedorov, Karpukhina, Simos, Tsitouras,
!    [Runge-Kutta Embedded Methods of Orders 8(7) for Use in Quadruple Precision Computations](https://www.researchgate.net/publication/363396601_Runge-Kutta_Embedded_Methods_of_Orders_87_for_Use_in_Quadruple_Precision_Computations),
!    Mathematics, Volume 10, Issue 18, 2022.
!  * [Rational coefficients](http://users.uoa.gr/~tsitourasc/t87.m)
!
!@note The floating point coefficients here were generated from the rational ones from the 2nd reference.

    module procedure rkk87

    real(wp),parameter :: a2 = 2.80032137795311131773988246233287895064681826798949202422973107165103409675642981593e-2_wp
    real(wp),parameter :: a3 = 1.17211542455798199481153380377682961936340701290456762390962847695113281163272045568e-1_wp
    real(wp),parameter :: a4 = 1.75817313683697299221730070566524442904511051935685143586444271542669921744908068353e-1_wp
    real(wp),parameter :: a5 = 4.12418774972264624675725499028219180367525003545123161749372304663713787609586013029e-1_wp
    real(wp),parameter :: a6 = 4.85788040125653655249646480463893555030416915791171976577551272195807115944917385569e-1_wp
    real(wp),parameter :: a7 = 1.61945436881864850373469095233292102203653467396439966236482648418795024659046870770e-1_wp
    real(wp),parameter :: a8 = 9.69999622636762726575177046252154115146102466697694310620259374331752600661643542686e-1_wp
    real(wp),parameter :: a9 = 9.67212187169379825484211792872972399851315982484116576998481473839941803252105313177e-1_wp
    real(wp),parameter :: a10 = 9.51733593009250954565256387258796094610870512261999333726263998154934269533352125669e-1_wp
    real(wp),parameter :: a11 = 9.68202771170154148019728796247765727264175032590143665416817858054805198295905065247e-1_wp

    real(wp),parameter :: b21 = 2.80032137795311131773988246233287895064681826798949202422973107165103409675642981593e-2_wp
    real(wp),parameter :: b31 = -1.2809147513766000602840952300143562083196457175706189219831636475705535103141265693e-1_wp
    real(wp),parameter :: b32 = 2.45303017593458205509562903379118613607705175036823097559603878184775233215504838718e-1_wp
    real(wp),parameter :: b41 = 4.39543284209243248054325176416311107261277629839212858966110678856674804362270170882e-2_wp
    real(wp),parameter :: b43 = 1.31862985262772974416297552924893332178383288951763857689833203657002441308681051264e-1_wp
    real(wp),parameter :: b51 = 3.37793581678770282454466907421349755434169654489904742664883921124336709641723527478e-1_wp
    real(wp),parameter :: b53 = -1.2272549349642134972638116912742760383098006898941228119302659681159830923493358553e0_wp
    real(wp),parameter :: b54 = 1.30188012825770783948507028288114558866194048988029573253542612869530190101485281039e0_wp
    real(wp),parameter :: b61 = 5.55707180742854953774795421040398153698425834244408078243992150711110988336315378908e-2_wp
    real(wp),parameter :: b64 = 2.51201707056335789089668842569664694678061738700922859182154307406956224053004373812e-1_wp
    real(wp),parameter :: b65 = 1.79015614995032370782498095790188986963634606811801882820741709651062551409519978967e-1_wp
    real(wp),parameter :: b71 = 6.68574347263499831237129311835436351646983012984884556347193166264996297617553477916e-2_wp
    real(wp),parameter :: b74 = 1.18464279145447863166879557612335382759089586154167432866512942838128687075481264338e-1_wp
    real(wp),parameter :: b75 = -4.9625776043451759740292422000972633501142622216483112027320727273129396510418976777e-2_wp
    real(wp),parameter :: b76 = 2.62494990535187638231690284383855963023506866196027758824536604264425581164307083575e-2_wp
    real(wp),parameter :: b81 = -2.8828852333227597022830828823050925772450620891223369172510088561729332090076138712e0_wp
    real(wp),parameter :: b84 = -6.2440208358504107644024024514118804244308557646280655968919670956442898657820780678e1_wp
    real(wp),parameter :: b85 = -1.8568693804204652514452785270528995249754499915069250253336433143578163151508120889e0_wp
    real(wp),parameter :: b86 = 3.63680272715359047470913994410370911642713568016445118939856192986207517504148311731e0_wp
    real(wp),parameter :: b87 = 6.45131598677305048496184230256252416917236861835199785338289493414147153759342124993e1_wp
    real(wp),parameter :: b91 = -2.8044177030062210609793681802747767950835025205736626228179150569096325843989251959e0_wp
    real(wp),parameter :: b94 = -6.0724923514423992066404703993882607779493861739697985128789442636875454172217394054e1_wp
    real(wp),parameter :: b95 = -1.8439925449308251307668788876911848442475862438554757484077098049668682602503860260e0_wp
    real(wp),parameter :: b96 = 3.58341131683436481398169260904220163311107056209428641899333714250541234007157548631e0_wp
    real(wp),parameter :: b97 = 6.27562314336350305301529075880419713175762252340352806429712952136086493920793363741e1_wp
    real(wp),parameter :: b98 = 9.03199061022739500562657637369442031043920598379602859131908802359360137694011980530e-4_wp
    real(wp),parameter :: b101 = -2.3867982748517262502656630171168908345533985055664349630981323816949508696589885131e0_wp
    real(wp),parameter :: b104 = -5.1627063109123959312754417147040225909854829131526049890819790306938208226330016825e1_wp
    real(wp),parameter :: b105 = -1.7536886231138950479955440820659955467185829948904221201507369859127252396249986529e0_wp
    real(wp),parameter :: b106 = 3.28049564858832505325605350268523413528007142886883523265902376116967563208653997102e0_wp
    real(wp),parameter :: b107 = 5.34320383180959473568498403628072891122624268061087757470806360734926014535143237753e1_wp
    real(wp),parameter :: b108 = -5.1871871749054814179692202838627648419279746926926280632135852213870053073114332763e-2_wp
    real(wp),parameter :: b109 = 5.86215051636139696546789708280127234169485458303572483088290904612000509294669736233e-2_wp
    real(wp),parameter :: b111 = -2.8321844833678190185701414338307927128819672795892717761836464805524068227918002961e0_wp
    real(wp),parameter :: b114 = -6.1331782300299804156579640111028383230359362036389968935124318342451005168121921975e1_wp
    real(wp),parameter :: b115 = -1.8486956144396106169929005145135397640092170913037834243348598217810735996138171447e0_wp
    real(wp),parameter :: b116 = 3.60245017376559155660368881029952545634233977218873589569763598469965203717427900757e0_wp
    real(wp),parameter :: b117 = 6.33778491267616171592057963772550300254984410272010439289876706968447695603368885389e1_wp
    real(wp),parameter :: b118 = -1.4415990912159637680387758372857367088383577073617671128542010599339480760810185327e-2_wp
    real(wp),parameter :: b119 = 1.79850128395124994545750132486994455015114960352855140561364324582618209513891071121e-2_wp
    real(wp),parameter :: b1110 = -3.0031531771736374212615868099161882544231529622685172138561311528231108610020573246e-3_wp
    real(wp),parameter :: b121 = -4.0538618880188790999955490264156166727875270359794492495465581102083340028760864442e0_wp
    real(wp),parameter :: b124 = -8.9573424620688719282465667174595500838253854698064448723443082397797973834679560630e1_wp
    real(wp),parameter :: b125 = -2.1741602779971559275874048256900391711850691989074830707636990301813342017567123130e0_wp
    real(wp),parameter :: b126 = 4.90277984566299792538628476237120023643679128602934227568174295156160923256501166825e0_wp
    real(wp),parameter :: b127 = 9.21308811959261498386202635372829188598935945642923099199446915685906551410074148629e1_wp
    real(wp),parameter :: b128 = -1.0615798848368367252894309217989593007445778391108496729169226136688468696047739072e2_wp
    real(wp),parameter :: b129 = 7.15244888413933960825621793700326950628446023762240029384782852503289700766263846501e2_wp
    real(wp),parameter :: b1210 = -5.8741641553458552194866804448440052315485267257872985945866582180724510668125792910e1_wp
    real(wp),parameter :: b1211 = -5.5057747263167612955576982434494393109010422446397491970269616705184666943368933333e2_wp
    real(wp),parameter :: b131 = -1.5426576117927492437318562966731433591238073288653428133646285261236434041035073113e1_wp
    real(wp),parameter :: b134 = -4.0842463342275502766697135673815027343155742825350936594898511656878424252559774704e2_wp
    real(wp),parameter :: b135 = -1.3703330210475051981990819170722577088815078067066750833460437452985770323423290388e1_wp
    real(wp),parameter :: b136 = 3.86861748380283839981840260187136546795985692447349074847554981151691738530542460349e1_wp
    real(wp),parameter :: b137 = 4.11064481162365098634151552138274804840423912301682160479008282157089539192861116480e2_wp
    real(wp),parameter :: b138 = -5.2625896422984832125592493280526052943397760170789326525937264023403370769397788590e3_wp
    real(wp),parameter :: b139 = 3.59120401468351311263204092337603786279414890192167416431431463233495380766373955611e4_wp
    real(wp),parameter :: b1310 = -2.9428565975170828789337889297586264433907248673913937095047691721267811559650558183e3_wp
    real(wp),parameter :: b1311 = -2.7717790023268800945373425815230530802580454445165684044343605511598675223912861429e4_wp

    real(wp),parameter :: c1 = 4.62710144559749507643438665696565266477463923324303941594264880671698479184226674443e-2_wp ! T87bvec
    real(wp),parameter :: c6 = 3.70024261130579576116287766947369066170603597694346848542347561458083670963945787412e-1_wp
    real(wp),parameter :: c7 = 2.59046906438325331026392759493100306344056209570957698029069680004143772202534147226e-1_wp
    real(wp),parameter :: c8 = 4.39526151941350658964709597497950710839169490812990167038362808341439048915871691690e2_wp
    real(wp),parameter :: c9 = 2.54699783254442554735577328130608901983844079828777855737604201012040528247471881063e-3_wp
    real(wp),parameter :: c10 = 3.71160265008765610052766362243075634874924739787417318444668226976710753590848050226e1_wp
    real(wp),parameter :: c11 = -4.6957326399687634573980555517438556136604535869046208947272277883667432727498150998e2_wp
    real(wp),parameter :: c12 = -6.8840856231271330609815727620309531232087818518011356960217496959521261426131545711e0_wp
    real(wp),parameter :: c13 = 1.37281997918834547346514047866805411030176899063475546305931321540062434963579604579e-1_wp

    real(wp),parameter :: e1 = 5.85295842711462028451365622251323842688965073989360183197454435266584210000889266151e-6_wp ! T87evec
    real(wp),parameter :: e6 = 9.39126876188594677381173972498082533162562941166681876536458073485843968573614935432e-5_wp
    real(wp),parameter :: e7 = -2.1858808218216614235946904150275890682710865952635557843966398644681448317233753310e-5_wp
    real(wp),parameter :: e8 = 3.82147409595844443746138596500996871329379548062006046484595468670810794523045735063e2_wp
    real(wp),parameter :: e9 = 7.57594006256575145190136681251802140063882902269863054941867158087541929651242123261e2_wp
    real(wp),parameter :: e10 = -5.0572248451354316398326863204362055571090849938158216870898444305578120258688628657e0_wp
    real(wp),parameter :: e11 = -1.1340723721698801541760405363985770050927554242983877738142540298190990694395498870e3_wp
    real(wp),parameter :: e12 = -6.2562494403371433261049312272180260925647282293121510390454736461360125268267236809e-1_wp
    real(wp),parameter :: e13 = 1.37281997918834547346514047866805411030176899063475546305931321540062434963579604579e-2_wp

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b104*f4+b105*f5+b106*f6+b107*f7+&
                            b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b114*f4+b115*f5+b116*f6+b117*f7+&
                            b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+h,    x+h*(b121*f1+b124*f4+b125*f5+b126*f6+b127*f7+&
                            b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+h,    x+h*(b131*f1+b134*f4+b135*f5+b136*f6+b137*f7+&
                            b138*f8+b139*f9+b1310*f10+b1311*f11),f13)

        xf = x+h*(c1*f1+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12+c13*f13)

        xerr = h*(e1*f1+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13)

    end associate

    end procedure rkk87
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's "most efficient" Runge-Kutta (16:9(8)) pair.
!
!### Reference
!  * J.H. Verner, "Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error",
!    SIAM Journal on Numerical Analysis, Vol. 15, Issue 4 (1978)
!  * [Jim Verner's Refuge for Runge-Kutta Pairs](https://www.sfu.ca/~jverner/)
!  * [Coefficients](https://www.sfu.ca/~jverner/RKV98.IIa.Efficient.000000349.081209.FLOAT6040OnWeb)

    module procedure rkv98e

    real(wp),parameter :: a2    = .3462e-1_wp
    real(wp),parameter :: a3    = .9702435063878044594828361677100617517633e-1_wp
    real(wp),parameter :: a4    = .1455365259581706689224254251565092627645_wp
    real(wp),parameter :: a5    = .561_wp
    real(wp),parameter :: a6    = .2290079115904850126662751771814700052182_wp
    real(wp),parameter :: a7    = .5449920884095149873337248228185299947818_wp
    real(wp),parameter :: a8    = .645_wp
    real(wp),parameter :: a9    = .4837500000000000000000000000000000000000_wp
    real(wp),parameter :: a10   = .6757e-1_wp
    real(wp),parameter :: a11   = .2500_wp
    real(wp),parameter :: a12   = .6590650618730998549405331618649220295334_wp
    real(wp),parameter :: a13   = .8206_wp
    real(wp),parameter :: a14   = .9012_wp

    real(wp),parameter :: b21   =  .3462e-1_wp
    real(wp),parameter :: b31   = -.389335438857287327017042687229284478532e-1_wp
    real(wp),parameter :: b32   =  .1359578945245091786499878854939346230295_wp
    real(wp),parameter :: b41   =  .3638413148954266723060635628912731569111e-1_wp
    real(wp),parameter :: b43   =  .1091523944686280016918190688673819470733_wp
    real(wp),parameter :: b51   =  2.025763914393969636805657604282571047511_wp
    real(wp),parameter :: b53   = -7.638023836496292020387602153091964592952_wp
    real(wp),parameter :: b54   =  6.173259922102322383581944548809393545442_wp
    real(wp),parameter :: b61   =  .5112275589406060872792270881648288397197e-1_wp
    real(wp),parameter :: b64   =  .1770823794555021537929910813839068684087_wp
    real(wp),parameter :: b65   =  .80277624092225014536138698108025283759e-3_wp
    real(wp),parameter :: b71   =  .1316006357975216279279871693164256985334_wp
    real(wp),parameter :: b74   = -.2957276252669636417685183174672273730699_wp
    real(wp),parameter :: b75   =  .878137803564295237421124704053886667082e-1_wp
    real(wp),parameter :: b76   =  .6213052975225274774321435005639430026100_wp
    real(wp),parameter :: b81   =  .7166666666666666666666666666666666666667e-1_wp
    real(wp),parameter :: b86   =  .3305533578915319409260346730051472207728_wp
    real(wp),parameter :: b87   =  .2427799754418013924072986603281861125606_wp
    real(wp),parameter :: b91   =  .7180664062500000000000000000000000000000e-1_wp
    real(wp),parameter :: b96   =  .3294380283228177160744825466257672816401_wp
    real(wp),parameter :: b97   =  .1165190029271822839255174533742327183599_wp
    real(wp),parameter :: b98   = -.3401367187500000000000000000000000000000e-1_wp
    real(wp),parameter :: b101  =  .4836757646340646986611287718844085773549e-1_wp
    real(wp),parameter :: b106  =  .3928989925676163974333190042057047002852e-1_wp
    real(wp),parameter :: b107  =  .1054740945890344608263649267140088017604_wp
    real(wp),parameter :: b108  = -.2143865284648312665982642293830533996214e-1_wp
    real(wp),parameter :: b109  = -.1041229174627194437759832813847147895623_wp
    real(wp),parameter :: b111  = -.2664561487201478635337289243849737340534e-1_wp
    real(wp),parameter :: b116  =  .3333333333333333333333333333333333333333e-1_wp
    real(wp),parameter :: b117  = -.1631072244872467239162704487554706387141_wp
    real(wp),parameter :: b118  =  .3396081684127761199487954930015522928244e-1_wp
    real(wp),parameter :: b119  =  .1572319413814626097110769806810024118077_wp
    real(wp),parameter :: b1110 =  .2152267478031879552303534778794770376960_wp
    real(wp),parameter :: b121  =  .3689009248708622334786359863227633989718e-1_wp
    real(wp),parameter :: b126  = -.1465181576725542928653609891758501156785_wp
    real(wp),parameter :: b127  =  .2242577768172024345345469822625833796001_wp
    real(wp),parameter :: b128  =  .2294405717066072637090897902753790803034e-1_wp
    real(wp),parameter :: b129  = -.35850052905728761357394424889330334334e-2_wp
    real(wp),parameter :: b1210 =  .8669223316444385506869203619044453906053e-1_wp
    real(wp),parameter :: b1211 =  .4383840651968337846196219974168630120572_wp
    real(wp),parameter :: b131  = -.4866012215113340846662212357570395295088_wp
    real(wp),parameter :: b136  = -6.304602650282852990657772792012007122988_wp
    real(wp),parameter :: b137  = -.281245618289472564778284183790118418111_wp
    real(wp),parameter :: b138  = -2.679019236219849057687906597489223155566_wp
    real(wp),parameter :: b139  =  .518815663924157511565311164615012522024_wp
    real(wp),parameter :: b1310 =  1.365353187603341710683633635235238678626_wp
    real(wp),parameter :: b1311 =  5.885091088503946585721274891680604830712_wp
    real(wp),parameter :: b1312 =  2.802808786272062889819965117517532194812_wp
    real(wp),parameter :: b141  =  .4185367457753471441471025246471931649633_wp
    real(wp),parameter :: b146  =  6.724547581906459363100870806514855026676_wp
    real(wp),parameter :: b147  = -.425444280164611790606983409697113064616_wp
    real(wp),parameter :: b148  =  3.343279153001265577811816947557982637749_wp
    real(wp),parameter :: b149  =  .617081663117537759528421117507709784737_wp
    real(wp),parameter :: b1410 = -.929966123939932833937749523988800852013_wp
    real(wp),parameter :: b1411 = -6.099948804751010722472962837945508844846_wp
    real(wp),parameter :: b1412 = -3.002206187889399044804158084895173690015_wp
    real(wp),parameter :: b1413 =  .2553202529443445472336424602988558373637_wp
    real(wp),parameter :: b151  = -.779374086122884664644623040843840506343_wp
    real(wp),parameter :: b156  = -13.93734253810777678786523664804936051203_wp
    real(wp),parameter :: b157  =  1.252048853379357320949735183924200895136_wp
    real(wp),parameter :: b158  = -14.69150040801686878191527989293072091588_wp
    real(wp),parameter :: b159  = -.494705058533141685655191992136962873577_wp
    real(wp),parameter :: b1510 =  2.242974909146236657906984549543692874755_wp
    real(wp),parameter :: b1511 =  13.36789380382864375813864978592679139881_wp
    real(wp),parameter :: b1512 =  14.39665048665068644512236935340272139005_wp
    real(wp),parameter :: b1513 = -.7975813331776800379127866056663258667437_wp
    real(wp),parameter :: b1514 =  .4409353709534277758753793068298041158235_wp
    real(wp),parameter :: b161  =  2.058051337466886442151242368989994043993_wp
    real(wp),parameter :: b166  =  22.35793772796803295519317565842520212899_wp
    real(wp),parameter :: b167  =  .90949810997556332745009198137971890783_wp
    real(wp),parameter :: b168  =  35.89110098240264104710550686568482456493_wp
    real(wp),parameter :: b169  = -3.442515027624453437985000403608480262211_wp
    real(wp),parameter :: b1610 = -4.865481358036368826566013387928704014496_wp
    real(wp),parameter :: b1611 = -18.90980381354342625688427480879773032857_wp
    real(wp),parameter :: b1612 = -34.26354448030451782929251177395134170515_wp
    real(wp),parameter :: b1613 =  1.264756521695642578827783499806516664686_wp

    real(wp),parameter :: c1    =  .1461197685842315252051541915018784713459e-1_wp
    real(wp),parameter :: c8    = -.3915211862331339089410228267288242030810_wp
    real(wp),parameter :: c9    =  .2310932500289506415909675644868993669908_wp
    real(wp),parameter :: c10   =  .1274766769992852382560589467488989175618_wp
    real(wp),parameter :: c11   =  .2246434176204157731566981937082069688984_wp
    real(wp),parameter :: c12   =  .5684352689748512932705226972873692126743_wp
    real(wp),parameter :: c13   =  .5825871557215827200814768021863420902155e-1_wp
    real(wp),parameter :: c14   =  .1364317403482215641609022744494239843327_wp
    real(wp),parameter :: c15   =  .3057013983082797397721005067920369646664e-1_wp

    real(wp),parameter :: d1    =  .1996996514886773085518508418098868756464e-1_wp
    real(wp),parameter :: d8    =  2.191499304949330054530747099310837524864_wp
    real(wp),parameter :: d9    =  .8857071848208438030833722031786358862953e-1_wp
    real(wp),parameter :: d10   =  .1140560234865965622484956605091432032674_wp
    real(wp),parameter :: d11   =  .2533163805345107065564577734569651977347_wp
    real(wp),parameter :: d12   = -2.056564386240941011158999594595981300493_wp
    real(wp),parameter :: d13   =  .3408096799013119935160094894224543812830_wp
    real(wp),parameter :: d16   =  .4834231373823958314376726739772871714902e-1_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10 = c10 - d10
    real(wp),parameter :: e11 = c11 - d11
    real(wp),parameter :: e12 = c12 - d12
    real(wp),parameter :: e13 = c13 - d13
    real(wp),parameter :: e14 = c14
    real(wp),parameter :: e15 = c15
    real(wp),parameter :: e16 =     - d16

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13), &
               f14 => me%funcs(:,14), &
               f15 => me%funcs(:,15), &
               f16 => me%funcs(:,16))

        call me%f(t,x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b106*f6+b107*f7+b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b116*f6+b117*f7+b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b126*f6+b127*f7+b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b131*f1+b136*f6+b137*f7+b138*f8+b139*f9+b1310*f10+b1311*f11+&
                               b1312*f12),f13)
        call me%f(t+a14*h,x+h*(b141*f1+b146*f6+b147*f7+b148*f8+b149*f9+b1410*f10+b1411*f11+&
                               b1412*f12+b1413*f13),f14)
        call me%f(t+h,    x+h*(b151*f1+b156*f6+b157*f7+b158*f8+b159*f9+b1510*f10+b1511*f11+&
                               b1512*f12+b1513*f13+b1514*f14),f15)
        call me%f(t+h,    x+h*(b161*f1+b166*f6+b167*f7+b168*f8+b169*f9+b1610*f10+b1611*f11+&
                               b1612*f12+b1613*f13),f16)

        xf = x + h*(c1*f1 + c8*f8 + c9*f9 + c10*f10 + c11*f11 + c12*f12 + c13*f13 + c14*f14 + c15*f15 )

        xerr =   h*(e1*f1 + e8*f8 + e9*f9 + e10*f10 + e11*f11 + e12*f12 + e13*f13 + e14*f14 + e15*f15 + e16*f16)

    end associate

    end procedure rkv98e
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's "most robust" Runge-Kutta (16:9(8)) pair.
!
!### Reference
!  * J.H. Verner, "Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error",
!    SIAM Journal on Numerical Analysis, Vol. 15, Issue 4 (1978)
!  * [Jim Verner's Refuge for Runge-Kutta Pairs](https://www.sfu.ca/~jverner/)
!  * [Coefficients](https://www.sfu.ca/~jverner/RKV98.IIa.Robust.000000351.081209.CoeffsOnlyFLOAT6040)

    module procedure rkv98r

    real(wp),parameter :: a2  =  .4e-1_wp
    real(wp),parameter :: a3  =  .9648736013787361245235039379666356743708e-1_wp
    real(wp),parameter :: a4  =  .1447310402068104186785255906949953511556_wp
    real(wp),parameter :: a5  =  .576_wp
    real(wp),parameter :: a6  =  .2272326564618766017153738192188229509142_wp
    real(wp),parameter :: a7  =  .5407673435381233982846261807811770490858_wp
    real(wp),parameter :: a8  =  .64_wp
    real(wp),parameter :: a9  =  .48_wp
    real(wp),parameter :: a10 =  .6754e-1_wp
    real(wp),parameter :: a11 =  .25_wp
    real(wp),parameter :: a12 =  .6770920153543242682384311058159603931192_wp
    real(wp),parameter :: a13 =  .8115_wp
    real(wp),parameter :: a14 =  .906_wp

    real(wp),parameter :: b21 =  .4e-1_wp
    real(wp),parameter :: b31 = -.198852731918229097650241511466089129345e-1_wp
    real(wp),parameter :: b32 =  .1163726333296965222173745449432724803716_wp
    real(wp),parameter :: b41 =  .3618276005170260466963139767374883778890e-1_wp
    real(wp),parameter :: b43 =  .1085482801551078140088941930212465133667_wp
    real(wp),parameter :: b51 =  2.272114264290177409193144938921415409241_wp
    real(wp),parameter :: b53 = -8.526886447976398578316416192982602292786_wp
    real(wp),parameter :: b54 =  6.830772183686221169123271254061186883545_wp
    real(wp),parameter :: b61 =  .5094385535389374394512668566783434123978e-1_wp
    real(wp),parameter :: b64 =  .1755865049809071110203693328749561646990_wp
    real(wp),parameter :: b65 =  .70229612707574674987780067603244497535e-3_wp
    real(wp),parameter :: b71 =  .1424783668683284782770955365543878809824_wp
    real(wp),parameter :: b74 = -.3541799434668684104094753917518523845155_wp
    real(wp),parameter :: b75 =  .7595315450295100889001534202778550159932e-1_wp
    real(wp),parameter :: b76 =  .6765157656337123215269906939508560510196_wp
    real(wp),parameter :: b81 =  .7111111111111111111111111111111111111111e-1_wp
    real(wp),parameter :: b86 =  .3279909287605898328568406057725491803016_wp
    real(wp),parameter :: b87 =  .2408979601282990560320482831163397085872_wp
    real(wp),parameter :: b91 =  .7125e-1_wp
    real(wp),parameter :: b96 =  .3268842451575245554847578757216915662785_wp
    real(wp),parameter :: b97 =  .1156157548424754445152421242783084337215_wp
    real(wp),parameter :: b98 = -.3375e-1_wp
    real(wp),parameter :: b101 =  .4822677322465810178387112087673611111111e-1_wp
    real(wp),parameter :: b106 =  .3948559980495400110769549704186108167677e-1_wp
    real(wp),parameter :: b107 =  .1058851161934658144373823566907778072121_wp
    real(wp),parameter :: b108 = -.2152006320474309346664428710937500000000e-1_wp
    real(wp),parameter :: b109 = -.1045374260183348238623046875000000000000_wp
    real(wp),parameter :: b111 = -.2609113435754923412210928689962011065179e-1_wp
    real(wp),parameter :: b116 =  .3333333333333333333333333333333333333333e-1_wp
    real(wp),parameter :: b117 = -.1652504006638105086724681598195267241410_wp
    real(wp),parameter :: b118 =  .3434664118368616658319419895678838776647e-1_wp
    real(wp),parameter :: b119 =  .1595758283215209043195814910843067811951_wp
    real(wp),parameter :: b1110 =  .2140857321828193385584684233447183324979_wp
    real(wp),parameter :: b121 = -.362842339625565859076509979091267105528e-1_wp
    real(wp),parameter :: b126 = -1.096167597427208807028761474420297770752_wp
    real(wp),parameter :: b127 =  .1826035504321331052308236240517254331348_wp
    real(wp),parameter :: b128 =  .708225444417068325613028685455625123741e-1_wp
    real(wp),parameter :: b129 = -.231364701848243126999929738482630407146e-1_wp
    real(wp),parameter :: b1210 =  .2711204726320932916455631550463654973432_wp
    real(wp),parameter :: b1211 =  1.308133749422980744437146904349994472286_wp
    real(wp),parameter :: b131 = -.5074635056416974879347823927726392374259_wp
    real(wp),parameter :: b136 = -6.631342198657237090355284142048733580937_wp
    real(wp),parameter :: b137 = -.252748010090880105270020973014860316405_wp
    real(wp),parameter :: b138 = -.4952612380036095562991116175550167835424_wp
    real(wp),parameter :: b139 =  .293252554525388690285739720360003594753_wp
    real(wp),parameter :: b1310 =  1.440108693768280908474851998204423941413_wp
    real(wp),parameter :: b1311 =  6.237934498647055877243623886838802127716_wp
    real(wp),parameter :: b1312 =  .7270192054526987638549835199880202544289_wp
    real(wp),parameter :: b141 =  .6130118256955931701496387847232542148725_wp
    real(wp),parameter :: b146 =  9.088803891640463313341034206647776279557_wp
    real(wp),parameter :: b147 = -.407378815629344868103315381138325162923_wp
    real(wp),parameter :: b148 =  1.790733389490374687043894756399015035977_wp
    real(wp),parameter :: b149 =  .714927166761755073724875250629602731782_wp
    real(wp),parameter :: b1410 = -1.438580857841722850237810322456327208949_wp
    real(wp),parameter :: b1411 = -8.263329312064740580595954649844133476994_wp
    real(wp),parameter :: b1412 = -1.537570570808865115231450725068826856201_wp
    real(wp),parameter :: b1413 =  .3453832827564871699090880801079644428793_wp
    real(wp),parameter :: b151 = -1.211697910343873872490625222495537087293_wp
    real(wp),parameter :: b156 = -19.05581871559595277753334676575234493500_wp
    real(wp),parameter :: b157 =  1.26306067538987510135943101851905310045_wp
    real(wp),parameter :: b158 = -6.913916969178458046793476128409110926069_wp
    real(wp),parameter :: b159 = -.676462266509498065300115641383621209887_wp
    real(wp),parameter :: b1510 =  3.367860445026607887090352785684064242560_wp
    real(wp),parameter :: b1511 =  18.00675164312590810020103216906571965203_wp
    real(wp),parameter :: b1512 =  6.838828926794279896350389904990814350968_wp
    real(wp),parameter :: b1513 = -1.031516451921950498420447675652291096155_wp
    real(wp),parameter :: b1514 =  .4129106232130622755368055554332539084021_wp
    real(wp),parameter :: b161 =  2.157389007494053627033175177985666660692_wp
    real(wp),parameter :: b166 =  23.80712219809580523172312179815279712750_wp
    real(wp),parameter :: b167 =  .88627792492165554903036801415266308369_wp
    real(wp),parameter :: b168 =  13.13913039759876381480201677314222971522_wp
    real(wp),parameter :: b169 = -2.604415709287714883747369630937415176632_wp
    real(wp),parameter :: b1610 = -5.193859949783872300189266203049579105962_wp
    real(wp),parameter :: b1611 = -20.41234071154150778768154893536134356354_wp
    real(wp),parameter :: b1612 = -12.30085625250572261314889445241581039623_wp
    real(wp),parameter :: b1613 =  1.521553095008539362178397458330791655267_wp

    real(wp),parameter :: c1 =  .1458885278405539719101539582255752917034e-1_wp
    real(wp),parameter :: c8 =  .2024197887889332650566666683195656097825e-2_wp
    real(wp),parameter :: c9 =  .2178047084569716646796256135839225745895_wp
    real(wp),parameter :: c10 =  .1274895340854389692868677968654808668201_wp
    real(wp),parameter :: c11 =  .2244617745463131861258531547137348031621_wp
    real(wp),parameter :: c12 =  .1787254491259903095100090833796054447157_wp
    real(wp),parameter :: c13 =  .7594344758096557172908303416513173076283e-1_wp
    real(wp),parameter :: c14 =  .1294845879197561516869001434704642286297_wp
    real(wp),parameter :: c15 =  .2947744761261941714007911131590716605202e-1_wp

    real(wp),parameter :: d1 =  .2034666655224434599707885098832906986649e-1_wp
    real(wp),parameter :: d8 =  1.069617650982700109541321983413338230042_wp
    real(wp),parameter :: d9 =  .7680834711303187278673130261850350530338e-1_wp
    real(wp),parameter :: d10 =  .1130778186885240437498706751119241126785_wp
    real(wp),parameter :: d11 =  .2552587357981962194892445789565762186511_wp
    real(wp),parameter :: d12 = -.9825898086919164036191607912120918904022_wp
    real(wp),parameter :: d13 =  .3981545824421514217762002137442675068982_wp
    real(wp),parameter :: d16 =  .4932600711506839027871318637915324696208e-1_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10  = c10  - d10
    real(wp),parameter :: e11  = c11  - d11
    real(wp),parameter :: e12  = c12  - d12
    real(wp),parameter :: e13  = c13  - d13
    real(wp),parameter :: e14  = c14
    real(wp),parameter :: e15  = c15
    real(wp),parameter :: e16  =      - d16

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13), &
               f14 => me%funcs(:,14), &
               f15 => me%funcs(:,15), &
               f16 => me%funcs(:,16))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b106*f6+b107*f7+b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b116*f6+b117*f7+b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b126*f6+b127*f7+b128*f8+b129*f9+b1210*f10+&
                               b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b131*f1+b136*f6+b137*f7+b138*f8+b139*f9+b1310*f10+&
                               b1311*f11+b1312*f12),f13)
        call me%f(t+a14*h,x+h*(b141*f1+b146*f6+b147*f7+b148*f8+b149*f9+b1410*f10+&
                               b1411*f11+b1412*f12+b1413*f13),f14)
        call me%f(t+h,    x+h*(b151*f1+b156*f6+b157*f7+b158*f8+b159*f9+b1510*f10+&
                               b1511*f11+b1512*f12+b1513*f13+b1514*f14),f15)
        call me%f(t+h,    x+h*(b161*f1+b166*f6+b167*f7+b168*f8+b169*f9+b1610*f10+&
                               b1611*f11+b1612*f12+b1613*f13),f16)

        xf = x+h*(c1*f1+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12+c13*f13+c14*f14+c15*f15)

        xerr = h*(e1*f1+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13+e14*f14+e15*f15+e16*f16)

    end associate

    end procedure rkv98r
!*****************************************************************************************

!*****************************************************************************************
!>
!  Tsitouras & Papakostas NEW6(4) Runge-Kutta method.
!
!### Reference
!  * C. Tsitouras and S. N. Papakostas, "Cheap Error Estimation for Runge-Kutta
!    methods", SIAM J. Sci. Comput. 20(1999) 2067-2088.
!  * [Rational coefficients](http://users.uoa.gr/~tsitourasc/publications.html)
!    (see [rktp64.m](http://users.uoa.gr/~tsitourasc/rktp64.m))

    module procedure rktp64

    real(wp),parameter :: b21 = 4.0_wp / 27.0_wp
    real(wp),parameter :: b31 = 1.0_wp / 18.0_wp
    real(wp),parameter :: b32 = 1.0_wp / 6.0_wp
    real(wp),parameter :: b41 = 66.0_wp / 343.0_wp
    real(wp),parameter :: b42 = -729.0_wp / 1372.0_wp
    real(wp),parameter :: b43 = 1053.0_wp / 1372.0_wp
    real(wp),parameter :: b51 = 13339.0_wp / 49152.0_wp
    real(wp),parameter :: b52 = -4617.0_wp / 16384.0_wp
    real(wp),parameter :: b53 = 5427.0_wp / 53248.0_wp
    real(wp),parameter :: b54 = 95207.0_wp / 159744.0_wp
    real(wp),parameter :: b61 = -6935.0_wp / 57122.0_wp
    real(wp),parameter :: b62 = 23085.0_wp / 48334.0_wp
    real(wp),parameter :: b63 = 33363360.0_wp / 273642941.0_wp
    real(wp),parameter :: b64 = 972160.0_wp / 118442467.0_wp
    real(wp),parameter :: b65 = 172687360.0_wp / 610434253.0_wp
    real(wp),parameter :: b71 = 611.0_wp / 1891.0_wp
    real(wp),parameter :: b72 = -4617.0_wp / 7564.0_wp
    real(wp),parameter :: b73 = 6041007.0_wp / 13176488.0_wp
    real(wp),parameter :: b74 = 12708836.0_wp / 22100117.0_wp
    real(wp),parameter :: b75 = -35840000.0_wp / 62461621.0_wp
    real(wp),parameter :: b76 = 6597591.0_wp / 7972456.0_wp
    real(wp),parameter :: a2  = 4.0_wp / 27.0_wp
    real(wp),parameter :: a3  = 2.0_wp / 9.0_wp
    real(wp),parameter :: a4  = 3.0_wp / 7.0_wp
    real(wp),parameter :: a5  = 11.0_wp / 16.0_wp
    real(wp),parameter :: a6  = 10.0_wp / 13.0_wp
    real(wp),parameter :: c1  = 131.0_wp / 1800.0_wp           ! 6th order formula
    real(wp),parameter :: c3  = 1121931.0_wp / 3902080.0_wp
    real(wp),parameter :: c4  = 319333.0_wp / 1682928.0_wp
    real(wp),parameter :: c5  = 262144.0_wp / 2477325.0_wp
    real(wp),parameter :: c6  = 4084223.0_wp / 15177600.0_wp
    real(wp),parameter :: c7  = 1891.0_wp / 25200.0_wp
    real(wp),parameter :: d1  = 2694253.0_wp / 26100360.0_wp    ! 4th order formula
    real(wp),parameter :: d3  = 83647323.0_wp / 535804360.0_wp
    real(wp),parameter :: d4  = 691202281.0_wp / 1789061040.0_wp
    real(wp),parameter :: d5  = -1275547648.0_wp / 10565208225.0_wp
    real(wp),parameter :: d6  = 2.0_wp / 5.0_wp
    real(wp),parameter :: d7  = 1891.0_wp / 25200.0_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7))

        call me%f(t,      x,f1)
        call me%f(t+a2*h, x+h*(b21*f1),f2)
        call me%f(t+a3*h, x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h, x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h, x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+a6*h, x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+h,    x+h*(b71*f1+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)

        xf = x + h*(c1*f1+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7)

        xerr =   h*(e1*f1+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7)

    end associate

    end procedure rktp64
!*****************************************************************************************

!*****************************************************************************************
!>
!  Fehlberg's 7(8) algorithm.
!
!### Reference
!  * E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order
!    Runge-Kutta Formulas with Stepsize Control",
!    [NASA TR R-2870](https://ntrs.nasa.gov/citations/19680027281).

    module procedure rkf78

    real(wp),parameter :: a1  = 2.0_wp/27.0_wp
    real(wp),parameter :: a2  = 1.0_wp/9.0_wp
    real(wp),parameter :: a3  = 1.0_wp/6.0_wp
    real(wp),parameter :: a4  = 5.0_wp/12.0_wp
    real(wp),parameter :: a5  = 1.0_wp/2.0_wp
    real(wp),parameter :: a6  = 5.0_wp/6.0_wp
    real(wp),parameter :: a7  = 1.0_wp/6.0_wp
    real(wp),parameter :: a8  = 2.0_wp/3.0_wp
    real(wp),parameter :: a9  = 1.0_wp/3.0_wp

    real(wp),parameter :: b10  = 2.0_wp/27.0_wp
    real(wp),parameter :: b20  = 1.0_wp/36.0_wp
    real(wp),parameter :: b21  = 1.0_wp/12.0_wp
    real(wp),parameter :: b30  = 1.0_wp/24.0_wp
    real(wp),parameter :: b32  = 1.0_wp/8.0_wp
    real(wp),parameter :: b40  = 5.0_wp/12.0_wp
    real(wp),parameter :: b42  = -25.0_wp/16.0_wp
    real(wp),parameter :: b43  = 25.0_wp/16.0_wp
    real(wp),parameter :: b50  = 1.0_wp/20.0_wp
    real(wp),parameter :: b53  = 1.0_wp/4.0_wp
    real(wp),parameter :: b54  = 1.0_wp/5.0_wp
    real(wp),parameter :: b60  = -25.0_wp/108.0_wp
    real(wp),parameter :: b63  = 125.0_wp/108.0_wp
    real(wp),parameter :: b64  = -65.0_wp/27.0_wp
    real(wp),parameter :: b65  = 125.0_wp/54.0_wp
    real(wp),parameter :: b70  = 31.0_wp/300.0_wp
    real(wp),parameter :: b74  = 61.0_wp/225.0_wp
    real(wp),parameter :: b75  = -2.0_wp/9.0_wp
    real(wp),parameter :: b76  = 13.0_wp/900.0_wp
    real(wp),parameter :: b80  = 2.0_wp
    real(wp),parameter :: b83  = -53.0_wp/6.0_wp
    real(wp),parameter :: b84  = 704.0_wp/45.0_wp
    real(wp),parameter :: b85  = -107.0_wp/9.0_wp
    real(wp),parameter :: b86  = 67.0_wp/90.0_wp
    real(wp),parameter :: b87  = 3.0_wp
    real(wp),parameter :: b90  = -91.0_wp/108.0_wp
    real(wp),parameter :: b93  = 23.0_wp/108.0_wp
    real(wp),parameter :: b94  = -976.0_wp/135.0_wp
    real(wp),parameter :: b95  = 311.0_wp/54.0_wp
    real(wp),parameter :: b96  = -19.0_wp/60.0_wp
    real(wp),parameter :: b97  = 17.0_wp/6.0_wp
    real(wp),parameter :: b98  = -1.0_wp/12.0_wp
    real(wp),parameter :: b100 = 2383.0_wp/4100.0_wp
    real(wp),parameter :: b103 = -341.0_wp/164.0_wp
    real(wp),parameter :: b104 = 4496.0_wp/1025.0_wp
    real(wp),parameter :: b105 = -301.0_wp/82.0_wp
    real(wp),parameter :: b106 = 2133.0_wp/4100.0_wp
    real(wp),parameter :: b107 = 45.0_wp/82.0_wp
    real(wp),parameter :: b108 = 45.0_wp/164.0_wp
    real(wp),parameter :: b109 = 18.0_wp/41.0_wp
    real(wp),parameter :: b110 = 3.0_wp/205.0_wp
    real(wp),parameter :: b115 = -6.0_wp/41.0_wp
    real(wp),parameter :: b116 = -3.0_wp/205.0_wp
    real(wp),parameter :: b117 = -3.0_wp/41.0_wp
    real(wp),parameter :: b118 = 3.0_wp/41.0_wp
    real(wp),parameter :: b119 = 6.0_wp/41.0_wp
    real(wp),parameter :: b120 = -1777.0_wp/4100.0_wp
    real(wp),parameter :: b123 = -341.0_wp/164.0_wp
    real(wp),parameter :: b124 = 4496.0_wp/1025.0_wp
    real(wp),parameter :: b125 = -289.0_wp/82.0_wp
    real(wp),parameter :: b126 = 2193.0_wp/4100.0_wp
    real(wp),parameter :: b127 = 51.0_wp/82.0_wp
    real(wp),parameter :: b128 = 33.0_wp/164.0_wp
    real(wp),parameter :: b129 = 12.0_wp/41.0_wp
    !real(wp),parameter :: b1211 = 1.0_wp

    real(wp),parameter :: c5  = 34.0_wp/105.0_wp
    real(wp),parameter :: c6  = 9.0_wp/35.0_wp
    real(wp),parameter :: c7  = 9.0_wp/35.0_wp
    real(wp),parameter :: c8  = 9.0_wp/280.0_wp
    real(wp),parameter :: c9  = 9.0_wp/280.0_wp
    real(wp),parameter :: c11 = 41.0_wp/840.0_wp
    real(wp),parameter :: c12 = 41.0_wp/840.0_wp

    associate (f0 => me%funcs(:,1), &
               f1 => me%funcs(:,2), &
               f2 => me%funcs(:,3), &
               f3 => me%funcs(:,4), &
               f4 => me%funcs(:,5), &
               f5 => me%funcs(:,6), &
               f6 => me%funcs(:,7), &
               f7 => me%funcs(:,8), &
               f8 => me%funcs(:,9), &
               f9  => me%funcs(:,10), &
               f10 => me%funcs(:,11), &
               f11 => me%funcs(:,12), &
               f12 => me%funcs(:,13))

        call me%f(t,x,f0)
        call me%f(t+h*a1,x+f0*b10*h,f1)
        call me%f(t+h*a2,x+(f0*b20+f1*b21)*h,f2)
        call me%f(t+h*a3,x+(f0*b30+f2*b32)*h,f3)
        call me%f(t+h*a4,x+(f0*b40+f2*b42+f3*b43)*h,f4)
        call me%f(t+h*a5,x+(f0*b50+f3*b53+f4*b54)*h,f5)
        call me%f(t+h*a6,x+(f0*b60+f3*b63+f4*b64+f5*b65)*h,f6)
        call me%f(t+h*a7,x+(f0*b70+f4*b74+f5*b75+f6*b76)*h,f7)
        call me%f(t+h*a8,x+(f0*b80+f3*b83+f4*b84+f5*b85+f6*b86+&
                    f7*b87)*h,f8)
        call me%f(t+h*a9,x+(f0*b90+f3*b93+f4*b94+f5*b95+f6*b96+&
                    f7*b97+f8*b98)*h,f9)
        call me%f(t+h,x+(f0*b100+f3*b103+f4*b104+f5*b105+&
                    f6*b106+f7*b107+f8*b108+f9*b109)*h,f10)
        call me%f(t,x+(f0*b110+f5*b115+f6*b116+f7*b117+f8*b118+&
                    f9*b119)*h,f11)
        call me%f(t+h,x+(f0*b120+f3*b123+f4*b124+f5*b125+f6*b126+&
                    f7*b127+f8*b128+f9*b129+f11)*h,f12)

        xf = x + h*(f5*c5+f6*c6+f7*c7+f8*c8+f9*c9+f11*c11+f12*c12)

        xerr = (41.0_wp/840.0_wp)*h*(f0+f10-f11-f12)
        ! is this negative ? does it matter ?

    end associate

    end procedure rkf78
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dormand-Prince RK8(7)13M method.
!
!### Reference
!  * P.J. Prince and J.R. Dormand, J. Comp. Appl. Math.,7, pp. 67-75, 1981
!  * [rksuite_90](http://www.netlib.org/ode/rksuite/)

    module procedure rkdp87

    real(wp),parameter :: b21   = 1.0_wp / 18.0_wp                     ! 5.55555555555555555555555555556e-2_wp
    real(wp),parameter :: b31   = 1.0_wp / 48.0_wp                     ! 2.08333333333333333333333333333e-2_wp
    real(wp),parameter :: b32   = 1.0_wp / 16.0_wp                     ! 6.25e-2_wp
    real(wp),parameter :: b41   = 1.0_wp / 32.0_wp                     ! 3.125e-2_wp
    real(wp),parameter :: b43   = 3.0_wp  / 32.0_wp                    ! 9.375e-2_wp
    real(wp),parameter :: b51   = 5.0_wp / 16.0_wp                     ! 3.125e-1_wp
    real(wp),parameter :: b53   = -75.0_wp / 64.0_wp                   ! -1.171875_wp
    real(wp),parameter :: b54   = 75.0_wp / 64.0_wp                    ! 1.171875_wp
    real(wp),parameter :: b61   = 3.0_wp / 80.0_wp                     ! 3.75e-2_wp
    real(wp),parameter :: b64   = 3.0_wp / 16.0_wp                     ! 1.875e-1_wp
    real(wp),parameter :: b65   = 3.0_wp / 20.0_wp                     ! 1.5e-1_wp
    real(wp),parameter :: b71   = 29443841.0_wp / 614563906.0_wp       ! 4.79101371111111111111111111111e-2_wp
    real(wp),parameter :: b74   = 77736538.0_wp / 692538347.0_wp       ! 1.12248712777777777777777777778e-1_wp
    real(wp),parameter :: b75   = -28693883.0_wp / 1125000000.0_wp     ! -2.55056737777777777777777777778e-2_wp
    real(wp),parameter :: b76   = 23124283.0_wp / 1800000000.0_wp      ! 1.28468238888888888888888888889e-2_wp
    real(wp),parameter :: b81   = 16016141.0_wp / 946692911.0_wp       ! 1.6917989787292281181431107136e-2_wp
    real(wp),parameter :: b84   = 61564180.0_wp / 158732637.0_wp       ! 3.87848278486043169526545744159e-1_wp
    real(wp),parameter :: b85   = 22789713.0_wp / 633445777.0_wp       ! 3.59773698515003278967008896348e-2_wp
    real(wp),parameter :: b86   = 545815736.0_wp / 2771057229.0_wp     ! 1.96970214215666060156715256072e-1_wp
    real(wp),parameter :: b87   = -180193667.0_wp / 1043307555.0_wp    ! -1.72713852340501838761392997002e-1_wp
    real(wp),parameter :: b91   = 39632708.0_wp / 573591083.0_wp       ! 6.90957533591923006485645489846e-2_wp
    real(wp),parameter :: b94   = -433636366.0_wp / 683701615.0_wp     ! -6.34247976728854151882807874972e-1_wp
    real(wp),parameter :: b95   = -421739975.0_wp / 2616292301.0_wp    ! -1.61197575224604080366876923982e-1_wp
    real(wp),parameter :: b96   = 100302831.0_wp / 723423059.0_wp      ! 1.38650309458825255419866950133e-1_wp
    real(wp),parameter :: b97   = 790204164.0_wp / 839813087.0_wp      ! 9.4092861403575626972423968413e-1_wp
    real(wp),parameter :: b98   = 800635310.0_wp / 3783071287.0_wp     ! 2.11636326481943981855372117132e-1_wp
    real(wp),parameter :: b101  = 246121993.0_wp / 1340847787.0_wp     ! 1.83556996839045385489806023537e-1_wp
    real(wp),parameter :: b104  = -37695042795.0_wp / 15268766246.0_wp ! -2.46876808431559245274431575997_wp
    real(wp),parameter :: b105  = -309121744.0_wp / 1061227803.0_wp    ! -2.91286887816300456388002572804e-1_wp
    real(wp),parameter :: b106  = -12992083.0_wp / 490766935.0_wp      ! -2.6473020233117375688439799466e-2_wp
    real(wp),parameter :: b107  = 6005943493.0_wp / 2108947869.0_wp    ! 2.84783876419280044916451825422_wp
    real(wp),parameter :: b108  = 393006217.0_wp / 1396673457.0_wp     ! 2.81387331469849792539403641827e-1_wp
    real(wp),parameter :: b109  = 123872331.0_wp / 1001029789.0_wp     ! 1.23744899863314657627030212664e-1_wp
    real(wp),parameter :: b111  = -1028468189.0_wp / 846180014.0_wp    ! -1.21542481739588805916051052503_wp
    real(wp),parameter :: b114  = 8478235783.0_wp / 508512852.0_wp     ! 1.66726086659457724322804132886e1_wp
    real(wp),parameter :: b115  = 1311729495.0_wp / 1432422823.0_wp    ! 9.15741828416817960595718650451e-1_wp
    real(wp),parameter :: b116  = -10304129995.0_wp / 1701304382.0_wp  ! -6.05660580435747094755450554309_wp
    real(wp),parameter :: b117  = -48777925059.0_wp / 3047939560.0_wp  ! -1.60035735941561781118417064101e1_wp
    real(wp),parameter :: b118  = 15336726248.0_wp / 1032824649.0_wp   ! 1.4849303086297662557545391898e1_wp
    real(wp),parameter :: b119  = -45442868181.0_wp / 3398467696.0_wp  ! -1.33715757352898493182930413962e1_wp
    real(wp),parameter :: b1110 = 3065993473.0_wp / 597172653.0_wp     ! 5.13418264817963793317325361166_wp
    real(wp),parameter :: b121  = 185892177.0_wp / 718116043.0_wp      ! 2.58860916438264283815730932232e-1_wp
    real(wp),parameter :: b124  = -3185094517.0_wp / 667107341.0_wp    ! -4.77448578548920511231011750971_wp
    real(wp),parameter :: b125  = -477755414.0_wp / 1098053517.0_wp    ! -4.3509301377703250944070041181e-1_wp
    real(wp),parameter :: b126  = -703635378.0_wp / 230739211.0_wp     ! -3.04948333207224150956051286631_wp
    real(wp),parameter :: b127  = 5731566787.0_wp / 1027545527.0_wp    ! 5.57792003993609911742367663447_wp
    real(wp),parameter :: b128  = 5232866602.0_wp / 850066563.0_wp     ! 6.15583158986104009733868912669_wp
    real(wp),parameter :: b129  = -4093664535.0_wp / 808688257.0_wp    ! -5.06210458673693837007740643391_wp
    real(wp),parameter :: b1210 = 3962137247.0_wp / 1805957418.0_wp    ! 2.19392617318067906127491429047_wp
    real(wp),parameter :: b1211 = 65686358.0_wp / 487910083.0_wp       ! 1.34627998659334941535726237887e-1_wp
    real(wp),parameter :: b131  = 403863854.0_wp / 491063109.0_wp      ! 8.22427599626507477963168204773e-1_wp
    real(wp),parameter :: b134  = -5068492393.0_wp / 434740067.0_wp    ! -1.16586732572776642839765530355e1_wp
    real(wp),parameter :: b135  = -411421997.0_wp / 543043805.0_wp     ! -7.57622116690936195881116154088e-1_wp
    real(wp),parameter :: b136  = 652783627.0_wp / 914296604.0_wp      ! 7.13973588159581527978269282765e-1_wp
    real(wp),parameter :: b137  = 11173962825.0_wp / 925320556.0_wp    ! 1.20757749868900567395661704486e1_wp
    real(wp),parameter :: b138  = -13158990841.0_wp / 6184727034.0_wp  ! -2.12765911392040265639082085897_wp
    real(wp),parameter :: b139  = 3936647629.0_wp / 1978049680.0_wp    ! 1.99016620704895541832807169835_wp
    real(wp),parameter :: b1310 = -160528059.0_wp / 685178525.0_wp     ! -2.34286471544040292660294691857e-1_wp
    real(wp),parameter :: b1311 = 248638103.0_wp / 1413531060.0_wp     ! 1.7589857770794226507310510589e-1_wp

    ! order 8
    real(wp),parameter :: c1  = 14005451.0_wp / 335480064.0_wp      ! 4.17474911415302462220859284685e-2_wp
    real(wp),parameter :: c6  = -59238493.0_wp / 1068277825.0_wp    ! -5.54523286112393089615218946547e-2_wp
    real(wp),parameter :: c7  = 181606767.0_wp / 758867731.0_wp     ! 2.39312807201180097046747354249e-1_wp
    real(wp),parameter :: c8  = 561292985.0_wp / 797845732.0_wp     ! 7.0351066940344302305804641089e-1_wp
    real(wp),parameter :: c9  = -1041891430.0_wp / 1371343529.0_wp  ! -7.59759613814460929884487677085e-1_wp
    real(wp),parameter :: c10 = 760417239.0_wp / 1151165299.0_wp    ! 6.60563030922286341461378594838e-1_wp
    real(wp),parameter :: c11 = 118820643.0_wp / 751138087.0_wp     ! 1.58187482510123335529614838601e-1_wp
    real(wp),parameter :: c12 = -528747749.0_wp / 2220607170.0_wp   ! -2.38109538752862804471863555306e-1_wp
    real(wp),parameter :: c13 = 1.0_wp / 4.0_wp                     ! 2.5e-1_wp

    ! order 7
    real(wp),parameter :: d1  = 13451932.0_wp / 455176623.0_wp      ! 2.9553213676353496981964883112e-2_wp
    real(wp),parameter :: d6  = -808719846.0_wp / 976000145.0_wp    ! -8.28606276487797039766805612689e-1_wp
    real(wp),parameter :: d7  = 1757004468.0_wp / 5645159321.0_wp   ! 3.11240900051118327929913751627e-1_wp
    real(wp),parameter :: d8  = 656045339.0_wp / 265891186.0_wp     ! 2.46734519059988698196468570407_wp
    real(wp),parameter :: d9  = -3867574721.0_wp / 1518517206.0_wp  ! -2.54694165184190873912738007542_wp
    real(wp),parameter :: d10 = 465885868.0_wp / 322736535.0_wp     ! 1.44354858367677524030187495069_wp
    real(wp),parameter :: d11 = 53011238.0_wp / 667516719.0_wp      ! 7.94155958811272872713019541622e-2_wp
    real(wp),parameter :: d12 = 2.0_wp / 45.0_wp                    ! 4.44444444444444444444444444445e-2_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10 = c10 - d10
    real(wp),parameter :: e11 = c11 - d11
    real(wp),parameter :: e12 = c12 - d12
    real(wp),parameter :: e13 = c13

    real(wp),parameter :: a2  = 1.0_wp / 18.0_wp                  ! 5.55555555555555555555555555556e-2_wp
    real(wp),parameter :: a3  = 1.0_wp / 12.0_wp                  ! 8.33333333333333333333333333334e-2_wp
    real(wp),parameter :: a4  = 1.0_wp / 8.0_wp                   ! 1.25e-1_wp
    real(wp),parameter :: a5  = 5.0_wp / 16.0_wp                  ! 3.125e-1_wp
    real(wp),parameter :: a6  = 3.0_wp / 8.0_wp                   ! 3.75e-1_wp
    real(wp),parameter :: a7  = 59.0_wp / 400.0_wp                ! 1.475e-1_wp
    real(wp),parameter :: a8  = 93.0_wp / 200.0_wp                ! 4.65e-1_wp
    real(wp),parameter :: a9  = 5490023248.0_wp / 9719169821.0_wp ! 5.64865451382259575398358501426e-1_wp
    real(wp),parameter :: a10 = 13.0_wp / 20.0_wp                 ! 6.5e-1_wp
    real(wp),parameter :: a11 = 1201146811.0_wp / 1299019798.0_wp ! 9.24656277640504446745013574318e-1_wp

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13))

        call me%f(t,      x,f1)
        call me%f(t+a2*h, x+h*(b21*f1),f2)
        call me%f(t+a3*h, x+h*(b31*f1  + b32*f2),f3)
        call me%f(t+a4*h, x+h*(b41*f1  + b43*f3),f4)
        call me%f(t+a5*h, x+h*(b51*f1  + b53*f3+b54*f4),f5)
        call me%f(t+a6*h, x+h*(b61*f1  + b64*f4  + b65 *f5),f6)
        call me%f(t+a7*h, x+h*(b71*f1  + b74*f4  + b75 *f5 + b76 *f6),f7)
        call me%f(t+a8*h, x+h*(b81*f1  + b84*f4  + b85 *f5 + b86 *f6 + b87 *f7),f8)
        call me%f(t+a9*h, x+h*(b91*f1  + b94*f4  + b95 *f5 + b96 *f6 + b97 *f7 + &
                               b98 *f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1 + b104*f4 + b105*f5 + b106*f6 + b107*f7 + &
                               b108*f8 + b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1 + b114*f4 + b115*f5 + b116*f6 + b117*f7 + &
                               b118*f8 + b119*f9 + b1110*f10),f11)
        call me%f(t+h,    x+h*(b121*f1 + b124*f4 + b125*f5 + b126*f6 + b127*f7 + &
                               b128*f8 + b129*f9 + b1210*f10 + b1211*f11),f12)
        call me%f(t+h,    x+h*(b131*f1 + b134*f4 + b135*f5 + b136*f6 + b137*f7 + &
                               b138*f8 + b139*f9 + b1310*f10 + b1311*f11),f13)

        xf = x + h*(c1*f1+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12+c13*f13)

        xerr = h*(e1*f1+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13)

    end associate

    end procedure rkdp87
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's 7(8) algorithm.
!
!### Reference
!  * [Mathematics Source Library](http://www.mymathlib.com/diffeq/embedded_runge_kutta/embedded_verner_7_8.html)

    module procedure rkv78

    real(wp),parameter :: c1  = 13.0_wp    / 288.0_wp
    real(wp),parameter :: c6  = 32.0_wp    / 125.0_wp
    real(wp),parameter :: c7  = 31213.0_wp / 144000.0_wp
    real(wp),parameter :: c8  = 2401.0_wp  / 12375.0_wp
    real(wp),parameter :: c9  = 1701.0_wp  / 14080.0_wp
    real(wp),parameter :: c10 = 2401.0_wp  / 19200.0_wp
    real(wp),parameter :: c11 = 19.0_wp    / 450.0_wp

    real(wp),parameter :: a2  = 1.0_wp / 4.0_wp
    real(wp),parameter :: a3  = 1.0_wp / 12.0_wp
    real(wp),parameter :: a4  = 1.0_wp / 8.0_wp
    real(wp),parameter :: a5  = 2.0_wp / 5.0_wp
    real(wp),parameter :: a6  = 1.0_wp / 2.0_wp
    real(wp),parameter :: a7  = 6.0_wp / 7.0_wp
    real(wp),parameter :: a8  = 1.0_wp / 7.0_wp
    real(wp),parameter :: a9  = 2.0_wp / 3.0_wp
    real(wp),parameter :: a10 = 2.0_wp / 7.0_wp
    real(wp),parameter :: a12 = 1.0_wp / 3.0_wp

    real(wp),parameter :: b31   =  5.0_wp       / 72.0_wp
    real(wp),parameter :: b32   =  1.0_wp       / 72.0_wp
    real(wp),parameter :: b41   =  1.0_wp       / 32.0_wp
    real(wp),parameter :: b43   =  3.0_wp       / 32.0_wp
    real(wp),parameter :: b51   =  106.0_wp     / 125.0_wp
    real(wp),parameter :: b53   = -408.0_wp     / 125.0_wp
    real(wp),parameter :: b54   =  352.0_wp     / 125.0_wp
    real(wp),parameter :: b61   =  1.0_wp       / 48.0_wp
    real(wp),parameter :: b64   =  8.0_wp       / 33.0_wp
    real(wp),parameter :: b65   =  125.0_wp     / 528.0_wp
    real(wp),parameter :: b71   = -13893.0_wp   / 26411.0_wp
    real(wp),parameter :: b74   =  39936.0_wp   / 26411.0_wp
    real(wp),parameter :: b75   = -64125.0_wp   / 26411.0_wp
    real(wp),parameter :: b76   =  60720.0_wp   / 26411.0_wp
    real(wp),parameter :: b81   =  37.0_wp      / 392.0_wp
    real(wp),parameter :: b85   =  1625.0_wp    / 9408.0_wp
    real(wp),parameter :: b86   = -2.0_wp       / 15.0_wp
    real(wp),parameter :: b87   =  61.0_wp      / 6720.0_wp
    real(wp),parameter :: b91   =  17176.0_wp   / 25515.0_wp
    real(wp),parameter :: b94   = -47104.0_wp   / 25515.0_wp
    real(wp),parameter :: b95   =  1325.0_wp    / 504.0_wp
    real(wp),parameter :: b96   = -41792.0_wp   / 25515.0_wp
    real(wp),parameter :: b97   =  20237.0_wp   / 145800.0_wp
    real(wp),parameter :: b98   =  4312.0_wp    / 6075.0_wp
    real(wp),parameter :: b101  = -23834.0_wp   / 180075.0_wp
    real(wp),parameter :: b104  = -77824.0_wp   / 1980825.0_wp
    real(wp),parameter :: b105  = -636635.0_wp  / 633864.0_wp
    real(wp),parameter :: b106  =  254048.0_wp  / 300125.0_wp
    real(wp),parameter :: b107  = -183.0_wp     / 7000.0_wp
    real(wp),parameter :: b108  =  8.0_wp       / 11.0_wp
    real(wp),parameter :: b109  = -324.0_wp     / 3773.0_wp
    real(wp),parameter :: b111  =  12733.0_wp   / 7600.0_wp
    real(wp),parameter :: b114  = -20032.0_wp   / 5225.0_wp
    real(wp),parameter :: b115  =  456485.0_wp  / 80256.0_wp
    real(wp),parameter :: b116  = -42599.0_wp   / 7125.0_wp
    real(wp),parameter :: b117  =  339227.0_wp  / 912000.0_wp
    real(wp),parameter :: b118  = -1029.0_wp    / 4180.0_wp
    real(wp),parameter :: b119  =  1701.0_wp    / 1408.0_wp
    real(wp),parameter :: b1110 =  5145.0_wp    / 2432.0_wp
    real(wp),parameter :: b121  = -27061.0_wp   / 204120.0_wp
    real(wp),parameter :: b124  =  40448.0_wp   / 280665.0_wp
    real(wp),parameter :: b125  = -1353775.0_wp / 1197504.0_wp
    real(wp),parameter :: b126  =  17662.0_wp   / 25515.0_wp
    real(wp),parameter :: b127  = -71687.0_wp   / 1166400.0_wp
    real(wp),parameter :: b128  =  98.0_wp      / 225.0_wp
    real(wp),parameter :: b129  =  1.0_wp       / 16.0_wp
    real(wp),parameter :: b1210 =  3773.0_wp    / 11664.0_wp
    real(wp),parameter :: b131  =  11203.0_wp   / 8680.0_wp
    real(wp),parameter :: b134  = -38144.0_wp   / 11935.0_wp
    real(wp),parameter :: b135  =  2354425.0_wp / 458304.0_wp
    real(wp),parameter :: b136  = -84046.0_wp   / 16275.0_wp
    real(wp),parameter :: b137  =  673309.0_wp  / 1636800.0_wp
    real(wp),parameter :: b138  =  4704.0_wp    / 8525.0_wp
    real(wp),parameter :: b139  =  9477.0_wp    / 10912.0_wp
    real(wp),parameter :: b1310 = -1029.0_wp    / 992.0_wp
    real(wp),parameter :: b1312 =  729.0_wp     / 341.0_wp

    real(wp),parameter :: e1  = -6600.0_wp   / 3168000.0_wp
    real(wp),parameter :: e6  = -135168.0_wp / 3168000.0_wp
    real(wp),parameter :: e7  = -14406.0_wp  / 3168000.0_wp
    real(wp),parameter :: e8  =  57624.0_wp  / 3168000.0_wp
    real(wp),parameter :: e9  =  54675.0_wp  / 3168000.0_wp
    real(wp),parameter :: e10 = -396165.0_wp / 3168000.0_wp
    real(wp),parameter :: e11 = -133760.0_wp / 3168000.0_wp
    real(wp),parameter :: e12 =  437400.0_wp / 3168000.0_wp
    real(wp),parameter :: e13 =  136400.0_wp / 3168000.0_wp

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13))

        call me%f(t,      x,f1)
        call me%f(t+a2*h, x+h*(a2*f1),f2)
        call me%f(t+a3*h, x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h, x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h, x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h, x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h, x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h, x+h*(b81*f1+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h, x+h*(b91*f1+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b104*f4+b105*f5+b106*f6+b107*f7+b108*f8+&
                               b109*f9),f10)
        call me%f(t+h,    x+h*(b111*f1+b114*f4+b115*f5+b116*f6+b117*f7+b118*f8+&
                               b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b124*f4+b125*f5+b126*f6+b127*f7+b128*f8+&
                               b129*f9+b1210*f10),f12)
        call me%f(t+h,    x+h*(b131*f1+b134*f4+b135*f5+b136*f6+b137*f7+b138*f8+&
                               b139*f9+b1310*f10+b1312*f12),f13)

        xf = x + h*(c1*f1+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11)

        xerr = h*(e1*f1+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13)

    end associate

    end procedure rkv78
!*****************************************************************************************

!*****************************************************************************************
!>
!  Verner's "Maple" (dverk78).
!
!### Reference
!  * [Coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK8/RKcoeff8c_2.pdf)

    module procedure dverk78

    real(wp),parameter :: a2  = 6.25000000000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: a3  = 1.05164319248826291079812206572769953051643192488262910798122065727699530516431924883e-1_wp
    real(wp),parameter :: a4  = 1.57746478873239436619718309859154929577464788732394366197183098591549295774647887324e-1_wp
    real(wp),parameter :: a5  = 3.90000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a6  = 4.66666666666666666666666666666666666666666666666666666666666666666666666666666666667e-1_wp
    real(wp),parameter :: a7  = 1.56000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a8  = 9.60000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a9  = 8.92267261582876569617601634626813346131102801365312428359450025561131899210951871087e-1_wp
    real(wp),parameter :: a10 = 9.16666666666666666666666666666666666666666666666666666666666666666666666666666666667e-1_wp
    real(wp),parameter :: a11 = 9.50000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp

    real(wp),parameter :: b21   = 6.25000000000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: b31   = 1.66880469042738433732284158786836826908241310145694196477771165333157001476779298640e-2_wp
    real(wp),parameter :: b32   = 8.84762723445524477065837906940862703608190614736934911503449491943838303687539950186e-2_wp
    real(wp),parameter :: b41   = 3.94366197183098591549295774647887323943661971830985915492957746478873239436619718310e-2_wp
    real(wp),parameter :: b43   = 1.18309859154929577464788732394366197183098591549295774647887323943661971830985915493e-1_wp
    real(wp),parameter :: b51   = 3.76656086575255102040816326530612244897959183673469387755102040816326530612244897959e-1_wp
    real(wp),parameter :: b53   = -1.4062762954400510204081632653061224489795918367346938775510204081632653061224489796e0_wp
    real(wp),parameter :: b54   = 1.41962020886479591836734693877551020408163265306122448979591836734693877551020408163e0_wp
    real(wp),parameter :: b61   = 4.78355492244381133270022158911047799936688825577714466603355492244381133270022158911e-2_wp
    real(wp),parameter :: b64   = 2.34464677432637812666112645898374622854105368760715435467810644685518346322874319640e-1_wp
    real(wp),parameter :: b65   = 1.84366440009590740673551804877187263818892415348179784538520472756710207016790131136e-1_wp
    real(wp),parameter :: b71   = 6.11454257142857142857142857142857142857142857142857142857142857142857142857142857143e-2_wp
    real(wp),parameter :: b74   = 1.12966312083735007492797513828957773063162531957288289301243638495099728858422070298e-1_wp
    real(wp),parameter :: b75   = -3.6523126216152081630500698710681045165713080391277981385292799325019115669575763968e-2_wp
    real(wp),parameter :: b76   = 1.84113884181313598519888991674375578168362627197039777983348751156336725254394079556e-2_wp
    real(wp),parameter :: b81   = -1.1290190344812819716873002143964435191540855755186649222833121580618541079059519076e0_wp
    real(wp),parameter :: b84   = -1.4000000000000000000000000000000000000000000000000000000000000000000000000000000000e0_wp
    real(wp),parameter :: b85   = -8.8802928704648604017469555116668654703824580186941491818867912879961932084706688647e0_wp
    real(wp),parameter :: b86   = 7.13094409677757261060034674920714023333661544255001489279410558828850531341195834447e0_wp
    real(wp),parameter :: b87   = 5.23836780816856976283390897685616875619992815166279921137599785776954200296466242780e0_wp
    real(wp),parameter :: b91   = -6.1529353321310240531758577059630882789056716837277145620689461335329804789378294402e-1_wp
    real(wp),parameter :: b94   = -4.9518730627307206123980419770218646265183975456624504456109720877843836846893436386e0_wp
    real(wp),parameter :: b95   = -5.1425568821223535028317455499886772909974359615422348128021564911695375075220376057e0_wp
    real(wp),parameter :: b96   = 4.40119257638430388133537282166699561228064232903290937377378080578921781089684014377e0_wp
    real(wp),parameter :: b97   = 7.19135317060505642064565222344470796214151883455022732703393848995802890190295274383e0_wp
    real(wp),parameter :: b98   = 9.44499265969278818394988712196051711534231335963244217175392212110442651632317187457e-3_wp
    real(wp),parameter :: b101  = -8.0487754950012651630611680571217899707181954913548370272898245110488470451176167195e-1_wp
    real(wp),parameter :: b104  = -4.5929131940948324444453893585077961341346223390010271981133375952914577811961889600e0_wp
    real(wp),parameter :: b105  = -6.5344796429871976349685844664568433689437520794427587560502836549097348190068258330e0_wp
    real(wp),parameter :: b106  = 5.43653721645848061969301952101216405498569628561398967433445354186586889917076802846e0_wp
    real(wp),parameter :: b107  = 7.41428376786239108774181233597389231941925460585428819853117348548056293090948270879e0_wp
    real(wp),parameter :: b108  = 1.62969198761492907780853867652646563859111189149676578185954495159286718271369038020e-2_wp
    real(wp),parameter :: b109  = -1.8180850948197735826159946407835863974001376137309207124952108889616530525944509467e-2_wp
    real(wp),parameter :: b111  = -1.0578581056602630497010123299247635673755482083807816057707586430609366121794742310e0_wp
    real(wp),parameter :: b114  = -2.4902747870066666124482517838284554072812073092990236244911099705917862954551978618e0_wp
    real(wp),parameter :: b115  = -8.3678653641187673586734269894057729098141483833620036736491918565911809345813949287e0_wp
    real(wp),parameter :: b116  = 6.76706391372510233295677955943576702708643348087685990605145394724880397320088890883e0_wp
    real(wp),parameter :: b117  = 6.10187575754750708348558966187981017887077884162790182554024225466956680755366458148e0_wp
    real(wp),parameter :: b118  = 2.26346926591268273851441242963489920479522697011540597935824660244066631601010686288e-2_wp
    real(wp),parameter :: b119  = 2.42115288970666032926841417691368020854458954458812827759456719481984746605406648566e-2_wp
    real(wp),parameter :: b1110 = -4.9787636043105826297506384222071115619706586609988170250163869647072076359128202386e-2_wp
    real(wp),parameter :: b121  = -2.2141388211521594228089507244618610021434452510302869093059447809375117564415940198e0_wp
    real(wp),parameter :: b124  = 4.56355804724646925572212018240182661579238223094527486516681347955670181620817635015e0_wp
    real(wp),parameter :: b125  = -1.6957036205517104297839745395366915857642603514682315438023365324576286639172315989e1_wp
    real(wp),parameter :: b126  = 1.31105217928641421165649661367441373070008285256157591924114546268073486094645142906e1_wp
    real(wp),parameter :: b127  = 2.65774224749870274100467824154490008698121216459786964112414318137911711432599035412e0_wp
    real(wp),parameter :: b128  = 2.74592934104436868684638073474053713826464018949259425746847719273619536662550065244e0_wp
    real(wp),parameter :: b129  = -1.1506549223705344448192637136257715209341884980498719927919640966806759385051562199e0_wp
    real(wp),parameter :: b1210 = 2.42302296121448799919236503348925313670760208048083006656846347016299049952947649144e0_wp
    real(wp),parameter :: b1211 = -4.1789444408283726338625504954661059040264279273698536826180777484478790720345919102e0_wp
    real(wp),parameter :: b131  = -1.5511316950026916120467765404925329008269305742185819600030671406615107733844590777e0_wp
    real(wp),parameter :: b134  = 5.46092722623409081818657767842421731964760251560133825662023106815177868484348803926e0_wp
    real(wp),parameter :: b135  = -1.1979220224305761165813607594501926859816750798269770560276587319873047174682994374e1_wp
    real(wp),parameter :: b136  = 9.35213287509985895980770600293832447863564422228774591854591472707049491657923415773e0_wp
    real(wp),parameter :: b137  = -2.6805151435786875122787398112144651989157843381870119309573397987877093382768663729e-1_wp
    real(wp),parameter :: b138  = -2.9111005885695086970119200745004024552717605363363761755145270197855298352241438457e-1_wp
    real(wp),parameter :: b139  = -5.7342608212997377743524954389733665898828809950550624616346021546025114366487668036e-1_wp
    real(wp),parameter :: b1310 = 8.49879473319296398230415986100741386767477221557113401924155562629859407659708957311e-1_wp

    real(wp),parameter :: c1  = 4.46288489653070217415613140881920071523054807312072680684782090651440586079475760217e-2_wp
    real(wp),parameter :: c6  = 3.53849672808056955006078714962021325996588002624194340841646462152138105292559006658e-1_wp
    real(wp),parameter :: c7  = 2.49214137824034653157413522433245474576906451975957952284766070502834548257871008029e-1_wp
    real(wp),parameter :: c8  = -5.2283384709802500404740481515533751847798124016472832324227905218125681954020960889e0_wp
    real(wp),parameter :: c9  = 5.15116549236861223325928760078465842871662817825128618469086436536651971244098141329e0_wp
    real(wp),parameter :: c10 = -9.5097928721643673005225329855637374707813853242252723834130793535769142379902016318e0_wp
    real(wp),parameter :: c11 = 1.00596341373431488638526005247320813694813862000943570400067314269066845716409977948e1_wp
    real(wp),parameter :: c12 = -1.2036094616454238602036053988308595036261658780444717005661665860383856284805907814e-1_wp

    real(wp),parameter :: d1  = 4.46938336974745394582507450063655006984853368676139647920294768627190144403875012741e-2_wp
    real(wp),parameter :: d6  = 3.54849078534857257151115499384161010948488673442599070019549247826104344436599329295e-1_wp
    real(wp),parameter :: d7  = 2.48972827546254539183557734416986167348192069231072959065041272432884166102268412123e-1_wp
    real(wp),parameter :: d8  = 2.66563173304884358413831674493708561666064947673721173917129352743937890031432599813e0_wp
    real(wp),parameter :: d9  = 3.92326394871838972220877485800106547858756064671589732933492272245546934442704232603e0_wp
    real(wp),parameter :: d10 = -5.6804277609880002002190441358794913525295674525389817820177718907811376030065962976e0_wp
    real(wp),parameter :: d13 = -5.5698366055781944192097144586617242171380875045541328036506435623541816671402726928e-1_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10  = c10  - d10
    real(wp),parameter :: e11  = c11
    real(wp),parameter :: e12  = c12
    real(wp),parameter :: e13  =      - d13

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13))

        call me%f(t,     x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b94*f4+b95*f5+b96*f6+b97*f7+&
                            b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b104*f4+b105*f5+b106*f6+&
                            b107*f7+b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b114*f4+b115*f5+b116*f6+&
                            b117*f7+b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+h,    x+h*(b121*f1+b124*f4+b125*f5+b126*f6+&
                            b127*f7+b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+h,    x+h*(b131*f1+b134*f4+b135*f5+b136*f6+&
                            b137*f7+b138*f8+b139*f9+b1310*f10),f13)

        xf = x+h*(c1*f1+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12)

        xerr = h*(e1*f1+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13)

    end associate

    end procedure dverk78
!*****************************************************************************************

!*****************************************************************************************
!>
!  Tsitouras & Papakostas NEW7(5) Runge-Kutta method.
!
!### Reference
!  * C. Tsitouras and S. N. Papakostas, "Cheap Error Estimation for Runge-Kutta
!    methods", SIAM J. Sci. Comput. 20(1999) 2067-2088.
!  * [Rational coefficients](http://users.uoa.gr/~tsitourasc/publications.html)
!    (see [rktp75.m](http://users.uoa.gr/~tsitourasc/rktp75.m))

    module procedure rktp75

    real(wp),parameter :: b21   = 1.0_wp / 18.0_wp
    real(wp),parameter :: b32   = 1.0_wp / 9.0_wp
    real(wp),parameter :: b41   = 1.0_wp / 24.0_wp
    real(wp),parameter :: b43   = 1.0_wp / 8.0_wp
    real(wp),parameter :: b51   = 2183971.0_wp / 4000000.0_wp
    real(wp),parameter :: b53   = -8340813.0_wp / 4000000.0_wp
    real(wp),parameter :: b54   = 3968421.0_wp / 2000000.0_wp
    real(wp),parameter :: b61   = 695768212.0_wp / 7463744411.0_wp
    real(wp),parameter :: b63   = -1803549175.0_wp / 7007942496.0_wp
    real(wp),parameter :: b64   = 3474507053.0_wp / 6790877290.0_wp
    real(wp),parameter :: b65   = 2188198899.0_wp / 15264927763.0_wp
    real(wp),parameter :: b71   = -11894934857.0_wp / 8390623634.0_wp
    real(wp),parameter :: b73   = 53094780276.0_wp / 9800512003.0_wp
    real(wp),parameter :: b74   = -8415376229.0_wp / 2277049503.0_wp
    real(wp),parameter :: b75   = -18647567697.0_wp / 10138317907.0_wp
    real(wp),parameter :: b76   = 27551494893.0_wp / 11905950217.0_wp
    real(wp),parameter :: b81   = 30828057951.0_wp / 7654644085.0_wp
    real(wp),parameter :: b83   = -4511704.0_wp / 324729.0_wp
    real(wp),parameter :: b84   = 16217851618.0_wp / 1651177175.0_wp
    real(wp),parameter :: b85   = 282768186839.0_wp / 40694064384.0_wp
    real(wp),parameter :: b86   = -104400780537.0_wp / 15869257619.0_wp
    real(wp),parameter :: b87   = 5409241639.0_wp / 9600177208.0_wp
    real(wp),parameter :: b91   = -133775720546.0_wp / 36753383835.0_wp
    real(wp),parameter :: b93   = 49608695511.0_wp / 4066590848.0_wp
    real(wp),parameter :: b94   = -59896475201.0_wp / 7901259813.0_wp
    real(wp),parameter :: b95   = -48035527651.0_wp / 5727379426.0_wp
    real(wp),parameter :: b96   = 86266718551.0_wp / 10188951048.0_wp
    real(wp),parameter :: b97   = -7751618114.0_wp / 23575802495.0_wp
    real(wp),parameter :: b98   = 2289274942.0_wp / 8464405725.0_wp

    real(wp),parameter :: a2  = 1.0_wp / 18.0_wp
    real(wp),parameter :: a3  = 1.0_wp / 9.0_wp
    real(wp),parameter :: a4  = 1.0_wp / 6.0_wp
    real(wp),parameter :: a5  = 89.0_wp / 200.0_wp
    real(wp),parameter :: a6  = 56482.0_wp / 115069.0_wp
    real(wp),parameter :: a7  = 74.0_wp / 95.0_wp
    real(wp),parameter :: a8  = 8.0_wp / 9.0_wp

    real(wp),parameter :: c1   = 597988726.0_wp / 12374436915.0_wp  ! 7th order formula
    real(wp),parameter :: c4   = 3138312158.0_wp / 11968408119.0_wp
    real(wp),parameter :: c5   = 480882843.0_wp / 7850665645.0_wp
    real(wp),parameter :: c6   = 988558885.0_wp / 3512253271.0_wp
    real(wp),parameter :: c7   = 5302636961.0_wp / 26425940286.0_wp
    real(wp),parameter :: c8   = 1259489433.0_wp / 12163586030.0_wp
    real(wp),parameter :: c9   = 1016647712.0_wp / 23899101975.0_wp

    real(wp),parameter :: d1   = 1421940313.0_wp / 46193547077.0_wp  ! 5th order formula
    real(wp),parameter :: d4   = 1943068601.0_wp / 5911217046.0_wp
    real(wp),parameter :: d5   = -3807140880.0_wp / 8205366359.0_wp
    real(wp),parameter :: d6   = 9377220888.0_wp / 11577671635.0_wp
    real(wp),parameter :: d7   = 586186883.0_wp / 5187186385.0_wp
    real(wp),parameter :: d8   = 1114095023.0_wp / 8014791121.0_wp
    !real(wp),parameter :: d9   = 1016647712.0_wp / 23899101975.0_wp ! == c9

    real(wp),parameter :: e1   = c1  - d1
    real(wp),parameter :: e4   = c4  - d4
    real(wp),parameter :: e5   = c5  - d5
    real(wp),parameter :: e6   = c6  - d6
    real(wp),parameter :: e7   = c7  - d7
    real(wp),parameter :: e8   = c8  - d8

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9))

        call me%f(t,      x,f1)
        call me%f(t+a2*h, x+h*(b21*f1),f2)
        call me%f(t+a3*h, x+h*(b32*f2),f3)
        call me%f(t+a4*h, x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h, x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h, x+h*(b61*f1+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h, x+h*(b71*f1+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h, x+h*(b81*f1+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t*h,    x+h*(b91*f1+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)

        xf = x + h*(c1*f1+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8+c9*f9)

        xerr = h*(e1*f1+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8)

    end associate

    end procedure rktp75
!*****************************************************************************************

!*****************************************************************************************
!>
!  Tanaka/Muramatsu/Yamashita 7th order (a.k.a. Tanyam7)
!
!### Reference
!  * M. Tanaka, S. Muramatsu and S. Yamashita,
!    "On the Optimization of Some Nine-Stage Seventh-order Runge-Kutta Method",
!    Information Processing Society of Japan, Vol. 33, No. 12 (1992) pages 1512-1526.
!  * [Rational coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK7/RKcoeff7d_4.pdf)
!  * [Tanyam7](https://github.com/SciML/OrdinaryDiffEq.jl/blob/f8f06c586c7a5616d7c251bce2d9e0fef0823eb4/src/tableaus/high_order_rk_tableaus.jl)
!
!@note The floating point coefficients here were generated from the rational ones from the 3rd reference.

    module procedure rktmy7

    real(wp),parameter :: a1 = 7.81664651011384679726388269101836941032920932418419855605785254026460056610810379655e-2_wp
    real(wp),parameter :: a2 = 1.17249697651707701958958240365275541154938139862762978340867788103969008491621556948e-1_wp
    real(wp),parameter :: a3 = 1.75874546477561552938437360547913311732407209794144467511301682155953512737432335422e-1_wp
    real(wp),parameter :: a4 = 4.98740110191398804083870799279316917774238190981016789826879390030550695876961641904e-1_wp
    real(wp),parameter :: a5 = 7.72121690184483890484305252570415231459539559725240465368497233385136757456140428136e-1_wp
    real(wp),parameter :: a6 = 9.91185669604776798407733864088711970429343190218936593687802104065965311344896218368e-1_wp
    real(wp),parameter :: a7 = 9.99501958209766146741222013447128336314037987005636927534919520974350847802956557173e-1_wp

    real(wp),parameter :: b21  = 7.81664651011384679726388269101836941032920932418419855605785254026460056610810379655e-2_wp
    real(wp),parameter :: b31  = 2.93124244129269254897395600913188852887345349656907445852169470259922521229053892370e-2_wp
    real(wp),parameter :: b32  = 8.79372732387807764692186802739566558662036048970722337556508410779767563687161677111e-2_wp
    real(wp),parameter :: b41  = 4.39686366193903882346093401369783279331018024485361168778254205389883781843580838556e-2_wp
    real(wp),parameter :: b43  = 1.31905909858171164703828020410934983799305407345608350633476261616965134553074251567e-1_wp
    real(wp),parameter :: b51  = 7.36183483773838301137467440432293617773745769566559352064158169276511636772818047624e-1_wp
    real(wp),parameter :: b53  = -2.8337999624233305489833449285093928914333838008174765752291037107738552931533021862e0_wp
    real(wp),parameter :: b54  = 2.59635658884089105192974828735641619143387622223193401299182493152789435225744578050e0_wp
    real(wp),parameter :: b61  = -1.2062819391370865506697620184830239019171754039710073663311602518203645149254198062e1_wp
    real(wp),parameter :: b63  = 4.82083810017524271975659667922264091479441542231216831647946440093325968906031102739e1_wp
    real(wp),parameter :: b64  = -3.8058630464634340873444069083787729069064369254270820186967532356616235833794046365e1_wp
    real(wp),parameter :: b65  = 2.68519054443726307306002772896197417175150863058445115085298809887242084990127458112e0_wp
    real(wp),parameter :: b71  = 1.05219572763201983263505214763646609723666799568263028237820996736734142311884342584e2_wp
    real(wp),parameter :: b73  = -4.1792888626241254423220090632342880293119367792924845836746399438762371771451577296e2_wp
    real(wp),parameter :: b74  = 3.32315550449933288661672145724565517116512854686499482568452142501473088581714784505e2_wp
    real(wp),parameter :: b75  = -1.9827591183572936803332556016677309852484292769162933988083092149281072338147267880e1_wp
    real(wp),parameter :: b76  = 1.21253990245498590876383571598269791392765963386781814296174940276352447040880996631e0_wp
    real(wp),parameter :: b81  = 1.14677557186317412777743962690232735459737977190570736656211705995064581461317681933e2_wp
    real(wp),parameter :: b83  = -4.5556121698960972374080326255401158427067549223697915892004464775745251977067045545e2_wp
    real(wp),parameter :: b84  = 3.62240955539231418592330010344238763220304805879492791461256324311641976192069398751e2_wp
    real(wp),parameter :: b85  = -2.1671904421828087451986108533232639587308358187016768553336461965905412178073526990e1_wp
    real(wp),parameter :: b86  = 1.31891320071378062521147328839749773179746572184191322809881380227032622170542838742e0_wp
    real(wp),parameter :: b87  = -4.8025566150346557548532221776442175423603809038769446508148646446010785455700778379e-3_wp
    real(wp),parameter :: b91  = 1.15213348705537687993035352351751211261221390205532219553682499922866411003708635464e2_wp
    real(wp),parameter :: b93  = -4.5769356568613233007695763193867218928825178792989579576839091929408867970634224781e2_wp
    real(wp),parameter :: b94  = 3.63936882188627343156610348429345771832876086443132263964527166628432626799709477320e2_wp
    real(wp),parameter :: b95  = -2.1776682078900292925704872499507823559172972633117930893260468010116663302476870140e1_wp
    real(wp),parameter :: b96  = 1.32506708878784683255232611830210096225232882939350575184015801179589744903384223122e0_wp
    real(wp),parameter :: b97  = -4.5181909867689836302434492066200005706653127046344820697054796285731244656302707584e-3_wp
    real(wp),parameter :: b98  = -5.3202693348599590527901201245120835437960233962812632873177926101911916720679884551e-4_wp
    real(wp),parameter :: b101 = 1.15189282458001936117920705886017744028740693856879151180637835667882356288187386531e2_wp
    real(wp),parameter :: b103 = -4.5759802227164297152415973359923257831016657330505652256715384079434141772791111922e2_wp
    real(wp),parameter :: b104 = 3.63861025631214796441099277168858079581214448605485290879259097757252663203289829524e2_wp
    real(wp),parameter :: b105 = -2.1772127540275557896757277605910244539671287153847928139832256704010642220109132278e1_wp
    real(wp),parameter :: b106 = 1.32488046450743160713676531081866548206715403286785238191455539176685996800150702540e0_wp
    real(wp),parameter :: b107 = -4.5057252106918316437819055679743451410766257174024328014118092080515035470892820222e-3_wp
    real(wp),parameter :: b108 = -5.3301659494291359595525498369189704335941061044130202397950934176800791138229684609e-4_wp

    real(wp),parameter :: c1 = 5.12601424974468664207334896467928192605824917480148880826975913119881823951586777925e-2_wp
    real(wp),parameter :: c4 = 2.75216384562126285009146859937694605711798961533809877651595696709236573932734806466e-1_wp
    real(wp),parameter :: c5 = 3.36966503407105433969810392957648333413899278562014739968412391495966278358775981408e-1_wp
    real(wp),parameter :: c6 = 1.89860722449065759713312922673232536516903910783954769517234286153030491598156826344e-1_wp
    real(wp),parameter :: c7 = 8.46109941851440210639662306074261643772685843568635085555647158979307742068711096848e0_wp
    real(wp),parameter :: c8 = -1.3015941672640542180639580507865346065998831988530342657791656834495604969014688205e2_wp
    real(wp),parameter :: c9 = 1.21845013554975275354886178352695475927358276806989281447140156789492750743174944785e2_wp

    real(wp),parameter :: d1  = 5.10024175903771829479213461329816195478958717593517724635264723499155306525942551495e-2_wp
    real(wp),parameter :: d4  = 2.76139295046665466902620362521695093962636297177632180749087744362855793802931918925e-1_wp
    real(wp),parameter :: d5  = 3.33788602069686023554754031916460354879718501346567668959954157311633875576165846369e-1_wp
    real(wp),parameter :: d6  = 2.01965311390813930477691388976007742664325756050214204738992402177761682278944674291e-1_wp
    real(wp),parameter :: d7  = 5.75507545904181112023147975998801190281904137375217710125886160168224573401735815006e0_wp
    real(wp),parameter :: d8  = -8.5617971085139353724114466889535156713873617800085942928170422377884412616327994845e1_wp
    real(wp),parameter :: d10 = 8.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e1_wp

    real(wp),parameter :: e1  = d1 - c1 ! -2.5772490706968347281214351381119971268661998866311561917111896207265174256442264298e-4_wp
    real(wp),parameter :: e4  = d4 - c4 ! 9.22910484539181893473502584000488250837335643822303097492047653619219870197112458608e-4_wp
    real(wp),parameter :: e5  = d5 - c5 ! -3.1779013374194104150563610411879785341807772154470710084582341843324027826101350395e-3_wp
    real(wp),parameter :: e6  = d6 - c6 ! 1.21045889417481707643784663027752061474218452662594352217581160247311906807878479474e-2_wp
    real(wp),parameter :: e7  = d7 - c7 ! -2.7060239594725909861651433007546045349078170619341737542976099881108316866697528184e0_wp
    real(wp),parameter :: e8  = d8 - c8 ! 4.45414456412660680822813381891183039461147020852174836497461459670716370738188872008e1_wp
    real(wp),parameter :: e9  =    - c9 ! -1.2184501355497527535488617835269547592735827680698928144714015678949275074317494479e2_wp
    real(wp),parameter :: e10 = d10     ! 8.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e1_wp

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10))

        call me%f(t+a1*h,x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+h,   x+h*(b81*f1+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+h,   x+h*(b91*f1+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+h,   x+h*(b101*f1+b103*f3+b104*f4+b105*f5+b106*f6+b107*f7+b108*f8),f10)

        xf = x+h*(c1*f1+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8+c9*f9)

        xerr = h*(e1*f1+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10)

    end associate

    end procedure rktmy7
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dormand-Prince 8(5)
!
!### References
!  * E. Hairer, S.P. Norsett and G. Wanner, Solving Ordinary
!    Differential Equations I. Nonstiff Problems. 2nd edition.
!    springer series in computational mathematics,
!    springer-verlag (1993)
!  * [Original DOP853 by Hairer](https://www.unige.ch/~hairer/prog/nonstiff/dop853.f)
!  * [Modernized DOP853](https://github.com/jacobwilliams/dop853)
!  * [DP8ConstantCache](https://github.com/SciML/OrdinaryDiffEq.jl/blob/master/src/tableaus/high_order_rk_tableaus.jl)

    module procedure rkdp85

    real(wp),parameter :: a7  = 1.0_wp / 4.0_wp
    real(wp),parameter :: a8  = 4.0_wp / 13.0_wp
    real(wp),parameter :: a9  = 127.0_wp / 195.0_wp
    real(wp),parameter :: a10 = 3.0_wp / 5.0_wp
    real(wp),parameter :: a11 = 6.0_wp / 7.0_wp
    real(wp),parameter :: a6  = 4.0_wp / 3.0_wp * a7
    real(wp),parameter :: a5  = (6.0_wp + sqrt(6.0_wp)) / 10.0_wp * a6
    real(wp),parameter :: a4  = (6.0_wp - sqrt(6.0_wp)) / 10.0_wp * a6
    real(wp),parameter :: a3  = 2.0_wp / 3.0_wp * a4
    real(wp),parameter :: a2  = 2.0_wp / 3.0_wp * a3

    real(wp),parameter :: b21   = 5.26001519587677318785587544488e-2_wp
    real(wp),parameter :: b31   = 1.97250569845378994544595329183e-2_wp
    real(wp),parameter :: b32   = 5.91751709536136983633785987549e-2_wp
    real(wp),parameter :: b41   = 2.95875854768068491816892993775e-2_wp
    real(wp),parameter :: b43   = 8.87627564304205475450678981324e-2_wp
    real(wp),parameter :: b51   = 2.41365134159266685502369798665e-1_wp
    real(wp),parameter :: b53   = -8.84549479328286085344864962717e-1_wp
    real(wp),parameter :: b54   = 9.24834003261792003115737966543e-1_wp
    real(wp),parameter :: b61   = 3.7037037037037037037037037037e-2_wp
    real(wp),parameter :: b64   = 1.70828608729473871279604482173e-1_wp
    real(wp),parameter :: b65   = 1.25467687566822425016691814123e-1_wp
    real(wp),parameter :: b71   = 3.7109375e-2_wp
    real(wp),parameter :: b74   = 1.70252211019544039314978060272e-1_wp
    real(wp),parameter :: b75   = 6.02165389804559606850219397283e-2_wp
    real(wp),parameter :: b76   = -1.7578125e-2_wp
    real(wp),parameter :: b81   = 3.70920001185047927108779319836e-2_wp
    real(wp),parameter :: b84   = 1.70383925712239993810214054705e-1_wp
    real(wp),parameter :: b85   = 1.07262030446373284651809199168e-1_wp
    real(wp),parameter :: b86   = -1.53194377486244017527936158236e-2_wp
    real(wp),parameter :: b87   = 8.27378916381402288758473766002e-3_wp
    real(wp),parameter :: b91   = 6.24110958716075717114429577812e-1_wp
    real(wp),parameter :: b94   = -3.36089262944694129406857109825_wp
    real(wp),parameter :: b95   = -8.68219346841726006818189891453e-1_wp
    real(wp),parameter :: b96   = 2.75920996994467083049415600797e1_wp
    real(wp),parameter :: b97   = 2.01540675504778934086186788979e1_wp
    real(wp),parameter :: b98   = -4.34898841810699588477366255144e1_wp
    real(wp),parameter :: b101  = 4.77662536438264365890433908527e-1_wp
    real(wp),parameter :: b104  = -2.48811461997166764192642586468e0_wp
    real(wp),parameter :: b105  = -5.90290826836842996371446475743e-1_wp
    real(wp),parameter :: b106  = 2.12300514481811942347288949897e1_wp
    real(wp),parameter :: b107  = 1.52792336328824235832596922938e1_wp
    real(wp),parameter :: b108  = -3.32882109689848629194453265587e1_wp
    real(wp),parameter :: b109  = -2.03312017085086261358222928593e-2_wp
    real(wp),parameter :: b111  = -9.3714243008598732571704021658e-1_wp
    real(wp),parameter :: b114  = 5.18637242884406370830023853209e0_wp
    real(wp),parameter :: b115  = 1.09143734899672957818500254654e0_wp
    real(wp),parameter :: b116  = -8.14978701074692612513997267357e0_wp
    real(wp),parameter :: b117  = -1.85200656599969598641566180701e1_wp
    real(wp),parameter :: b118  = 2.27394870993505042818970056734e1_wp
    real(wp),parameter :: b119  = 2.49360555267965238987089396762e0_wp
    real(wp),parameter :: b1110 = -3.0467644718982195003823669022e0_wp
    real(wp),parameter :: b121  = 2.27331014751653820792359768449e0_wp
    real(wp),parameter :: b124  = -1.05344954667372501984066689879e1_wp
    real(wp),parameter :: b125  = -2.00087205822486249909675718444e0_wp
    real(wp),parameter :: b126  = -1.79589318631187989172765950534e1_wp
    real(wp),parameter :: b127  = 2.79488845294199600508499808837e1_wp
    real(wp),parameter :: b128  = -2.85899827713502369474065508674e0_wp
    real(wp),parameter :: b129  = -8.87285693353062954433549289258e0_wp
    real(wp),parameter :: b1210 = 1.23605671757943030647266201528e1_wp
    real(wp),parameter :: b1211 = 6.43392746015763530355970484046e-1_wp

    real(wp),parameter :: c1  = 5.42937341165687622380535766363e-2_wp
    real(wp),parameter :: c6  = 4.45031289275240888144113950566_wp
    real(wp),parameter :: c7  = 1.89151789931450038304281599044_wp
    real(wp),parameter :: c8  = -5.8012039600105847814672114227_wp
    real(wp),parameter :: c9  = 3.1116436695781989440891606237e-1_wp
    real(wp),parameter :: c10 = -1.52160949662516078556178806805e-1_wp
    real(wp),parameter :: c11 = 2.01365400804030348374776537501e-1_wp
    real(wp),parameter :: c12 = 4.47106157277725905176885569043e-2_wp

    real(wp),parameter :: e1  = 0.1312004499419488073250102996e-01_wp
    real(wp),parameter :: e6  = -0.1225156446376204440720569753e+01_wp
    real(wp),parameter :: e7  = -0.4957589496572501915214079952_wp
    real(wp),parameter :: e8  = 0.1664377182454986536961530415e+01_wp
    real(wp),parameter :: e9  = -0.3503288487499736816886487290_wp
    real(wp),parameter :: e10 = 0.3341791187130174790297318841_wp
    real(wp),parameter :: e11 = 0.8192320648511571246570742613e-01_wp
    real(wp),parameter :: e12 = -0.2235530786388629525884427845e-01_wp

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b104*f4+b105*f5+b106*f6+b107*f7+&
                               b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b114*f4+b115*f5+b116*f6+b117*f7+&
                               b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t,      x+h*(b121*f1+b124*f4+b125*f5+b126*f6+b127*f7+&
                               b128*f8+b129*f9+b1210*f10+b1211*f11),f12)

        xf = x+h*(c1*f1+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12)

        xerr = h*(e1*f1+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12)

    end associate

    end procedure rkdp85
!*****************************************************************************************

!*****************************************************************************************
!>
!  Tsitouras & Papakostas NEW8(6): 12-stage, 8th and 6th order Runge-Kutta method.
!
!### Reference
!  * C. Tsitouras and S. N. Papakostas, "Cheap Error Estimation for Runge-Kutta
!    methods", SIAM J. Sci. Comput. 20(1999) 2067-2088.
!  * [Matlab version](https://www.mathworks.com/matlabcentral/fileexchange/2911-ode86)
!  * [Rational coefficients](http://users.uoa.gr/~tsitourasc/publications.html)
!    (see [rktp86.m](http://users.uoa.gr/~tsitourasc/rktp86.m))

    module procedure rktp86

    real(wp),parameter :: b21   = 9.0_wp /142.0_wp                     ! 0.06338028169014085_wp
    real(wp),parameter :: b31   = 178422123.0_wp /9178574137.0_wp      ! 0.0194389804273365_wp
    real(wp),parameter :: b32   = 685501333.0_wp /8224473205.0_wp      ! 0.0833489654490278_wp
    real(wp),parameter :: b41   = 12257.0_wp /317988.0_wp              ! 0.03854547970363662_wp
    real(wp),parameter :: b43   = 12257.0_wp /105996.0_wp              ! 0.1156364391109098_wp
    real(wp),parameter :: b51   = 2584949729.0_wp /6554704252.0_wp     ! 0.394365577701125_wp
    real(wp),parameter :: b53   = -9163901916.0_wp /6184003973.0_wp    ! -1.481871932167337_wp
    real(wp),parameter :: b54   = 26222057794.0_wp /17776421907.0_wp   ! 1.475103253691018_wp
    real(wp),parameter :: b61   = 4418011.0_wp /96055225.0_wp          ! 0.04599448910769821_wp
    real(wp),parameter :: b64   = 2947922107.0_wp /12687381736.0_wp    ! 0.2323507062639547_wp
    real(wp),parameter :: b65   = 3229973413.0_wp /17234960414.0_wp    ! 0.1874082292858813_wp
    real(wp),parameter :: b71   = 2875139539.0_wp /47877267651.0_wp    ! 0.06005228953244051_wp
    real(wp),parameter :: b74   = 2702377211.0_wp /24084535832.0_wp    ! 0.1122038319463678_wp
    real(wp),parameter :: b75   = -135707089.0_wp /4042230341.0_wp     ! -0.03357232951906142_wp
    real(wp),parameter :: b76   = 299874140.0_wp /17933325691.0_wp     ! 0.01672161344565858_wp
    real(wp),parameter :: b81   = -7872176137.0_wp /5003514694.0_wp    ! -1.573329273208686_wp
    real(wp),parameter :: b84   = -35136108789.0_wp /26684798878.0_wp  ! -1.316708773022366_wp
    real(wp),parameter :: b85   = -114433184681.0_wp /9760995895.0_wp  ! -11.72351529618177_wp
    real(wp),parameter :: b86   = 299204996517.0_wp /32851421233.0_wp  ! 9.10782502817387_wp
    real(wp),parameter :: b87   = 254.0_wp /39.0_wp                    ! 6.512820512820513_wp
    real(wp),parameter :: b91   = -3559950777.0_wp /7399971898.0_wp    ! -0.4810762562439125_wp
    real(wp),parameter :: b94   = -29299291531.0_wp /4405504148.0_wp   ! -6.650610360746391_wp
    real(wp),parameter :: b95   = -42434013379.0_wp /9366905709.0_wp   ! -4.530206099782573_wp
    real(wp),parameter :: b96   = 20642871700.0_wp /5300635453.0_wp    ! 3.894414525020157_wp
    real(wp),parameter :: b97   = 12951197050.0_wp /1499985011.0_wp    ! 8.63421764552553_wp
    real(wp),parameter :: b98   = 59527523.0_wp /6331620793.0_wp       ! 0.0094016247886815_wp
    real(wp),parameter :: b101  = -8196723582.0_wp /10570795981.0_wp   !  -0.7754121446230568_wp
    real(wp),parameter :: b104  = -46181454005.0_wp /5775132776.0_wp   !  -7.996604718235832_wp
    real(wp),parameter :: b105  = -196277106011.0_wp /29179424052.0_wp !  -6.726558607230183_wp
    real(wp),parameter :: b106  = 63575135343.0_wp /11491868333.0_wp   !  5.532184454327406_wp
    real(wp),parameter :: b107  = 120535663067.0_wp /11060780187.0_wp  !  10.89757332024991_wp
    real(wp),parameter :: b108  = 195434294.0_wp /9727139945.0_wp      !  0.02009165028004539_wp
    real(wp),parameter :: b109  = -617468037.0_wp /15757346105.0_wp    !  -0.03918604268037686_wp
    real(wp),parameter :: b111  = -6373809055.0_wp /5357779452.0_wp    ! -1.189636324544999_wp
    real(wp),parameter :: b114  = -150772749657.0_wp /21151088080.0_wp ! -7.128368483301215_wp
    real(wp),parameter :: b115  = -58076657383.0_wp /6089469394.0_wp   ! -9.53722789710108_wp
    real(wp),parameter :: b116  = 9252721190.0_wp /1221566797.0_wp     ! 7.57447010898087_wp
    real(wp),parameter :: b117  = 132381309631.0_wp /11748965576.0_wp  ! 11.26748638207092_wp
    real(wp),parameter :: b118  = 704633904.0_wp /13813696331.0_wp     ! 0.05100980122305832_wp
    real(wp),parameter :: b119  = 656417033.0_wp /8185349658.0_wp      ! 0.0801941346950826_wp
    real(wp),parameter :: b1110 = -1669806516.0_wp /10555289849.0_wp   ! -0.1581961783984735_wp
    real(wp),parameter :: b121  = -2726346953.0_wp /6954959789.0_wp    ! -0.3920003904712727_wp
    real(wp),parameter :: b124  = 24906446731.0_wp /6359105161.0_wp    ! 3.916659042493857_wp
    real(wp),parameter :: b125  = -65277767625.0_wp /23298960463.0_wp  ! -2.801745928908056_wp
    real(wp),parameter :: b126  = 39128152317.0_wp /16028215273.0_wp   ! 2.441204566481742_wp
    real(wp),parameter :: b127  = -40638357893.0_wp /16804059016.0_wp  ! -2.418365577882472_wp
    real(wp),parameter :: b128  = -7437361171.0_wp /21911114743.0_wp   ! -0.3394332629003293_wp
    real(wp),parameter :: b129  = 1040125706.0_wp /5334949109.0_wp     ! 0.1949645038310337_wp
    real(wp),parameter :: b1210 = -1129865134.0_wp /5812907645.0_wp    ! -0.1943717676250815_wp
    real(wp),parameter :: b1211 = 6253441118.0_wp /10543852725.0_wp    ! 0.5930888149805792_wp

    real(wp),parameter :: a2  = 9.0_wp / 142.0_wp                      ! 0.06338028169014085_wp
    real(wp),parameter :: a3  = 24514.0_wp / 238491.0_wp               ! 0.1027879458763643_wp
    real(wp),parameter :: a4  = 12257.0_wp / 79497.0_wp                ! 0.1541819188145465_wp
    real(wp),parameter :: a5  = 50.0_wp / 129.0_wp                     ! 0.3875968992248062_wp
    real(wp),parameter :: a6  = 34.0_wp / 73.0_wp                      ! 0.4657534246575343_wp
    real(wp),parameter :: a7  = 23.0_wp / 148.0_wp                     ! 0.1554054054054054_wp
    real(wp),parameter :: a8  = 142.0_wp / 141.0_wp                    ! 1.00709219858156_wp
    real(wp),parameter :: a9  = 14183175345.0_wp / 16188232343.0_wp    ! 0.876141078561489_wp
    real(wp),parameter :: a10 = 83.0_wp / 91.0_wp                      ! 0.912087912087912_wp
    real(wp),parameter :: a11 = 143.0_wp / 149.0_wp                    ! 0.959731543624161_wp

    real(wp),parameter :: c1   = 438853193.0_wp / 9881496838.0_wp      ! 0.04441161093250152_wp
    real(wp),parameter :: c6   = 11093525429.0_wp / 31342013414.0_wp   ! 0.3539506311373312_wp
    real(wp),parameter :: c7   = 481311443.0_wp / 1936695762.0_wp      ! 0.2485219684184965_wp
    real(wp),parameter :: c8   = -3375294558.0_wp / 10145424253.0_wp   ! -0.3326913171720666_wp
    real(wp),parameter :: c9   = 9830993862.0_wp / 5116981057.0_wp     ! 1.921248828652836_wp
    real(wp),parameter :: c10  = -138630849943.0_wp / 50747474617.0_wp ! -2.731778300088252_wp
    real(wp),parameter :: c11  = 7152278206.0_wp / 5104393345.0_wp     ! 1.401200440989917_wp
    real(wp),parameter :: c12  = 5118195927.0_wp / 53798651926.0_wp    ! 0.0951361371292365_wp

    real(wp),parameter :: d1   = 289283091.0_wp / 6008696510.0_wp      ! -0.003732456673269437_wp
    real(wp),parameter :: d6   = 3034152487.0_wp / 7913336319.0_wp     ! -0.02947203216019219_wp
    real(wp),parameter :: d7   = 7170564158.0_wp / 30263027435.0_wp    ! 0.01158056612815422_wp
    real(wp),parameter :: d8   = 7206303747.0_wp / 16758195910.0_wp    ! -0.7627079959184843_wp
    real(wp),parameter :: d9   = -1059739258.0_wp / 8472387467.0_wp    ! 2.046330367018225_wp
    real(wp),parameter :: d10  = 16534129531.0_wp / 11550853505.0_wp   ! -4.163198889384351_wp
    real(wp),parameter :: d11  = -3.0_wp / 2.0_wp                      ! 2.901200440989918_wp

    real(wp),parameter :: e1   = c1  - d1
    real(wp),parameter :: e6   = c6  - d6
    real(wp),parameter :: e7   = c7  - d7
    real(wp),parameter :: e8   = c8  - d8
    real(wp),parameter :: e9   = c9  - d9
    real(wp),parameter :: e10  = c10 - d10
    real(wp),parameter :: e11  = c11 - d11

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12))

        call me%f(t,      x,f1)
        call me%f(t+a2*h, x+h*(b21*f1),f2)
        call me%f(t+a3*h, x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h, x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h, x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h, x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h, x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h, x+h*(b81*f1+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h, x+h*(b91*f1+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b104*f4+b105*f5+b106*f6+b107*f7+b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b114*f4+b115*f5+b116*f6+b117*f7+b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t*h,    x+h*(b121*f1+b124*f4+b125*f5+b126*f6+b127*f7+b128*f8+b129*f9+b1210*f10+b1211*f11),f12)

        xf = x + h*(c1*f1+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12)

        xerr = h*(e1*f1+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11)

    end associate

    end procedure rktp86
!*****************************************************************************************

!*****************************************************************************************
!>
!  Enright Verner 8(7)
!
!### Reference
!  * W. H. Enright, "The Relative Efficiency of Alternative Defect Control Schemes
!    for High-Order Continuous Runge-Kutta Formulas",
!    SIAM Journal on Numerical Analysis, Vol. 30, No. 5. (Oct., 1993), pp. 1419-1445.

    module procedure rkev87

    real(wp),parameter :: a2 = 5.56000000000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: a3 = 1.02577772963604852686308492201039861351819757365684575389948006932409012131715771231e-1_wp
    real(wp),parameter :: a4 = 1.53866659445407279029462738301559792027729636048526863084922010398613518197573656846e-1_wp
    real(wp),parameter :: a5 = 3.84600000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a6 = 4.61500000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a7 = 1.53800000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a8 = 8.57100000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a9 = 9.50522279549898543264080746155892215372083974054663210568629277143998134232923935826e-1_wp
    real(wp),parameter :: a10 = 7.22200000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a11 = 9.37500000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp

    real(wp),parameter :: b21 = 5.56000000000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: b31 = 7.95367668507191489463769912616233707314310387997043291210079553116635514850736965522e-3_wp
    real(wp),parameter :: b32 = 9.46240962785329377916707930748775242786766534857141424778472114012426569832084015753e-2_wp
    real(wp),parameter :: b41 = 3.84666648613518197573656845753899480069324090121317157712305025996533795493934142114e-2_wp
    real(wp),parameter :: b43 = 1.15399994584055459272097053726169844020797227036395147313691507798960138648180242634e-1_wp
    real(wp),parameter :: b51 = 3.84391795249995742463414694148930014496173869183551650147262740706965192166135782038e-1_wp
    real(wp),parameter :: b53 = -1.4413754967533895577013957427459670216869972115118037596261793190397908645560378797e0_wp
    real(wp),parameter :: b54 = 1.44158370150339381523798104859703700719082334232825210947891657833282567238990209761e0_wp
    real(wp),parameter :: b61 = 4.61679927225245821615407392990543100173216372016370799947874024989156378037664283266e-2_wp
    real(wp),parameter :: b64 = 2.30766671478580119613626743161835287145952544968066036170902600446677938235230248901e-1_wp
    real(wp),parameter :: b65 = 1.84565335798895298224832517539110402836725817830296883834309997054406423961003322773e-1_wp
    real(wp),parameter :: b71 = 5.98340656981684945502396474355291084616704020934744687941026109734450878397200860064e-2_wp
    real(wp),parameter :: b74 = 1.11070988365806962422446377950496848082610655711386443707870631101667161815495770356e-1_wp
    real(wp),parameter :: b75 = -3.4214310915191931585007005185985812353861433634302098328491649125465037270146208540e-2_wp
    real(wp),parameter :: b76 = 1.71092568512164746123209797999598558095803758294411858265184070503527876149303521772e-2_wp
    real(wp),parameter :: b81 = -5.3795007752787302261474715462946928642562492201550692194017067742239449282090711696e-1_wp
    real(wp),parameter :: b84 = -6.9376482130983202468485892167839192961199723756828602894007094394246115740310891917e0_wp
    real(wp),parameter :: b85 = -4.6624538209733341287003258854737026552374536190626584633882834158844336356509780740e0_wp
    real(wp),parameter :: b86 = 3.99515211159952739816366225688709123778305091676102567472916353273143970250297438270e0_wp
    real(wp),parameter :: b87 = 9.0_wp
    real(wp),parameter :: b91 = -1.6324274407986590789038658936937338307876068928047009736797164524613917396815405352e0_wp
    real(wp),parameter :: b94 = -1.0827155649128676300106788103992512405570868288293796026024992785621740125352159452e1_wp
    real(wp),parameter :: b95 = -1.2412770216529556426254179980877757507654516931216971035750731127725169216761599855e1_wp
    real(wp),parameter :: b96 = 9.72736897958029723333620569705898647991679730938513196401603844136574931065998854059e0_wp
    real(wp),parameter :: b97 = 1.61993509145171614713050347911446705109706857016123854149928694302903467232127746865e1_wp
    real(wp),parameter :: b98 = -1.0384430809066835611232576348376103150240692462738613298483822870379681784453944854e-1_wp
    real(wp),parameter :: b101 = 4.37969506182387858999383959550753086479857907336724342215547180958125611435049537720e-1_wp
    real(wp),parameter :: b104 = 3.93953188040780609055484179804938065477461949555863432014746148406805874085095872275e0_wp
    real(wp),parameter :: b105 = 2.86077034685671490951716275075256969722509129259119025387241581075503457766785108234e0_wp
    real(wp),parameter :: b106 = -1.7743107088675159424540294493688191717730257147970606451633813362863675504027539094e0_wp
    real(wp),parameter :: b107 = -4.8953905177642011781075138200194429107217851579522239829952612873605231471994530162e0_wp
    real(wp),parameter :: b108 = 2.13024858819919749451668262421798635017462468974251987239988291318720261714005089224e-1_wp
    real(wp),parameter :: b109 = -5.9395365635111487961513501386239991002220291711516275316770143453048494065657506474e-2_wp
    real(wp),parameter :: b111 = -1.4741971530084051096783684495220401762953156744718612221415953506898842409410764394e0_wp
    real(wp),parameter :: b114 = -1.0994004568773831856436455451565667915965033574293371252746458409074590501516335418e1_wp
    real(wp),parameter :: b115 = -1.1347103595550558407033481028962542986773842505250429938230098582640399568841966955e1_wp
    real(wp),parameter :: b116 = 8.95698732805847648932199138433047040420587034365927857662263761085553916622223226594e0_wp
    real(wp),parameter :: b117 = 1.58937788726768491199870127047315726092390614749240868516326055420934681415123903105e1_wp
    real(wp),parameter :: b118 = -9.8752574205231410895382886685088801840992767797867610760540167194860976207728613471e-2_wp
    real(wp),parameter :: b119 = 4.88850458495220344526104177023129683434865807354905894219062956659502007705067202486e-3_wp
    real(wp),parameter :: b1110 = -4.0968137822510287105773140969344294040959548433844633187412729158670403045658227312e-3_wp
    real(wp),parameter :: b121 = -2.6300593327636131335587362170718842216114396527353516008958114998659179172525086061e0_wp
    real(wp),parameter :: b124 = -9.1742180511635168771961994716051325461520269343576437631118099519514590372371556623e0_wp
    real(wp),parameter :: b125 = -1.9181392627635586309607750914002594644548673509102310272301267608284090817942461230e1_wp
    real(wp),parameter :: b126 = 1.46425586936966487212372354845143649157868001572824958943201356999170406197491766176e1_wp
    real(wp),parameter :: b127 = 1.75293194641808432064264792329037561372134510518264860784855289743130207876560918202e1_wp
    real(wp),parameter :: b128 = -3.7191756017725564567993148040646183259986035050825855493662910931190949416983788302e-1_wp
    real(wp),parameter :: b129 = -7.0099615383151451311120261802995845194696552530300984200157512415971100577511657140e-1_wp
    real(wp),parameter :: b1210 = 5.10160166123419814594373796567007747993680583372732066087135277690995003233586244844e-2_wp
    real(wp),parameter :: b1211 = 8.35689551081652570030668604041209869059346704560318853832715091573927364648452890489e-1_wp
    real(wp),parameter :: b131 = 2.15760322561474994411124388359746470315016189993958868463227784863178137907983232935e-1_wp
    real(wp),parameter :: b134 = 8.34514732667823366560737998962040201361659061600363822667194131727201870364735969189e0_wp
    real(wp),parameter :: b135 = 2.18566238546511021703127210271410825362935054094744855941151130759498192966438607636e0_wp
    real(wp),parameter :: b136 = -1.6872364802758613307171751762566037387388156104034081739315755043398409544172643181e0_wp
    real(wp),parameter :: b137 = -8.7118979009844764699498732005474293440649791305528456441763927818947185541425600579e0_wp
    real(wp),parameter :: b138 = 2.44414593429891254272241354773491646097520831036964631654327459178700700237090814351e-2_wp
    real(wp),parameter :: b139 = 8.46378799450539085004935994965690483699123804486085782773911124365252659256135973281e-2_wp
    real(wp),parameter :: b1310 = 5.43485007267475889689554161135858132263172930458903122118464018149985401390772696043e-1_wp

    real(wp),parameter :: c1 = 4.39177036443990264767390168245612938204165409329946310979955485867750206154532516854e-2_wp
    real(wp),parameter :: c6 = 3.51024625301198059003711947067197522471081955372891881652834064609672316461737698276e-1_wp
    real(wp),parameter :: c7 = 2.46142826354923978927496217611357561294520506174049687891933955470328003987048218576e-1_wp
    real(wp),parameter :: c8 = 9.00324493052912782292485743663086167980166191074574843765667123724307486893548248673e-1_wp
    real(wp),parameter :: c9 = 4.54941872725474687205992745995487359662269870465745627778562103310520198596208865435e0_wp
    real(wp),parameter :: c10 = 4.80250151923706131251435406595354490562463095633152715752232119918349986859035476853e-3_wp
    real(wp),parameter :: c11 = -4.7410543521200460045854095885960988103202068879671787514735082136428459116707380269e0_wp
    real(wp),parameter :: c12 = -3.5457652500737177548746515059093087677430164120112009787806583305262240211772839944e-1_wp

    real(wp),parameter :: d1 = 4.33210538122143048906481448990561225500794142960488184912531716454290716526337909125e-2_wp
    real(wp),parameter :: d6 = 3.38299786239324428483566710567298525845838208179552347378195505093356569277661366407e-1_wp
    real(wp),parameter :: d7 = 2.48479842819169689968723882111796748697722834015831878169551896051150521744166052704e-1_wp
    real(wp),parameter :: d8 = 2.23788296726381915717640915256391015765072451399129220461423268948426689142758885858e-1_wp
    real(wp),parameter :: d9 = -4.0201628625024655295655283841782384327701940439930130919929624261829397741432866016e-2_wp
    real(wp),parameter :: d10 = 1.23292316936198306185965165238573766500855475391670951932862563960817353162350588140e-1_wp
    real(wp),parameter :: d13 = 6.30203320917360100491104657686662049681335571576969144866432185626491927618621819940e-2_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10  = c10  - d10
    real(wp),parameter :: e11  = c11
    real(wp),parameter :: e12  = c12
    real(wp),parameter :: e13  =      - d13

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13))

        call me%f(t+h,    x,f1)
        call me%f(t+a2*h, x+h*(b21*f1),f2)
        call me%f(t+a3*h, x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h, x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h, x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h, x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h, x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h, x+h*(b81*f1+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h, x+h*(b91*f1+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b104*f4+b105*f5+b106*f6+b107*f7+b108*f8+&
                               b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b114*f4+b115*f5+b116*f6+b117*f7+b118*f8+&
                               b119*f9+b1110*f10),f11)
        call me%f(t+h,    x+h*(b121*f1+b124*f4+b125*f5+b126*f6+b127*f7+b128*f8+&
                               b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+h,    x+h*(b131*f1+b134*f4+b135*f5+b136*f6+b137*f7+b138*f8+&
                               b139*f9+b1310*f10),f13)

        xf = x+h*(c1*f1+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12)

        xerr = h*(e1*f1+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13)

    end associate

    end procedure rkev87
!*****************************************************************************************

!*****************************************************************************************
!>
!  Fehlberg 8(9) method.
!
!### Reference
!  * E. Fehlberg, "Classical Fifth-, Sixth-, Seventh-, and Eighth-Order
!    Runge-Kutta Formulas with Stepsize Control",
!    [NASA TR R-2870](https://ntrs.nasa.gov/citations/19680027281).

    module procedure rkf89

    real(wp),parameter :: a1  = 0.44368940376498183109599404281370_wp
    real(wp),parameter :: a2  = 0.66553410564747274664399106422055_wp
    real(wp),parameter :: a3  = 0.99830115847120911996598659633083_wp
    real(wp),parameter :: a4  = 0.3155_wp
    real(wp),parameter :: a5  = 0.50544100948169068626516126737384_wp
    real(wp),parameter :: a6  = 0.17142857142857142857142857142857_wp
    real(wp),parameter :: a7  = 0.82857142857142857142857142857143_wp
    real(wp),parameter :: a8  = 0.66543966121011562534953769255586_wp
    real(wp),parameter :: a9  = 0.24878317968062652069722274560771_wp
    real(wp),parameter :: a10 = 0.1090_wp
    real(wp),parameter :: a11 = 0.8910_wp
    real(wp),parameter :: a12 = 0.3995_wp
    real(wp),parameter :: a13 = 0.6005_wp
    real(wp),parameter :: a14 = 1.0_wp
    real(wp),parameter :: a16 = 1.0_wp

    real(wp),parameter :: b1    = 0.44368940376498183109599404281370_wp
    real(wp),parameter :: b20   = 0.16638352641186818666099776605514_wp
    real(wp),parameter :: b21   = 0.49915057923560455998299329816541_wp
    real(wp),parameter :: b30   = 0.24957528961780227999149664908271_wp
    real(wp),parameter :: b32   = 0.74872586885340683997448994724812_wp
    real(wp),parameter :: b40   = 0.20661891163400602426556710393185_wp
    real(wp),parameter :: b42   = 0.17707880377986347040380997288319_wp
    real(wp),parameter :: b43   = -0.68197715413869494669377076815048e-1_wp
    real(wp),parameter :: b50   = 0.10927823152666408227903890926157_wp
    real(wp),parameter :: b53   = 0.40215962642367995421990563690087e-2_wp
    real(wp),parameter :: b54   = 0.39214118169078980444392330174325_wp
    real(wp),parameter :: b60   = 0.98899281409164665304844765434355e-1_wp
    real(wp),parameter :: b63   = 0.35138370227963966951204487356703e-2_wp
    real(wp),parameter :: b64   = 0.12476099983160016621520625872489_wp
    real(wp),parameter :: b65   = -0.55745546834989799643742901466348e-1_wp
    real(wp),parameter :: b70   = -0.36806865286242203724153101080691_wp
    real(wp),parameter :: b74   = -0.22273897469476007645024020944166e+1_wp
    real(wp),parameter :: b75   = 0.13742908256702910729565691245744e+1_wp
    real(wp),parameter :: b76   = 0.20497390027111603002159354092206e+1_wp
    real(wp),parameter :: b80   = 0.45467962641347150077351950603349e-1_wp
    real(wp),parameter :: b85   = 0.32542131701589147114677469648853_wp
    real(wp),parameter :: b86   = 0.28476660138527908888182420573687_wp
    real(wp),parameter :: b87   = 0.97837801675979152435868397271099e-2_wp
    real(wp),parameter :: b90   = 0.60842071062622057051094145205182e-1_wp
    real(wp),parameter :: b95   = -0.21184565744037007526325275251206e-1_wp
    real(wp),parameter :: b96   = 0.19596557266170831957464490662983_wp
    real(wp),parameter :: b97   = -0.42742640364817603675144835342899e-2_wp
    real(wp),parameter :: b98   = 0.17434365736814911965323452558189e-1_wp
    real(wp),parameter :: b100  = 0.54059783296931917365785724111182e-1_wp
    real(wp),parameter :: b106  = 0.11029825597828926530283127648228_wp
    real(wp),parameter :: b107  = -0.12565008520072556414147763782250e-2_wp
    real(wp),parameter :: b108  = 0.36790043477581460136384043566339e-2_wp
    real(wp),parameter :: b109  = -0.57780542770972073040840628571866e-1_wp
    real(wp),parameter :: b110  = 0.12732477068667114646645181799160_wp
    real(wp),parameter :: b117  = 0.11448805006396105323658875721817_wp
    real(wp),parameter :: b118  = 0.28773020709697992776202201849198_wp
    real(wp),parameter :: b119  = 0.50945379459611363153735885079465_wp
    real(wp),parameter :: b1110 = -0.14799682244372575900242144449640_wp
    real(wp),parameter :: b120  = -0.36526793876616740535848544394333e-2_wp
    real(wp),parameter :: b125  = 0.81629896012318919777819421247030e-1_wp
    real(wp),parameter :: b126  = -0.38607735635693506490517694343215_wp
    real(wp),parameter :: b127  = 0.30862242924605106450474166025206e-1_wp
    real(wp),parameter :: b128  = -0.58077254528320602815829374733518e-1_wp
    real(wp),parameter :: b129  = 0.33598659328884971493143451362322_wp
    real(wp),parameter :: b1210 = 0.41066880401949958613549622786417_wp
    real(wp),parameter :: b1211 = -0.11840245972355985520633156154536e-1_wp
    real(wp),parameter :: b130  =  -0.12375357921245143254979096135669e+1_wp
    real(wp),parameter :: b135  =  -0.24430768551354785358734861366763e+2_wp
    real(wp),parameter :: b136  = 0.54779568932778656050436528991173_wp
    real(wp),parameter :: b137  =  -0.44413863533413246374959896569346e+1_wp
    real(wp),parameter :: b138  = 0.10013104813713266094792617851022e+2_wp
    real(wp),parameter :: b139  =  -0.14995773102051758447170985073142e+2_wp
    real(wp),parameter :: b1310 = 0.58946948523217013620824539651427e+1_wp
    real(wp),parameter :: b1311 = 0.17380377503428984877616857440542e+1_wp
    real(wp),parameter :: b1312 = 0.27512330693166730263758622860276e+2_wp
    real(wp),parameter :: b140  = -0.35260859388334522700502958875588_wp
    real(wp),parameter :: b145  = -0.18396103144848270375044198988231_wp
    real(wp),parameter :: b146  = -0.65570189449741645138006879985251_wp
    real(wp),parameter :: b147  = -0.39086144880439863435025520241310_wp
    real(wp),parameter :: b148  = 0.26794646712850022936584423271209_wp
    real(wp),parameter :: b149  = -0.10383022991382490865769858507427e+1_wp
    real(wp),parameter :: b1410 = 0.16672327324258671664727346168501e+1_wp
    real(wp),parameter :: b1411 = 0.49551925855315977067732967071441_wp
    real(wp),parameter :: b1412 = 0.11394001132397063228586738141784e+1_wp
    real(wp),parameter :: b1413 = 0.51336696424658613688199097191534e-1_wp
    real(wp),parameter :: b150  = 0.10464847340614810391873002406755e-2_wp
    real(wp),parameter :: b158  = -0.67163886844990282237778446178020e-2_wp
    real(wp),parameter :: b159  = 0.81828762189425021265330065248999e-2_wp
    real(wp),parameter :: b1510 = -0.42640342864483347277142138087561e-2_wp
    real(wp),parameter :: b1511 = 0.28009029474168936545976331153703e-3_wp
    real(wp),parameter :: b1512 = -0.87835333876238676639057813145633e-2_wp
    real(wp),parameter :: b1513 = 0.10254505110825558084217769664009e-1_wp
    real(wp),parameter :: b160  = -0.13536550786174067080442168889966e+1_wp
    real(wp),parameter :: b165  = -0.18396103144848270375044198988231_wp
    real(wp),parameter :: b166  = -0.65570189449741645138006879985251_wp
    real(wp),parameter :: b167  = -0.39086144880439863435025520241310_wp
    real(wp),parameter :: b168  = 0.27466285581299925758962207732989_wp
    real(wp),parameter :: b169  = -0.10464851753571915887035188572676e+1_wp
    real(wp),parameter :: b1610 = 0.16714967667123155012004488306588e+1_wp
    real(wp),parameter :: b1611 = 0.49523916825841808131186990740287_wp
    real(wp),parameter :: b1612 = 0.11481836466273301905225795954930e+1_wp
    real(wp),parameter :: b1613 = 0.41082191313833055603981327527525e-1_wp
    real(wp),parameter :: b1615 =  1.0_wp

    real(wp),parameter :: c0  = 0.32256083500216249913612900960247e-1_wp
    real(wp),parameter :: c8  = 0.25983725283715403018887023171963_wp
    real(wp),parameter :: c9  = 0.92847805996577027788063714302190e-1_wp
    real(wp),parameter :: c10 = 0.16452339514764342891647731842800_wp
    real(wp),parameter :: c11 = 0.17665951637860074367084298397547_wp
    real(wp),parameter :: c12 = 0.23920102320352759374108933320941_wp
    real(wp),parameter :: c13 = 0.39484274604202853746752118829325e-2_wp
    real(wp),parameter :: c14 = 0.30726495475860640406368305522124e-1_wp

    associate (f0 => me%funcs(:,1), &
               f1 => me%funcs(:,2), &
               f2 => me%funcs(:,3), &
               f3 => me%funcs(:,4), &
               f4 => me%funcs(:,5), &
               f5 => me%funcs(:,6), &
               f6 => me%funcs(:,7), &
               f7 => me%funcs(:,8), &
               f8 => me%funcs(:,9), &
               f9  => me%funcs(:,10), &
               f10 => me%funcs(:,11), &
               f11 => me%funcs(:,12), &
               f12 => me%funcs(:,13), &
               f13 => me%funcs(:,14), &
               f14 => me%funcs(:,15), &
               f15 => me%funcs(:,16), &
               f16 => me%funcs(:,17))

        call me%f(t,x,f0)
        call me%f(t+h*a1,x+f0*b1*h,f1)
        call me%f(t+h*a2,x+(f0*b20+f1*b21)*h,f2)
        call me%f(t+h*a3,x+(f0*b30+f2*b32)*h,f3)
        call me%f(t+h*a4,x+(f0*b40+f2*b42+f3*b43)*h,f4)
        call me%f(t+h*a5,x+(f0*b50+f3*b53+f4*b54)*h,f5)
        call me%f(t+h*a6,x+(f0*b60+f3*b63+f4*b64+f5*b65)*h,f6)
        call me%f(t+h*a7,x+(f0*b70+f4*b74+f5*b75+f6*b76)*h,f7)
        call me%f(t+h*a8,x+(f0*b80+f5*b85+f6*b86+f7*b87)*h,f8)
        call me%f(t+h*a9,x+(f0*b90+f5*b95+f6*b96+f7*b97+f8*b98)*h,f9)
        call me%f(t+h*a10,x+(f0*b100+f6*b106+f7*b107+f8*b108+&
                f9*b109)*h,f10)
        call me%f(t+h*a11,x+(f0*b110+f7*b117+f8*b118+f9*b119+&
                f10*b1110)*h,f11)
        call me%f(t+h*a12,x+(f0*b120+f5*b125+f6*b126+f7*b127+&
                f8*b128+f9*b129+f10*b1210+f11*b1211)*h,f12)
        call me%f(t+h*a13,x+(f0*b130+f5*b135+f6*b136+f7*b137+&
                f8*b138+f9*b139+f10*b1310+f11*b1311+f12*b1312)*h,f13)
        call me%f(t+h*a14,x+(f0*b140+f5*b145+f6*b146+f7*b147+f8*b148+&
                f9*b149+f10*b1410+f11*b1411+f12*b1412+f13*b1413)*h,f14)
        call me%f(t,x+(f0*b150+f8*b158+f9*b159+f10*b1510+f11*b1511+&
                f12*b1512+f13*b1513)*h,f15)
        call me%f(t+h*a16,x+(f0*b160+f5*b165+f6*b166+f7*b167+f8*b168+&
                f9*b169+f10*b1610+f11*b1611+f12*b1612+f13*b1613+&
                f15*b1615)*h,f16)

        xf = x+h*(f0*c0+f8*c8+f9*c9+f10*c10+f11*c11+f12*c12+f13*c13+f14*c14)

        xerr = c14*h*(f0+f14-f15-f16)

    end associate

    end procedure rkf89
!*****************************************************************************************

!*****************************************************************************************
!>
!  Runge Kutta Verner 8(9)
!
!### Reference
!  * J. H. Verner, "Explicit Runge-Kutta Methods with Estimates of the
!    Local Truncation Error", SIAM Journal on Numerical Analysis,
!   15(4), 772-790, 1978.

    module procedure rkv89

    real(wp),parameter :: s6 = sqrt(6.0_wp)

    real(wp),parameter :: a2  = 1.0_wp/12.0_wp
    real(wp),parameter :: a3  = 1.0_wp/9.0_wp
    real(wp),parameter :: a4  = 1.0_wp/6.0_wp
    real(wp),parameter :: a5  = 2.0_wp*(1.0_wp+s6)/15.0_wp
    real(wp),parameter :: a6  = (6.0_wp+s6)/15.0_wp
    real(wp),parameter :: a7  = (6.0_wp-s6)/15.0_wp
    real(wp),parameter :: a8  = 2.0_wp/3.0_wp
    real(wp),parameter :: a9  = 1.0_wp/2.0_wp
    real(wp),parameter :: a10 = 1.0_wp/3.0_wp
    real(wp),parameter :: a11 = 1.0_wp/4.0_wp
    real(wp),parameter :: a12 = 4.0_wp/3.0_wp
    real(wp),parameter :: a13 = 5.0_wp/6.0_wp
    real(wp),parameter :: a15 = 1.0_wp/6.0_wp

    real(wp),parameter :: b31   = 1.0_wp/27.0_wp
    real(wp),parameter :: b32   = 2.0_wp/27.0_wp
    real(wp),parameter :: b41   = 1.0_wp/24.0_wp
    real(wp),parameter :: b43   = 3.0_wp/24.0_wp
    real(wp),parameter :: b51   = (4.0_wp+94.0_wp*s6)/375.0_wp
    real(wp),parameter :: b53   = -(282.0_wp+252.0_wp*s6)/375.0_wp
    real(wp),parameter :: b54   = (328.0_wp+208.0_wp*s6)/375.0_wp
    real(wp),parameter :: b61   = (9.0_wp-s6)/150.0_wp
    real(wp),parameter :: b64   = (312.0_wp+32.0_wp*s6)/1425.0_wp
    real(wp),parameter :: b65   = (69.0_wp+29.0_wp*s6)/570.0_wp
    real(wp),parameter :: b71   = (927.0_wp-347.0_wp*s6)/1250.0_wp
    real(wp),parameter :: b74   = (-16248.0_wp+7328.0_wp*s6)/9375.0_wp
    real(wp),parameter :: b75   = (-489.0_wp+179.0_wp*s6)/3750.0_wp
    real(wp),parameter :: b76   = (14268.0_wp-5798.0_wp*s6)/9375.0_wp
    real(wp),parameter :: b81   = 4.0_wp/54.0_wp
    real(wp),parameter :: b86   = (16.0_wp-s6)/54.0_wp
    real(wp),parameter :: b87   = (16.0_wp+s6)/54.0_wp
    real(wp),parameter :: b91   = 38.0_wp/512.0_wp
    real(wp),parameter :: b96   = (118.0_wp-23.0_wp*s6)/512.0_wp
    real(wp),parameter :: b97   = (118.0_wp+23.0_wp*s6)/512.0_wp
    real(wp),parameter :: b98   = -18.0_wp/512.0_wp
    real(wp),parameter :: b101  = 11.0_wp/144.0_wp
    real(wp),parameter :: b106  = (266.0_wp-s6)/864.0_wp
    real(wp),parameter :: b107  = (266.0_wp+s6)/864.0_wp
    real(wp),parameter :: b108  = -1.0_wp/16.0_wp
    real(wp),parameter :: b109  = -8.0_wp/27.0_wp
    real(wp),parameter :: b111  = (5034.0_wp-271.0_wp*s6)/61440.0_wp
    real(wp),parameter :: b117  = (7859.0_wp-1626.0_wp*s6)/10240.0_wp
    real(wp),parameter :: b118  = (-2232.0_wp+813.0_wp*s6)/20480.0_wp
    real(wp),parameter :: b119  = (-594.0_wp+271.0_wp*s6)/960.0_wp
    real(wp),parameter :: b1110 = (657.0_wp-813.0_wp*s6)/5120.0_wp
    real(wp),parameter :: b121  = (5996.0_wp-3794.0_wp*s6)/405.0_wp
    real(wp),parameter :: b126  = (-4342.0_wp-338.0_wp*s6)/9.0_wp
    real(wp),parameter :: b127  = (154922.0_wp-40458.0_wp*s6)/135.0_wp
    real(wp),parameter :: b128  = (-4176.0_wp+3794.0_wp*s6)/45.0_wp
    real(wp),parameter :: b129  = (-340864.0_wp+242816.0_wp*s6)/405.0_wp
    real(wp),parameter :: b1210 = (26304.0_wp-15176.0_wp*s6)/45.0_wp
    real(wp),parameter :: b1211 = -26624.0_wp/81.0_wp
    real(wp),parameter :: b131  = (3793.0_wp+2168.0_wp*s6)/103680.0_wp
    real(wp),parameter :: b136  = (4042.0_wp+2263.0_wp*s6)/13824.0_wp
    real(wp),parameter :: b137  = (-231278.0_wp+40717.0_wp*s6)/69120.0_wp
    real(wp),parameter :: b138  = (7947.0_wp-2168.0_wp*s6)/11520.0_wp
    real(wp),parameter :: b139  = (1048.0_wp-542.0_wp*s6)/405.0_wp
    real(wp),parameter :: b1310 = (-1383.0_wp+542.0_wp*s6)/720.0_wp
    real(wp),parameter :: b1311 = 2624.0_wp/1053.0_wp
    real(wp),parameter :: b1312 = 3.0_wp/1664.0_wp
    real(wp),parameter :: b141  = -137.0_wp/1296.0_wp
    real(wp),parameter :: b146  = (5642.0_wp-337.0_wp*s6)/864.0_wp
    real(wp),parameter :: b147  = (5642.0_wp+337.0_wp*s6)/864.0_wp
    real(wp),parameter :: b148  = -299.0_wp/48.0_wp
    real(wp),parameter :: b149  = 184.0_wp/81.0_wp
    real(wp),parameter :: b1410 = -44.0_wp/9.0_wp
    real(wp),parameter :: b1411 = -5120.0_wp/1053.0_wp
    real(wp),parameter :: b1412 = -11.0_wp/468.0_wp
    real(wp),parameter :: b1413 = 16.0_wp/9.0_wp
    real(wp),parameter :: b151  = (33617.0_wp-2168.0_wp*s6)/518400.0_wp
    real(wp),parameter :: b156  = (-3846.0_wp+31.0_wp*s6)/13824.0_wp
    real(wp),parameter :: b157  = (155338.0_wp-52807.0_wp*s6)/345600.0_wp
    real(wp),parameter :: b158  = (-12537.0_wp+2168.0_wp*s6)/57600.0_wp
    real(wp),parameter :: b159  = (92.0_wp+542.0_wp*s6)/2025.0_wp
    real(wp),parameter :: b1510 = (-1797.0_wp-542.0_wp*s6)/3600.0_wp
    real(wp),parameter :: b1511 = 320.0_wp/567.0_wp
    real(wp),parameter :: b1512 = -1.0_wp/1920.0_wp
    real(wp),parameter :: b1513 = 4.0_wp/105.0_wp
    real(wp),parameter :: b161  = (-36487.0_wp-30352.0_wp*s6)/279600.0_wp
    real(wp),parameter :: b166  = (-29666.0_wp-4499.0_wp*s6)/7456.0_wp
    real(wp),parameter :: b167  = (2779182.0_wp-615973.0_wp*s6)/186400.0_wp
    real(wp),parameter :: b168  = (-94329.0_wp+91056.0_wp*s6)/93200.0_wp
    real(wp),parameter :: b169  = (-232192.0_wp+121408.0_wp*s6)/17475.0_wp
    real(wp),parameter :: b1610 = (101226.0_wp-22764.0_wp*s6)/5825.0_wp
    real(wp),parameter :: b1611 = -169984.0_wp/9087.0_wp
    real(wp),parameter :: b1612 = -87.0_wp/30290.0_wp
    real(wp),parameter :: b1613 = 492.0_wp/1165.0_wp
    real(wp),parameter :: b1615 = 1260.0_wp/233.0_wp

    real(wp),parameter :: c1  = 103.0_wp/1680.0_wp
    real(wp),parameter :: c8  = -27.0_wp/140.0_wp
    real(wp),parameter :: c9  = 76.0_wp/105.0_wp
    real(wp),parameter :: c10 = -201.0_wp/280.0_wp
    real(wp),parameter :: c11 = 1024.0_wp/1365.0_wp
    real(wp),parameter :: c12 = 3.0_wp/7280.0_wp
    real(wp),parameter :: c13 = 12.0_wp/35.0_wp
    real(wp),parameter :: c14 = 9.0_wp/280.0_wp

    real(wp),parameter :: e1  = -1911.0_wp/109200.0_wp
    real(wp),parameter :: e8  = 34398.0_wp/109200.0_wp
    real(wp),parameter :: e9  = -61152.0_wp/109200.0_wp
    real(wp),parameter :: e10 = 114660.0_wp/109200.0_wp
    real(wp),parameter :: e11 = -114688.0_wp/109200.0_wp
    real(wp),parameter :: e12 = -63.0_wp/109200.0_wp
    real(wp),parameter :: e13 = -13104.0_wp/109200.0_wp
    real(wp),parameter :: e14 = -3510.0_wp/109200.0_wp
    real(wp),parameter :: e15 = 39312.0_wp/109200.0_wp
    real(wp),parameter :: e16 = 6058.0_wp/109200.0_wp

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13), &
               f14 => me%funcs(:,14), &
               f15 => me%funcs(:,15), &
               f16 => me%funcs(:,16))

        call me%f(t,x,f1)
        call me%f(t+a2*h,x+h*(a2*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b106*f6+b107*f7+b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b117*f7+b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b126*f6+b127*f7+b128*f8+b129*f9+&
                    b1210*f10+b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b131*f1+b136*f6+b137*f7+b138*f8+b139*f9+&
                    b1310*f10+b1311*f11+b1312*f12),f13)
        call me%f(t+h,x+h*(b141*f1+b146*f6+b147*f7+b148*f8+b149*f9+b1410*f10+&
                    b1411*f11+b1412*f12+b1413*f13),f14)
        call me%f(t+a15*h,x+h*(b151*f1+b156*f6+b157*f7+b158*f8+b159*f9+b1510*f10+&
                    b1511*f11+b1512*f12+b1513*f13),f15)
        call me%f(t+h,x+h*(b161*f1+b166*f6+b167*f7+b168*f8+b169*f9+b1610*f10+&
                    b1611*f11+b1612*f12+b1613*f13+b1615*f15),f16)

        xf = x+h*(c1*f1+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12+c13*f13+c14*f14)

        xerr = h*(e1*f1+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13+e14*f14+e15*f15+e16*f16)

    end associate

    end procedure rkv89
!*****************************************************************************************

!*****************************************************************************************
!>
!  Tsitouras combined order 9(8) Runge-Kutta scheme A.
!
!### Reference
!  * Ch. Tsitouras, "Optimized explicit Runge-Kutta pairs of order 9(8)",
!    Applied Numerical Mathematics, 38 (2001) 123-134.
!  * [Coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK9/RKcoeff9b_1.pdf)

    module procedure rkt98a

    real(wp),parameter :: a2  = .2040816326530612244897959183673469387755102040816326530612244897959183673469387755102e-1_wp
    real(wp),parameter :: a3  = .8813293914998103008637915939241511230858121177369086646032265578682773616755710635430e-1_wp
    real(wp),parameter :: a4  = .1321994087249715451295687390886226684628718176605362996904839836802416042513356595314_wp
    real(wp),parameter :: a5  = .4285714285714285714285714285714285714285714285714285714285714285714285714285714285714_wp
    real(wp),parameter :: a6  = .5364755392243287681395100999813264375851395225813758811703296868095847858703057159708_wp
    real(wp),parameter :: a7  = .2254292226804331366223946619234354671767652393233860235915750750951771188915990459340_wp
    real(wp),parameter :: a8  = .6349206349206349206349206349206349206349206349206349206349206349206349206349206349206_wp
    real(wp),parameter :: a9  = .4761904761904761904761904761904761904761904761904761904761904761904761904761904761905_wp
    real(wp),parameter :: a10 = 1.055555555555555555555555555555555555555555555555555555555555555555555555555555555556_wp
    real(wp),parameter :: a11 = .7777777777777777777777777777777777777777777777777777777777777777777777777777777777778_wp
    real(wp),parameter :: a12 = .1474169624260994760478401587101660956180420042646161803585052367559161448324216358210_wp
    real(wp),parameter :: a13 = .9375_wp
    real(wp),parameter :: a14 = .975_wp

    real(wp),parameter :: b21   = .2040816326530612244897959183673469387755102040816326530612244897959183673469387755102e-1_wp
    real(wp),parameter :: b31   = -.1021687274487683147769723606535303770211896504959152824720881616958241197178852529434_wp
    real(wp),parameter :: b32   = .1903016665987493448633515200459454893297708622696061489324108174826518558854423592977_wp
    real(wp),parameter :: b41   = .3304985218124288628239218477215566711571795441513407492262099592006040106283391488286e-1_wp
    real(wp),parameter :: b43   = .9914955654372865884717655431646700134715386324540222476786298776018120318850174464858e-1_wp
    real(wp),parameter :: b51   = .9439263832171291883556670672811977550046744677129090342789268623214132543830715838461_wp
    real(wp),parameter :: b53   = -3.630115063093482037077114300073110289641773922761700400575782997067611129391489310214_wp
    real(wp),parameter :: b54   = 3.114760108447781420150018661363341106065670883620219937725427563317626446436989154939_wp
    real(wp),parameter :: b61   = .2056923328616261711782848556408566522300362983753802984037858427390204762957826067614e-1_wp
    real(wp),parameter :: b64   = .2604824183197402987552340723087785980849282464799395547519553603229784733524946656053_wp
    real(wp),parameter :: b65   = .2554238876184258522664475421084621742772076462638982965779957422127042648882327896893_wp
    real(wp),parameter :: b71   = .4318483705109263028835194255758669000001323804220689783671769343298616516403726355252e-1_wp
    real(wp),parameter :: b74   = .1769284839807689645520842945002498666673573450203833044058383710943493957801668363127_wp
    real(wp),parameter :: b75   = .7715414621876983794315795767271715967327448157389024021256382790156907385043923034755e-2_wp
    real(wp),parameter :: b76   = -.2399512973305442012357370901672805457932791896593202672237372222315349437648976965990e-2_wp
    real(wp),parameter :: b81   = .7054673721340388007054673721340388007054673721340388007054673721340388007054673721340e-1_wp
    real(wp),parameter :: b86   = .2389860715558522381270320269011306632810238539566724845073599194488366776462554669041_wp
    real(wp),parameter :: b87   = .3253878261513788024373418708061003772833500437505585560570139782583943629181184308032_wp
    real(wp),parameter :: b91   = .7068452380952380952380952380952380952380952380952380952380952380952380952380952380952e-1_wp
    real(wp),parameter :: b96   = .1146981694865827822571846471014964620252663838260725205469644176957364421846758824211_wp
    real(wp),parameter :: b97   = .3242899257515124558380534481365987760699717114120227175482736775423587959105622128170_wp
    real(wp),parameter :: b98   = -.3348214285714285714285714285714285714285714285714285714285714285714285714285714285714e-1_wp
    real(wp),parameter :: b101  = .3828039568865740740740740740740740740740740740740740740740740740740740740740740740741_wp
    real(wp),parameter :: b106  = -26.28104570763124221540974003468081947028606261133704083478349807677632293296095495655_wp
    real(wp),parameter :: b107  = -1.748424643449004698170506878899427443294184302243206078796748836803923980619291957032_wp
    real(wp),parameter :: b108  = 8.7107466796875_wp
    real(wp),parameter :: b109  = 19.99147527006172839506172839506172839506172839506172839506172839506172839506172839506_wp
    real(wp),parameter :: b111  = .5775788905284689042365660435243263588459243112336590023001577005230527581251090920939e-1_wp
    real(wp),parameter :: b116  = .76_wp
    real(wp),parameter :: b117  = .4018553287915560050245043972333732847203199147775749045996022855946715894896836169645_wp
    real(wp),parameter :: b118  = .6626589626062798854099077293240416432854142660801015341291381713518812647417195998478e-1_wp
    real(wp),parameter :: b119  = -.5150696954976787388652658281684056884304318728657759205657506333922693112630203285629_wp
    real(wp),parameter :: b1110 = .6968359170425632653891831427973381274755878134602740100996538387882097264431620182023e-2_wp
    real(wp),parameter :: b121  = .6994103035906395715686489597582150022600838999784034656313893569317997374351308591645e-1_wp
    real(wp),parameter :: b126  = -.1206166060565470490794750098024428160239891793653556676806869026855595723414769418541_wp
    real(wp),parameter :: b127  = .1326963017314458683010911057621835709637828919128544554925159376674464191887834137030_wp
    real(wp),parameter :: b128  = .2031532538077627922059108299487811850373691592058111305356962697560067377203999264697_wp
    real(wp),parameter :: b129  = -.7548640175023889467177231137902524280970115417595113975094208653922765875893450604889e-1_wp
    real(wp),parameter :: b1210 = .3789103686439011295643140618205140110558445127918206929994555147217682251832650190596e-2_wp
    real(wp),parameter :: b1211 = -.6605971935182620916042249241335724188598654843850115173121147228314743697169599255572e-1_wp
    real(wp),parameter :: b131  = .6878836016002817993752515099013030041019173079702440987327706356852301558017737433194_wp
    real(wp),parameter :: b136  = -3.658808527366451517196946897056976218504208192589081007775064092270089158621094491595_wp
    real(wp),parameter :: b137  = 5.845503378564669890314102776059891925747966958611333131112559317698644813283634646256_wp
    real(wp),parameter :: b138  = 5.411412752004688858844863696709300082070631013287790747126942423884017249319289413984_wp
    real(wp),parameter :: b139  = -1.582095988269673568992894109712396974043574141076643374159775089066374202615336309544_wp
    real(wp),parameter :: b1310 = .9945248865706299667788846635419585634800918858472131439118016352130452805983178532975e-1_wp
    real(wp),parameter :: b1311 = -1.329593487105234404098771549303922926244482088345956610970650644487259700977298059950_wp
    real(wp),parameter :: b1312 = -4.536254218085344054923493892951394749476260046442408298457962714965473684250800727799_wp
    real(wp),parameter :: b141  = .8075602959459265637003617465986992156750424242709811245269878174761587489984770070790_wp
    real(wp),parameter :: b146  = -7.835750139557992772198249358151231578755580407654683121798887409512496231362858547948_wp
    real(wp),parameter :: b147  = 6.174200025941229770123184423146954177664197990793732802584504944027199598185299195828_wp
    real(wp),parameter :: b148  = 7.111311052783799148723822283846058963658585479679811643355604800924488305671070483858_wp
    real(wp),parameter :: b149  = 1.126413982652184268358560776506518924191271296281999997786352511455840556312961264408_wp
    real(wp),parameter :: b1410 = .1138838501589788637237573235699770016403876336825059128989257049544545698305937427348_wp
    real(wp),parameter :: b1411 = -1.434774678607913207528871735039082294103761283429995473644501054462770211806990144762_wp
    real(wp),parameter :: b1412 = -5.064314742122892898464863809545209351959607473464534474081890709755483907888716564210_wp
    real(wp),parameter :: b1413 = -.2352964719331973643770165093268505801053566015981841162709660510739142793983643698681e-1_wp
    real(wp),parameter :: b151  = .8101674645254722523232876762637530370910583262175105775933018316393086542069855754015_wp
    real(wp),parameter :: b156  = -11.92841245677955898272791316012006359130670880301152754325896491950568584401093883276_wp
    real(wp),parameter :: b157  = 5.364312898490360554625912421385267920693980736534035810795373384314434259846321616737_wp
    real(wp),parameter :: b158  = 8.059954459759511795757228635532122933456240951117021423047033797467500489813129448587_wp
    real(wp),parameter :: b159  = 4.618063439001721483758160286714057772615332315540952842312029459375380582946251202099_wp
    real(wp),parameter :: b1510 = .1128601659181894505188796244218558774756228116520111404880300475877949274532116052613_wp
    real(wp),parameter :: b1511 = -1.315207706945002903017387913180865178001333874822374203640801449667142626811434784416_wp
    real(wp),parameter :: b1512 = -4.689892399903609125554297739238443813841792924181292032608257980483530407060762654113_wp
    real(wp),parameter :: b1513 = .6721374554329681653036240358001936392819624140428230276841754345869866338681915837036e-2_wp
    real(wp),parameter :: b1514 = -.3856723862141420733690607213568689457521916318676624500458592507392990272144509263344e-1_wp
    real(wp),parameter :: b161  = .7081867006810967019347772392865651712376321476396730568850847595024088622008058350820_wp
    real(wp),parameter :: b166  = -11.35112404359189673946143658241126665471485966274613155801176412146453899273199024599_wp
    real(wp),parameter :: b167  = 4.439836852763619566651386303128039432194840880353792439043045399385549660055953598204_wp
    real(wp),parameter :: b168  = 7.160528966855846780330688332694210230082237210593958603950542885510151860241609708315_wp
    real(wp),parameter :: b169  = 4.944752289066181389539519911976094691889352151362799884071179481037363375554630307368_wp
    real(wp),parameter :: b1610 = .9819535610867987425927035019405551701168012197193983668338466050076531884883427929320e-1_wp
    real(wp),parameter :: b1611 = -1.023835739158951766446520340031695827655784562575140398451539531509310691718294156648_wp
    real(wp),parameter :: b1612 = -3.935185699031579433419722235418354733643010565663948821548988100702690247383937415641_wp
    real(wp),parameter :: b1613 = -.2186260374551646198917238469418771662941856089051019740795631241695933765898516172461e-1_wp
    real(wp),parameter :: b1614 = -.1949207994747991139879059472346010977266916004643284521298911984273980740862674825940e-1_wp

    real(wp),parameter :: c1  = .4153556008805959168895898921867201921459132555672551671643920848892579905394945027033e-1_wp
    real(wp),parameter :: c8  = -.4252280874169805589327845543703559909471287347438636516704849004566228405918914211737_wp
    real(wp),parameter :: c9  = .4911269629417608841870622286295630148000839362387724323042417629187646342598357379221_wp
    real(wp),parameter :: c10 = .4541782417588474254373676812638040365345450854512311037987257309434313503153553995250_wp
    real(wp),parameter :: c11 = 1.006032649094428065183781872196896559026233580812897577690817173850511464781870422371_wp
    real(wp),parameter :: c12 = .2396980714287725918428584296767300637015778481901015943085981245718488414129283930348_wp
    real(wp),parameter :: c13 = -4.455491297731407466229188603532881977410859597554953958343743900169006073451667494575_wp
    real(wp),parameter :: c14 = 9.288978775706101824452250592593505318906312756833224753343846524732597445318251252011_wp
    real(wp),parameter :: c15 = -6.410061645100351588399537404906702274594586970014904598917670494111219851867862508616_wp
    real(wp),parameter :: c16 = .7692307692307692307692307692307692307692307692307692307692307692307692307692307692308_wp

    real(wp),parameter :: d1  = .3999250134104311539086473895855274325427452646000198100229842429713750574052674683348e-1_wp
    real(wp),parameter :: d8  = .4011325788584545964379727040448630014695811775770738360963521200500104385050340643944e-1_wp
    real(wp),parameter :: d9  = .3925127095625233142417946847596602750245079278580359243007567443309335209235627547386_wp
    real(wp),parameter :: d10 = -1.051197094764485752167499211132792259622314942028765288915320450078532078636064340575_wp
    real(wp),parameter :: d11 = -.2505083298864024791053208779904093804217869849340728193482139880406203405457882973484_wp
    real(wp),parameter :: d12 = .2465428288939690633785173766277520619424499467256191777438442341744800285916120244248_wp
    real(wp),parameter :: d13 = 6.839965112378725948727677787215511467384166222758213034610265848044786251996046243769_wp
    real(wp),parameter :: d14 = -16.02665175464198790087906253807353043847748558382750862377249679396395516268962930751_wp
    real(wp),parameter :: d15 = 10.0_wp
    real(wp),parameter :: d16 = .7692307692307692307692307692307692307692307692307692307692307692307692307692307692308_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10 = c10 - d10
    real(wp),parameter :: e11 = c11 - d11
    real(wp),parameter :: e12 = c12 - d12
    real(wp),parameter :: e13 = c13 - d13
    real(wp),parameter :: e14 = c14 - d14
    real(wp),parameter :: e15 = c15 - d15
    real(wp),parameter :: e16 = c16 - d16

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13), &
               f14 => me%funcs(:,14), &
               f15 => me%funcs(:,15), &
               f16 => me%funcs(:,16))

        call me%f(t,x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b106*f6+b107*f7+b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b116*f6+b117*f7+b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b126*f6+b127*f7+b128*f8+b129*f9+b1210*f10+&
                               b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b131*f1+b136*f6+b137*f7+b138*f8+b139*f9+b1310*f10+&
                               b1311*f11+b1312*f12),f13)
        call me%f(t+a14*h,x+h*(b141*f1+b146*f6+b147*f7+b148*f8+b149*f9+b1410*f10+&
                               b1411*f11+b1412*f12+b1413*f13),f14)
        call me%f(t+h,x+h*(b151*f1+b156*f6+b157*f7+b158*f8+b159*f9+b1510*f10+b1511*f11+&
                               b1512*f12+b1513*f13+b1514*f14),f15)
        call me%f(t+h,x+h*(b161*f1+b166*f6+b167*f7+b168*f8+b169*f9+b1610*f10+b1611*f11+&
                               b1612*f12+b1613*f13+b1614*f14),f16)

        xf = x+h*(c1*f1+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16)

        xerr = h*(e1*f1+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13+e14*f14+e15*f15+e16*f16)

    end associate

    end procedure rkt98a
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sharp 9(8)
!
!### Reference
!  * Philip W. Sharp, "High order explicit Runge-Kutta pairs for ephemerides of the Solar System and the Moon",
!    Journal of Applied Mathematics & Decision Sciences, 4(2), 183-192 (2000),
!  * https://github.com/SciML/DiffEqDevTools.jl/blob/master/src/ode_tableaus.jl
!
!@note The floating point coefficients here were generated from the rational ones from the 2nd reference.

    module procedure rks98

    real(wp),parameter :: a2  = 2.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: a3  = 9.66220283800537781651565568869531890039129854640901213159807057025513165826690659288e-2_wp
    real(wp),parameter :: a4  = 1.44933042570080667247734835330429783505869478196135181973971058553826974874003598893e-1_wp
    real(wp),parameter :: a5  = 3.11111111111111111111111111111111111111111111111111111111111111111111111111111111111e-1_wp
    real(wp),parameter :: a6  = 3.51498773299780208885007017507765081905783415195317477342800010797639951702224305104e-1_wp
    real(wp),parameter :: a7  = 1.47701226700219791114992982492234918094216584804682522657199989202360048297775694896e-1_wp
    real(wp),parameter :: a8  = 4.16000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a9  = 3.12000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a10 = 1.05000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: a11 = 5.87002096436058700209643605870020964360587002096436058700209643605870020964360587002e-1_wp
    real(wp),parameter :: a12 = 6.83818495973559057197009883791356922634313147538567178788853669325746483233925079300e-1_wp
    real(wp),parameter :: a13 = 8.79003558718861209964412811387900355871886120996441281138790035587188612099644128114e-1_wp
    real(wp),parameter :: a14 = 9.16000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp

    real(wp),parameter :: b21   = 2.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: b31   = -1.3677338082684416521484257374903677512719677066274978006646616817250635870287502595e-1_wp
    real(wp),parameter :: b32   = 2.33395409206897943379999130635989964131109756126839901382446873875057675285544091877e-1_wp
    real(wp),parameter :: b41   = 3.62332606425201668119337088326074458764673695490337954934927646384567437185008997233e-2_wp
    real(wp),parameter :: b43   = 1.08699781927560500435801126497822337629402108647101386480478293915370231155502699170e-1_wp
    real(wp),parameter :: b51   = 1.93101200330205381129632393375249728863798118161035750556803447904498441444456645259e-1_wp
    real(wp),parameter :: b53   = -6.4770998560312388589569152951216097849732130945128175207869408225841451367147303820e-1_wp
    real(wp),parameter :: b54   = 7.65719896384029615877170247248022360744634302401357112633001745465027183338127504049e-1_wp
    real(wp),parameter :: b61   = 4.77439440833806232500188714529354943400951843563751259289256424427344686758530074484e-2_wp
    real(wp),parameter :: b64   = 1.96932175795354826650179476488669790168656892603887237546677770684075728220363810345e-1_wp
    real(wp),parameter :: b65   = 1.06822653421044758984808669566159797397031338235055113867196597670829754806007487311e-1_wp
    real(wp),parameter :: b71   = 5.35653520740100302121103332775188587103176832181869303040831220432325640950516524391e-2_wp
    real(wp),parameter :: b74   = 1.20625588947728919056244246899316675798038285391744744331826391017415590396850177402e-1_wp
    real(wp),parameter :: b75   = -6.7750756225149701080673518583857606163917837693866673698932599085495183512419090449e-2_wp
    real(wp),parameter :: b76   = 4.12610419036305429273119208992569897497784538886175217202230752272070773182929555031e-2_wp
    real(wp),parameter :: b81   = 4.62222222222222222222222222222222222222222222222222222222222222222222222222222222222e-2_wp
    real(wp),parameter :: b86   = 1.56583674083394386420831384025620810581726829112411811849222219222877791193826581916e-1_wp
    real(wp),parameter :: b87   = 2.13194103694383391356946393752156967196050948665365965928555558554899986583951195862e-1_wp
    real(wp),parameter :: b91   = 4.63125000000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: b96   = 7.51502406476090389349073807809004819189545346828427154623710864742465169193996381623e-2_wp
    real(wp),parameter :: b97   = 2.12474759352390961065092619219099518081045465317157284537628913525753483080600361838e-1_wp
    real(wp),parameter :: b98   = -2.1937500000000000000000000000000000000000000000000000000000000000000000000000000000e-2_wp
    real(wp),parameter :: b101  = 4.71370805167570855745159483211372150835054795000175063898322887854066734358040684850e-2_wp
    real(wp),parameter :: b106  = 2.09088767177245700619002712371830037495557938247264604576291946948488435273048771651e-1_wp
    real(wp),parameter :: b107  = 1.00170199845336072001441249674316806102355300084726572203232579503666601315305977833e-1_wp
    real(wp),parameter :: b108  = -4.1746535935132619792943174258604390602569938027379993697699660376037253597563110535e-2_wp
    real(wp),parameter :: b109  = -2.0964951160420623840201673610867966807884877980462868947165715486152445642659570743e-1_wp
    real(wp),parameter :: b111  = 9.62506254804954308501791978456184013295337530291135748904639061412584879741027523020e-2_wp
    real(wp),parameter :: b116  = -3.2889823940483436676586225555059183413399004899899513082733172247062419719368801086e0_wp
    real(wp),parameter :: b117  = 1.58120592819482939980472941802825082360990331356773185743485106618875616830910113284e-1_wp
    real(wp),parameter :: b118  = 1.38472174110756949611966846419084262312633249684929784010383708600379765979254809143e0_wp
    real(wp),parameter :: b119  = 2.38689153107685450091794555753665319888363091085120276623574076954818022830367973855e0_wp
    real(wp),parameter :: b1110 = -1.5000000000000000000000000000000000000000000000000000000000000000000000000000000000e-1_wp
    real(wp),parameter :: b121  = -4.1407181780829333826951338677923830282957788550316015757490123469403067795923041646e-1_wp
    real(wp),parameter :: b126  = 2.79362329285576366638884698134000509392713675687003189893845038096478481618920557795e0_wp
    real(wp),parameter :: b127  = -3.9420148894844370214218382106635188284591903158605020970611433089765204398968285433e0_wp
    real(wp),parameter :: b128  = -1.8050833493703232778218655623536664178736717202844850014226094729833601601906484738e0_wp
    real(wp),parameter :: b129  = -8.2411175299103457871280197945625000732690195817565658068273621079502124203745163579e-2_wp
    real(wp),parameter :: b1210 = 3.76724363839931467375143823390080271058360595326428259518770377058250406454692000855e0_wp
    real(wp),parameter :: b1211 = 3.66532796680637812441222026292597668018700554869965598789627155511871004748252089945e-1_wp
    real(wp),parameter :: b131  = 2.57118803244482153847626462807856681730780184786345438393372906485058163365316011894e-1_wp
    real(wp),parameter :: b136  = 3.22367498544415690028414563287198307551779583958713835991093894064695605463707833852e0_wp
    real(wp),parameter :: b137  = 3.95644809715768343577645732188103093651324923567481376817610794743631132119165042871e0_wp
    real(wp),parameter :: b138  = 1.90269610866406552179805731106267938194469931412329236667880503143330311252532589886e0_wp
    real(wp),parameter :: b139  = -5.2391055401717914982701824261828721839527015786103641427597273886467794920713496236e0_wp
    real(wp),parameter :: b1310 = -2.6899481313272987041207973348420457663574444684767403019558990390373975528574973750e0_wp
    real(wp),parameter :: b1311 = -1.6326510182468950685438619936368781848861272140333241607130497353632849427158674841e0_wp
    real(wp),parameter :: b1312 = 1.10077025395445846919296783742614641536163480794527995340824137263302194802498793285e0_wp
    real(wp),parameter :: b141  = 4.08583369788787674770196821641246834353608695949114599877074475305978651235285469124e0_wp
    real(wp),parameter :: b146  = -2.2271040812705511368267126642675337195508778742206509440033684681620630123740836983e1_wp
    real(wp),parameter :: b147  = 2.14948200377041255968873743962135304520531006690258805591495319234107676337073934764e1_wp
    real(wp),parameter :: b148  = -1.2574667847036709064117981533059871307009863478980024060295805266513615915043138735e1_wp
    real(wp),parameter :: b149  = 2.54025651049869277589813851770135045700861797774341329374584242615333576062269117780e1_wp
    real(wp),parameter :: b1410 = -2.4865730830721493170953431486117421078941351149011810402690627604263543458446989032e1_wp
    real(wp),parameter :: b1411 = 1.82479073749516085147180202069097519978378918397970543443181892233252604987304701845e1_wp
    real(wp),parameter :: b1412 = -9.4291787815393124480125158715016538695138573450523345127986553789149287859337960255e0_wp
    real(wp),parameter :: b1413 = 8.25492056472487433062307536805028087460591469502464576121882769983546032147130645121e-1_wp
    real(wp),parameter :: b151  = -2.0602212042704652199244807912284589395147988173289204904243925620741467390453149503e0_wp
    real(wp),parameter :: b156  = -2.3055564596876178326198929677498175792624263330011017577312032231421737767030559692e0_wp
    real(wp),parameter :: b157  = -1.8379647671170325920306936223300488315443940076512542750939032591666020672221106554e1_wp
    real(wp),parameter :: b158  = 9.80063543421329918398827905876150569685353384114341904048681595591515045410441012888e-1_wp
    real(wp),parameter :: b159  = 6.44589158658994320826489839207649833491520903216082715269469597755156737618623074351e0_wp
    real(wp),parameter :: b1510 = 1.69192749452203241565519919050371391213914448119219462829620457478758678378197682618e1_wp
    real(wp),parameter :: b1511 = -2.0988423575394962681126735718162381051185885412147351718156202906941587529583151802e0_wp
    real(wp),parameter :: b1512 = 1.32027878990329475096946854788916481157084451169381843791872814006897180501111028426e0_wp
    real(wp),parameter :: b1513 = 1.29768753187020321933985725604863066333871547241993002451338834678114639100987388885e-1_wp
    real(wp),parameter :: b1514 = 4.89900743459928848448110776111870354430304809243733908347583718104632373992549620659e-2_wp
    real(wp),parameter :: b161  = -1.2158274587431439669291103556715073078571331121887250332403319959071088711190342783e-1_wp
    real(wp),parameter :: b166  = -1.3009191687212087953937489554001604506419126739405937395616109697496506450739882872e1_wp
    real(wp),parameter :: b167  = -1.0172115280165551555631444224047088737887421449467885159610413901710577571727155577e1_wp
    real(wp),parameter :: b168  = -9.1376040597684553208500728397293269920338377345009007549988772949417238792208104990e0_wp
    real(wp),parameter :: b169  = 2.20271190656799492286306009834337099718633750373570345566819174880669811072575443296e1_wp
    real(wp),parameter :: b1610 = 6.06122452890837591335001768945892462299766158185132002339883804559430763723850436801e0_wp
    real(wp),parameter :: b1611 = 9.49756499938619826937987552815088909338119306846899551533354213035888156520270940573e0_wp
    real(wp),parameter :: b1612 = -4.7967103151651448708172582989203880579450635055043503690520505457817795045514591312e0_wp
    real(wp),parameter :: b1613 = 6.51295494211030686568681751222035336828933052420596087187186975501127983652453403424e-1_wp

    real(wp),parameter :: c1  = 2.96275053023277765257147159271391983963018401625397769109581065606146193195367573656e-2_wp
    real(wp),parameter :: c8  = 7.80838132003046860057542718319735163112389396834314055277737601781852412778315993413e-2_wp
    real(wp),parameter :: c9  = 2.07620199958252939386280504658544994598465755307530373969021587633956273026455132933e-1_wp
    real(wp),parameter :: c10 = 1.69689395533794848088161847963068489566652421567777796351282799972583366756393555465e-1_wp
    real(wp),parameter :: c11 = 1.51705614234815145844037646207281350814128482549399473556449026004053183774844489538e-1_wp
    real(wp),parameter :: c12 = 1.44245948285115701641251053349684808091393327368972824014235232832191510849137389325e-1_wp
    real(wp),parameter :: c13 = 1.66761530218905580484076634450811227110090575704792045335148871651256233134779966440e-1_wp
    real(wp),parameter :: c14 = 1.92534285621243594762130197569399694644610012365298208475543983102390200679811369880e-2_wp
    real(wp),parameter :: c15 = 3.30125647043589625485103058545564456472676564190264834875762168569205517930399726031e-2_wp

    real(wp),parameter :: d1  = 2.98603077327271866543455374899994274225585056921637520652930266673233165901602953218e-2_wp
    real(wp),parameter :: d8  = 6.57755975911182224511842125160987804765191941684313617933281553333569953681029565846e-2_wp
    real(wp),parameter :: d9  = 2.13956779658119140898029891022265645728716698807857092033794693592749028494691587115e-1_wp
    real(wp),parameter :: d10 = 1.68721578796992869389907340999574096253253203866862151079585793545930377244267574219e-1_wp
    real(wp),parameter :: d11 = 1.71768360736786458177293485795001832452979570477507671772344473084105210367066149609e-1_wp
    real(wp),parameter :: d12 = 1.25492680790662476798676554483056276508051861757938232565244001424475853044387070116e-1_wp
    real(wp),parameter :: d13 = 1.89263652417235144013217502395309419983557239873730317305944966458721798394126740078e-1_wp
    real(wp),parameter :: d16 = 3.51610422763585016173454752986945211743637253555094213844648898933374204971976269572e-2_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10 = c10 - d10
    real(wp),parameter :: e11 = c11 - d11
    real(wp),parameter :: e12 = c12 - d12
    real(wp),parameter :: e13 = c13 - d13
    real(wp),parameter :: e14 = c14
    real(wp),parameter :: e15 = c15
    real(wp),parameter :: e16 =     - d16

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13), &
               f14 => me%funcs(:,14), &
               f15 => me%funcs(:,15), &
               f16 => me%funcs(:,16))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b106*f6+b107*f7+b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b116*f6+b117*f7+b118*f8+b119*f9+&
                               b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b126*f6+b127*f7+b128*f8+b129*f9+&
                               b1210*f10+b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b131*f1+b136*f6+b137*f7+b138*f8+b139*f9+&
                               b1310*f10+b1311*f11+b1312*f12),f13)
        call me%f(t+a14*h,x+h*(b141*f1+b146*f6+b147*f7+b148*f8+b149*f9+&
                               b1410*f10+b1411*f11+b1412*f12+b1413*f13),f14)
        call me%f(t+h,    x+h*(b151*f1+b156*f6+b157*f7+b158*f8+b159*f9+&
                               b1510*f10+b1511*f11+b1512*f12+b1513*f13+b1514*f14),f15)
        call me%f(t+h,    x+h*(b161*f1+b166*f6+b167*f7+b168*f8+b169*f9+&
                               b1610*f10+b1611*f11+b1612*f12+b1613*f13),f16)

        xf = x+h*(c1*f1+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12+c13*f13+c14*f14+c15*f15)

        xerr = h*(e1*f1+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+e13*f13+e14*f14+e15*f15+e16*f16)

    end associate

    end procedure rks98
!*****************************************************************************************

!*****************************************************************************************
!>
!  Feagin's RK8(10) method -- a 10th-order method with an embedded 8th-order method.
!
!### Reference
!  * T. Feagin, "[A Tenth-Order Runge-Kutta Method with Error Estimate](http://sce.uhcl.edu/feagin/courses/rk10.pdf)",
!    [coefficient file](http://sce.uhcl.edu/rungekutta/rk108.txt)

    module procedure rkf108

    real(wp),parameter :: a1  = 0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: a2  = 0.539357840802981787532485197881302436857273449701009015505500_wp
    real(wp),parameter :: a3  = 0.809036761204472681298727796821953655285910174551513523258250_wp
    real(wp),parameter :: a4  = 0.309036761204472681298727796821953655285910174551513523258250_wp
    real(wp),parameter :: a5  = 0.981074190219795268254879548310562080489056746118724882027805_wp
    real(wp),parameter :: a6  = 0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a7  = 0.354017365856802376329264185948796742115824053807373968324184_wp
    real(wp),parameter :: a8  = 0.882527661964732346425501486979669075182867844268052119663791_wp
    real(wp),parameter :: a9  = 0.642615758240322548157075497020439535959501736363212695909875_wp
    real(wp),parameter :: a10 = 0.357384241759677451842924502979560464040498263636787304090125_wp
    real(wp),parameter :: a11 = 0.117472338035267653574498513020330924817132155731947880336209_wp
    real(wp),parameter :: a12 = 0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a13 = 0.309036761204472681298727796821953655285910174551513523258250_wp
    real(wp),parameter :: a14 = 0.539357840802981787532485197881302436857273449701009015505500_wp
    real(wp),parameter :: a15 = 0.100000000000000000000000000000000000000000000000000000000000_wp

    real(wp),parameter :: c0  = 0.0333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: c1  = 0.0250000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c2  = 0.0333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: c4  = 0.0500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c6  = 0.0400000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c8  = 0.189237478148923490158306404106012326238162346948625830327194_wp
    real(wp),parameter :: c9  = 0.277429188517743176508360262560654340428504319718040836339472_wp
    real(wp),parameter :: c10 = 0.277429188517743176508360262560654340428504319718040836339472_wp
    real(wp),parameter :: c11 = 0.189237478148923490158306404106012326238162346948625830327194_wp
    real(wp),parameter :: c12 = -0.0400000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c13 = -0.0500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c14 = -0.0333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: c15 = -0.0250000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c16 = 0.0333333333333333333333333333333333333333333333333333333333333_wp

    real(wp),parameter :: b10 = 0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b20 = -0.915176561375291440520015019275342154318951387664369720564660_wp
    real(wp),parameter :: b21 = 1.45453440217827322805250021715664459117622483736537873607016_wp
    real(wp),parameter :: b30 = 0.202259190301118170324681949205488413821477543637878380814562_wp
    real(wp),parameter :: b32 = 0.606777570903354510974045847616465241464432630913635142443687_wp
    real(wp),parameter :: b40 = 0.184024714708643575149100693471120664216774047979591417844635_wp
    real(wp),parameter :: b42 = 0.197966831227192369068141770510388793370637287463360401555746_wp
    real(wp),parameter :: b43 = -0.0729547847313632629185146671595558023015011608914382961421311_wp
    real(wp),parameter :: b50 = 0.0879007340206681337319777094132125475918886824944548534041378_wp
    real(wp),parameter :: b53 = 0.410459702520260645318174895920453426088035325902848695210406_wp
    real(wp),parameter :: b54 = 0.482713753678866489204726942976896106809132737721421333413261_wp
    real(wp),parameter :: b60 = 0.0859700504902460302188480225945808401411132615636600222593880_wp
    real(wp),parameter :: b63 = 0.330885963040722183948884057658753173648240154838402033448632_wp
    real(wp),parameter :: b64 = 0.489662957309450192844507011135898201178015478433790097210790_wp
    real(wp),parameter :: b65 = -0.0731856375070850736789057580558988816340355615025188195854775_wp
    real(wp),parameter :: b70 = 0.120930449125333720660378854927668953958938996999703678812621_wp
    real(wp),parameter :: b74 = 0.260124675758295622809007617838335174368108756484693361887839_wp
    real(wp),parameter :: b75 = 0.0325402621549091330158899334391231259332716675992700000776101_wp
    real(wp),parameter :: b76 = -0.0595780211817361001560122202563305121444953672762930724538856_wp
    real(wp),parameter :: b80 = 0.110854379580391483508936171010218441909425780168656559807038_wp
    real(wp),parameter :: b85 = -0.0605761488255005587620924953655516875526344415354339234619466_wp
    real(wp),parameter :: b86 = 0.321763705601778390100898799049878904081404368603077129251110_wp
    real(wp),parameter :: b87 = 0.510485725608063031577759012285123416744672137031752354067590_wp
    real(wp),parameter :: b90 = 0.112054414752879004829715002761802363003717611158172229329393_wp
    real(wp),parameter :: b95 = -0.144942775902865915672349828340980777181668499748506838876185_wp
    real(wp),parameter :: b96 = -0.333269719096256706589705211415746871709467423992115497968724_wp
    real(wp),parameter :: b97 = 0.499269229556880061353316843969978567860276816592673201240332_wp
    real(wp),parameter :: b98 = 0.509504608929686104236098690045386253986643232352989602185060_wp
    real(wp),parameter :: b100 = 0.113976783964185986138004186736901163890724752541486831640341_wp
    real(wp),parameter :: b105 = -0.0768813364203356938586214289120895270821349023390922987406384_wp
    real(wp),parameter :: b106 = 0.239527360324390649107711455271882373019741311201004119339563_wp
    real(wp),parameter :: b107 = 0.397774662368094639047830462488952104564716416343454639902613_wp
    real(wp),parameter :: b108 = 0.0107558956873607455550609147441477450257136782823280838547024_wp
    real(wp),parameter :: b109 = -0.327769124164018874147061087350233395378262992392394071906457_wp
    real(wp),parameter :: b110 = 0.0798314528280196046351426864486400322758737630423413945356284_wp
    real(wp),parameter :: b115 = -0.0520329686800603076514949887612959068721311443881683526937298_wp
    real(wp),parameter :: b116 = -0.0576954146168548881732784355283433509066159287152968723021864_wp
    real(wp),parameter :: b117 = 0.194781915712104164976306262147382871156142921354409364738090_wp
    real(wp),parameter :: b118 = 0.145384923188325069727524825977071194859203467568236523866582_wp
    real(wp),parameter :: b119 = -0.0782942710351670777553986729725692447252077047239160551335016_wp
    real(wp),parameter :: b1110 = -0.114503299361098912184303164290554670970133218405658122674674_wp
    real(wp),parameter :: b120 = 0.985115610164857280120041500306517278413646677314195559520529_wp
    real(wp),parameter :: b123 = 0.330885963040722183948884057658753173648240154838402033448632_wp
    real(wp),parameter :: b124 = 0.489662957309450192844507011135898201178015478433790097210790_wp
    real(wp),parameter :: b125 = -1.37896486574843567582112720930751902353904327148559471526397_wp
    real(wp),parameter :: b126 = -0.861164195027635666673916999665534573351026060987427093314412_wp
    real(wp),parameter :: b127 = 5.78428813637537220022999785486578436006872789689499172601856_wp
    real(wp),parameter :: b128 = 3.28807761985103566890460615937314805477268252903342356581925_wp
    real(wp),parameter :: b129 = -2.38633905093136384013422325215527866148401465975954104585807_wp
    real(wp),parameter :: b1210 = -3.25479342483643918654589367587788726747711504674780680269911_wp
    real(wp),parameter :: b1211 = -2.16343541686422982353954211300054820889678036420109999154887_wp
    real(wp),parameter :: b130 = 0.895080295771632891049613132336585138148156279241561345991710_wp
    real(wp),parameter :: b132 = 0.197966831227192369068141770510388793370637287463360401555746_wp
    real(wp),parameter :: b133 = -0.0729547847313632629185146671595558023015011608914382961421311_wp
    real(wp),parameter :: b135 = -0.851236239662007619739049371445966793289359722875702227166105_wp
    real(wp),parameter :: b136 = 0.398320112318533301719718614174373643336480918103773904231856_wp
    real(wp),parameter :: b137 = 3.63937263181035606029412920047090044132027387893977804176229_wp
    real(wp),parameter :: b138 = 1.54822877039830322365301663075174564919981736348973496313065_wp
    real(wp),parameter :: b139 = -2.12221714704053716026062427460427261025318461146260124401561_wp
    real(wp),parameter :: b1310 = -1.58350398545326172713384349625753212757269188934434237975291_wp
    real(wp),parameter :: b1311 = -1.71561608285936264922031819751349098912615880827551992973034_wp
    real(wp),parameter :: b1312 = -0.0244036405750127452135415444412216875465593598370910566069132_wp
    real(wp),parameter :: b140 = -0.915176561375291440520015019275342154318951387664369720564660_wp
    real(wp),parameter :: b141 = 1.45453440217827322805250021715664459117622483736537873607016_wp
    real(wp),parameter :: b144 = -0.777333643644968233538931228575302137803351053629547286334469_wp
    real(wp),parameter :: b146 = -0.0910895662155176069593203555807484200111889091770101799647985_wp
    real(wp),parameter :: b1412 = 0.0910895662155176069593203555807484200111889091770101799647985_wp
    real(wp),parameter :: b1413 = 0.777333643644968233538931228575302137803351053629547286334469_wp
    real(wp),parameter :: b150 = 0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b152 = -0.157178665799771163367058998273128921867183754126709419409654_wp
    real(wp),parameter :: b1514 = 0.157178665799771163367058998273128921867183754126709419409654_wp
    real(wp),parameter :: b160 = 0.181781300700095283888472062582262379650443831463199521664945_wp
    real(wp),parameter :: b161 = 0.675000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b162 = 0.342758159847189839942220553413850871742338734703958919937260_wp
    real(wp),parameter :: b164 = 0.259111214548322744512977076191767379267783684543182428778156_wp
    real(wp),parameter :: b165 = -0.358278966717952089048961276721979397739750634673268802484271_wp
    real(wp),parameter :: b166 = -1.04594895940883306095050068756409905131588123172378489286080_wp
    real(wp),parameter :: b167 = 0.930327845415626983292300564432428777137601651182965794680397_wp
    real(wp),parameter :: b168 = 1.77950959431708102446142106794824453926275743243327790536000_wp
    real(wp),parameter :: b169 = 0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b1610 = -0.282547569539044081612477785222287276408489375976211189952877_wp
    real(wp),parameter :: b1611 = -0.159327350119972549169261984373485859278031542127551931461821_wp
    real(wp),parameter :: b1612 = -0.145515894647001510860991961081084111308650130578626404945571_wp
    real(wp),parameter :: b1613 = -0.259111214548322744512977076191767379267783684543182428778156_wp
    real(wp),parameter :: b1614 = -0.342758159847189839942220553413850871742338734703958919937260_wp
    real(wp),parameter :: b1615 = -0.675000000000000000000000000000000000000000000000000000000000_wp

    associate (f0 => me%funcs(:,1), &
               f1 => me%funcs(:,2), &
               f2 => me%funcs(:,3), &
               f3 => me%funcs(:,4), &
               f4 => me%funcs(:,5), &
               f5 => me%funcs(:,6), &
               f6 => me%funcs(:,7), &
               f7 => me%funcs(:,8), &
               f8 => me%funcs(:,9), &
               f9  => me%funcs(:,10), &
               f10 => me%funcs(:,11), &
               f11 => me%funcs(:,12), &
               f12 => me%funcs(:,13), &
               f13 => me%funcs(:,14), &
               f14 => me%funcs(:,15), &
               f15 => me%funcs(:,16), &
               f16 => me%funcs(:,17))

        call me%f(t,x,f0)
        call me%f(t+a1*h,x+h*(b10*f0),f1)
        call me%f(t+a2*h,x+h*(b20*f0+b21*f1),f2)
        call me%f(t+a3*h,x+h*(b30*f0+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b40*f0+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b50*f0+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b60*f0+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b70*f0+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b80*f0+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b90*f0+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b100*f0+b105*f5+b106*f6+b107*f7+b108*f8+&
                    b109*f9),f10)
        call me%f(t+a11*h,x+h*(b110*f0+b115*f5+b116*f6+b117*f7+b118*f8+&
                    b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b120*f0+b123*f3+b124*f4+b125*f5+b126*f6+&
                    b127*f7+b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b130*f0+b132*f2+b133*f3+b135*f5+b136*f6+&
                    b137*f7+b138*f8+b139*f9+b1310*f10+b1311*f11+&
                    b1312*f12),f13)
        call me%f(t+a14*h,x+h*(b140*f0+b141*f1+b144*f4+b146*f6+b1412*f12+&
                    b1413*f13),f14)
        call me%f(t+a15*h,x+h*(b150*f0+b152*f2+b1514*f14),f15)
        call me%f(t+h,x+h*(b160*f0+b161*f1+b162*f2+b164*f4+b165*f5+&
                    b166*f6+b167*f7+b168*f8+b169*f9+b1610*f10+b1611*f11+&
                    b1612*f12+b1613*f13+b1614*f14+b1615*f15),f16)

        xf = x+h*(c0*f0+c1*f1+c2*f2+c4*f4+c6*f6+c8*f8+c9*f9+&
                c10*f10+c11*f11+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16)

        xerr = (1.0_wp/360.0_wp)*h*(f1-f15)

    end associate

    end procedure rkf108
!*****************************************************************************************

!*****************************************************************************************
!>
!  A modification of Curtis' order 10 Runge-Kutta scheme with an order 8 embedded scheme
!
!### Reference
!  * A.R.Curtis, "High-order Explicit Runge-Kutta Formulae, Their uses, and Limitations",
!    J. Inst. Maths Applics (1975) 16, 35-55.
!  * [Rational coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10a(8)_2.pdf)

    module procedure rkc108

    real(wp),parameter :: a2 = .1452518960316150517617548528770033320314511251329947060838468741983976455607179673401_wp
    real(wp),parameter :: a3 = .1452518960316150517617548528770033320314511251329947060838468741983976455607179673401_wp
    real(wp),parameter :: a4 = .2178778440474225776426322793155049980471766876994920591257703112975964683410769510101_wp
    real(wp),parameter :: a5 = .5446946101185564441065806982887624951179417192487301478144257782439911708526923775252_wp
    real(wp),parameter :: a6 = .6536335321422677329278968379465149941415300630984761773773109338927894050232308530303_wp
    real(wp),parameter :: a7 = .2746594919905254008808021630247618520892150865127407293922085868737635475402543533498_wp
    real(wp),parameter :: a8 = .7735775201106609448405825008093973718589542913426807556412662673054607938029043386501_wp
    real(wp),parameter :: a9 = .5801831400829957086304368756070480288942157185070105667309497004790955953521782539876_wp
    real(wp),parameter :: a10 = .1174723380352676535744985130203309248171321557319478803362088220814723414805867429383_wp
    real(wp),parameter :: a11 = .3573842417596774518429245029795604640404982636367873040901247917361510345429002009092_wp
    real(wp),parameter :: a12 = .6426157582403225481570754970204395359595017363632126959098752082638489654570997990908_wp
    real(wp),parameter :: a13 = .1174723380352676535744985130203309248171321557319478803362088220814723414805867429383_wp
    real(wp),parameter :: a14 = .8825276619647323464255014869796690751828678442680521196637911779185276585194132570617_wp
    real(wp),parameter :: a15 = .3573842417596774518429245029795604640404982636367873040901247917361510345429002009092_wp
    real(wp),parameter :: a16 = .6426157582403225481570754970204395359595017363632126959098752082638489654570997990908_wp
    real(wp),parameter :: a17 = .8825276619647323464255014869796690751828678442680521196637911779185276585194132570617_wp
    real(wp),parameter :: a18 = 1
    real(wp),parameter :: a19 = .3510848126232741617357001972386587771203155818540433925049309664694280078895463510848_wp
    real(wp),parameter :: a20 = .6157407407407407407407407407407407407407407407407407407407407407407407407407407407407_wp
    real(wp),parameter :: a21 = 1

    real(wp),parameter :: b21 = .1452518960316150517617548528770033320314511251329947060838468741983976455607179673401_wp
    real(wp),parameter :: b31 = .7262594801580752588087742643850166601572556256649735304192343709919882278035898367003e-1_wp
    real(wp),parameter :: b32 = .7262594801580752588087742643850166601572556256649735304192343709919882278035898367003e-1_wp
    real(wp),parameter :: b41 = .5446946101185564441065806982887624951179417192487301478144257782439911708526923775252e-1_wp
    real(wp),parameter :: b42 = 0
    real(wp),parameter :: b43 = .1634083830355669332319742094866287485353825157746190443443277334731973512558077132576_wp
    real(wp),parameter :: b51 = .5446946101185564441065806982887624951179417192487301478144257782439911708526923775252_wp
    real(wp),parameter :: b52 = 0
    real(wp),parameter :: b53 = -2.042604787944586665399677618582859356692281447182738054304096668414966890697596415720_wp
    real(wp),parameter :: b54 = 2.042604787944586665399677618582859356692281447182738054304096668414966890697596415720_wp
    real(wp),parameter :: b61 = .6536335321422677329278968379465149941415300630984761773773109338927894050232308530303e-1_wp
    real(wp),parameter :: b62 = 0
    real(wp),parameter :: b63 = 0
    real(wp),parameter :: b64 = .3268167660711338664639484189732574970707650315492380886886554669463947025116154265151_wp
    real(wp),parameter :: b65 = .2614534128569070931711587351786059976566120252393904709509243735571157620092923412121_wp
    real(wp),parameter :: b71 = .8233707757482716585173454344310125296066814318521742241762319051772963627695955263034e-1_wp
    real(wp),parameter :: b72 = 0
    real(wp),parameter :: b73 = 0
    real(wp),parameter :: b74 = .2119171963202803561687843468555305553175658807629274312902985594840086570224567152664_wp
    real(wp),parameter :: b75 = -.3997343508054218311577932550061320162379840049816347807630118786107674477850206579628e-1_wp
    real(wp),parameter :: b76 = .2037865317596006197606259822674324543477946306275935376058802473310199901934015124941e-1_wp
    real(wp),parameter :: b81 = .8595305779007343831562027786771081909543936570474230618236291858949564375587825985001e-1_wp
    real(wp),parameter :: b82 = 0
    real(wp),parameter :: b83 = 0
    real(wp),parameter :: b84 = 0
    real(wp),parameter :: b85 = 0
    real(wp),parameter :: b86 = .2911769478058850960337179621761553399856026049598393013981874594942289837064329700000_wp
    real(wp),parameter :: b87 = .3964475145147024104912442607655312127779123206780991480607158892217361663405931088001_wp
    real(wp),parameter :: b91 = .8612093485606967549983047372292119178898514571588438099912534616486575243508895957628e-1_wp
    real(wp),parameter :: b92 = 0
    real(wp),parameter :: b93 = 0
    real(wp),parameter :: b94 = 0
    real(wp),parameter :: b95 = 0
    real(wp),parameter :: b96 = .1397464826824442089036313891001189801074425314582326737716288563521183595455090268480_wp
    real(wp),parameter :: b97 = .3951098495815674599900526056001284215294125840404176924334653987770478924197803010468_wp
    real(wp),parameter :: b98 = -.4079412703708563576307759281612056453162454270752418047326990081493640904820003348350e-1_wp
    real(wp),parameter :: b101 = .7233144422337948077616348229119326315582930871089020733092900891206129381937795204778e-1_wp
    real(wp),parameter :: b102 = 0
    real(wp),parameter :: b103 = 0
    real(wp),parameter :: b104 = 0
    real(wp),parameter :: b105 = 0
    real(wp),parameter :: b106 = .2200276284689998102140972735735070061373242800181187459951219347361114857342828430157_wp
    real(wp),parameter :: b107 = .8789533425436734013369780264792573637952226487753296416823846876217040795688489371334e-1_wp
    real(wp),parameter :: b108 = -.4445383996260350863990674880611108986832860648196030000580004690002268108984238641730e-1_wp
    real(wp),parameter :: b109 = -.2183282289488754689095532966861839909872150913926337371522805434288481649401165594213_wp
    real(wp),parameter :: b111 = .8947100936731114228785441966773836169071038390882857211057269158522704971585365845223e-1_wp
    real(wp),parameter :: b112 = 0
    real(wp),parameter :: b113 = 0
    real(wp),parameter :: b114 = 0
    real(wp),parameter :: b115 = 0
    real(wp),parameter :: b116 = .3946008170285561860741397654755022300929434262701385530048127140223687993778661654316_wp
    real(wp),parameter :: b117 = .3443011367963333487713764986067104675654371857504670290688086760696354596195596354011_wp
    real(wp),parameter :: b118 = -.7946682664292661290694938113119430997053815140863772328764150866582492425892231395780e-1_wp
    real(wp),parameter :: b119 = -.3915218947895966123834967996391962853380545808840091268064277812752553499114569444180_wp
    real(wp),parameter :: b1110 = 0
    real(wp),parameter :: b121 = .3210006877963209212945282736072241886741425314298532400216927262619488479186214523312e-1_wp
    real(wp),parameter :: b122 = 0
    real(wp),parameter :: b123 = 0
    real(wp),parameter :: b124 = 0
    real(wp),parameter :: b125 = 0
    real(wp),parameter :: b126 = 0
    real(wp),parameter :: b127 = 0
    real(wp),parameter :: b128 = -.1846375997512050141835163881753227910996323204749769226655464078048769505209525299752e-3_wp
    real(wp),parameter :: b129 = .1560894025313219860759149162557283383430181475726228517203663063649626288079337909898_wp
    real(wp),parameter :: b1210 = .1934496857654560252749984220385188727138526287670744309970093278715606577140084022992_wp
    real(wp),parameter :: b1211 = .2611612387636636496908928477536452288263163392010050661129958478089356710938164130987_wp
    real(wp),parameter :: b131 = .4423749328524996327035388417792688154433173133294892285295756457561276315648477233732e-1_wp
    real(wp),parameter :: b132 = 0
    real(wp),parameter :: b133 = 0
    real(wp),parameter :: b134 = 0
    real(wp),parameter :: b135 = 0
    real(wp),parameter :: b136 = 0
    real(wp),parameter :: b137 = 0
    real(wp),parameter :: b138 = .4640774434539039636406222168781981616534115643208114455689698789119941732444857047798e-2_wp
    real(wp),parameter :: b139 = .4704660282615136532130927218172390570903230981414159347904277946537920001824903276586e-1_wp
    real(wp),parameter :: b1310 = .8620749948011488160369445167416002799205317397013619044391270706339561700281526529703e-1_wp
    real(wp),parameter :: b1311 = -.2607983024682138093233254079066687623148682426317395111719299641390118652802949600035e-1_wp
    real(wp),parameter :: b1312 = -.3858020174396621532493277639159499581333235076531298977820093139813399390137768850940e-1_wp
    real(wp),parameter :: b141 = .2318046717429411567006043539613275607940758021709332569729352990777336390158311630529e-1_wp
    real(wp),parameter :: b142 = 0
    real(wp),parameter :: b143 = 0
    real(wp),parameter :: b144 = 0
    real(wp),parameter :: b145 = 0
    real(wp),parameter :: b146 = 0
    real(wp),parameter :: b147 = 0
    real(wp),parameter :: b148 = .3197856784116367067302124322582100058864027838197120089129330601737324659881765852593_wp
    real(wp),parameter :: b149 = .5933233331841898686063939886797828376866051205773280426848164018120869674204443797948_wp
    real(wp),parameter :: b1410 = -1.937519548878479314706815782408229952008442222624773168771865465659822020582450444783_wp
    real(wp),parameter :: b1411 = .1803950557030502357344063195737827904476240180662764468232042537858892203518134072359_wp
    real(wp),parameter :: b1412 = -.4554014298857220726863505256926549022316460712353658688873150702827663762861750674926_wp
    real(wp),parameter :: b1413 = 2.158764106255762807077594619172645539322916635447781333204724468181634037726021280742_wp
    real(wp),parameter :: b151 = .2624364325798105891527733985858552391723553030719144065844544880498188553839263944447e-1_wp
    real(wp),parameter :: b152 = 0
    real(wp),parameter :: b153 = 0
    real(wp),parameter :: b154 = 0
    real(wp),parameter :: b155 = 0
    real(wp),parameter :: b156 = 0
    real(wp),parameter :: b157 = 0
    real(wp),parameter :: b158 = .4863139423867266106526843913609225996253073727381961544415263239431571586043622332760e-1_wp
    real(wp),parameter :: b159 = .4274382538346478867636942429421724367591866585774144180215122660980822123988151132213e-1_wp
    real(wp),parameter :: b1510 = -.4862259869465547771298976981868643277396586803130813159599600102115609499827986711663_wp
    real(wp),parameter :: b1511 = .1326047194917652331781527125743684254490968718259563958293167893998110899691451568372_wp
    real(wp),parameter :: b1512 = -.9402962152946515651634831658142934852383791641671387741034606371378082209616938685225e-1_wp
    real(wp),parameter :: b1513 = .6993864679941022534190304512277131176659196396138275832136258135631963192299339871223_wp
    real(wp),parameter :: b1514 = -.1197020013028860976492784934312243036670658451195397948726104511062042521592125912599e-1_wp
    real(wp),parameter :: b161 = .5568066641536216461090823068917803436066365804361903532125349474551476120813558125830e-1_wp
    real(wp),parameter :: b162 = 0
    real(wp),parameter :: b163 = 0
    real(wp),parameter :: b164 = 0
    real(wp),parameter :: b165 = 0
    real(wp),parameter :: b166 = 0
    real(wp),parameter :: b167 = 0
    real(wp),parameter :: b168 = -.4324853319508358432896036654421685136736530810118924113940744870078036705505610668088_wp
    real(wp),parameter :: b169 = -.9979726994172038714656907882931844552238093285811791155499130927685987422432191170216_wp
    real(wp),parameter :: b1610 = 2.707893755718926115778725270396739994070337972517006747100005607751792006959604868323_wp
    real(wp),parameter :: b1611 = -1.024823023512132929313567156576969954855232272749038347671818195935585095295127839150_wp
    real(wp),parameter :: b1612 = 1.334565206642246959252239602313589265188981560552694580059808406200559397799055652161_wp
    real(wp),parameter :: b1613 = -2.587748998830690939658228913150922979184368065866213469477796089200252812362701917187_wp
    real(wp),parameter :: b1614 = .8992773696348355846430438306111181223414632598285854300924423251352733205187087732678e-1_wp
    real(wp),parameter :: b1615 = 1.497578446211167333777988534023066333042434967475357134513165331964695787890042760189_wp
    real(wp),parameter :: b171 = -.8434891199686377639125188391985671318383858641413517143104162188088468627447515172982e-3_wp
    real(wp),parameter :: b172 = 0
    real(wp),parameter :: b173 = 0
    real(wp),parameter :: b174 = 0
    real(wp),parameter :: b175 = 0
    real(wp),parameter :: b176 = 0
    real(wp),parameter :: b177 = 0
    real(wp),parameter :: b178 = .7602144218856081893754106886111596435015500427480120290148318740899211421773423234728_wp
    real(wp),parameter :: b179 = 1.769083927820959377467464871522349066447068428702073590698445112684989184432409492025_wp
    real(wp),parameter :: b1710 = -4.499239797622297101452915424261016593995695456495268863455643396071539024609271033574_wp
    real(wp),parameter :: b1711 = 1.490558190212043468817221563278239942209691100326719140478588601720867838040211450448_wp
    real(wp),parameter :: b1712 = -2.552203480132132516997563217309689292804518121743365818482497611667126218719069737195_wp
    real(wp),parameter :: b1713 = 4.795167551528575994217413424533259845001657006088189480440731104737960266616292993321_wp
    real(wp),parameter :: b1714 = -.9161854401769482236671414092387917470686251714192236693920061138984202381209109248553e-1_wp
    real(wp),parameter :: b1715 = -1.525735678746850818217653470352135651821164556169070505816135230784807058389577753184_wp
    real(wp),parameter :: b1716 = .7371445601564892133467497107205798584829803038168267854389817508169123996459113657504_wp
    real(wp),parameter :: b181 = .1017366974111576638766809656369828971944080018220332809259398740674738807023371082700_wp
    real(wp),parameter :: b182 = 0
    real(wp),parameter :: b183 = 0
    real(wp),parameter :: b184 = 0
    real(wp),parameter :: b185 = 0
    real(wp),parameter :: b186 = 0
    real(wp),parameter :: b187 = 0
    real(wp),parameter :: b188 = -1.696217553209432810711666838709742166182992092906177246174096517233561845662947862824_wp
    real(wp),parameter :: b189 = -3.825235846211624254528740857512255693551264719132875740261231165548583482101116676418_wp
    real(wp),parameter :: b1810 = 9.754768979885866648856431516333641627109105703674164986615824197909762854575668793816_wp
    real(wp),parameter :: b1811 = -2.520767789227152291196336314591227486393143379933686189126240710041836742414125694941_wp
    real(wp),parameter :: b1812 = 5.472417145227780046950992000565734793413395536531652419585004300790370984185945495978_wp
    real(wp),parameter :: b1813 = -9.781098113458736121002383874108051372067873053264954833376114258940736444388841687929_wp
    real(wp),parameter :: b1814 = .3189152692455334369024560213486753019540464785641163242047782111839399471147176681561_wp
    real(wp),parameter :: b1815 = 3.447227036527756718156475010324322155277035924051392880570525223655410460762027138915_wp
    real(wp),parameter :: b1816 = -.6051983612219277832241707671295607127814820499715293613761402732652780120810041653591_wp
    real(wp),parameter :: b1817 = .3334525350307787459202631378414806560287636505658634784117511174230383993073398823363_wp
    real(wp),parameter :: b191 = -.1012987737478284424676828232882617689682012456457322189102956361570156443805900941944_wp
    real(wp),parameter :: b192 = 0
    real(wp),parameter :: b193 = 0
    real(wp),parameter :: b194 = 0
    real(wp),parameter :: b195 = 0
    real(wp),parameter :: b196 = -.2409389328948775401304659380663043147167897928467308244359962659633933617326533285822e-1_wp
    real(wp),parameter :: b197 = -.6679880790275182076676283582867036095782150170801495251932447614617249253864579543857_wp
    real(wp),parameter :: b198 = 1.600262798493100648047998296908183265688507618079976446601985464092263571149154964705_wp
    real(wp),parameter :: b199 = 3.706958893826695766827011000213884379914407774639901049574259778345288538246990591819_wp
    real(wp),parameter :: b1910 = -8.581755560147929325446798534254342948628755672447282004336563881429983605741487870996_wp
    real(wp),parameter :: b1911 = .5607314974300953986559644699099897253584501767603091982484141468619493221310582281877e-1_wp
    real(wp),parameter :: b1912 = -4.547761497422899514520768375507009011918601407646237921467449197008085790456674001879_wp
    real(wp),parameter :: b1913 = 9.255775439941294621826928846245618922061242300726600002589630404152665447900428712156_wp
    real(wp),parameter :: b1914 = -.3450876657451631707159097079770789925142348071643902737346329921538351794816584861003_wp
    real(wp),parameter :: b1915 = 0
    real(wp),parameter :: b1916 = 0
    real(wp),parameter :: b1917 = 0
    real(wp),parameter :: b1918 = 0
    real(wp),parameter :: b201 = .3826909723812638609001259641818040193828105314579492422836388985468479567237561247336e-1_wp
    real(wp),parameter :: b202 = 0
    real(wp),parameter :: b203 = 0
    real(wp),parameter :: b204 = 0
    real(wp),parameter :: b205 = 0
    real(wp),parameter :: b206 = .7786978965202527814624406274393101840018332461648638653990700950184871893714491273096_wp
    real(wp),parameter :: b207 = .4859454140913448249612202172501868752761599132465501266008866131088163955018926230543_wp
    real(wp),parameter :: b208 = 1.814925350154666364151014269029611427420766367555858499108920245656959783343309816408_wp
    real(wp),parameter :: b209 = 4.551165245704657956889158854062833952834232753889932986749613143631480116805870313264_wp
    real(wp),parameter :: b2010 = -7.173770670344544101351160462586215092596352548535380880420409450623251883641801862305_wp
    real(wp),parameter :: b2011 = -.3943009017000923237232456850787591816773705728833192412204243696911216045268772747196_wp
    real(wp),parameter :: b2012 = -6.036544185898100312430357626685382432626027303329497026597513524312479466987506315664_wp
    real(wp),parameter :: b2013 = 7.338904299721887701527380004651998686389416058019429466200740313593568240326087171554_wp
    real(wp),parameter :: b2014 = -.4143158595971836110248598960027762194900538872022960061452263646470675916118824501965_wp
    real(wp),parameter :: b2015 = 0
    real(wp),parameter :: b2016 = 0
    real(wp),parameter :: b2017 = 0
    real(wp),parameter :: b2018 = 0
    real(wp),parameter :: b2019 = -.3732349451502749258108621577582478607301443393311959731632798508493352335121760204375_wp
    real(wp),parameter :: b211 = .2162339046022045866878628785550588026780578552494608097931198882276791962912244674840e-1_wp
    real(wp),parameter :: b212 = 0
    real(wp),parameter :: b213 = 0
    real(wp),parameter :: b214 = 0
    real(wp),parameter :: b215 = 0
    real(wp),parameter :: b216 = .4611834700744369218866370212060318930941322187829670117414118166503940620998117275429_wp
    real(wp),parameter :: b217 = .1940797759547798743610542713744618433967649025379792966207862125676964319674160574624_wp
    real(wp),parameter :: b218 = .7041001229739959807963554405302474570280838416767002383409508232534658577705201658489_wp
    real(wp),parameter :: b219 = 2.877431096792763528910415905652149398266490601780194388811216042455337979365709745445_wp
    real(wp),parameter :: b2110 = 0
    real(wp),parameter :: b2111 = -.4332742088749107411735902392606181444105337491234912425673655059805456011404518143074_wp
    real(wp),parameter :: b2112 = -2.234178753588834452567105459024473991729105867012210449973082203886376638514123583334_wp
    real(wp),parameter :: b2113 = .2235678086885984010238782832657956960650576194069632574873732156942360146780276407657_wp
    real(wp),parameter :: b2114 = .1293532338308457711442786069651741293532338308457711442786069651741293532338308457711_wp
    real(wp),parameter :: b2115 = 0
    real(wp),parameter :: b2116 = 0
    real(wp),parameter :: b2117 = 0
    real(wp),parameter :: b2118 = 0
    real(wp),parameter :: b2119 = .1418136968194278394808045812385429206355105705182818920178205766092934777719870449624_wp
    real(wp),parameter :: b2120 = -1.085699633131323582531514699802817081967439754938101617737029931360398856861850276906_wp

    real(wp),parameter :: c1 = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1_wp
    real(wp),parameter :: c2 = 0
    real(wp),parameter :: c3 = 0
    real(wp),parameter :: c4 = 0
    real(wp),parameter :: c5 = 0
    real(wp),parameter :: c6 = 0
    real(wp),parameter :: c7 = 0
    real(wp),parameter :: c8 = 0
    real(wp),parameter :: c9 = 0
    real(wp),parameter :: c10 = 0
    real(wp),parameter :: c11 = 0
    real(wp),parameter :: c12 = .1387145942588715882541801312803271702142521598590204181697361204933422401935856968980_wp
    real(wp),parameter :: c13 = .1892374781489234901583064041060123262381623469486258303271944256799821862794952728707_wp
    real(wp),parameter :: c14 = .9461873907446174507915320205300616311908117347431291516359721283999109313974763643533e-1_wp
    real(wp),parameter :: c15 = .2774291885177431765083602625606543404285043197180408363394722409866844803871713937960_wp
    real(wp),parameter :: c16 = .1387145942588715882541801312803271702142521598590204181697361204933422401935856968980_wp
    real(wp),parameter :: c17 = .9461873907446174507915320205300616311908117347431291516359721283999109313974763643533e-1_wp
    real(wp),parameter :: c18 = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1_wp
    real(wp),parameter :: c19 = 0
    real(wp),parameter :: c20 = 0
    real(wp),parameter :: c21 = 0

    real(wp),parameter :: d1 = .3339829895931337572271945815422988633728883413227543303554098429202731077409488318421e-1_wp
    real(wp),parameter :: d2 = 0
    real(wp),parameter :: d3 = 0
    real(wp),parameter :: d4 = 0
    real(wp),parameter :: d5 = 0
    real(wp),parameter :: d6 = 0
    real(wp),parameter :: d7 = 0
    real(wp),parameter :: d8 = 0
    real(wp),parameter :: d9 = .5024509803921568627450980392156862745098039215686274509803921568627450980392156862745e-1_wp
    real(wp),parameter :: d10 = -.1423859191318858946753152353981644782061337055184060977838998119673893661279423564924_wp
    real(wp),parameter :: d11 = .2126013199429258434998789109063801828540550730541648287733608913970804891935883227446_wp
    real(wp),parameter :: d12 = .3254854965632843133622967470840062095221514741629108993207015882688341071771214986692_wp
    real(wp),parameter :: d13 = .3312629399585921325051759834368530020703933747412008281573498964803312629399585921325_wp
    real(wp),parameter :: d14 = .1887845809230650005639203350759631573314744764356665687950807917985316096487827551997_wp
    real(wp),parameter :: d15 = 0
    real(wp),parameter :: d16 = 0
    real(wp),parameter :: d17 = 0
    real(wp),parameter :: d18 = 0
    real(wp),parameter :: d19 = .6159811094287144604404508847679200761569154839698701533406779753675999506388462962070e-1_wp
    real(wp),parameter :: d20 = -.9440109660594088037957791636147830082275023098021053120999763935859315478744850552959e-1_wp
    real(wp),parameter :: d21 = .3341117040855897708234682470384970584684876341854831047975628586614323631403861184369e-1_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e2  = c2  - d2
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10 = c10 - d10
    real(wp),parameter :: e11 = c11 - d11
    real(wp),parameter :: e12 = c12 - d12
    real(wp),parameter :: e13 = c13 - d13
    real(wp),parameter :: e14 = c14 - d14
    real(wp),parameter :: e15 = c15 - d15
    real(wp),parameter :: e16 = c16 - d16
    real(wp),parameter :: e17 = c17 - d17
    real(wp),parameter :: e18 = c18 - d18
    real(wp),parameter :: e19 = c19 - d19
    real(wp),parameter :: e20 = c20 - d20
    real(wp),parameter :: e21 = c21 - d21

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13), &
               f14 => me%funcs(:,14), &
               f15 => me%funcs(:,15), &
               f16 => me%funcs(:,16), &
               f17 => me%funcs(:,17), &
               f18 => me%funcs(:,18), &
               f19 => me%funcs(:,19), &
               f20 => me%funcs(:,20), &
               f21 => me%funcs(:,21))

        call me%f(t+h,x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b82*f2+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b92*f2+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b102*f2+b103*f3+b104*f4+b105*f5+b106*f6+b107*f7+&
                            b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b112*f2+b113*f3+b114*f4+b115*f5+b116*f6+b117*f7+&
                            b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b122*f2+b123*f3+b124*f4+b125*f5+b126*f6+b127*f7+&
                            b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b131*f1+b132*f2+b133*f3+b134*f4+b135*f5+b136*f6+b137*f7+&
                            b138*f8+b139*f9+b1310*f10+b1311*f11+b1312*f12),f13)
        call me%f(t+a14*h,x+h*(b141*f1+b142*f2+b143*f3+b144*f4+b145*f5+b146*f6+b147*f7+&
                            b148*f8+b149*f9+b1410*f10+b1411*f11+b1412*f12+b1413*f13),f14)
        call me%f(t+a15*h,x+h*(b151*f1+b152*f2+b153*f3+b154*f4+b155*f5+b156*f6+b157*f7+&
                            b158*f8+b159*f9+b1510*f10+b1511*f11+b1512*f12+b1513*f13+&
                            b1514*f14),f15)
        call me%f(t+a16*h,x+h*(b161*f1+b162*f2+b163*f3+b164*f4+b165*f5+b166*f6+b167*f7+&
                            b168*f8+b169*f9+b1610*f10+b1611*f11+b1612*f12+b1613*f13+&
                            b1614*f14+b1615*f15),f16)
        call me%f(t+a17*h,x+h*(b171*f1+b172*f2+b173*f3+b174*f4+b175*f5+b176*f6+b177*f7+&
                            b178*f8+b179*f9+b1710*f10+b1711*f11+b1712*f12+b1713*f13+&
                            b1714*f14+b1715*f15+b1716*f16),f17)
        call me%f(t+a18*h,x+h*(b181*f1+b182*f2+b183*f3+b184*f4+b185*f5+b186*f6+b187*f7+&
                            b188*f8+b189*f9+b1810*f10+b1811*f11+b1812*f12+b1813*f13+&
                            b1814*f14+b1815*f15+b1816*f16+b1817*f17),f18)
        call me%f(t+a19*h,x+h*(b191*f1+b192*f2+b193*f3+b194*f4+b195*f5+b196*f6+b197*f7+&
                            b198*f8+b199*f9+b1910*f10+b1911*f11+b1912*f12+b1913*f13+&
                            b1914*f14+b1915*f15+b1916*f16+b1917*f17+b1918*f18),f19)
        call me%f(t+a20*h,x+h*(b201*f1+b202*f2+b203*f3+b204*f4+b205*f5+b206*f6+b207*f7+&
                            b208*f8+b209*f9+b2010*f10+b2011*f11+b2012*f12+b2013*f13+&
                            b2014*f14+b2015*f15+b2016*f16+b2017*f17+b2018*f18+b2019*f19),f20)
        call me%f(t+a21*h,x+h*(b211*f1+b212*f2+b213*f3+b214*f4+b215*f5+b216*f6+b217*f7+&
                            b218*f8+b219*f9+b2110*f10+b2111*f11+b2112*f12+b2113*f13+&
                            b2114*f14+b2115*f15+b2116*f16+b2117*f17+b2118*f18+b2119*&
                            f19+b2120*f20),f21)

        xf = x+h*(c1*f1+c2*f2+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+&
                c11*f11+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16+c17*f17+c18*f18+&
                c19*f19+c20*f20+c21*f21)

        xerr = h*(e1*f1+e2*f2+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+&
                e11*f11+e12*f12+e13*f13+e14*f14+e15*f15+e16*f16+e17*f17+e18*f18+&
                e19*f19+e20*f20+e21*f21)

    end associate

    end procedure rkc108
!*****************************************************************************************

!*****************************************************************************************
!>
!  Baker 10(9) method.
!
!### Reference
!  * Tom Baker, University of Teeside, "RK-Aid software: RK10921M".
!    [Coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10c_1.pdf)

    module procedure rkb109

    real(wp),parameter :: a2  = .2232129192123735665527132860096509572180113541116006004801774831897914255866976731834_wp
    real(wp),parameter :: a3  = .3348193788185603498290699290144764358270170311674009007202662247846871383800465097751_wp
    real(wp),parameter :: a4  = .5022290682278405247436048935217146537405255467511013510803993371770307075700697646627_wp
    real(wp),parameter :: a5  = .1176948756548443596524604744040430749724814058786280779261159420482531581524139111558_wp
    real(wp),parameter :: a6  = .6425923677604462086219057105219705193480713109262431864517703346402040737649802005250_wp
    real(wp),parameter :: a7  = .181826565311210497_wp
    real(wp),parameter :: a8  = .4341610334077954337462856744299909865481150413847250582144329859884319977831972426112_wp
    real(wp),parameter :: a9  = .7122335424182418496009012272374473135415214567547251841147975874613222331084545787886_wp
    real(wp),parameter :: a10 = .1894400309592476109308204841234726559856268372193596062465558450034952083309882121192_wp
    real(wp),parameter :: a11 = .4943921538691045853333333333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a12 = .6403413702953303301268868454436798781085931758963391614454137298592697657577711636312_wp
    real(wp),parameter :: a13 = .741588230803656878_wp
    real(wp),parameter :: a14 = .382911329862520855_wp
    real(wp),parameter :: a15 = .107157755822487322_wp
    real(wp),parameter :: a16 = .875691376241245711_wp
    real(wp),parameter :: a17 = .964069299370187816_wp
    real(wp),parameter :: a18 = .281729610523717922_wp
    real(wp),parameter :: a19 = .631145492177176748_wp
    real(wp),parameter :: a20 = .973803039377034607_wp

    real(wp),parameter :: b21 = .2232129192123735665527132860096509572180113541116006004801774831897914255866976731834_wp
    real(wp),parameter :: b31 = .8370484470464008745726748225361910895675425779185022518006655619617178459501162744378e-1_wp
    real(wp),parameter :: b32 = .2511145341139202623718024467608573268702627733755506755401996685885153537850348823314_wp
    real(wp),parameter :: b41 = .1255572670569601311859012233804286634351313866877753377700998342942576768925174411657_wp
    real(wp),parameter :: b43 = .3766718011708803935577036701412859903053941600633260133102995028827730306775523234970_wp
    real(wp),parameter :: b51 = .8645012631860283043095827588443080342589988168209518203030156172654040053726457363040e-1_wp
    real(wp),parameter :: b53 = .5236243791730284147504959465381950604441622286481307774905930763200738501679225239292e-1_wp
    real(wp),parameter :: b54 = -.2111768858106131225354739613420723449783469866828018185324492731029462740164291486751e-1_wp
    real(wp),parameter :: b61 = -.2639257721366917655020719039008537814312675653905044467083967854693494521584367325145e-1_wp
    real(wp),parameter :: b64 = .3321586951309942525210016931443520748342752809336835485070926725007327816351124496012_wp
    real(wp),parameter :: b65 = .3368262498431211326511112077677038226569227865316100826155173406864062373457114241753_wp
    real(wp),parameter :: b71 = .4214481430128891924590207161639307925193144098983055869572587559196013842334308634337e-1_wp
    real(wp),parameter :: b75 = .1395091008993089985837339976891755843541033762656629420011519949721775047502000865276_wp
    real(wp),parameter :: b76 = .1726501106125791703639306944313363939651827445064993031221294358623568264568271289984e-3_wp
    real(wp),parameter :: b81 = .1913223517063317850561878098464623128757461736096094048517648270556161690028960671678_wp
    real(wp),parameter :: b85 = -.6549476653170516368338607340883454058106339891146285998978549868910881698171713079544_wp
    real(wp),parameter :: b86 = .1755876323742502725238471568246235150972597859217101866902108545773998108136024750683e-1_wp
    real(wp),parameter :: b87 = .8802275837810902582715738829894117279732768782975732345915020603661640175161122358910_wp
    real(wp),parameter :: b91 = .5819519506251173970176101222224200992345451475572200041442763494089781743712740682029e-1_wp
    real(wp),parameter :: b96 = .1647908558976784212367930935320516004070340774715814512653291320117464889726338996259_wp
    real(wp),parameter :: b97 = .2562758618556534830885936612274458691344630164445656645146446091037555911278059725788_wp
    real(wp),parameter :: b98 = .2329716296023982055737534602557078340765698480828560679203962114049223355708872997637_wp
    real(wp),parameter :: b101 = .6609808969751805186432114363691927852621669632797821827493854746683225501745446560252e-1_wp
    real(wp),parameter :: b106 = .4993798303607936315852000088060192246997357096577012552548527217200430403080630524414e-1_wp
    real(wp),parameter :: b107 = .1536216094454890707317822140262379594167730891287251273659790090180235843883503581847_wp
    real(wp),parameter :: b108 = -.5414130381994316899529893137177555939937633423733000345305400547297139497578508331128e-1_wp
    real(wp),parameter :: b109 = -.2607634739989570582850394304851094502796018496578386146679297818039354012983783360089e-1_wp
    real(wp),parameter :: b111 = .5977281413067455622264838600535046960155059547948827360548395666199015131925709118733e-1_wp
    real(wp),parameter :: b118 = .1656649460691327819663527933964608353435263909407310250007028968559622442637053740488_wp
    real(wp),parameter :: b119 = -.1270467077219322361480177163110966635379558631317440513037999854071349138024611223099e-2_wp
    real(wp),parameter :: b1110 = .2702248607465165695058123310946329950236359055444314752401844796694522868883954793203_wp
    real(wp),parameter :: b121 = .5714646346927726543617292683289117330020257083064077657640682066128476784510949443172e-1_wp
    real(wp),parameter :: b129 = .1706647333279140031893907789070497592665517669717541507696626200796144815409471923569e-1_wp
    real(wp),parameter :: b1210 = .2853767151165243146804641383356191775396044057521502467001987179873305495674085297549_wp
    real(wp),parameter :: b1211 = .2807517183767373496913107023844645513421310226163727230918419292026930001911584202089_wp
    real(wp),parameter :: b131 = .5713690770667411245888754904665701448619979494564888202366599688508571201235483609660e-1_wp
    real(wp),parameter :: b139 = .8662146075905090149738288126241361791710263955784173047391548833627524888649185757901e-1_wp
    real(wp),parameter :: b1310 = .2854329582389984323634810972884390808440109342965556282630983683305184760043391707767_wp
    real(wp),parameter :: b1311 = .2797633923959254653770985744082321055561046971559230633717557229339774187704436927395_wp
    real(wp),parameter :: b1312 = .3263351170300796630314989799425818119658193404403069586756442351414314432637044280815e-1_wp
    real(wp),parameter :: b141 = .5792154217958510846632846508298421800727675597066414244319766319648769731268367294802e-1_wp
    real(wp),parameter :: b149 = -.1168598764394854203108670771355968386139754171212251879659829443523204385801129100146_wp
    real(wp),parameter :: b1410 = .2791270024477207066381851603610170508129292534731657465038373367821947287200953472810_wp
    real(wp),parameter :: b1411 = .5523482312158897030618178268050450731680950103624425704577340542168627357243106364454e-1_wp
    real(wp),parameter :: b1412 = .3313821420230448163586175200236513161592367149984143944502756986720155621978125806284e-1_wp
    real(wp),parameter :: b1413 = .7434962435080700826430991700872593086103623514130960252814696908475018275512156807813e-1_wp
    real(wp),parameter :: b151 = .5424349131266048906033564077410202025889503282676783408313456162421294081865905389248e-1_wp
    real(wp),parameter :: b159 = .2555744036484019017045493997261865838594175905367142496400185786785191604031397243218e-1_wp
    real(wp),parameter :: b1510 = .1065318948685044129706483445621162791089602467505369099004177225283624953902959013228_wp
    real(wp),parameter :: b1511 = .1494739176822617600298548318084537294306548147599349812604239821749187453609545343658_wp
    real(wp),parameter :: b1512 = -.7541418150821686840498891896239704084250667222662001839556328840119093341810397465136e-1_wp
    real(wp),parameter :: b1513 = .1183662716974706173695161845106353658054818835708868187585164205844835807880512638150e-2_wp
    real(wp),parameter :: b1514 = -.154418469614537368_wp
    real(wp),parameter :: b161 = -.3130073355625165399871063999971367682939799979170688711523373232868473829204122305445e-1_wp
    real(wp),parameter :: b169 = -.4393041016682237597697452523147415408933247130414835042917723489831082747083637622576_wp
    real(wp),parameter :: b1610 = -.2971894957150042139937620479069507765251785957115322181603409909075380714667735160813_wp
    real(wp),parameter :: b1611 = .1262980195755398451585609791771777378208362550125140868265029793320289593624402155837_wp
    real(wp),parameter :: b1612 = -.1352235663389501393036813956309918309028161801076417057198393439731017100387729485379e-1_wp
    real(wp),parameter :: b1613 = .6695605432945227295340251006073274395173466715429726933128280272846122961086155806635_wp
    real(wp),parameter :: b1614 = .398926850454282556_wp
    real(wp),parameter :: b1615 = .462222650490275222_wp
    real(wp),parameter :: b171 = .4254608456033344945622257854593736284131499735662276010271273463216872625702497770866e-1_wp
    real(wp),parameter :: b179 = -.1821182296512653493976221650402057582117417514147759470622260585285199838506585092559_wp
    real(wp),parameter :: b1710 = .1693488568788993457692776359468516487310835997985566768493328170663073221084671579068_wp
    real(wp),parameter :: b1711 = -.2812932311380423097505532359005039285971305910469541561543126052387690343228512975388e-2_wp
    real(wp),parameter :: b1712 = .4367359500435692592480947784835035112081291953371300785838045346974343617081011595250_wp
    real(wp),parameter :: b1713 = .9649623253501190961728232014648167901724237884864786722867265364689659268793239634042e-2_wp
    real(wp),parameter :: b1714 = .1827541189026459104590716992261562936435133038172843022420870171219768684076842803782_wp
    real(wp),parameter :: b1715 = .8163992629788446963341645082372801271393454416361108721409416282456258720143539344240e-1_wp
    real(wp),parameter :: b1716 = .2263259013959999629673163223583858004580131789671757969088706528737681492423808136362_wp
    real(wp),parameter :: b181 = -.6162786199706763955052568808743247751165729542444646544961722501628952288996567177280e-1_wp
    real(wp),parameter :: b189 = .2630819816424975496840706912577088603030984130761613342651746722632362416592729338065e-1_wp
    real(wp),parameter :: b1810 = -.6563819353303824389148917307603424726525093364227772615481314402265589817962663323807_wp
    real(wp),parameter :: b1811 = -.4996969738639946601661934398570589720521936732243676171454742217015940903119194293034_wp
    real(wp),parameter :: b1812 = -.7123015625596798250364056005999154454586809433543859495632863935368185749281956104361e-1_wp
    real(wp),parameter :: b1813 = .1945887064918253295880838290230513152388881501421435771759995471289411971821731631008_wp
    real(wp),parameter :: b1814 = .7502438829684451516052952921664213906375154664124515160998662893556068266252405773932_wp
    real(wp),parameter :: b1815 = .6778525665536439867432200838066450041417881571145524785032579156045245492715999778952_wp
    real(wp),parameter :: b1816 = -.1046057647455310786652080097009039359249473747848876107492887123406055065619567484248_wp
    real(wp),parameter :: b1817 = .2627894853849749889545315434384080663867415921515384464319901932333376180798673115539e-1_wp
    real(wp),parameter :: b191 = .3796521632874261201186203436841209430559695655768644122346210733272382403437031805294e-1_wp
    real(wp),parameter :: b199 = -.2594056416058702079475722192315214454017994671360402614116801042085357278273610676615_wp
    real(wp),parameter :: b1910 = .6023088303591701879982724849459819049819908747310071139535368513067871259409896536125_wp
    real(wp),parameter :: b1911 = .6119387226913542078978281038711879095845432514371708320113689988639552031445275756169_wp
    real(wp),parameter :: b1912 = -.8557593281472216275173535332798856826091404689475696899046474643882545069592793853088_wp
    real(wp),parameter :: b1913 = 1.016014044479150879819181655935044451622873646716581194764492422092270586935018600847_wp
    real(wp),parameter :: b1914 = .4174635465318728144170599907012271233953357864579855698710555478932012900105054820193_wp
    real(wp),parameter :: b1915 = -.3142068748294011550388298801522627398568106507178091710884531834553635025761917127745e-2_wp
    real(wp),parameter :: b1916 = -.3113859536805386540643212859648902182769708893494791481799285293081168071417215791223_wp
    real(wp),parameter :: b1917 = .7475323182002611486490630088842356498085791607202230555314272159767158496743797941577e-1_wp
    real(wp),parameter :: b1918 = -.6996051078512155679294752334324570751847195000321862661699180193471489380787256603442_wp
    real(wp),parameter :: b201 = .6354247858004162057833209973904648790423789887640110847597201944718467473335513915457e-1_wp
    real(wp),parameter :: b209 = .1335709427214855386564319367816597150656631989678512452281616080142611136536750546091_wp
    real(wp),parameter :: b2010 = .3491939710685402856556370077059343902437022368395394555525545826876029022896488736825_wp
    real(wp),parameter :: b2011 = -.1284299219085288104043584513269729325317465049017238832115899338924554224999281557779_wp
    real(wp),parameter :: b2012 = .5820143454609020234416665664957204705761540086346612646427774305741185500271652747165_wp
    real(wp),parameter :: b2013 = -.7032295749061967719609757517752619736970110102217160910226632819061710474763494562046e-1_wp
    real(wp),parameter :: b2014 = .4108572358925710179476579003323188290642373757595878674150744283920916547885055706399_wp
    real(wp),parameter :: b2015 = -.2510654779239962466269627907503763505655390551209676232633905831376503530736867220512e-1_wp
    real(wp),parameter :: b2016 = .1223166254706007785807460575472930729088496060945489634909944334258869465746143128690_wp
    real(wp),parameter :: b2017 = .4712374572639742166850752217533416565767093458666637954607896698933218387654817641078e-1_wp
    real(wp),parameter :: b2018 = -.2035773172368828280351383367281992122282011746329188150975183458730713499142844139626_wp
    real(wp),parameter :: b2019 = -.3073795611150731392306884484695711542343125736903452146138998032605691134742962145162_wp
    real(wp),parameter :: b211 = .1189924137880676460852663367282125023346761041034714988771622602337896848257576311721_wp
    real(wp),parameter :: b219 = 1.567065389218647018255284948985537561396036403483411695550932448674279128530361565244_wp
    real(wp),parameter :: b2110 = 1.062742213102872753657741254837796588064400661866785937264672961506012193460062190291_wp
    real(wp),parameter :: b2111 = -.7051696022721412819771965938099389611228930171636494130149946862206926743626267591091_wp
    real(wp),parameter :: b2112 = .9907825344767892999287937933800084396119601938859814415204705439463102044422773242413_wp
    real(wp),parameter :: b2113 = .2507980962931107922266871779870423357231706420300879746540980799200371363012966217081_wp
    real(wp),parameter :: b2114 = 1.972591757153629542226038011120263810650524345125694323573613060238829589123920677616_wp
    real(wp),parameter :: b2115 = -.3262628387941532893880299499584174546880007304318604226454131079685364715701264155684_wp
    real(wp),parameter :: b2116 = -.7195961905900189101450624277119790470362703727002334732811111315403329989928142118605_wp
    real(wp),parameter :: b2117 = .5802031535291661528133435810283051023192180079596550160030993013203826047471987749923_wp
    real(wp),parameter :: b2118 = -1.503104563178041688806298917914932279967077722444117887454135595533548208509663052219_wp
    real(wp),parameter :: b2119 = -2.082917406986733532241895178393857460921252779279252368270349439324592549474879051363_wp
    real(wp),parameter :: b2120 = -.2061249557411945026346720362780411363644917364359743227780446952519376385207652951444_wp

    real(wp),parameter :: c1 = .3074440935793207675883517938984149786411995438646211366049687901929864736726360396701e-1_wp
    real(wp),parameter :: c12 = .4447886124335158463373029422601699728444313324875120673928896655863659883998986985001_wp
    real(wp),parameter :: c13 = -.7714472898672551027236427412075689933241937304249667830486089002567926592070578328099e-1_wp
    real(wp),parameter :: c14 = .1895969236666960081782143661904924178003392273539755384045555627359634132139704039829_wp
    real(wp),parameter :: c15 = .1696794742660741325357704245298007854728500376148224262050684467711119689665579207267_wp
    real(wp),parameter :: c16 = .3107725327623489988155072749206969138501610275590999188855016092183529050525521871700_wp
    real(wp),parameter :: c17 = -1.052173416444910172223599400939351222118848154980691993896969290833000795415995203779_wp
    real(wp),parameter :: c18 = .1098353169596310609348506603862869818226618225743191040640071521542528953059954349550_wp
    real(wp),parameter :: c19 = -.1334572232187388442842841059278351057902328711498735111118805829160094435254500793407_wp
    real(wp),parameter :: c20 = 1.154026328574134208065260150796827370454610532524502151460002781580501491609214462821_wp
    real(wp),parameter :: c21 = -.1466682293699578048454932174861727128676735353276311367588113332911578050533016457222_wp

    real(wp),parameter :: d1 = .3069935889151056339784198488152917237594701959064431113839566366053496143740989806580e-1_wp
    real(wp),parameter :: d12 = .4501037370765134863000400416812391552231563060566733706571606677666390164501322407860_wp
    real(wp),parameter :: d13 = -.6514757933632400121177566572810133388907129562810217527103583351082868427677317818171e-1_wp
    real(wp),parameter :: d14 = .1916476309177840038043731680068457542540807928558914659725060136615214781805707317846_wp
    real(wp),parameter :: d15 = .1699093037015861786030207592719576231809404332930270025745281206261397945871067257969_wp
    real(wp),parameter :: d16 = .2831822563446878149008449047225439184357774091336334653238641676396106848038134341285_wp
    real(wp),parameter :: d17 = -.7938273328067640993878535581971311480056976834520424653124835645605083935506842817243_wp
    real(wp),parameter :: d18 = .1085805168491408111596918257428235537727685528809234551377846603272747280143198404809_wp
    real(wp),parameter :: d19 = -.1442902422009048623094173280045594665841429401482516755746485620108546739010869727583_wp
    real(wp),parameter :: d20 = .8718101111217405681350791198631736702436128801469513492598526664004710882551915616216_wp
    real(wp),parameter :: d21 = -.102667760558970463391845252240320899007371474729348103905924_wp

    real(wp),parameter :: e1  = c1   - d1
    real(wp),parameter :: e12 = c12  - d12
    real(wp),parameter :: e13 = c13  - d13
    real(wp),parameter :: e14 = c14  - d14
    real(wp),parameter :: e15 = c15  - d15
    real(wp),parameter :: e16 = c16  - d16
    real(wp),parameter :: e17 = c17  - d17
    real(wp),parameter :: e18 = c18  - d18
    real(wp),parameter :: e19 = c19  - d19
    real(wp),parameter :: e20 = c20  - d20
    real(wp),parameter :: e21 = c21  - d21

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13), &
               f14 => me%funcs(:,14), &
               f15 => me%funcs(:,15), &
               f16 => me%funcs(:,16), &
               f17 => me%funcs(:,17), &
               f18 => me%funcs(:,18), &
               f19 => me%funcs(:,19), &
               f20 => me%funcs(:,20), &
               f21 => me%funcs(:,21))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b106*f6+b107*f7+b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b131*f1+b139*f9+b1310*f10+b1311*f11+b1312*f12),f13)
        call me%f(t+a14*h,x+h*(b141*f1+b149*f9+b1410*f10+b1411*f11+b1412*f12+b1413*f13),f14)
        call me%f(t+a15*h,x+h*(b151*f1+b159*f9+b1510*f10+b1511*f11+b1512*f12+b1513*f13+&
                            b1514*f14),f15)
        call me%f(t+a16*h,x+h*(b161*f1+b169*f9+b1610*f10+b1611*f11+b1612*f12+b1613*f13+&
                            b1614*f14+b1615*f15),f16)
        call me%f(t+a17*h,x+h*(b171*f1+b179*f9+b1710*f10+b1711*f11+b1712*f12+b1713*f13+&
                            b1714*f14+b1715*f15+b1716*f16),f17)
        call me%f(t+a18*h,x+h*(b181*f1+b189*f9+b1810*f10+b1811*f11+b1812*f12+b1813*f13+&
                            b1814*f14+b1815*f15+b1816*f16+b1817*f17),f18)
        call me%f(t+a19*h,x+h*(b191*f1+b199*f9+b1910*f10+b1911*f11+b1912*f12+b1913*f13+&
                            b1914*f14+b1915*f15+b1916*f16+b1917*f17+b1918*f18),f19)
        call me%f(t+a20*h,x+h*(b201*f1+b209*f9+b2010*f10+b2011*f11+b2012*f12+b2013*f13+&
                            b2014*f14+b2015*f15+b2016*f16+b2017*f17+b2018*f18+b2019*f19),f20)
        call me%f(t+h,    x+h*(b211*f1+b219*f9+b2110*f10+b2111*f11+b2112*f12+b2113*f13+&
                            b2114*f14+b2115*f15+b2116*f16+b2117*f17+b2118*f18+b2119*f19+&
                            b2120*f20),f21)

        xf = x+h*(c1*f1+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16+&
                c17*f17+c18*f18+c19*f19+c20*f20+c21*f21)

        xerr = h*(e1*f1+e12*f12+e13*f13+e14*f14+e15*f15+e16*f16+&
                e17*f17+e18*f18+e19*f19+e20*f20+e21*f21)

    end associate

    end procedure rkb109
!*****************************************************************************************

!*****************************************************************************************
!>
!  Feagin's RK12(10) method -- a 12th-order method with an embedded 10th-order method.
!
!### Reference
!  * [coefficient file](http://sce.uhcl.edu/rungekutta/rk1210.txt)

    module procedure rkf1210

    real(wp),parameter :: a0  = 0
    real(wp),parameter :: a1  = 0.200000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: a2  = 0.555555555555555555555555555555555555555555555555555555555556_wp
    real(wp),parameter :: a3  = 0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a4  = 0.333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a6  = 0.671835709170513812712245661002797570438953420568682550710222_wp
    real(wp),parameter :: a7  = 0.288724941110620201935458488967024976908118598341806976469674_wp
    real(wp),parameter :: a8  = 0.562500000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: a9  = 0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a10 = 0.947695431179199287562380162101836721649589325892740646458322_wp
    real(wp),parameter :: a11 = 0.0548112876863802643887753674810754475842153612931128785028369_wp
    real(wp),parameter :: a12 = 0.0848880518607165350639838930162674302064148175640019542045934_wp
    real(wp),parameter :: a13 = 0.265575603264642893098114059045616835297201264164077621448665_wp
    real(wp),parameter :: a14 = 0.500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: a15 = 0.734424396735357106901885940954383164702798735835922378551335_wp
    real(wp),parameter :: a16 = 0.915111948139283464936016106983732569793585182435998045795407_wp
    real(wp),parameter :: a17 = 0.947695431179199287562380162101836721649589325892740646458322_wp
    real(wp),parameter :: a18 = 0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a19 = 0.288724941110620201935458488967024976908118598341806976469674_wp
    real(wp),parameter :: a20 = 0.671835709170513812712245661002797570438953420568682550710222_wp
    real(wp),parameter :: a21 = 0.333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a22 = 0.555555555555555555555555555555555555555555555555555555555556_wp
    real(wp),parameter :: a23 = 0.200000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: a24 = 1

    real(wp),parameter :: c0  = 0.0238095238095238095238095238095238095238095238095238095238095_wp
    real(wp),parameter :: c1  = 0.0234375000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c2  = 0.0312500000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c3  = 0
    real(wp),parameter :: c4  = 0.0416666666666666666666666666666666666666666666666666666666667_wp
    real(wp),parameter :: c5  = 0
    real(wp),parameter :: c6  = 0.0500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c7  = 0.0500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c8  = 0
    real(wp),parameter :: c9  = 0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c10 =   0.0714285714285714285714285714285714285714285714285714285714286_wp
    real(wp),parameter :: c11 =   0
    real(wp),parameter :: c12 =   0.138413023680782974005350203145033146748813640089941234591267_wp
    real(wp),parameter :: c13 =   0.215872690604931311708935511140681138965472074195773051123019_wp
    real(wp),parameter :: c14 =   0.243809523809523809523809523809523809523809523809523809523810_wp
    real(wp),parameter :: c15 =   0.215872690604931311708935511140681138965472074195773051123019_wp
    real(wp),parameter :: c16 =   0.138413023680782974005350203145033146748813640089941234591267_wp
    real(wp),parameter :: c17 =  -0.0714285714285714285714285714285714285714285714285714285714286_wp
    real(wp),parameter :: c18 =  -0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c19 =  -0.0500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c20 =  -0.0500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c21 =  -0.0416666666666666666666666666666666666666666666666666666666667_wp
    real(wp),parameter :: c22 =  -0.0312500000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c23 =  -0.0234375000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c24 =   0.0238095238095238095238095238095238095238095238095238095238095_wp

    real(wp),parameter :: b10   =    0.200000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b20   =   -0.216049382716049382716049382716049382716049382716049382716049_wp
    real(wp),parameter :: b21   =    0.771604938271604938271604938271604938271604938271604938271605_wp
    real(wp),parameter :: b30   =    0.208333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: b31   =    0
    real(wp),parameter :: b32   =    0.625000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b40   =    0.193333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: b41   =    0
    real(wp),parameter :: b42   =    0.220000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b43   =   -0.0800000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b50   =    0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b51   =    0
    real(wp),parameter :: b52   =    0
    real(wp),parameter :: b53   =    0.400000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b54   =    0.500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b60   =    0.103364471650010477570395435690481791543342708330349879244197_wp
    real(wp),parameter :: b61   =    0
    real(wp),parameter :: b62   =    0
    real(wp),parameter :: b63   =    0.124053094528946761061581889237115328211074784955180298044074_wp
    real(wp),parameter :: b64   =    0.483171167561032899288836480451962508724109257517289177302380_wp
    real(wp),parameter :: b65   =   -0.0387530245694763252085681443767620580395733302341368038804290_wp
    real(wp),parameter :: b70   =    0.124038261431833324081904585980175168140024670698633612292480_wp
    real(wp),parameter :: b71   =    0
    real(wp),parameter :: b72   =    0
    real(wp),parameter :: b73   =    0
    real(wp),parameter :: b74   =    0.217050632197958486317846256953159942875916353757734167684657_wp
    real(wp),parameter :: b75   =    0.0137455792075966759812907801835048190594443990939408530842918_wp
    real(wp),parameter :: b76   =   -0.0661095317267682844455831341498149531672668252085016565917546_wp
    real(wp),parameter :: b80   =    0.0914774894856882983144991846980432197088832099976660100090486_wp
    real(wp),parameter :: b81   =    0
    real(wp),parameter :: b82   =    0
    real(wp),parameter :: b83   =    0
    real(wp),parameter :: b84   =    0
    real(wp),parameter :: b85   =   -0.00544348523717469689965754944144838611346156873847009178068318_wp
    real(wp),parameter :: b86   =    0.0680716801688453518578515120895103863112751730758794372203952_wp
    real(wp),parameter :: b87   =    0.408394315582641046727306852653894780093303185664924644551239_wp
    real(wp),parameter :: b90   =    0.0890013652502551018954509355423841780143232697403434118692699_wp
    real(wp),parameter :: b91   =    0
    real(wp),parameter :: b92   =    0
    real(wp),parameter :: b93   =    0
    real(wp),parameter :: b94   =    0
    real(wp),parameter :: b95   =    0.00499528226645532360197793408420692800405891149406814091955810_wp
    real(wp),parameter :: b96   =    0.397918238819828997341739603001347156083435060931424970826304_wp
    real(wp),parameter :: b97   =    0.427930210752576611068192608300897981558240730580396406312359_wp
    real(wp),parameter :: b98   =   -0.0865117637557827005740277475955029103267246394128995965941585_wp
    real(wp),parameter :: b100  =    0.0695087624134907543112693906409809822706021061685544615255758_wp
    real(wp),parameter :: b101  =    0
    real(wp),parameter :: b102  =    0
    real(wp),parameter :: b103  =    0
    real(wp),parameter :: b104  =    0
    real(wp),parameter :: b105  =    0.129146941900176461970759579482746551122871751501482634045487_wp
    real(wp),parameter :: b106  =    1.53073638102311295076342566143214939031177504112433874313011_wp
    real(wp),parameter :: b107  =    0.577874761129140052546751349454576715334892100418571882718036_wp
    real(wp),parameter :: b108  =   -0.951294772321088980532340837388859453930924498799228648050949_wp
    real(wp),parameter :: b109  =   -0.408276642965631951497484981519757463459627174520978426909934_wp
    real(wp),parameter :: b110  =    0.0444861403295135866269453507092463581620165501018684152933313_wp
    real(wp),parameter :: b111  =    0
    real(wp),parameter :: b112  =    0
    real(wp),parameter :: b113  =    0
    real(wp),parameter :: b114  =    0
    real(wp),parameter :: b115  =   -0.00380476867056961731984232686574547203016331563626856065717964_wp
    real(wp),parameter :: b116  =    0.0106955064029624200721262602809059154469206077644957399593972_wp
    real(wp),parameter :: b117  =    0.0209616244499904333296674205928919920806734650660039898074652_wp
    real(wp),parameter :: b118  =   -0.0233146023259321786648561431551978077665337818756053603898847_wp
    real(wp),parameter :: b119  =    0.00263265981064536974369934736325334761174975280887405725010964_wp
    real(wp),parameter :: b1110 =    0.00315472768977025060103545855572111407955208306374459723959783_wp
    real(wp),parameter :: b120  =    0.0194588815119755475588801096525317761242073762016273186231215_wp
    real(wp),parameter :: b121  =    0
    real(wp),parameter :: b122  =    0
    real(wp),parameter :: b123  =    0
    real(wp),parameter :: b124  =    0
    real(wp),parameter :: b125  =    0
    real(wp),parameter :: b126  =    0
    real(wp),parameter :: b127  =    0
    real(wp),parameter :: b128  =    0.0000678512949171812509306121653452367476194364781259165332321534_wp
    real(wp),parameter :: b129  =   -0.0000429795859049273623271005330230162343568863387724883603675550_wp
    real(wp),parameter :: b1210 =    0.0000176358982260285155407485928953302139937553442829975734148981_wp
    real(wp),parameter :: b1211 =    0.0653866627415027051009595231385181033549511358787382098351924_wp
    real(wp),parameter :: b130  =    0.206836835664277105916828174798272361078909196043446411598231_wp
    real(wp),parameter :: b131  =    0
    real(wp),parameter :: b132  =    0
    real(wp),parameter :: b133  =    0
    real(wp),parameter :: b134  =    0
    real(wp),parameter :: b135  =    0
    real(wp),parameter :: b136  =    0
    real(wp),parameter :: b137  =    0
    real(wp),parameter :: b138  =    0.0166796067104156472828045866664696450306326505094792505215514_wp
    real(wp),parameter :: b139  =   -0.00879501563200710214457024178249986591130234990219959208704979_wp
    real(wp),parameter :: b1310 =    0.00346675455362463910824462315246379209427513654098596403637231_wp
    real(wp),parameter :: b1311 =   -0.861264460105717678161432562258351242030270498966891201799225_wp
    real(wp),parameter :: b1312 =    0.908651882074050281096239478469262145034957129939256789178785_wp
    real(wp),parameter :: b140  =    0.0203926084654484010091511314676925686038504449562413004562382_wp
    real(wp),parameter :: b141  =    0
    real(wp),parameter :: b142  =    0
    real(wp),parameter :: b143  =    0
    real(wp),parameter :: b144  =    0
    real(wp),parameter :: b145  =    0
    real(wp),parameter :: b146  =    0
    real(wp),parameter :: b147  =    0
    real(wp),parameter :: b148  =    0.0869469392016685948675400555583947505833954460930940959577347_wp
    real(wp),parameter :: b149  =   -0.0191649630410149842286436611791405053287170076602337673587681_wp
    real(wp),parameter :: b1410 =    0.00655629159493663287364871573244244516034828755253746024098838_wp
    real(wp),parameter :: b1411 =    0.0987476128127434780903798528674033899738924968006632201445462_wp
    real(wp),parameter :: b1412 =    0.00535364695524996055083260173615567408717110247274021056118319_wp
    real(wp),parameter :: b1413 =    0.301167864010967916837091303817051676920059229784957479998077_wp
    real(wp),parameter :: b150  =    0.228410433917778099547115412893004398779136994596948545722283_wp
    real(wp),parameter :: b151  =    0
    real(wp),parameter :: b152  =    0
    real(wp),parameter :: b153  =    0
    real(wp),parameter :: b154  =    0
    real(wp),parameter :: b155  =    0
    real(wp),parameter :: b156  =    0
    real(wp),parameter :: b157  =    0
    real(wp),parameter :: b158  =   -0.498707400793025250635016567442511512138603770959682292383042_wp
    real(wp),parameter :: b159  =    0.134841168335724478552596703792570104791700727205981058201689_wp
    real(wp),parameter :: b1510 =   -0.0387458244055834158439904226924029230935161059142806805674360_wp
    real(wp),parameter :: b1511 =   -1.27473257473474844240388430824908952380979292713250350199641_wp
    real(wp),parameter :: b1512 =    1.43916364462877165201184452437038081875299303577911839630524_wp
    real(wp),parameter :: b1513 =   -0.214007467967990254219503540827349569639028092344812795499026_wp
    real(wp),parameter :: b1514 =    0.958202417754430239892724139109781371059908874605153648768037_wp
    real(wp),parameter :: b160  =    2.00222477655974203614249646012506747121440306225711721209798_wp
    real(wp),parameter :: b161  =    0
    real(wp),parameter :: b162  =    0
    real(wp),parameter :: b163  =    0
    real(wp),parameter :: b164  =    0
    real(wp),parameter :: b165  =    0
    real(wp),parameter :: b166  =    0
    real(wp),parameter :: b167  =    0
    real(wp),parameter :: b168  =    2.06701809961524912091954656438138595825411859673341600679555_wp
    real(wp),parameter :: b169  =    0.623978136086139541957471279831494466155292316167021080663140_wp
    real(wp),parameter :: b1610 =   -0.0462283685500311430283203554129062069391947101880112723185773_wp
    real(wp),parameter :: b1611 =   -8.84973288362649614860075246727118949286604835457092701094630_wp
    real(wp),parameter :: b1612 =    7.74257707850855976227437225791835589560188590785037197433615_wp
    real(wp),parameter :: b1613 =   -0.588358519250869210993353314127711745644125882130941202896436_wp
    real(wp),parameter :: b1614 =   -1.10683733362380649395704708016953056176195769617014899442903_wp
    real(wp),parameter :: b1615 =   -0.929529037579203999778397238291233214220788057511899747507074_wp
    real(wp),parameter :: b170  =    3.13789533412073442934451608989888796808161259330322100268310_wp
    real(wp),parameter :: b171  =    0
    real(wp),parameter :: b172  =    0
    real(wp),parameter :: b173  =    0
    real(wp),parameter :: b174  =    0
    real(wp),parameter :: b175  =    0.129146941900176461970759579482746551122871751501482634045487_wp
    real(wp),parameter :: b176  =    1.53073638102311295076342566143214939031177504112433874313011_wp
    real(wp),parameter :: b177  =    0.577874761129140052546751349454576715334892100418571882718036_wp
    real(wp),parameter :: b178  =    5.42088263055126683050056840891857421941300558851862156403363_wp
    real(wp),parameter :: b179  =    0.231546926034829304872663800877643660904880180835945693836936_wp
    real(wp),parameter :: b1710 =    0.0759292995578913560162301311785251873561801342333194895292058_wp
    real(wp),parameter :: b1711 =  -12.3729973380186513287414553402595806591349822617535905976253_wp
    real(wp),parameter :: b1712 =    9.85455883464769543935957209317369202080367765721777101906955_wp
    real(wp),parameter :: b1713 =    0.0859111431370436529579357709052367772889980495122329601159540_wp
    real(wp),parameter :: b1714 =   -5.65242752862643921117182090081762761180392602644189218673969_wp
    real(wp),parameter :: b1715 =   -1.94300935242819610883833776782364287728724899124166920477873_wp
    real(wp),parameter :: b1716 =   -0.128352601849404542018428714319344620742146491335612353559923_wp
    real(wp),parameter :: b180  =    1.38360054432196014878538118298167716825163268489922519995564_wp
    real(wp),parameter :: b181  =    0
    real(wp),parameter :: b182  =    0
    real(wp),parameter :: b183  =    0
    real(wp),parameter :: b184  =    0
    real(wp),parameter :: b185  =    0.00499528226645532360197793408420692800405891149406814091955810_wp
    real(wp),parameter :: b186  =    0.397918238819828997341739603001347156083435060931424970826304_wp
    real(wp),parameter :: b187  =    0.427930210752576611068192608300897981558240730580396406312359_wp
    real(wp),parameter :: b188  =   -1.30299107424475770916551439123047573342071475998399645982146_wp
    real(wp),parameter :: b189  =    0.661292278669377029097112528107513072734573412294008071500699_wp
    real(wp),parameter :: b1810 =   -0.144559774306954349765969393688703463900585822441545655530145_wp
    real(wp),parameter :: b1811 =   -6.96576034731798203467853867461083919356792248105919255460819_wp
    real(wp),parameter :: b1812 =    6.65808543235991748353408295542210450632193197576935120716437_wp
    real(wp),parameter :: b1813 =   -1.66997375108841486404695805725510845049807969199236227575796_wp
    real(wp),parameter :: b1814 =    2.06413702318035263832289040301832647130604651223986452170089_wp
    real(wp),parameter :: b1815 =   -0.674743962644306471862958129570837723192079875998405058648892_wp
    real(wp),parameter :: b1816 =   -0.00115618834794939500490703608435907610059605754935305582045729_wp
    real(wp),parameter :: b1817 =   -0.00544057908677007389319819914241631024660726585015012485938593_wp
    real(wp),parameter :: b190  =    0.951236297048287669474637975894973552166903378983475425758226_wp
    real(wp),parameter :: b191  =    0
    real(wp),parameter :: b192  =    0
    real(wp),parameter :: b193  =    0
    real(wp),parameter :: b194  =    0.217050632197958486317846256953159942875916353757734167684657_wp
    real(wp),parameter :: b195  =    0.0137455792075966759812907801835048190594443990939408530842918_wp
    real(wp),parameter :: b196  =   -0.0661095317267682844455831341498149531672668252085016565917546_wp
    real(wp),parameter :: b197  =    0
    real(wp),parameter :: b198  =    0.152281696736414447136604697040747131921486432699422112099617_wp
    real(wp),parameter :: b199  =   -0.337741018357599840802300793133998004354643424457539667670080_wp
    real(wp),parameter :: b1910 =   -0.0192825981633995781534949199286824400469353110630787982121133_wp
    real(wp),parameter :: b1911 =   -3.68259269696866809932409015535499603576312120746888880201882_wp
    real(wp),parameter :: b1912 =    3.16197870406982063541533528419683854018352080342887002331312_wp
    real(wp),parameter :: b1913 =   -0.370462522106885290716991856022051125477943482284080569177386_wp
    real(wp),parameter :: b1914 =   -0.0514974200365440434996434456698127984941168616474316871020314_wp
    real(wp),parameter :: b1915 =   -0.000829625532120152946787043541792848416659382675202720677536554_wp
    real(wp),parameter :: b1916 =    0.00000279801041419278598986586589070027583961355402640879503213503_wp
    real(wp),parameter :: b1917 =    0.0418603916412360287969841020776788461794119440689356178942252_wp
    real(wp),parameter :: b1918 =    0.279084255090877355915660874555379649966282167560126269290222_wp
    real(wp),parameter :: b200  =    0.103364471650010477570395435690481791543342708330349879244197_wp
    real(wp),parameter :: b201  =    0
    real(wp),parameter :: b202  =    0
    real(wp),parameter :: b203  =    0.124053094528946761061581889237115328211074784955180298044074_wp
    real(wp),parameter :: b204  =    0.483171167561032899288836480451962508724109257517289177302380_wp
    real(wp),parameter :: b205  =   -0.0387530245694763252085681443767620580395733302341368038804290_wp
    real(wp),parameter :: b206  =    0
    real(wp),parameter :: b207  =   -0.438313820361122420391059788940960176420682836652600698580091_wp
    real(wp),parameter :: b208  =    0
    real(wp),parameter :: b209  =   -0.218636633721676647685111485017151199362509373698288330593486_wp
    real(wp),parameter :: b2010 =   -0.0312334764394719229981634995206440349766174759626578122323015_wp
    real(wp),parameter :: b2011 =    0
    real(wp),parameter :: b2012 =    0
    real(wp),parameter :: b2013 =    0
    real(wp),parameter :: b2014 =    0
    real(wp),parameter :: b2015 =    0
    real(wp),parameter :: b2016 =    0
    real(wp),parameter :: b2017 =    0.0312334764394719229981634995206440349766174759626578122323015_wp
    real(wp),parameter :: b2018 =    0.218636633721676647685111485017151199362509373698288330593486_wp
    real(wp),parameter :: b2019 =    0.438313820361122420391059788940960176420682836652600698580091_wp
    real(wp),parameter :: b210  =    0.193333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: b211  =    0
    real(wp),parameter :: b212  =    0.220000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b213  =   -0.0800000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b214  =    0
    real(wp),parameter :: b215  =    0
    real(wp),parameter :: b216  =    0.0984256130499315928152900286856048243348202521491288575952143_wp
    real(wp),parameter :: b217  =   -0.196410889223054653446526504390100417677539095340135532418849_wp
    real(wp),parameter :: b218  =    0
    real(wp),parameter :: b219  =    0.436457930493068729391826122587949137609670676712525034763317_wp
    real(wp),parameter :: b2110 =    0.0652613721675721098560370939805555698350543810708414716730270_wp
    real(wp),parameter :: b2111 =    0
    real(wp),parameter :: b2112 =    0
    real(wp),parameter :: b2113 =    0
    real(wp),parameter :: b2114 =    0
    real(wp),parameter :: b2115 =    0
    real(wp),parameter :: b2116 =    0
    real(wp),parameter :: b2117 =   -0.0652613721675721098560370939805555698350543810708414716730270_wp
    real(wp),parameter :: b2118 =   -0.436457930493068729391826122587949137609670676712525034763317_wp
    real(wp),parameter :: b2119 =    0.196410889223054653446526504390100417677539095340135532418849_wp
    real(wp),parameter :: b2120 =   -0.0984256130499315928152900286856048243348202521491288575952143_wp
    real(wp),parameter :: b220  =   -0.216049382716049382716049382716049382716049382716049382716049_wp
    real(wp),parameter :: b221  =    0.771604938271604938271604938271604938271604938271604938271605_wp
    real(wp),parameter :: b222  =    0
    real(wp),parameter :: b223  =    0
    real(wp),parameter :: b224  =   -0.666666666666666666666666666666666666666666666666666666666667_wp
    real(wp),parameter :: b225  =    0
    real(wp),parameter :: b226  =   -0.390696469295978451446999802258495981249099665294395945559163_wp
    real(wp),parameter :: b227  =    0
    real(wp),parameter :: b228  =    0
    real(wp),parameter :: b229  =    0
    real(wp),parameter :: b2210 =    0
    real(wp),parameter :: b2211 =    0
    real(wp),parameter :: b2212 =    0
    real(wp),parameter :: b2213 =    0
    real(wp),parameter :: b2214 =    0
    real(wp),parameter :: b2215 =    0
    real(wp),parameter :: b2216 =    0
    real(wp),parameter :: b2217 =    0
    real(wp),parameter :: b2218 =    0
    real(wp),parameter :: b2219 =    0
    real(wp),parameter :: b2220 =    0.390696469295978451446999802258495981249099665294395945559163_wp
    real(wp),parameter :: b2221 =    0.666666666666666666666666666666666666666666666666666666666667_wp
    real(wp),parameter :: b230  =    0.200000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b231  =    0
    real(wp),parameter :: b232  =   -0.164609053497942386831275720164609053497942386831275720164609_wp
    real(wp),parameter :: b233  =    0
    real(wp),parameter :: b234  =    0
    real(wp),parameter :: b235  =    0
    real(wp),parameter :: b236  =    0
    real(wp),parameter :: b237  =    0
    real(wp),parameter :: b238  =    0
    real(wp),parameter :: b239  =    0
    real(wp),parameter :: b2310 =    0
    real(wp),parameter :: b2311 =    0
    real(wp),parameter :: b2312 =    0
    real(wp),parameter :: b2313 =    0
    real(wp),parameter :: b2314 =    0
    real(wp),parameter :: b2315 =    0
    real(wp),parameter :: b2316 =    0
    real(wp),parameter :: b2317 =    0
    real(wp),parameter :: b2318 =    0
    real(wp),parameter :: b2319 =    0
    real(wp),parameter :: b2320 =    0
    real(wp),parameter :: b2321 =    0
    real(wp),parameter :: b2322 =    0.164609053497942386831275720164609053497942386831275720164609_wp
    real(wp),parameter :: b240  =    1.47178724881110408452949550989023611293535315518571691939396_wp
    real(wp),parameter :: b241  =    0.787500000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b242  =    0.421296296296296296296296296296296296296296296296296296296296_wp
    real(wp),parameter :: b243  =    0
    real(wp),parameter :: b244  =    0.291666666666666666666666666666666666666666666666666666666667_wp
    real(wp),parameter :: b245  =    0
    real(wp),parameter :: b246  =    0.348600717628329563206854421629657569274689947367847465753757_wp
    real(wp),parameter :: b247  =    0.229499544768994849582890233710555447073823569666506700662510_wp
    real(wp),parameter :: b248  =    5.79046485790481979159831978177003471098279506036722411333192_wp
    real(wp),parameter :: b249  =    0.418587511856506868874073759426596207226461447604248151080016_wp
    real(wp),parameter :: b2410 =    0.307039880222474002649653817490106690389251482313213999386651_wp
    real(wp),parameter :: b2411 =   -4.68700905350603332214256344683853248065574415794742040470287_wp
    real(wp),parameter :: b2412 =    3.13571665593802262152038152399873856554395436199962915429076_wp
    real(wp),parameter :: b2413 =    1.40134829710965720817510506275620441055845017313930508348898_wp
    real(wp),parameter :: b2414 =   -5.52931101439499023629010306005764336421276055777658156400910_wp
    real(wp),parameter :: b2415 =   -0.853138235508063349309546894974784906188927508039552519557498_wp
    real(wp),parameter :: b2416 =    0.103575780373610140411804607167772795518293914458500175573749_wp
    real(wp),parameter :: b2417 =   -0.140474416950600941142546901202132534870665923700034957196546_wp
    real(wp),parameter :: b2418 =   -0.418587511856506868874073759426596207226461447604248151080016_wp
    real(wp),parameter :: b2419 =   -0.229499544768994849582890233710555447073823569666506700662510_wp
    real(wp),parameter :: b2420 =   -0.348600717628329563206854421629657569274689947367847465753757_wp
    real(wp),parameter :: b2421 =   -0.291666666666666666666666666666666666666666666666666666666667_wp
    real(wp),parameter :: b2422 =   -0.421296296296296296296296296296296296296296296296296296296296_wp
    real(wp),parameter :: b2423 =   -0.787500000000000000000000000000000000000000000000000000000000_wp

    associate (f0 => me%funcs(:,1), &
               f1 => me%funcs(:,2), &
               f2 => me%funcs(:,3), &
               f3 => me%funcs(:,4), &
               f4 => me%funcs(:,5), &
               f5 => me%funcs(:,6), &
               f6 => me%funcs(:,7), &
               f7 => me%funcs(:,8), &
               f8 => me%funcs(:,9), &
               f9  => me%funcs(:,10), &
               f10 => me%funcs(:,11), &
               f11 => me%funcs(:,12), &
               f12 => me%funcs(:,13), &
               f13 => me%funcs(:,14), &
               f14 => me%funcs(:,15), &
               f15 => me%funcs(:,16), &
               f16 => me%funcs(:,17), &
               f17 => me%funcs(:,18), &
               f18 => me%funcs(:,19), &
               f19 => me%funcs(:,20), &
               f20 => me%funcs(:,21), &
               f21 => me%funcs(:,22), &
               f22 => me%funcs(:,23), &
               f23 => me%funcs(:,24), &
               f24 => me%funcs(:,25))

        call me%f(t+a0*h,  x,f0)
        call me%f(t+a1*h,  x+h*(b10*f0),f1)
        call me%f(t+a2*h,  x+h*(b20*f0  + b21 *f1),f2)
        call me%f(t+a3*h,  x+h*(b30*f0  + b31 *f1 + b32*f2),f3)
        call me%f(t+a4*h,  x+h*(b40*f0  + b41 *f1 + b42*f2  + b43*f3),f4)
        call me%f(t+h,     x+h*(b50*f0  + b51 *f1 + b52*f2  + b53*f3  + &
                                b54*f4),f5)
        call me%f(t+a6*h,  x+h*(b60*f0  + b61 *f1 + b62*f2  + b63*f3  + &
                                b64*f4  + b65*f5),f6)
        call me%f(t+a7*h,  x+h*(b70*f0  + b71 *f1 + b72*f2  + b73*f3  + &
                                b74*f4  + b75*f5  + b76*f6),f7)
        call me%f(t+a8*h,  x+h*(b80*f0  + b81 *f1 + b82*f2  + b83*f3  + &
                                b84*f4  + b85*f5  + b86*f6  + b87*f7),f8)
        call me%f(t+a9*h,  x+h*(b90*f0  + b91 *f1 + b92*f2  + b93*f3  + &
                                b94*f4  + b95*f5  + b96*f6  + b97*f7  + &
                                b98*f8),f9)
        call me%f(t+a10*h, x+h*(b100*f0 + b101*f1 + b102*f2 + b103*f3 + &
                                b104*f4 + b105*f5 + b106*f6 + b107*f7 + &
                                b108*f8 + b109*f9),f10)
        call me%f(t+a11*h, x+h*(b110*f0 + b111*f1 + b112*f2 + b113*f3 + &
                                b114*f4 + b115*f5 + b116*f6 + b117*f7 + &
                                b118*f8 + b119*f9 + b1110*f10),f11)
        call me%f(t+a12*h, x+h*(b120*f0 + b121*f1 + b122*f2 + b123*f3 + &
                                b124*f4 + b125*f5 + b126*f6 + b127*f7 + &
                                b128*f8 + b129*f9 + b1210*f10 + b1211*f11),f12)
        call me%f(t+a13*h, x+h*(b130*f0 + b131*f1 + b132*f2 + b133*f3 + &
                                b134*f4 + b135*f5 + b136*f6 + b137*f7 + &
                                b138*f8 + b139*f9 + b1310*f10 + b1311*f11 + &
                                b1312*f12),f13)
        call me%f(t+a14*h, x+h*(b140*f0 + b141*f1 + b142*f2 + b143*f3 + &
                                b144*f4 + b145*f5 + b146*f6 + b147*f7 + &
                                b148*f8 + b149*f9 + b1410*f10 + b1411*f11 + &
                                b1412*f12 + b1413*f13),f14)
        call me%f(t+a15*h, x+h*(b150*f0 + b151*f1 + b152*f2 + b153*f3 + &
                                b154*f4 + b155*f5 + b156*f6 + b157*f7 + &
                                b158*f8 + b159*f9 + b1510*f10 + b1511*f11 + &
                                b1512*f12 + b1513*f13 + b1514*f14),f15)
        call me%f(t+a16*h, x+h*(b160*f0 + b161*f1 + b162*f2 + b163*f3 + &
                                b164*f4 + b165*f5 + b166*f6 + b167*f7 + &
                                b168*f8 + b169*f9 + b1610*f10 + b1611*f11 + &
                                b1612*f12 + b1613*f13 + b1614*f14 + b1615*f15),f16)
        call me%f(t+a17*h, x+h*(b170*f0 + b171*f1 + b172*f2 + b173*f3 + &
                                b174*f4 + b175*f5 + b176*f6 + b177*f7 + &
                                b178*f8 + b179*f9 + b1710*f10 + b1711*f11 + &
                                b1712*f12 + b1713*f13 + b1714*f14 + b1715*f15 + &
                                b1716*f16),f17)
        call me%f(t+a18*h, x+h*(b180*f0 + b181*f1 + b182*f2 + b183*f3 + &
                                b184*f4 + b185*f5 + b186*f6 + b187*f7 + &
                                b188*f8 + b189*f9 + b1810*f10 + b1811*f11 + &
                                b1812*f12 + b1813*f13 + b1814*f14 + b1815*f15 + &
                                b1816*f16 + b1817*f17),f18)
        call me%f(t+a19*h, x+h*(b190*f0 + b191*f1 + b192*f2 + b193*f3 + &
                                b194*f4 + b195*f5 + b196*f6 + b197*f7 + &
                                b198*f8 + b199*f9 + b1910*f10 + b1911*f11 + &
                                b1912*f12 + b1913*f13 + b1914*f14 + b1915*f15 + &
                                b1916*f16 + b1917*f17 + b1918*f18),f19)
        call me%f(t+a20*h, x+h*(b200*f0 + b201*f1 + b202*f2 + b203*f3 + &
                                b204*f4 + b205*f5 + b206*f6 + b207*f7 + &
                                b208*f8 + b209*f9 + b2010*f10 + b2011*f11 + &
                                b2012*f12 + b2013*f13 + b2014*f14 + b2015*f15 + &
                                b2016*f16 + b2017*f17 + b2018*f18 + b2019*f19),f20)
        call me%f(t+a21*h, x+h*(b210*f0 + b211*f1 + b212*f2 + b213*f3 + &
                                b214*f4 + b215*f5 + b216*f6 + b217*f7 + &
                                b218*f8 + b219*f9 + b2110*f10 + b2111*f11 + &
                                b2112*f12 + b2113*f13 + b2114*f14 + b2115*f15 + &
                                b2116*f16 + b2117*f17 + b2118*f18 + b2119*f19 + &
                                b2120*f20),f21)
        call me%f(t+a22*h, x+h*(b220*f0 + b221*f1 + b222*f2 + b223*f3 + &
                                b224*f4 + b225*f5 + b226*f6 + b227*f7 + &
                                b228*f8 + b229*f9 + b2210*f10 + b2211*f11 + &
                                b2212*f12 + b2213*f13 + b2214*f14 + b2215*f15 + &
                                b2216*f16 + b2217*f17 + b2218*f18 + b2219*f19 + &
                                b2220*f20 + b2221*f21),f22)
        call me%f(t+a23*h, x+h*(b230*f0 + b231*f1 + b232*f2 + b233*f3 + &
                                b234*f4 + b235*f5 + b236*f6 + b237*f7 + &
                                b238*f8 + b239*f9 + b2310*f10 + b2311*f11 + &
                                b2312*f12 + b2313*f13 + b2314*f14 + b2315*f15 + &
                                b2316*f16 + b2317*f17 + b2318*f18 + b2319*f19 + &
                                b2320*f20 + b2321*f21 + b2322*f22),f23)
        call me%f(t+a24*h, x+h*(b240*f0 + b241*f1 + b242*f2 + b243*f3 + &
                                b244*f4 + b245*f5 + b246*f6 + b247*f7 + &
                                b248*f8 + b249*f9 + b2410*f10 + b2411*f11 + &
                                b2412*f12 + b2413*f13 + b2414*f14 + b2415*f15 + &
                                b2416*f16 + b2417*f17 + b2418*f18 + b2419*f19 + &
                                b2420*f20 + b2421*f21 + b2422*f22 + b2423*f23),f24)

        xf = x+h*(  c0*f0   + &
                    c1*f1   + &
                    c2*f2   + &
                    c3*f3   + &
                    c4*f4   + &
                    c5*f5   + &
                    c6*f6   + &
                    c7*f7   + &
                    c8*f8   + &
                    c9*f9   + &
                    c10*f10 + &
                    c11*f11 + &
                    c12*f12 + &
                    c13*f13 + &
                    c14*f14 + &
                    c15*f15 + &
                    c16*f16 + &
                    c17*f17 + &
                    c18*f18 + &
                    c19*f19 + &
                    c20*f20 + &
                    c21*f21 + &
                    c22*f22 + &
                    c23*f23 + &
                    c24*f24 )

        xerr = (49.0_wp/640.0_wp)*h*(f1-f23)

    end associate

    end procedure rkf1210
!*****************************************************************************************

!*****************************************************************************************
!>
!  Feagin's RK14(12) - a 14th-order method with an embedded 12th-order method.
!
!### Reference
!  * [coefficient file](http://sce.uhcl.edu/rungekutta/rk1412.txt)

    module procedure rkf1412

    real(wp),parameter :: a0  = 0
    real(wp),parameter :: a1  = 0.111111111111111111111111111111111111111111111111111111111111_wp
    real(wp),parameter :: a2  = 0.555555555555555555555555555555555555555555555555555555555556_wp
    real(wp),parameter :: a3  = 0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a4  = 0.333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a6  = 0.669986979272772921764683785505998513938845229638460353285142_wp
    real(wp),parameter :: a7  = 0.297068384213818357389584716808219413223332094698915687379168_wp
    real(wp),parameter :: a8  = 0.727272727272727272727272727272727272727272727272727272727273_wp
    real(wp),parameter :: a9  = 0.140152799042188765276187487966946717629806463082532936287323_wp
    real(wp),parameter :: a10 = 0.700701039770150737151099854830749337941407049265546408969222_wp
    real(wp),parameter :: a11 = 0.363636363636363636363636363636363636363636363636363636363636_wp
    real(wp),parameter :: a12 = 0.263157894736842105263157894736842105263157894736842105263158_wp
    real(wp),parameter :: a13 = 0.0392172246650270859125196642501208648863714315266128052078483_wp
    real(wp),parameter :: a14 = 0.812917502928376762983393159278036506189612372617238550774312_wp
    real(wp),parameter :: a15 = 0.166666666666666666666666666666666666666666666666666666666667_wp
    real(wp),parameter :: a16 = 0.900000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: a17 = 0.0641299257451966923312771193896682809481096651615083225402924_wp
    real(wp),parameter :: a18 = 0.204149909283428848927744634301023405027149505241333751628870_wp
    real(wp),parameter :: a19 = 0.395350391048760565615671369827324372352227297456659450554577_wp
    real(wp),parameter :: a20 = 0.604649608951239434384328630172675627647772702543340549445423_wp
    real(wp),parameter :: a21 = 0.795850090716571151072255365698976594972850494758666248371130_wp
    real(wp),parameter :: a22 = 0.935870074254803307668722880610331719051890334838491677459708_wp
    real(wp),parameter :: a23 = 0.166666666666666666666666666666666666666666666666666666666667_wp
    real(wp),parameter :: a24 = 0.812917502928376762983393159278036506189612372617238550774312_wp
    real(wp),parameter :: a25 = 0.0392172246650270859125196642501208648863714315266128052078483_wp
    real(wp),parameter :: a26 = 0.363636363636363636363636363636363636363636363636363636363636_wp
    real(wp),parameter :: a27 = 0.700701039770150737151099854830749337941407049265546408969222_wp
    real(wp),parameter :: a28 = 0.140152799042188765276187487966946717629806463082532936287323_wp
    real(wp),parameter :: a29 = 0.297068384213818357389584716808219413223332094698915687379168_wp
    real(wp),parameter :: a30 = 0.669986979272772921764683785505998513938845229638460353285142_wp
    real(wp),parameter :: a31 = 0.333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a32 = 0.555555555555555555555555555555555555555555555555555555555556_wp
    real(wp),parameter :: a33 = 0.111111111111111111111111111111111111111111111111111111111111_wp

    real(wp),parameter :: c0  = 0.0178571428571428571428571428571428571428571428571428571428571_wp
    real(wp),parameter :: c1  = 0.00585937500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c2  = 0.0117187500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c3  = 0
    real(wp),parameter :: c4  = 0.0175781250000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c5  = 0
    real(wp),parameter :: c6  = 0.0234375000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c7  = 0.0292968750000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c8  = 0
    real(wp),parameter :: c9  = 0.0351562500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c10 = 0.0410156250000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c11 = 0.0468750000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c12 = 0
    real(wp),parameter :: c13 = 0.0527343750000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c14 = 0.0585937500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c15 = 0.0644531250000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c16 = 0
    real(wp),parameter :: c17 = 0.105352113571753019691496032887878162227673083080523884041670_wp
    real(wp),parameter :: c18 = 0.170561346241752182382120338553874085887555487802790804737501_wp
    real(wp),parameter :: c19 = 0.206229397329351940783526485701104894741914286259542454077972_wp
    real(wp),parameter :: c20 = 0.206229397329351940783526485701104894741914286259542454077972_wp
    real(wp),parameter :: c21 = 0.170561346241752182382120338553874085887555487802790804737501_wp
    real(wp),parameter :: c22 = 0.105352113571753019691496032887878162227673083080523884041670_wp
    real(wp),parameter :: c23 = -0.0644531250000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c24 = -0.0585937500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c25 = -0.0527343750000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c26 = -0.0468750000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c27 = -0.0410156250000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c28 = -0.0351562500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c29 = -0.0292968750000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c30 = -0.0234375000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c31 = -0.0175781250000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c32 = -0.0117187500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c33 = -0.00585937500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c34 = 0.0178571428571428571428571428571428571428571428571428571428571_wp

    real(wp),parameter :: b10   =     0.111111111111111111111111111111111111111111111111111111111111_wp
    real(wp),parameter :: b20   =    -0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: b21   =     1.38888888888888888888888888888888888888888888888888888888889_wp
    real(wp),parameter :: b30   =     0.208333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: b31   =     0
    real(wp),parameter :: b32   =     0.625000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b40   =     0.193333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: b41   =     0
    real(wp),parameter :: b42   =     0.220000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b43   =    -0.0800000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b50   =     0.100000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b51   =     0
    real(wp),parameter :: b52   =     0
    real(wp),parameter :: b53   =     0.400000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b54   =     0.500000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b60   =     0.103484561636679776672993546511910344499744798201971316606663_wp
    real(wp),parameter :: b61   =     0
    real(wp),parameter :: b62   =     0
    real(wp),parameter :: b63   =     0.122068887306407222589644082868962077139592714834162134741275_wp
    real(wp),parameter :: b64   =     0.482574490331246622475134780125688112865919023850168049679402_wp
    real(wp),parameter :: b65   =    -0.0381409600015606999730886240005620205664113072478411477421970_wp
    real(wp),parameter :: b70   =     0.124380526654094412881516420868799316268491466359671423163289_wp
    real(wp),parameter :: b71   =     0
    real(wp),parameter :: b72   =     0
    real(wp),parameter :: b73   =     0
    real(wp),parameter :: b74   =     0.226120282197584301422238662979202901196752320742633143965145_wp
    real(wp),parameter :: b75   =     0.0137885887618080880607695837016477814530969417491493385363543_wp
    real(wp),parameter :: b76   =    -0.0672210133996684449749399507414305856950086341525382182856200_wp
    real(wp),parameter :: b80   =     0.0936919065659673815530885456083005933866349695217750085655603_wp
    real(wp),parameter :: b81   =     0
    real(wp),parameter :: b82   =     0
    real(wp),parameter :: b83   =     0
    real(wp),parameter :: b84   =     0
    real(wp),parameter :: b85   =    -0.00613406843450510987229498995641664735620914507128858871007099_wp
    real(wp),parameter :: b86   =     0.216019825625503063708860097659866573490979433278117320188668_wp
    real(wp),parameter :: b87   =     0.423695063515761937337619073960976753205867469544123532683116_wp
    real(wp),parameter :: b90   =     0.0838479812409052664616968791372814085980533139224911131069335_wp
    real(wp),parameter :: b91   =     0
    real(wp),parameter :: b92   =     0
    real(wp),parameter :: b93   =     0
    real(wp),parameter :: b94   =     0
    real(wp),parameter :: b95   =    -0.0117949367100973814319755056031295775367961960590736150777613_wp
    real(wp),parameter :: b96   =    -0.247299020568812652339473838743194598325992840353340132697498_wp
    real(wp),parameter :: b97   =     0.0978080858367729012259313014081291665503740655476733940756599_wp
    real(wp),parameter :: b98   =     0.217590689243420631360008651767860318344168120024782176879989_wp
    real(wp),parameter :: b100   =     0.0615255359769428227954562389614314714333423969064821107453940_wp
    real(wp),parameter :: b101   =     0
    real(wp),parameter :: b102   =     0
    real(wp),parameter :: b103   =     0
    real(wp),parameter :: b104   =     0
    real(wp),parameter :: b105   =     0.00592232780324503308042990005798046524738389560444257136834990_wp
    real(wp),parameter :: b106   =     0.470326159963841112217224303205894113455362530746108825010848_wp
    real(wp),parameter :: b107   =     0.299688863848679000853981837096192399136831121671781279184194_wp
    real(wp),parameter :: b108   =    -0.247656877593994914689992276329810825853958069263947095548189_wp
    real(wp),parameter :: b109   =     0.110895029771437682893999851839061714522445173600678718208625_wp
    real(wp),parameter :: b110   =     0.0419700073362782579861792864787277787213483656543104611245994_wp
    real(wp),parameter :: b111   =     0
    real(wp),parameter :: b112   =     0
    real(wp),parameter :: b113   =     0
    real(wp),parameter :: b114   =     0
    real(wp),parameter :: b115   =    -0.00317987696266205093901912847692712407988609169703103952205634_wp
    real(wp),parameter :: b116   =     0.806397714906192077260821711520379506393543111567419750119748_wp
    real(wp),parameter :: b117   =     0.0975983126412388979093522850684288851314672048003054550357187_wp
    real(wp),parameter :: b118   =     0.778575578158398909027512446452927238999763460594181964958853_wp
    real(wp),parameter :: b119   =     0.204890423831599428189499202098105603312029235081420653574829_wp
    real(wp),parameter :: b1110  =   -1.56261579627468188307070943950527825211462892236424360892806_wp
    real(wp),parameter :: b120   =     0.0437726782233730163574465242495339811688214967071614123256973_wp
    real(wp),parameter :: b121   =     0
    real(wp),parameter :: b122   =     0
    real(wp),parameter :: b123   =     0
    real(wp),parameter :: b124   =     0
    real(wp),parameter :: b125   =     0
    real(wp),parameter :: b126   =     0
    real(wp),parameter :: b127   =     0
    real(wp),parameter :: b128   =     0.00624365027520195208794358628580933625281631216903095917201250_wp
    real(wp),parameter :: b129   =     0.200043097109577314994435165469647856829066232218264969608768_wp
    real(wp),parameter :: b1210  =   -0.00805328367804983036823857162048902911923392887337029314844206_wp
    real(wp),parameter :: b1211  =    0.0211517528067396521915711903523399601316877825157550573051221_wp
    real(wp),parameter :: b130   =     0.0283499250363514563095023591920717312247137654896477097768495_wp
    real(wp),parameter :: b131   =     0
    real(wp),parameter :: b132   =     0
    real(wp),parameter :: b133   =     0
    real(wp),parameter :: b134   =     0
    real(wp),parameter :: b135   =     0
    real(wp),parameter :: b136   =     0
    real(wp),parameter :: b137   =     0
    real(wp),parameter :: b138   =     0.00249163204855817407538949148805995149459884653585417680098222_wp
    real(wp),parameter :: b139   =     0.0230138787854593149638399846373742768772087122638142234223658_wp
    real(wp),parameter :: b1310  =   -0.00322155956692977098724476092467120878189463604760620461043308_wp
    real(wp),parameter :: b1311  =    0.00988442549447664668946335414487885256040819982786014648129297_wp
    real(wp),parameter :: b1312  =   -0.0213010771328887351384307642875927384886634565429572466632092_wp
    real(wp),parameter :: b140   =     0.343511894290243001049432234735147943083353174980701426268122_wp
    real(wp),parameter :: b141   =     0
    real(wp),parameter :: b142   =     0
    real(wp),parameter :: b143   =     0
    real(wp),parameter :: b144   =     0
    real(wp),parameter :: b145   =     0
    real(wp),parameter :: b146   =     0
    real(wp),parameter :: b147   =     0
    real(wp),parameter :: b148   =     0.210451912023627385609097011999010655788807405225626700040882_wp
    real(wp),parameter :: b149   =     1.03427452057230411936482926828825709938667999698324740166559_wp
    real(wp),parameter :: b1410  =    0.00600303645864422487051240448206640574939078092406156945568306_wp
    real(wp),parameter :: b1411  =    0.855938125099619537578012106002407728915062652616416005816477_wp
    real(wp),parameter :: b1412  =   -0.977235005036766810872264852372525633013107656892839677696022_wp
    real(wp),parameter :: b1413  =   -0.660026980479294694616225013856327693720573981219974874776419_wp
    real(wp),parameter :: b150   =    -0.0143574001672168069538206399935076366657755954378399880691949_wp
    real(wp),parameter :: b151   =     0
    real(wp),parameter :: b152   =     0
    real(wp),parameter :: b153   =     0
    real(wp),parameter :: b154   =     0
    real(wp),parameter :: b155   =     0
    real(wp),parameter :: b156   =     0
    real(wp),parameter :: b157   =     0
    real(wp),parameter :: b158   =    -0.0366253270049039970293685796848974791733119081733552207318285_wp
    real(wp),parameter :: b159   =     0.0350254975636213681976849406979846524346789082471103574920148_wp
    real(wp),parameter :: b1510  =    0.0360946016362113508931786658758335239823689929864237671348749_wp
    real(wp),parameter :: b1511  =   -0.0265219967553681106351595946834601923649627012457464284442911_wp
    real(wp),parameter :: b1512  =    0.0445699011305698119638911537508839908104336323082226770910408_wp
    real(wp),parameter :: b1513  =    0.124343093331358243286225595741786448038973408895106741855721_wp
    real(wp),parameter :: b1514  =    0.00413829693239480694403512496204335960426192908674476033832967_wp
    real(wp),parameter :: b160   =     0.356032404425120290975609116398089176264106222379748802654822_wp
    real(wp),parameter :: b161   =     0
    real(wp),parameter :: b162   =     0
    real(wp),parameter :: b163   =     0
    real(wp),parameter :: b164   =     0
    real(wp),parameter :: b165   =     0
    real(wp),parameter :: b166   =     0
    real(wp),parameter :: b167   =     0
    real(wp),parameter :: b168   =    -0.450192758947562595966821779075956175110645100214763601190349_wp
    real(wp),parameter :: b169   =     0.430527907083710898626656292808782917793030154094709462877146_wp
    real(wp),parameter :: b1610  =    0.511973029011022237668556960394071692077125787030651386389972_wp
    real(wp),parameter :: b1611  =    0.908303638886404260390159124638110213997496214819904630546596_wp
    real(wp),parameter :: b1612  =   -1.23921093371933931757372469151534028854413889248605726186520_wp
    real(wp),parameter :: b1613  =   -0.649048661671761465141672348879062553905402831967191097656668_wp
    real(wp),parameter :: b1614  =    0.251708904586819292210480529948970541404887852931447491219418_wp
    real(wp),parameter :: b1615  =    0.779906470345586398810756795282334476023540593411550187024263_wp
    real(wp),parameter :: b170   =     0.0130935687406513066406881206418834980127470438213192487844956_wp
    real(wp),parameter :: b171   =     0
    real(wp),parameter :: b172   =     0
    real(wp),parameter :: b173   =     0
    real(wp),parameter :: b174   =     0
    real(wp),parameter :: b175   =     0
    real(wp),parameter :: b176   =     0
    real(wp),parameter :: b177   =     0
    real(wp),parameter :: b178   =     0
    real(wp),parameter :: b179   =     0
    real(wp),parameter :: b1710  =    0
    real(wp),parameter :: b1711  =    0
    real(wp),parameter :: b1712  =   -0.0000932053067985113945908461962767108237858631509684667142124826_wp
    real(wp),parameter :: b1713  =    0.0505374334262299359640090443138590726770942344716122381702746_wp
    real(wp),parameter :: b1714  =    8.04470341944487979109579109610197797641311868930865361048975e-7_wp
    real(wp),parameter :: b1715  =    0.000591726029494171190528755742777717259844340971924321528178248_wp
    real(wp),parameter :: b1716  =   -4.01614722154557337064691684906375587732264247950093804676867e-7_wp
    real(wp),parameter :: b180   =     0.0207926484466053012541944544000765652167255206144373407979758_wp
    real(wp),parameter :: b181   =     0
    real(wp),parameter :: b182   =     0
    real(wp),parameter :: b183   =     0
    real(wp),parameter :: b184   =     0
    real(wp),parameter :: b185   =     0
    real(wp),parameter :: b186   =     0
    real(wp),parameter :: b187   =     0
    real(wp),parameter :: b188   =     0
    real(wp),parameter :: b189   =     0
    real(wp),parameter :: b1810  =    0
    real(wp),parameter :: b1811  =    0
    real(wp),parameter :: b1812  =    0.000582695918800085915101902697837284108951406103029871570103075_wp
    real(wp),parameter :: b1813  =   -0.00801700732358815939083342186525852746640558465919633524655451_wp
    real(wp),parameter :: b1814  =    4.03847643847136940375170821743560570484117290330895506618968e-6_wp
    real(wp),parameter :: b1815  =    0.0854609998055506144225056114567535602510114622033622491802597_wp
    real(wp),parameter :: b1816  =   -2.04486480935804242706707569691004307904442837552677456232848e-6_wp
    real(wp),parameter :: b1817  =    0.105328578824431893399799402979093997354240904235172843146582_wp
    real(wp),parameter :: b190   =     1.40153449795736021415446247355771306718486452917597731683689_wp
    real(wp),parameter :: b191   =     0
    real(wp),parameter :: b192   =     0
    real(wp),parameter :: b193   =     0
    real(wp),parameter :: b194   =     0
    real(wp),parameter :: b195   =     0
    real(wp),parameter :: b196   =     0
    real(wp),parameter :: b197   =     0
    real(wp),parameter :: b198   =     0
    real(wp),parameter :: b199   =     0
    real(wp),parameter :: b1910  =    0
    real(wp),parameter :: b1911  =    0
    real(wp),parameter :: b1912  =   -0.230252000984221261616272410367415621261130298274455611733277_wp
    real(wp),parameter :: b1913  =   -7.21106840466912905659582237106874247165856493509961561958267_wp
    real(wp),parameter :: b1914  =    0.00372901560694836335236995327852132340217759566678662385552634_wp
    real(wp),parameter :: b1915  =   -4.71415495727125020678778179392224757011323373221820091641216_wp
    real(wp),parameter :: b1916  =   -0.00176367657545349242053841995032797673574903886695600132759652_wp
    real(wp),parameter :: b1917  =    7.64130548038698765563029310880237651185173367813936997648198_wp
    real(wp),parameter :: b1918  =    3.50602043659751834989896082949744710968212949893375368243588_wp
    real(wp),parameter :: b200   =     11.9514650694120686799372385830716401674473610826553517297976_wp
    real(wp),parameter :: b201   =     0
    real(wp),parameter :: b202   =     0
    real(wp),parameter :: b203   =     0
    real(wp),parameter :: b204   =     0
    real(wp),parameter :: b205   =     0
    real(wp),parameter :: b206   =     0
    real(wp),parameter :: b207   =     0
    real(wp),parameter :: b208   =     0
    real(wp),parameter :: b209   =     0
    real(wp),parameter :: b2010  =    0
    real(wp),parameter :: b2011  =    0
    real(wp),parameter :: b2012  =    7.79480932108175968783516700231764388220284279598980948538579_wp
    real(wp),parameter :: b2013  =   -56.4501393867325792523560991120904281440468100061340556540132_wp
    real(wp),parameter :: b2014  =    0.0912376306930644901344530449290276645709607450403673704844997_wp
    real(wp),parameter :: b2015  =   -12.7336279925434886201945524309199275038162717529918963305155_wp
    real(wp),parameter :: b2016  =    -0.0396895921904719712313542810939736674712383070433147873009352_wp
    real(wp),parameter :: b2017  =    54.4392141883570886996225765155307791861438378423305337073797_wp
    real(wp),parameter :: b2018  =    -3.64411637921569236846406990361350645806721478409266709351203_wp
    real(wp),parameter :: b2019  =    -0.804503249910509910899030787958579499315694913210787878260459_wp
    real(wp),parameter :: b210   =     -148.809426507100488427838868268647625561930612082148597076690_wp
    real(wp),parameter :: b211   =     0
    real(wp),parameter :: b212   =     0
    real(wp),parameter :: b213   =     0
    real(wp),parameter :: b214   =     0
    real(wp),parameter :: b215   =     0
    real(wp),parameter :: b216   =     0
    real(wp),parameter :: b217   =     0
    real(wp),parameter :: b218   =     0
    real(wp),parameter :: b219   =     0
    real(wp),parameter :: b2110  =    0
    real(wp),parameter :: b2111  =    0
    real(wp),parameter :: b2112  =    -91.7295278291256484357935662402321623495228729036354276506427_wp
    real(wp),parameter :: b2113  =    707.656144971598359834575719286335716154821128966649565194286_wp
    real(wp),parameter :: b2114  =    -1.10563611857482440905296961311590930801338308942637769555540_wp
    real(wp),parameter :: b2115  =    176.134591883811372587859898076055660406999516762301689616841_wp
    real(wp),parameter :: b2116  =    0.491384824214880662268898345164454557416884631402764792538746_wp
    real(wp),parameter :: b2117  =    -684.278000449814944358237535610895081956077167893600278300805_wp
    real(wp),parameter :: b2118  =    27.9910604998398258984224332124380407446002518400668657974589_wp
    real(wp),parameter :: b2119  =    13.1939710030282333443670964371153238435064159623744975073252_wp
    real(wp),parameter :: b2120  =    1.25128781283980445450114974148056006317268830077396406361417_wp
    real(wp),parameter :: b220   =     -9.67307946948196763644126118433219395839951408571877262880482_wp
    real(wp),parameter :: b221   =     0
    real(wp),parameter :: b222   =     0
    real(wp),parameter :: b223   =     0
    real(wp),parameter :: b224   =     0
    real(wp),parameter :: b225   =     0
    real(wp),parameter :: b226   =     0
    real(wp),parameter :: b227   =     0
    real(wp),parameter :: b228   =     0
    real(wp),parameter :: b229   =     0
    real(wp),parameter :: b2210  =    0
    real(wp),parameter :: b2211  =    0
    real(wp),parameter :: b2212  =    -4.46990150858505531443846227701960360497830681408751431146712_wp
    real(wp),parameter :: b2213  =    45.5127128690952681968241950400052751178905907817398483534845_wp
    real(wp),parameter :: b2214  =    -0.0713085086183826912791492024438246129930559805352394367050813_wp
    real(wp),parameter :: b2215  =    11.2273614068412741582590624479939384207826800776794485051540_wp
    real(wp),parameter :: b2216  =    0.126244376717622724516237912909138809361786889819105426371393_wp
    real(wp),parameter :: b2217  =    -43.5439339549483313605810624907242107623814304467621407753424_wp
    real(wp),parameter :: b2218  =    0.787174307543058978398792994996550902064546091443233850464377_wp
    real(wp),parameter :: b2219  =    0.532264696744684215669300708603886690785395776821503851830821_wp
    real(wp),parameter :: b2220  =    0.422422733996325326010225127471388772575086538809603346825334_wp
    real(wp),parameter :: b2221  =    0.0859131249503067107308438031499859443441115056294154956487671_wp
    real(wp),parameter :: b230   =     -10.0664032447054702403396606900426891472202824757968765569183_wp
    real(wp),parameter :: b231   =     0
    real(wp),parameter :: b232   =     0
    real(wp),parameter :: b233   =     0
    real(wp),parameter :: b234   =     0
    real(wp),parameter :: b235   =     0
    real(wp),parameter :: b236   =     0
    real(wp),parameter :: b237   =     0
    real(wp),parameter :: b238   =    -0.0366253270049039970293685796848974791733119081733552207318285_wp
    real(wp),parameter :: b239   =    0.0350254975636213681976849406979846524346789082471103574920148_wp
    real(wp),parameter :: b2310  =    0.0360946016362113508931786658758335239823689929864237671348749_wp
    real(wp),parameter :: b2311  =    -0.0265219967553681106351595946834601923649627012457464284442911_wp
    real(wp),parameter :: b2312  =    -6.27088972181464143590553149478871603839356122957396018530209_wp
    real(wp),parameter :: b2313  =    48.2079237442562989090702103008195063923492593141636117832993_wp
    real(wp),parameter :: b2314  =    -0.0694471689136165640882395180583732834557754169149088630301342_wp
    real(wp),parameter :: b2315  =    12.6810690204850295698341370913609807066108483811412127009785_wp
    real(wp),parameter :: b2316  =    0.0119671168968323754838161435501011294100927813964199613229864_wp
    real(wp),parameter :: b2317  =    -46.7249764992482408003358268242662695593201321659795608950429_wp
    real(wp),parameter :: b2318  =    1.33029613326626711314710039298216591399033511191227101321435_wp
    real(wp),parameter :: b2319  =    1.00766787503398298353438903619926657771162717793661719708370_wp
    real(wp),parameter :: b2320  =    0.0209512051933665091664122388475480702892770753864487241177616_wp
    real(wp),parameter :: b2321  =    0.0210134706331264177317735424331396407424412188443757490871603_wp
    real(wp),parameter :: b2322  =    0.00952196014417121794175101542454575907376360233658356240547761_wp
    real(wp),parameter :: b240   =     -409.478081677743708772589097409370357624424341606752069725341_wp
    real(wp),parameter :: b241   =     0
    real(wp),parameter :: b242   =     0
    real(wp),parameter :: b243   =     0
    real(wp),parameter :: b244   =     0
    real(wp),parameter :: b245   =     0
    real(wp),parameter :: b246   =     0
    real(wp),parameter :: b247   =     0
    real(wp),parameter :: b248   =    0.210451912023627385609097011999010655788807405225626700040882_wp
    real(wp),parameter :: b249   =    1.03427452057230411936482926828825709938667999698324740166559_wp
    real(wp),parameter :: b2410  =    0.00600303645864422487051240448206640574939078092406156945568306_wp
    real(wp),parameter :: b2411  =    0.855938125099619537578012106002407728915062652616416005816477_wp
    real(wp),parameter :: b2412  =    -250.516998547447860492777657729316130386584050420782075966990_wp
    real(wp),parameter :: b2413  =    1946.42466652388427766053750328264758595829850895761428240231_wp
    real(wp),parameter :: b2414  =    -3.04503882102310365506105809086860882786950544097602101685174_wp
    real(wp),parameter :: b2415  =    490.626379528281713521208265299168083841598542274061671576230_wp
    real(wp),parameter :: b2416  =    1.56647589531270907115484067013597445739595615245966775329993_wp
    real(wp),parameter :: b2417  =    -1881.97428994011173362217267377035870619215906638453056643641_wp
    real(wp),parameter :: b2418  =    75.2592224724847175278837713643303149821620618914245864351135_wp
    real(wp),parameter :: b2419  =    34.5734356980331067622434344736554689696728644793551014989002_wp
    real(wp),parameter :: b2420  =    3.21147679440968961435417361847073755169022966748891627882572_wp
    real(wp),parameter :: b2421  =    -0.460408041738414391307201404237058848867245095265382820823055_wp
    real(wp),parameter :: b2422  =    -0.0870718339841810522431884137957986245724252047388936572215438_wp
    real(wp),parameter :: b2423  =    -7.39351814158303067567016952195521063999185773249132944724553_wp
    real(wp),parameter :: b250   =     3.43347475853550878921093496257596781120623891072008459930197_wp
    real(wp),parameter :: b251   =     0
    real(wp),parameter :: b252   =     0
    real(wp),parameter :: b253   =     0
    real(wp),parameter :: b254   =     0
    real(wp),parameter :: b255   =     0
    real(wp),parameter :: b256   =     0
    real(wp),parameter :: b257   =     0
    real(wp),parameter :: b258   =     0.00249163204855817407538949148805995149459884653585417680098222_wp
    real(wp),parameter :: b259   =     0.0230138787854593149638399846373742768772087122638142234223658_wp
    real(wp),parameter :: b2510  =    -0.00322155956692977098724476092467120878189463604760620461043308_wp
    real(wp),parameter :: b2511  =    0.00988442549447664668946335414487885256040819982786014648129297_wp
    real(wp),parameter :: b2512  =    2.16252799377922507788307841904757354045759225335732707916530_wp
    real(wp),parameter :: b2513  =    -16.2699864546457421328065640660139489006987552040228852402716_wp
    real(wp),parameter :: b2514  =    -0.128534502120524552843583417470935010538029037542654506231743_wp
    real(wp),parameter :: b2515  =    -8.98915042666504253089307820833379330486511746063552853023189_wp
    real(wp),parameter :: b2516  =    -0.00348595363232025333387080201851013650192401767250513765000963_wp
    real(wp),parameter :: b2517  =    15.7936194113339807536235187388695574135853387025139738341334_wp
    real(wp),parameter :: b2518  =    -0.574403330914095065628165482017335820148383663195675408024658_wp
    real(wp),parameter :: b2519  =    -0.345602039021393296692722496608124982535237228827655306030152_wp
    real(wp),parameter :: b2520  =    -0.00662241490206585091731619991383757781133067992707418687587487_wp
    real(wp),parameter :: b2521  =    -0.00777788129242204164032546458607364309759347209626759111946150_wp
    real(wp),parameter :: b2522  =    -0.00356084192402274913338827232697437364675240818791706587952939_wp
    real(wp),parameter :: b2523  =    4.79282506449930799649797749629840189457296934139359048988332_wp
    real(wp),parameter :: b2524  =    0.153725464873068577844576387402512082757034273069877432944621_wp
    real(wp),parameter :: b260   =     32.3038520871985442326994734440031535091364975047784630088983_wp
    real(wp),parameter :: b261   =     0
    real(wp),parameter :: b262   =     0
    real(wp),parameter :: b263   =     0
    real(wp),parameter :: b264   =     0
    real(wp),parameter :: b265   =    -0.00317987696266205093901912847692712407988609169703103952205634_wp
    real(wp),parameter :: b266   =    0.806397714906192077260821711520379506393543111567419750119748_wp
    real(wp),parameter :: b267   =    0.0975983126412388979093522850684288851314672048003054550357187_wp
    real(wp),parameter :: b268   =    0.778575578158398909027512446452927238999763460594181964958853_wp
    real(wp),parameter :: b269   =    0.204890423831599428189499202098105603312029235081420653574829_wp
    real(wp),parameter :: b2610  =    -1.56261579627468188307070943950527825211462892236424360892806_wp
    real(wp),parameter :: b2611  =    0
    real(wp),parameter :: b2612  =    16.3429891882310570648504243973927174708753353504154550405647_wp
    real(wp),parameter :: b2613  =    -154.544555293543621230730189631471036399316683669609116705323_wp
    real(wp),parameter :: b2614  =    1.56971088703334872692034283417621761466263593582497085955201_wp
    real(wp),parameter :: b2615  =    3.27685545087248131321429817269900731165522404974733504794135_wp
    real(wp),parameter :: b2616  =    -0.0503489245193653176348040727199783626534081095691632396802451_wp
    real(wp),parameter :: b2617  =    153.321151858041665070593767885914694011224363102594556731397_wp
    real(wp),parameter :: b2618  =    7.17568186327720495846766484814784143567826308034865369443637_wp
    real(wp),parameter :: b2619  =    -2.94036748675300481945917659896930989215320594380777597403592_wp
    real(wp),parameter :: b2620  =    -0.0665845946076803144470749676022628870281920493197256887985612_wp
    real(wp),parameter :: b2621  =    -0.0462346054990843661229248668562217261176966514016859284197145_wp
    real(wp),parameter :: b2622  =    -0.0204198733585679401539388228617269778848579774821581777675337_wp
    real(wp),parameter :: b2623  =    -53.3523106438735850515953441165998107974045090495791591218714_wp
    real(wp),parameter :: b2624  =    -1.35548714715078654978732186705996404017554501614191325114947_wp
    real(wp),parameter :: b2625  =    -1.57196275801232751882901735171459249177687219114442583461866_wp
    real(wp),parameter :: b270   =     -16.6451467486341512872031294403931758764560371130818978459405_wp
    real(wp),parameter :: b271   =     0
    real(wp),parameter :: b272   =     0
    real(wp),parameter :: b273   =     0
    real(wp),parameter :: b274   =     0
    real(wp),parameter :: b275   =     0.00592232780324503308042990005798046524738389560444257136834990_wp
    real(wp),parameter :: b276   =     0.470326159963841112217224303205894113455362530746108825010848_wp
    real(wp),parameter :: b277   =     0.299688863848679000853981837096192399136831121671781279184194_wp
    real(wp),parameter :: b278   =    -0.247656877593994914689992276329810825853958069263947095548189_wp
    real(wp),parameter :: b279   =     0.110895029771437682893999851839061714522445173600678718208625_wp
    real(wp),parameter :: b2710  =    0
    real(wp),parameter :: b2711  =    -0.491719043846229147070666628704194097678081907210673044988866_wp
    real(wp),parameter :: b2712  =    -11.4743154427289496968389492564352536350842454130853175250727_wp
    real(wp),parameter :: b2713  =    80.2593166576230272541702485886484400152793366623589989106256_wp
    real(wp),parameter :: b2714  =    -0.384132303980042847625312526759029103746926841342088219165648_wp
    real(wp),parameter :: b2715  =    7.28147667468107583471326950926136115767612581862877764249646_wp
    real(wp),parameter :: b2716  =    -0.132699384612248379510571708176035274836827341616751884314074_wp
    real(wp),parameter :: b2717  =    -81.0799832525730726674679289752255240006070716633632990308935_wp
    real(wp),parameter :: b2718  =    -1.25037492835620639521768185656179119962253747492403205797494_wp
    real(wp),parameter :: b2719  =    2.59263594969543681023776379504377324994226447359296887778718_wp
    real(wp),parameter :: b2720  =    -0.301440298346404539830163997260526875264431537275641495291993_wp
    real(wp),parameter :: b2721  =    0.221384460789832337451706451572773791695246839057318414301020_wp
    real(wp),parameter :: b2722  =    0.0827577274771892931955989870974693152996276435429809890551210_wp
    real(wp),parameter :: b2723  =    18.9960662040611520464672450037243263998175161412237156872211_wp
    real(wp),parameter :: b2724  =    0.269231946409639685623468015128334167460051910348912845121977_wp
    real(wp),parameter :: b2725  =    1.62674827447066537462989364929628933988125029284183680279020_wp
    real(wp),parameter :: b2726  =    0.491719043846229147070666628704194097678081907210673044988866_wp
    real(wp),parameter :: b280   =     0.0838479812409052664616968791372814085980533139224911131069335_wp
    real(wp),parameter :: b281   =    0
    real(wp),parameter :: b282   =    0
    real(wp),parameter :: b283   =    0
    real(wp),parameter :: b284   =    0
    real(wp),parameter :: b285   =    -0.0117949367100973814319755056031295775367961960590736150777613_wp
    real(wp),parameter :: b286   =    -0.247299020568812652339473838743194598325992840353340132697498_wp
    real(wp),parameter :: b287   =    0.0978080858367729012259313014081291665503740655476733940756599_wp
    real(wp),parameter :: b288   =    0.217590689243420631360008651767860318344168120024782176879989_wp
    real(wp),parameter :: b289   =    0
    real(wp),parameter :: b2810  =    0.137585606763325224865659632196787746647447222975084865975440_wp
    real(wp),parameter :: b2811  =    0.0439870229715046685058790092341545026046103890294261359042581_wp
    real(wp),parameter :: b2812  =    0
    real(wp),parameter :: b2813  =    -0.513700813768193341957004456618630303738757363641964030086972_wp
    real(wp),parameter :: b2814  =    0.826355691151315508644211308399153458701423158616168576922372_wp
    real(wp),parameter :: b2815  =    25.7018139719811832625873882972519939511136556341960074626615_wp
    real(wp),parameter :: b2816  =    0
    real(wp),parameter :: b2817  =    0
    real(wp),parameter :: b2818  =    0
    real(wp),parameter :: b2819  =    0
    real(wp),parameter :: b2820  =    0
    real(wp),parameter :: b2821  =    0
    real(wp),parameter :: b2822  =    0
    real(wp),parameter :: b2823  =    -25.7018139719811832625873882972519939511136556341960074626615_wp
    real(wp),parameter :: b2824  =    -0.826355691151315508644211308399153458701423158616168576922372_wp
    real(wp),parameter :: b2825  =    0.513700813768193341957004456618630303738757363641964030086972_wp
    real(wp),parameter :: b2826  =    -0.0439870229715046685058790092341545026046103890294261359042581_wp
    real(wp),parameter :: b2827  =    -0.137585606763325224865659632196787746647447222975084865975440_wp
    real(wp),parameter :: b290   =     0.124380526654094412881516420868799316268491466359671423163289_wp
    real(wp),parameter :: b291   =    0
    real(wp),parameter :: b292   =    0
    real(wp),parameter :: b293   =    0
    real(wp),parameter :: b294   =    0.226120282197584301422238662979202901196752320742633143965145_wp
    real(wp),parameter :: b295   =    0.0137885887618080880607695837016477814530969417491493385363543_wp
    real(wp),parameter :: b296   =    -0.0672210133996684449749399507414305856950086341525382182856200_wp
    real(wp),parameter :: b297   =    0
    real(wp),parameter :: b298   =    0
    real(wp),parameter :: b299   =    -0.856238975085428354755349769879501772112121597411563802855067_wp
    real(wp),parameter :: b2910  =    -1.96337522866858908928262850028093813988180440518267404553576_wp
    real(wp),parameter :: b2911  =    -0.232332822724119401237246257308921847250108199230419994978218_wp
    real(wp),parameter :: b2912  =    0
    real(wp),parameter :: b2913  =    4.30660719086453349461668936876562947772432562053478092626764_wp
    real(wp),parameter :: b2914  =    -2.92722963249465482659787911202390446687687394950633612630592_wp
    real(wp),parameter :: b2915  =    -82.3131666397858944454492334105458707735761966428138676971041_wp
    real(wp),parameter :: b2916  =    0
    real(wp),parameter :: b2917  =    0
    real(wp),parameter :: b2918  =    0
    real(wp),parameter :: b2919  =    0
    real(wp),parameter :: b2920  =    0
    real(wp),parameter :: b2921  =    0
    real(wp),parameter :: b2922  =    0
    real(wp),parameter :: b2923  =    82.3131666397858944454492334105458707735761966428138676971041_wp
    real(wp),parameter :: b2924  =    2.92722963249465482659787911202390446687687394950633612630592_wp
    real(wp),parameter :: b2925  =    -4.30660719086453349461668936876562947772432562053478092626764_wp
    real(wp),parameter :: b2926  =    0.232332822724119401237246257308921847250108199230419994978218_wp
    real(wp),parameter :: b2927  =    1.96337522866858908928262850028093813988180440518267404553576_wp
    real(wp),parameter :: b2928  =    0.856238975085428354755349769879501772112121597411563802855067_wp
    real(wp),parameter :: b300   =     0.103484561636679776672993546511910344499744798201971316606663_wp
    real(wp),parameter :: b301   =    0
    real(wp),parameter :: b302   =    0
    real(wp),parameter :: b303   =    0.122068887306407222589644082868962077139592714834162134741275_wp
    real(wp),parameter :: b304   =    0.482574490331246622475134780125688112865919023850168049679402_wp
    real(wp),parameter :: b305   =    -0.0381409600015606999730886240005620205664113072478411477421970_wp
    real(wp),parameter :: b306   =    0
    real(wp),parameter :: b307   =    -0.550499525310802324138388507020508177411414311000037561712836_wp
    real(wp),parameter :: b308   =    0
    real(wp),parameter :: b309   =    -0.711915811585189227887648262043794387578291882406745570495765_wp
    real(wp),parameter :: b3010  =    -0.584129605671551340432988730158480872095335329645227595707052_wp
    real(wp),parameter :: b3011  =    0
    real(wp),parameter :: b3012  =    0
    real(wp),parameter :: b3013  =    2.11046308125864932128717300046622750300375054278936987850718_wp
    real(wp),parameter :: b3014  =    -0.0837494736739572135525742023001037992695260175335123517729291_wp
    real(wp),parameter :: b3015  =    5.10021499072320914075295969043344113107545060862804249161191_wp
    real(wp),parameter :: b3016  =    0
    real(wp),parameter :: b3017  =    0
    real(wp),parameter :: b3018  =    0
    real(wp),parameter :: b3019  =    0
    real(wp),parameter :: b3020  =    0
    real(wp),parameter :: b3021  =    0
    real(wp),parameter :: b3022  =    0
    real(wp),parameter :: b3023  =    -5.10021499072320914075295969043344113107545060862804249161191_wp
    real(wp),parameter :: b3024  =    0.0837494736739572135525742023001037992695260175335123517729291_wp
    real(wp),parameter :: b3025  =    -2.11046308125864932128717300046622750300375054278936987850718_wp
    real(wp),parameter :: b3026  =    0
    real(wp),parameter :: b3027  =    0.584129605671551340432988730158480872095335329645227595707052_wp
    real(wp),parameter :: b3028  =    0.711915811585189227887648262043794387578291882406745570495765_wp
    real(wp),parameter :: b3029  =    0.550499525310802324138388507020508177411414311000037561712836_wp
    real(wp),parameter :: b310   =     0.193333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: b311   =    0
    real(wp),parameter :: b312   =    0.220000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b313   =    -0.0800000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b314   =    0
    real(wp),parameter :: b315   =    0
    real(wp),parameter :: b316   =    0.109993425580724703919462404865068340845119058295846426463652_wp
    real(wp),parameter :: b317   =    -0.254297048076270161384068506997153122141835626976703920846242_wp
    real(wp),parameter :: b318   =    0
    real(wp),parameter :: b319   =    0.865570777116694254343770343821098281832847401233011859346737_wp
    real(wp),parameter :: b3110  =    3.32416449114093083106799552786572018336860092936986407160200_wp
    real(wp),parameter :: b3111  =    0
    real(wp),parameter :: b3112  =    0
    real(wp),parameter :: b3113  =    -12.0102223315977933882352385148661841260301942633996815127277_wp
    real(wp),parameter :: b3114  =    0.476601466242493239430442776862061899602963782003580209476163_wp
    real(wp),parameter :: b3115  =    -29.0243011221036390525802623213654099596251221332470910692353_wp
    real(wp),parameter :: b3116  =    0
    real(wp),parameter :: b3117  =    0
    real(wp),parameter :: b3118  =    0
    real(wp),parameter :: b3119  =    0
    real(wp),parameter :: b3120  =    0
    real(wp),parameter :: b3121  =    0
    real(wp),parameter :: b3122  =    0
    real(wp),parameter :: b3123  =    29.0243011221036390525802623213654099596251221332470910692353_wp
    real(wp),parameter :: b3124  =    -0.476601466242493239430442776862061899602963782003580209476163_wp
    real(wp),parameter :: b3125  =    12.0102223315977933882352385148661841260301942633996815127277_wp
    real(wp),parameter :: b3126  =    0
    real(wp),parameter :: b3127  =    -3.32416449114093083106799552786572018336860092936986407160200_wp
    real(wp),parameter :: b3128  =    -0.865570777116694254343770343821098281832847401233011859346737_wp
    real(wp),parameter :: b3129  =    0.254297048076270161384068506997153122141835626976703920846242_wp
    real(wp),parameter :: b3130  =    -0.109993425580724703919462404865068340845119058295846426463652_wp
    real(wp),parameter :: b320   =     -0.833333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: b321   =    1.38888888888888888888888888888888888888888888888888888888889_wp
    real(wp),parameter :: b322   =    0
    real(wp),parameter :: b323   =    0
    real(wp),parameter :: b324   =    -0.750000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b325   =    0
    real(wp),parameter :: b326   =    -0.492529543718026304422682049114021320200214681580657784719074_wp
    real(wp),parameter :: b327   =    0
    real(wp),parameter :: b328   =    0
    real(wp),parameter :: b329   =    0
    real(wp),parameter :: b3210  =    0
    real(wp),parameter :: b3211  =    0
    real(wp),parameter :: b3212  =    0
    real(wp),parameter :: b3213  =    0
    real(wp),parameter :: b3214  =    0
    real(wp),parameter :: b3215  =    0
    real(wp),parameter :: b3216  =    0
    real(wp),parameter :: b3217  =    0
    real(wp),parameter :: b3218  =    0
    real(wp),parameter :: b3219  =    0
    real(wp),parameter :: b3220  =    0
    real(wp),parameter :: b3221  =    0
    real(wp),parameter :: b3222  =    0
    real(wp),parameter :: b3223  =    0
    real(wp),parameter :: b3224  =    0
    real(wp),parameter :: b3225  =    0
    real(wp),parameter :: b3226  =    0
    real(wp),parameter :: b3227  =    0
    real(wp),parameter :: b3228  =    0
    real(wp),parameter :: b3229  =    0
    real(wp),parameter :: b3230  =    0.492529543718026304422682049114021320200214681580657784719074_wp
    real(wp),parameter :: b3231  =    0.750000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b330  =    0.111111111111111111111111111111111111111111111111111111111111_wp
    real(wp),parameter :: b331  =    0
    real(wp),parameter :: b332  =   -0.222222222222222222222222222222222222222222222222222222222222_wp
    real(wp),parameter :: b333  =    0
    real(wp),parameter :: b334  =    0
    real(wp),parameter :: b335  =    0
    real(wp),parameter :: b336  =    0
    real(wp),parameter :: b337  =    0
    real(wp),parameter :: b338  =    0
    real(wp),parameter :: b339  =    0
    real(wp),parameter :: b3310  =    0
    real(wp),parameter :: b3311  =    0
    real(wp),parameter :: b3312  =    0
    real(wp),parameter :: b3313  =    0
    real(wp),parameter :: b3314  =    0
    real(wp),parameter :: b3315  =    0
    real(wp),parameter :: b3316  =    0
    real(wp),parameter :: b3317  =    0
    real(wp),parameter :: b3318  =    0
    real(wp),parameter :: b3319  =    0
    real(wp),parameter :: b3320  =    0
    real(wp),parameter :: b3321  =    0
    real(wp),parameter :: b3322  =    0
    real(wp),parameter :: b3323  =    0
    real(wp),parameter :: b3324  =    0
    real(wp),parameter :: b3325  =    0
    real(wp),parameter :: b3326  =    0
    real(wp),parameter :: b3327  =    0
    real(wp),parameter :: b3328  =    0
    real(wp),parameter :: b3329  =    0
    real(wp),parameter :: b3330  =    0
    real(wp),parameter :: b3331  =    0
    real(wp),parameter :: b3332  =    0.222222222222222222222222222222222222222222222222222222222222_wp
    real(wp),parameter :: b340  =    0.285835140388971558796088842163836414852927537894596466840753_wp
    real(wp),parameter :: b341  =    0.291666666666666666666666666666666666666666666666666666666667_wp
    real(wp),parameter :: b342  =    0.218750000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b343  =    0
    real(wp),parameter :: b344  =    0.164062500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b345  =    0
    real(wp),parameter :: b346  =    0.218194354945556658327188241581352107093288824322187941141516_wp
    real(wp),parameter :: b347  =    0.180392898478697766863635221946775437719620053641849228562435_wp
    real(wp),parameter :: b348  =    0
    real(wp),parameter :: b349  =    0.205713839404845018859120755122929542277570094982808905393991_wp
    real(wp),parameter :: b3410  =    0.242715791581770239970282927959446515762745971386670541948576_wp
    real(wp),parameter :: b3411  =    0.246465780813629305833609291181891407799228103869305705137021_wp
    real(wp),parameter :: b3412  =   -3.44991940790890824979834154601622662060370460614931644223924_wp
    real(wp),parameter :: b3413  =    0.228875562160036081760729060738458584294220372552740218459295_wp
    real(wp),parameter :: b3414  =    0.283290599702151415321527419056733335978436595493855789831434_wp
    real(wp),parameter :: b3415  =    3.21085125837766640960131490544236787005557320332238705967955_wp
    real(wp),parameter :: b3416  =   -0.223538777364845699920233756214162507964125230083674032084065_wp
    real(wp),parameter :: b3417  =   -0.707121157204419073518727286207487212130091231955206160635271_wp
    real(wp),parameter :: b3418  =    3.21123345150287080408174729202856500893260034443022374267639_wp
    real(wp),parameter :: b3419  =    1.40954348309669766030414474301123175769045945573548986335553_wp
    real(wp),parameter :: b3420  =   -0.151362053443742613121602276742518111090963026203676055891793_wp
    real(wp),parameter :: b3421  =    0.372350574527014276454724080214619984397121028202148298716575_wp
    real(wp),parameter :: b3422  =    0.252978746406361336722199907762141285915775728129414319261111_wp
    real(wp),parameter :: b3423  =   -3.21085125837766640960131490544236787005557320332238705967955_wp
    real(wp),parameter :: b3424  =   -0.283290599702151415321527419056733335978436595493855789831434_wp
    real(wp),parameter :: b3425  =   -0.228875562160036081760729060738458584294220372552740218459295_wp
    real(wp),parameter :: b3426  =   -0.246465780813629305833609291181891407799228103869305705137021_wp
    real(wp),parameter :: b3427  =   -0.242715791581770239970282927959446515762745971386670541948576_wp
    real(wp),parameter :: b3428  =   -0.205713839404845018859120755122929542277570094982808905393991_wp
    real(wp),parameter :: b3429  =   -0.180392898478697766863635221946775437719620053641849228562435_wp
    real(wp),parameter :: b3430  =   -0.218194354945556658327188241581352107093288824322187941141516_wp
    real(wp),parameter :: b3431  =   -0.164062500000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b3432  =   -0.218750000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: b3433  =   -0.291666666666666666666666666666666666666666666666666666666667_wp

    associate (f0 => me%funcs(:,1), &
               f1 => me%funcs(:,2), &
               f2 => me%funcs(:,3), &
               f3 => me%funcs(:,4), &
               f4 => me%funcs(:,5), &
               f5 => me%funcs(:,6), &
               f6 => me%funcs(:,7), &
               f7 => me%funcs(:,8), &
               f8 => me%funcs(:,9), &
               f9  => me%funcs(:,10), &
               f10 => me%funcs(:,11), &
               f11 => me%funcs(:,12), &
               f12 => me%funcs(:,13), &
               f13 => me%funcs(:,14), &
               f14 => me%funcs(:,15), &
               f15 => me%funcs(:,16), &
               f16 => me%funcs(:,17), &
               f17 => me%funcs(:,18), &
               f18 => me%funcs(:,19), &
               f19 => me%funcs(:,20), &
               f20 => me%funcs(:,21), &
               f21 => me%funcs(:,22), &
               f22 => me%funcs(:,23), &
               f23 => me%funcs(:,24), &
               f24 => me%funcs(:,25), &
               f25 => me%funcs(:,26), &
               f26 => me%funcs(:,27), &
               f27 => me%funcs(:,28), &
               f28 => me%funcs(:,29), &
               f29 => me%funcs(:,30), &
               f30 => me%funcs(:,31), &
               f31 => me%funcs(:,32), &
               f32 => me%funcs(:,33), &
               f33 => me%funcs(:,34), &
               f34 => me%funcs(:,35))

        call me%f(t+a0*h,  x,f0)
        call me%f(t+a1*h,  x+h*(b10*f0),f1)
        call me%f(t+a2*h,  x+h*(b20*f0  + b21 *f1),f2)
        call me%f(t+a3*h,  x+h*(b30*f0  + b31 *f1 + b32*f2),f3)
        call me%f(t+a4*h,  x+h*(b40*f0  + b41 *f1 + b42*f2  + b43*f3),f4)
        call me%f(t+h,     x+h*(b50*f0  + b51 *f1 + b52*f2  + b53*f3  + &
                                b54*f4),f5)
        call me%f(t+a6*h,  x+h*(b60*f0  + b61 *f1 + b62*f2  + b63*f3  + &
                                b64*f4  + b65*f5),f6)
        call me%f(t+a7*h,  x+h*(b70*f0  + b71 *f1 + b72*f2  + b73*f3  + &
                                b74*f4  + b75*f5  + b76*f6),f7)
        call me%f(t+a8*h,  x+h*(b80*f0  + b81 *f1 + b82*f2  + b83*f3  + &
                                b84*f4  + b85*f5  + b86*f6  + b87*f7),f8)
        call me%f(t+a9*h,  x+h*(b90*f0  + b91 *f1 + b92*f2  + b93*f3  + &
                                b94*f4  + b95*f5  + b96*f6  + b97*f7  + &
                                b98*f8),f9)
        call me%f(t+a10*h, x+h*(b100*f0 + b101*f1 + b102*f2 + b103*f3 + &
                                b104*f4 + b105*f5 + b106*f6 + b107*f7 + &
                                b108*f8 + b109*f9),f10)
        call me%f(t+a11*h, x+h*(b110*f0 + b111*f1 + b112*f2 + b113*f3 + &
                                b114*f4 + b115*f5 + b116*f6 + b117*f7 + &
                                b118*f8 + b119*f9 + b1110*f10),f11)
        call me%f(t+a12*h, x+h*(b120*f0 + b121*f1 + b122*f2 + b123*f3 + &
                                b124*f4 + b125*f5 + b126*f6 + b127*f7 + &
                                b128*f8 + b129*f9 + b1210*f10 + b1211*f11),f12)
        call me%f(t+a13*h, x+h*(b130*f0 + b131*f1 + b132*f2 + b133*f3 + &
                                b134*f4 + b135*f5 + b136*f6 + b137*f7 + &
                                b138*f8 + b139*f9 + b1310*f10 + b1311*f11 + &
                                b1312*f12),f13)
        call me%f(t+a14*h, x+h*(b140*f0 + b141*f1 + b142*f2 + b143*f3 + &
                                b144*f4 + b145*f5 + b146*f6 + b147*f7 + &
                                b148*f8 + b149*f9 + b1410*f10 + b1411*f11 + &
                                b1412*f12 + b1413*f13),f14)
        call me%f(t+a15*h, x+h*(b150*f0 + b151*f1 + b152*f2 + b153*f3 + &
                                b154*f4 + b155*f5 + b156*f6 + b157*f7 + &
                                b158*f8 + b159*f9 + b1510*f10 + b1511*f11 + &
                                b1512*f12 + b1513*f13 + b1514*f14),f15)
        call me%f(t+a16*h, x+h*(b160*f0 + b161*f1 + b162*f2 + b163*f3 + &
                                b164*f4 + b165*f5 + b166*f6 + b167*f7 + &
                                b168*f8 + b169*f9 + b1610*f10 + b1611*f11 + &
                                b1612*f12 + b1613*f13 + b1614*f14 + b1615*f15),f16)
        call me%f(t+a17*h, x+h*(b170*f0 + b171*f1 + b172*f2 + b173*f3 + &
                                b174*f4 + b175*f5 + b176*f6 + b177*f7 + &
                                b178*f8 + b179*f9 + b1710*f10 + b1711*f11 + &
                                b1712*f12 + b1713*f13 + b1714*f14 + b1715*f15 + &
                                b1716*f16),f17)
        call me%f(t+a18*h, x+h*(b180*f0 + b181*f1 + b182*f2 + b183*f3 + &
                                b184*f4 + b185*f5 + b186*f6 + b187*f7 + &
                                b188*f8 + b189*f9 + b1810*f10 + b1811*f11 + &
                                b1812*f12 + b1813*f13 + b1814*f14 + b1815*f15 + &
                                b1816*f16 + b1817*f17),f18)
        call me%f(t+a19*h, x+h*(b190*f0 + b191*f1 + b192*f2 + b193*f3 + &
                                b194*f4 + b195*f5 + b196*f6 + b197*f7 + &
                                b198*f8 + b199*f9 + b1910*f10 + b1911*f11 + &
                                b1912*f12 + b1913*f13 + b1914*f14 + b1915*f15 + &
                                b1916*f16 + b1917*f17 + b1918*f18),f19)
        call me%f(t+a20*h, x+h*(b200*f0 + b201*f1 + b202*f2 + b203*f3 + &
                                b204*f4 + b205*f5 + b206*f6 + b207*f7 + &
                                b208*f8 + b209*f9 + b2010*f10 + b2011*f11 + &
                                b2012*f12 + b2013*f13 + b2014*f14 + b2015*f15 + &
                                b2016*f16 + b2017*f17 + b2018*f18 + b2019*f19),f20)
        call me%f(t+a21*h, x+h*(b210*f0 + b211*f1 + b212*f2 + b213*f3 + &
                                b214*f4 + b215*f5 + b216*f6 + b217*f7 + &
                                b218*f8 + b219*f9 + b2110*f10 + b2111*f11 + &
                                b2112*f12 + b2113*f13 + b2114*f14 + b2115*f15 + &
                                b2116*f16 + b2117*f17 + b2118*f18 + b2119*f19 + &
                                b2120*f20),f21)
        call me%f(t+a22*h, x+h*(b220*f0 + b221*f1 + b222*f2 + b223*f3 + &
                                b224*f4 + b225*f5 + b226*f6 + b227*f7 + &
                                b228*f8 + b229*f9 + b2210*f10 + b2211*f11 + &
                                b2212*f12 + b2213*f13 + b2214*f14 + b2215*f15 + &
                                b2216*f16 + b2217*f17 + b2218*f18 + b2219*f19 + &
                                b2220*f20 + b2221*f21),f22)
        call me%f(t+a23*h, x+h*(b230*f0 + b231*f1 + b232*f2 + b233*f3 + &
                                b234*f4 + b235*f5 + b236*f6 + b237*f7 + &
                                b238*f8 + b239*f9 + b2310*f10 + b2311*f11 + &
                                b2312*f12 + b2313*f13 + b2314*f14 + b2315*f15 + &
                                b2316*f16 + b2317*f17 + b2318*f18 + b2319*f19 + &
                                b2320*f20 + b2321*f21 + b2322*f22),f23)
        call me%f(t+a24*h, x+h*(b240*f0 + b241*f1 + b242*f2 + b243*f3 + &
                                b244*f4 + b245*f5 + b246*f6 + b247*f7 + &
                                b248*f8 + b249*f9 + b2410*f10 + b2411*f11 + &
                                b2412*f12 + b2413*f13 + b2414*f14 + b2415*f15 + &
                                b2416*f16 + b2417*f17 + b2418*f18 + b2419*f19 + &
                                b2420*f20 + b2421*f21 + b2422*f22 + b2423*f23),f24)
        call me%f(t+a25*h, x+h*(b250*f0 + b251*f1 + b252*f2 + b253*f3 + &
                                b254*f4 + b255*f5 + b256*f6 + b257*f7 + &
                                b258*f8 + b259*f9 + b2510*f10 + b2511*f11 + &
                                b2512*f12 + b2513*f13 + b2514*f14 + b2515*f15 + &
                                b2516*f16 + b2517*f17 + b2518*f18 + b2519*f19 + &
                                b2520*f20 + b2521*f21 + b2522*f22 + b2523*f23 + &
                                b2524*f24),f25)
        call me%f(t+a26*h, x+h*(b260*f0 + b261*f1 + b262*f2 + b263*f3 + &
                                b264*f4 + b265*f5 + b266*f6 + b267*f7 + &
                                b268*f8 + b269*f9 + b2610*f10 + b2611*f11 + &
                                b2612*f12 + b2613*f13 + b2614*f14 + b2615*f15 + &
                                b2616*f16 + b2617*f17 + b2618*f18 + b2619*f19 + &
                                b2620*f20 + b2621*f21 + b2622*f22 + b2623*f23 + &
                                b2624*f24 + b2625*f25),f26)
        call me%f(t+a27*h, x+h*(b270*f0 + b271*f1 + b272*f2 + b273*f3 + &
                                b274*f4 + b275*f5 + b276*f6 + b277*f7 + &
                                b278*f8 + b279*f9 + b2710*f10 + b2711*f11 + &
                                b2712*f12 + b2713*f13 + b2714*f14 + b2715*f15 + &
                                b2716*f16 + b2717*f17 + b2718*f18 + b2719*f19 + &
                                b2720*f20 + b2721*f21 + b2722*f22 + b2723*f23 + &
                                b2724*f24 + b2725*f25 + b2726*f26),f27)
        call me%f(t+a28*h, x+h*(b280*f0 + b281*f1 + b282*f2 + b283*f3 + &
                                b284*f4 + b285*f5 + b286*f6 + b287*f7 + &
                                b288*f8 + b289*f9 + b2810*f10 + b2811*f11 + &
                                b2812*f12 + b2813*f13 + b2814*f14 + b2815*f15 + &
                                b2816*f16 + b2817*f17 + b2818*f18 + b2819*f19 + &
                                b2820*f20 + b2821*f21 + b2822*f22 + b2823*f23 + &
                                b2824*f24 + b2825*f25 + b2826*f26 + b2827*f27),f28)
        call me%f(t+a29*h, x+h*(b290*f0 + b291*f1 + b292*f2 + b293*f3 + &
                                b294*f4 + b295*f5 + b296*f6 + b297*f7 + &
                                b298*f8 + b299*f9 + b2910*f10 + b2911*f11 + &
                                b2912*f12 + b2913*f13 + b2914*f14 + b2915*f15 + &
                                b2916*f16 + b2917*f17 + b2918*f18 + b2919*f19 + &
                                b2920*f20 + b2921*f21 + b2922*f22 + b2923*f23 + &
                                b2924*f24 + b2925*f25 + b2926*f26 + b2927*f27 + &
                                b2928*f28),f29)
        call me%f(t+a30*h, x+h*(b300*f0 + b301*f1 + b302*f2 + b303*f3 + &
                                b304*f4 + b305*f5 + b306*f6 + b307*f7 + &
                                b308*f8 + b309*f9 + b3010*f10 + b3011*f11 + &
                                b3012*f12 + b3013*f13 + b3014*f14 + b3015*f15 + &
                                b3016*f16 + b3017*f17 + b3018*f18 + b3019*f19 + &
                                b3020*f20 + b3021*f21 + b3022*f22 + b3023*f23 + &
                                b3024*f24 + b3025*f25 + b3026*f26 + b3027*f27 + &
                                b3028*f28 + b3029*f29),f30)
        call me%f(t+a31*h, x+h*(b310*f0 + b311*f1 + b312*f2 + b313*f3 + &
                                b314*f4 + b315*f5 + b316*f6 + b317*f7 + &
                                b318*f8 + b319*f9 + b3110*f10 + b3111*f11 + &
                                b3112*f12 + b3113*f13 + b3114*f14 + b3115*f15 + &
                                b3116*f16 + b3117*f17 + b3118*f18 + b3119*f19 + &
                                b3120*f20 + b3121*f21 + b3122*f22 + b3123*f23 + &
                                b3124*f24 + b3125*f25 + b3126*f26 + b3127*f27 + &
                                b3128*f28 + b3129*f29 + b3130*f30),f31)
        call me%f(t+a32*h, x+h*(b320*f0 + b321*f1 + b322*f2 + b323*f3 + &
                                b324*f4 + b325*f5 + b326*f6 + b327*f7 + &
                                b328*f8 + b329*f9 + b3210*f10 + b3211*f11 + &
                                b3212*f12 + b3213*f13 + b3214*f14 + b3215*f15 + &
                                b3216*f16 + b3217*f17 + b3218*f18 + b3219*f19 + &
                                b3220*f20 + b3221*f21 + b3222*f22 + b3223*f23 + &
                                b3224*f24 + b3225*f25 + b3226*f26 + b3227*f27 + &
                                b3228*f28 + b3229*f29 + b3230*f30 + b3231*f31),f32)
        call me%f(t+a33*h, x+h*(b330*f0 + b331*f1 + b332*f2 + b333*f3 + &
                                b334*f4 + b335*f5 + b336*f6 + b337*f7 + &
                                b338*f8 + b339*f9 + b3310*f10 + b3311*f11 + &
                                b3312*f12 + b3313*f13 + b3314*f14 + b3315*f15 + &
                                b3316*f16 + b3317*f17 + b3318*f18 + b3319*f19 + &
                                b3320*f20 + b3321*f21 + b3322*f22 + b3323*f23 + &
                                b3324*f24 + b3325*f25 + b3326*f26 + b3327*f27 + &
                                b3328*f28 + b3329*f29 + b3330*f30 + b3331*f31 + &
                                b3332*f32),f33)
        call me%f(t+h,     x+h*(b340*f0 + b341*f1 + b342*f2 + b343*f3 + &
                                b344*f4 + b345*f5 + b346*f6 + b347*f7 + &
                                b348*f8 + b349*f9 + b3410*f10 + b3411*f11 + &
                                b3412*f12 + b3413*f13 + b3414*f14 + b3415*f15 + &
                                b3416*f16 + b3417*f17 + b3418*f18 + b3419*f19 + &
                                b3420*f20 + b3421*f21 + b3422*f22 + b3423*f23 + &
                                b3424*f24 + b3425*f25 + b3426*f26 + b3427*f27 + &
                                b3428*f28 + b3429*f29 + b3430*f30 + b3431*f31 + &
                                b3432*f32 + b3433*f33),f34)

        xf = x+h*(  c0*f0   + &
                    c1*f1   + &
                    c2*f2   + &
                    c3*f3   + &
                    c4*f4   + &
                    c5*f5   + &
                    c6*f6   + &
                    c7*f7   + &
                    c8*f8   + &
                    c9*f9   + &
                    c10*f10 + &
                    c11*f11 + &
                    c12*f12 + &
                    c13*f13 + &
                    c14*f14 + &
                    c15*f15 + &
                    c16*f16 + &
                    c17*f17 + &
                    c18*f18 + &
                    c19*f19 + &
                    c20*f20 + &
                    c21*f21 + &
                    c22*f22 + &
                    c23*f23 + &
                    c24*f24 + &
                    c25*f25 + &
                    c26*f26 + &
                    c27*f27 + &
                    c28*f28 + &
                    c29*f29 + &
                    c30*f30 + &
                    c31*f31 + &
                    c32*f32 + &
                    c33*f33 + &
                    c34*f34 )

        xerr = (1.0_wp/1000.0_wp)*h*(f1-f33)

    end associate

    end procedure rkf1412
!*****************************************************************************************

!*****************************************************************************************
!>
!  A modification of Hiroshi Ono's order 12 Runge-Kutta scheme with an order 9 embedded scheme.
!
!### Reference
!  * Hiroshi Ono, On the 25 stage 12th order explicit Runge-Kutta method,
!    Transactions of the Japan Society for Industrial and applied Mathematics,
!    Vol. 6, No. 3, (2006) pages 177 to 186.
!  * [Rational coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK12/RKcoeff12h(9)_1.pdf)

    module procedure rko129

    real(wp),parameter :: a2 = .4351851851851851851851851851851851851851851851851851851851851851851851851851851851852_wp
    real(wp),parameter :: a3 = .4429824561403508771929824561403508771929824561403508771929824561403508771929824561404_wp
    real(wp),parameter :: a4 = .6644736842105263157894736842105263157894736842105263157894736842105263157894736842105_wp
    real(wp),parameter :: a5 = .1069403994175161223216143124609943831911795298522987310172664863740378614520490950697_wp
    real(wp),parameter :: a6 = .1644736842105263157894736842105263157894736842105263157894736842105263157894736842105_wp
    real(wp),parameter :: a7 = .5843251088534107402031930333817126269956458635703918722786647314949201741654571843251_wp
    real(wp),parameter :: a8 = .6382358235823582358235823582358235823582358235823582358235823582358235823582358235824e-1_wp
    real(wp),parameter :: a9 = .2_wp
    real(wp),parameter :: a10 = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a11 = .9446116054065563496847375720237340591364440847545279548637407620203329280139199114088_wp
    real(wp),parameter :: a12 = .5179584680428461035835615460185207449756996828425551866437683588600807071775316523299e-1_wp
    real(wp),parameter :: a13 = .8488805186071653506398389301626743020641481756400195420459339398355773991365476236893e-1_wp
    real(wp),parameter :: a14 = .2655756032646428930981140590456168352972012641640776214486652703185222349414361456016_wp
    real(wp),parameter :: a15 = .5_wp
    real(wp),parameter :: a16 = .7344243967353571069018859409543831647027987358359223785513347296814777650585638543984_wp
    real(wp),parameter :: a17 = .9151119481392834649360161069837325697935851824359980457954066060164422600863452376311_wp
    real(wp),parameter :: a18 = .9446116054065563496847375720237340591364440847545279548637407620203329280139199114088_wp
    real(wp),parameter :: a19 = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333_wp
    real(wp),parameter :: a20 = .2_wp
    real(wp),parameter :: a21 = .5843251088534107402031930333817126269956458635703918722786647314949201741654571843251_wp
    real(wp),parameter :: a22 = .1644736842105263157894736842105263157894736842105263157894736842105263157894736842105_wp
    real(wp),parameter :: a23 = .4429824561403508771929824561403508771929824561403508771929824561403508771929824561404_wp
    real(wp),parameter :: a24 = .4351851851851851851851851851851851851851851851851851851851851851851851851851851851852_wp
    real(wp),parameter :: a26 = .4970267001007476028032930885363848318550815967070143979104342007017543082576391973337_wp
    real(wp),parameter :: a27 = .8043478260869565217391304347826086956521739130434782608695652173913043478260869565217_wp
    real(wp),parameter :: a28 = .8717948717948717948717948717948717948717948717948717948717948717948717948717948717949_wp

    real(wp),parameter :: b21 = .4351851851851851851851851851851851851851851851851851851851851851851851851851851851852_wp
    real(wp),parameter :: b31 = .2175227402212137286104398734798923400326123258875071216675507357419304139407870179368_wp
    real(wp),parameter :: b32 = .2254597159191371485825425826604585371603701302528437555254317203984204632521954382036_wp
    real(wp),parameter :: b41 = .1661184210526315789473684210526315789473684210526315789473684210526315789473684210526_wp
    real(wp),parameter :: b42 = 0
    real(wp),parameter :: b43 = .4983552631578947368421052631578947368421052631578947368421052631578947368421052631579_wp
    real(wp),parameter :: b51 = .8681163193918508865080629713483674348617570159685880673878476405659784672427434203840e-1_wp
    real(wp),parameter :: b52 = 0
    real(wp),parameter :: b53 = .3456981948164897000739691452644149026866903085576649936607708344667115529230058847943e-1_wp
    real(wp),parameter :: b54 = -.1444105200331793633658889920028385056366520260032657508759536112923114056452583544813e-1_wp
    real(wp),parameter :: b61 = .3850951504524952575244726520324404860898823844950484490138343541601000605579674908780e-1_wp
    real(wp),parameter :: b62 = 0
    real(wp),parameter :: b63 = 0
    real(wp),parameter :: b64 = .9889604363651382462812798900870548282288802143553903392828423536490053424426345577709e-4_wp
    real(wp),parameter :: b65 = .1258652731216402762123982910182735616976625577395859318541619645591514091994326716669_wp
    real(wp),parameter :: b71 = .5247404461891304721708365639844356354386624470737997647397329116437616959059854701778_wp
    real(wp),parameter :: b72 = 0
    real(wp),parameter :: b73 = 0
    real(wp),parameter :: b74 = .7610651429965941990931946190296373660533922130208900643092560661920652146268815424885e-1_wp
    real(wp),parameter :: b75 = -2.135538596825204678401302463372339475781988628449683891486081018353295131152366999963_wp
    real(wp),parameter :: b76 = 2.119016745189825526524339470866652730733632823644186992594087231585247087949150559862_wp
    real(wp),parameter :: b81 = .3572122856624484555412656871410309454657254402954350529356994466234212462699695879099e-1_wp
    real(wp),parameter :: b82 = 0
    real(wp),parameter :: b83 = 0
    real(wp),parameter :: b84 = 0
    real(wp),parameter :: b85 = .4596205641509305161863483311675360886293026632692961086374314175410528945129931429031e-1_wp
    real(wp),parameter :: b86 = -.1800016957713219812376781513215261066354472582005125759344515744644655947135707070493e-1_wp
    real(wp),parameter :: b87 = .1404669540301245333646491248782654898654978218139650184903068535815036288843799818693e-3_wp
    real(wp),parameter :: b91 = .1888176809184108413680419706237505704281864015973418413153299364663524364341955734913e-1_wp
    real(wp),parameter :: b92 = 0
    real(wp),parameter :: b93 = 0
    real(wp),parameter :: b94 = 0
    real(wp),parameter :: b95 = 0
    real(wp),parameter :: b96 = .8381198329740939383758467268684128106638290469173528975872944409244863365387981657651e-1_wp
    real(wp),parameter :: b97 = .9031585241436450659339626877786556727134244569192456020771509670971435444573825121328e-5_wp
    real(wp),parameter :: b98 = .9729721702550808557495179062390587533407132090396133365371679075124515126725605224924e-1_wp
    real(wp),parameter :: b101 = -.4080067703469846387001599007585094416645285179656696882006645696781274941696191338570e-1_wp
    real(wp),parameter :: b102 = 0
    real(wp),parameter :: b103 = 0
    real(wp),parameter :: b104 = 0
    real(wp),parameter :: b105 = 0
    real(wp),parameter :: b106 = -.6005539308646711394173617417282635784603029781694346062080009050806164919590318269156_wp
    real(wp),parameter :: b107 = .1585222658367901823092810289995660365558277114283363946284113044047141269675809741481e-2_wp
    real(wp),parameter :: b108 = .3026658851734428583799288358366834830670493200050903757742452562250879940584583879948_wp
    real(wp),parameter :: b109 = .6704368334008921764176894190107687125274815661799611686408713261126274393811928758984_wp
    real(wp),parameter :: b111 = 6.344326927733666873696263058980778771328558370914152084380122688263116113189631874466_wp
    real(wp),parameter :: b112 = 0
    real(wp),parameter :: b113 = 0
    real(wp),parameter :: b114 = 0
    real(wp),parameter :: b115 = 0
    real(wp),parameter :: b116 = 0
    real(wp),parameter :: b117 = 1.975263319684766813850955500728188106557880392559127720135666889111189136253813975729_wp
    real(wp),parameter :: b118 = -13.82822337504897849061527131653164721900473695395858210726164571154383357718661916508_wp
    real(wp),parameter :: b119 = 14.82423926991066347292103944504014293055447712090877566105421769783705551827529470930_wp
    real(wp),parameter :: b1110 = -8.370994536873562320168249116193728530299734845668945403444620801647194262518201483006_wp
    real(wp),parameter :: b121 = -.9910783781470375735018424706990342484478259963427121190616343432259493286903501070750e-1_wp
    real(wp),parameter :: b122 = 0
    real(wp),parameter :: b123 = 0
    real(wp),parameter :: b124 = 0
    real(wp),parameter :: b125 = 0
    real(wp),parameter :: b126 = -1.046319958132641377892211373672215159710875883372049858922517534871447317958500644882_wp
    real(wp),parameter :: b127 = -.4662578801256029967826301701912208605727776822301916928051722653634090636623902615221e-3_wp
    real(wp),parameter :: b128 = .4243518408197806737520879584209906999975132950474577957967837636354467847566731274758_wp
    real(wp),parameter :: b129 = .8350767448298191909277817680995907515365477233544915021061336714983604939879298121048_wp
    real(wp),parameter :: b1210 = -.6204379562043795620437956204379562043795620437956204379562043795620437956204379562044e-1_wp
    real(wp),parameter :: b1211 = .3051106025934401220442410373760488176964149504195270785659801678108314263920671243326e-3_wp
    real(wp),parameter :: b131 = .1731635805893703684448829089734942070445562382116745420114094326722705203174775140332e-1_wp
    real(wp),parameter :: b132 = 0
    real(wp),parameter :: b133 = 0
    real(wp),parameter :: b134 = 0
    real(wp),parameter :: b135 = 0
    real(wp),parameter :: b136 = 0
    real(wp),parameter :: b137 = 0
    real(wp),parameter :: b138 = 0
    real(wp),parameter :: b139 = .8690027159880925811795093256330430220184719404902985010762237927304116502630885425916e-3_wp
    real(wp),parameter :: b1310 = -.9198044746158460357290146202828770811642343120551843217420771155968903404182588708338e-4_wp
    real(wp),parameter :: b1311 = .1833594777349928706130754092439075395240645268074916590065144720493485758190652981058e-6_wp
    real(wp),parameter :: b1312 = .6679448817377525524901838117990401028051762116902291244289142812068791591710992924480e-1_wp
    real(wp),parameter :: b141 = .1497702285787817250603134142375551527881374815295532946701659971640303704469333256722e-1_wp
    real(wp),parameter :: b142 = 0
    real(wp),parameter :: b143 = 0
    real(wp),parameter :: b144 = 0
    real(wp),parameter :: b145 = 0
    real(wp),parameter :: b146 = 0
    real(wp),parameter :: b147 = 0
    real(wp),parameter :: b148 = 0
    real(wp),parameter :: b149 = .1299053198125687130806056561360945267563487713700848629960379786258472215118648828622_wp
    real(wp),parameter :: b1410 = .4444252090193421415880477613428889941728587016143899156135432416265141679149665438354e-2_wp
    real(wp),parameter :: b1411 = -.2185380330434085597487818850115824768660155325792369439283645364044817815712239904765e-5_wp
    real(wp),parameter :: b1412 = .6235770138025566019351747531432668514647398310437277172705875380249580436468314885283e-1_wp
    real(wp),parameter :: b1413 = .5389349250407735998767659637686133399860483467584655047185578940287507515886082812093e-1_wp
    real(wp),parameter :: b151 = .1956886038861791180930007029621719896471785620210633245797118818081004119873004962774_wp
    real(wp),parameter :: b152 = 0
    real(wp),parameter :: b153 = 0
    real(wp),parameter :: b154 = 0
    real(wp),parameter :: b155 = 0
    real(wp),parameter :: b156 = 0
    real(wp),parameter :: b157 = 0
    real(wp),parameter :: b158 = 0
    real(wp),parameter :: b159 = 1.132878041352190960246053936692850204846641933665191609652843227265803479941175715645_wp
    real(wp),parameter :: b1510 = .9754686357770750130309903534180324613492470735016258589778436907188769670738143543989_wp
    real(wp),parameter :: b1511 = .4607633555393391875625052489818876130584589587566695260680043303699764560744838206512e-3_wp
    real(wp),parameter :: b1512 = -.4833419929869440350491944569648386990066349467457051943375195864364996227787596872297_wp
    real(wp),parameter :: b1513 = .2837255261146423361717978265203569959020643630988946562818139270625464117810614651744_wp
    real(wp),parameter :: b1514 = -1.604879577498682731680210867877554840351555444499826924680761144749197624460666828087_wp
    real(wp),parameter :: b161 = -.6628581660952109109950216641584960475848453023977418371554567406302394052049642993363_wp
    real(wp),parameter :: b162 = 0
    real(wp),parameter :: b163 = 0
    real(wp),parameter :: b164 = 0
    real(wp),parameter :: b165 = 0
    real(wp),parameter :: b166 = 0
    real(wp),parameter :: b167 = 0
    real(wp),parameter :: b168 = 0
    real(wp),parameter :: b169 = -5.301219753823165631138506283750981613687493380976143446419698251739201243786272433269_wp
    real(wp),parameter :: b1610 = -5.493744530005151771950832780438489668352801158207831962667055860430416893664761589330_wp
    real(wp),parameter :: b1611 = .6448107716343851659685115053410322161204197332247306929859781244843876115320017653542e-2_wp
    real(wp),parameter :: b1612 = 2.226911096857986220827171118161838780361655511952909244621862410327793805385281866476_wp
    real(wp),parameter :: b1613 = -.8260094546883369994121948528168152647220748853112540011094140493170779782006783468229_wp
    real(wp),parameter :: b1614 = 9.736973734199539440165209948098457305642565726332220268712929325968631361580234056276_wp
    real(wp),parameter :: b1615 = 1.047923362573352907746375340805459350884588027111516805638308114257144242834404582751_wp
    real(wp),parameter :: b171 = 9.451896878619703179615895147719645120861585484250220485534769741441232402797149004402_wp
    real(wp),parameter :: b172 = 0
    real(wp),parameter :: b173 = 0
    real(wp),parameter :: b174 = 0
    real(wp),parameter :: b175 = 0
    real(wp),parameter :: b176 = 0
    real(wp),parameter :: b177 = 0
    real(wp),parameter :: b178 = 0
    real(wp),parameter :: b179 = 74.07837676951819392708315469048815926453257791083353986411194358002277837032529121349_wp
    real(wp),parameter :: b1710 = 80.08971633421252928663933464653279322909492698803384552922422585968931770071180510457_wp
    real(wp),parameter :: b1711 = -.1241702484260160323471778040612055855620194337087768444994146945534716737633295008233_wp
    real(wp),parameter :: b1712 = -32.04108125365225923900659816438329678215226279318337501510105142207810131441066638158_wp
    real(wp),parameter :: b1713 = 15.51919421000708125838831543558938623465692911684926717234376165413330242668122537342_wp
    real(wp),parameter :: b1714 = -136.4444237346563024309541869705527862321862115141467999073204884062451779712319267938_wp
    real(wp),parameter :: b1715 = -11.36109896858298168444633695622136642878717316645309025294277498926997344146696251530_wp
    real(wp),parameter :: b1716 = 1.746701961099335199963616081872403749335232589961167014444435282876535760443759733215_wp
    real(wp),parameter :: b181 = 1.059086740089530572084686602455620391446019980075022020682291957643879144548290466243_wp
    real(wp),parameter :: b182 = 0
    real(wp),parameter :: b183 = 0
    real(wp),parameter :: b184 = 0
    real(wp),parameter :: b185 = 0
    real(wp),parameter :: b186 = 0
    real(wp),parameter :: b187 = 1.975263319684766813850955500728188106557880392559127720135666889111189136253813975729_wp
    real(wp),parameter :: b188 = -13.82822337504897849061527131653164721900473695395858210726164571154383357718661916508_wp
    real(wp),parameter :: b189 = -26.72676722061492102175215840759265382960786498230207636418025793583833107088763635551_wp
    real(wp),parameter :: b1810 = -53.49798986553787152824955918158085193982743760760492627590775850823540448031315562431_wp
    real(wp),parameter :: b1811 = .7179812487812968675250181644918743516316099586081748670146732996043783143974128607587e-1_wp
    real(wp),parameter :: b1812 = 18.05559723664965139501142471462605880255358278669785143935434349294122071206578011157_wp
    real(wp),parameter :: b1813 = -8.765163819232982113792324392383185125858432915348234208512916467037776994687248464210_wp
    real(wp),parameter :: b1814 = 76.87522358555575651988028075659565616302805855797883225625477685501643635288672453023_wp
    real(wp),parameter :: b1815 = 6.541007260506904941709851688460080123863085658071943624909334991972860537373864701434_wp
    real(wp),parameter :: b1816 = -.8461837296739539683436767973109554382501241857205608248520628152174372445350689400749_wp
    real(wp),parameter :: b1817 = .3096334815052354314802658810823658907325235844531318754050068324709258105543338930304e-1_wp
    real(wp),parameter :: b191 = -.9776673260040047976710026469624134014469714141641816253185153958371390334954366206102e-1_wp
    real(wp),parameter :: b192 = 0
    real(wp),parameter :: b193 = 0
    real(wp),parameter :: b194 = 0
    real(wp),parameter :: b195 = 0
    real(wp),parameter :: b196 = -.6005539308646711394173617417282635784603029781694346062080009050806164919590318269156_wp
    real(wp),parameter :: b197 = .1585222658367901823092810289995660365558277114283363946284113044047141269675809741481e-2_wp
    real(wp),parameter :: b198 = .3026658851734428583799288358366834830670493200050903757742452562250879940584583879948_wp
    real(wp),parameter :: b199 = .4285248952296136656501504776128834912424841890282487604085652387338305228046948965999_wp
    real(wp),parameter :: b1910 = -.2494936049956423479331713062821820760626388763661853215936426940474748269791918051828e-1_wp
    real(wp),parameter :: b1911 = .2011342149874705950638155217740310555133571768677965217932407885434178872356469888792e-1_wp
    real(wp),parameter :: b1912 = .1265809640150575117740069284741749647793461966567090378701137188708013785394651655531_wp
    real(wp),parameter :: b1913 = -.7625505511021937468044360160133448561674800889340550320120549846220783051800508863794e-2_wp
    real(wp),parameter :: b1914 = .2257683215130023640661938534278959810874704491725768321513002364066193853427895981087_wp
    real(wp),parameter :: b1915 = -.2745492489167708483879226183906024865124767290578864277433486220883760763766017543321e-1_wp
    real(wp),parameter :: b1916 = .7783761260627937338705544952370637303437265442198231246348540202794422441986844495010e-2_wp
    real(wp),parameter :: b1917 = .1680830476266175444138000910160343423421707378563384687234091879582812751321589095322e-4_wp
    real(wp),parameter :: b1918 = -.2135549195295375067495229039525877007339581782873806009604806379884885927885993014667e-1_wp
    real(wp),parameter :: b201 = .4918722777774214635454694806095471334564471116346374258920200467475606642975663633169e-1_wp
    real(wp),parameter :: b202 = 0
    real(wp),parameter :: b203 = 0
    real(wp),parameter :: b204 = 0
    real(wp),parameter :: b205 = 0
    real(wp),parameter :: b206 = .8381198329740939383758467268684128106638290469173528975872944409244863365387981657651e-1_wp
    real(wp),parameter :: b207 = .9031585241436450659339626877786556727134244569192456020771509670971435444573825121328e-5_wp
    real(wp),parameter :: b208 = .9729721702550808557495179062390587533407132090396133365371679075124515126725605224924e-1_wp
    real(wp),parameter :: b209 = -.7780700010332922911083090560520956762627768671838912755501179322945892692289272010345e-1_wp
    real(wp),parameter :: b2010 = .2322816011687352870985186646759446719321113455516507860765563208248738749278233611786_wp
    real(wp),parameter :: b2011 = -.1069713522971685698540328329852038445947544770807964964519088069861375575524890552230e-1_wp
    real(wp),parameter :: b2012 = -.1273357855008010164642613326804825893357151796749893640146110897003842493364412495521_wp
    real(wp),parameter :: b2013 = .1195596306045770335798408097988199315656849089154404333183111529360783290389646594900_wp
    real(wp),parameter :: b2014 = .1268882175226586102719033232628398791540785498489425981873111782477341389728096676737_wp
    real(wp),parameter :: b2015 = .2044989775051124744376278118609406952965235173824130879345603271983640081799591002045e-1_wp
    real(wp),parameter :: b2016 = -.3909856506984242164691157929655038995858933888334921275260056794557284704385445582377e-2_wp
    real(wp),parameter :: b2017 = -.9117410119332842109125953032188840693633119317454401715027062207348542639139803659651e-5_wp
    real(wp),parameter :: b2018 = .1126093151077456154585623418197514536662349793737574398626301058923382946852045056475e-1_wp
    real(wp),parameter :: b2019 = -.3209868434922071245903287586373535845929558438862699119277778588606558307508436673461_wp
    real(wp),parameter :: b211 = .5247404461891304721708365639844356354386624470737997647397329116437616959059854701778_wp
    real(wp),parameter :: b212 = 0
    real(wp),parameter :: b213 = 0
    real(wp),parameter :: b214 = .7610651429965941990931946190296373660533922130208900643092560661920652146268815424885e-1_wp
    real(wp),parameter :: b215 = -2.135538596825204678401302463372339475781988628449683891486081018353295131152366999963_wp
    real(wp),parameter :: b216 = 2.119016745189825526524339470866652730733632823644186992594087231585247087949150559862_wp
    real(wp),parameter :: b217 = 0
    real(wp),parameter :: b218 = 0
    real(wp),parameter :: b219 = -.6068669292751351984511141135476525247489731841233590933089036540462466843434512014723_wp
    real(wp),parameter :: b2110 = -.6975023816048118989159659127206041759560640169551283418432821453966450135301718130674_wp
    real(wp),parameter :: b2111 = -.2539552135383387231958801508125376392283214913791972115937570999107334296436460657368e-1_wp
    real(wp),parameter :: b2112 = 0
    real(wp),parameter :: b2113 = 0
    real(wp),parameter :: b2114 = 0
    real(wp),parameter :: b2115 = 0
    real(wp),parameter :: b2116 = 0
    real(wp),parameter :: b2117 = 0
    real(wp),parameter :: b2118 = .2539552135383387231958801508125376392283214913791972115937570999107334296436460657368e-1_wp
    real(wp),parameter :: b2119 = .6975023816048118989159659127206041759560640169551283418432821453966450135301718130674_wp
    real(wp),parameter :: b2120 = .6068669292751351984511141135476525247489731841233590933089036540462466843434512014723_wp
    real(wp),parameter :: b221 = .3850951504524952575244726520324404860898823844950484490138343541601000605579674908780e-1_wp
    real(wp),parameter :: b222 = 0
    real(wp),parameter :: b223 = 0
    real(wp),parameter :: b224 = .9889604363651382462812798900870548282288802143553903392828423536490053424426345577709e-4_wp
    real(wp),parameter :: b225 = .1258652731216402762123982910182735616976625577395859318541619645591514091994326716669_wp
    real(wp),parameter :: b226 = 0
    real(wp),parameter :: b227 = .1148546657824708136802241649776446744025406100582655444423341878411637297467351936831_wp
    real(wp),parameter :: b228 = 0
    real(wp),parameter :: b229 = -.1299246462925958527077352114079744607750386324042347411360003531859939834460301561195_wp
    real(wp),parameter :: b2210 = -.3664591598580916638621456622089859363144583004356560608229075984337239720284817960193_wp
    real(wp),parameter :: b2211 = 0
    real(wp),parameter :: b2212 = 0
    real(wp),parameter :: b2213 = 0
    real(wp),parameter :: b2214 = 0
    real(wp),parameter :: b2215 = 0
    real(wp),parameter :: b2216 = 0
    real(wp),parameter :: b2217 = 0
    real(wp),parameter :: b2218 = 0
    real(wp),parameter :: b2219 = .3664591598580916638621456622089859363144583004356560608229075984337239720284817960193_wp
    real(wp),parameter :: b2220 = .1299246462925958527077352114079744607750386324042347411360003531859939834460301561195_wp
    real(wp),parameter :: b2221 = -.1148546657824708136802241649776446744025406100582655444423341878411637297467351936831_wp
    real(wp),parameter :: b231 = .2175227402212137286104398734798923400326123258875071216675507357419304139407870179368_wp
    real(wp),parameter :: b232 = .2254597159191371485825425826604585371603701302528437555254317203984204632521954382036_wp
    real(wp),parameter :: b233 = 0
    real(wp),parameter :: b234 = 0
    real(wp),parameter :: b235 = 0
    real(wp),parameter :: b236 = -.7003676470588235294117647058823529411764705882352941176470588235294117647058823529412_wp
    real(wp),parameter :: b237 = -.3841432262252079340469853296205751191707035527677198880717241749265840367609120672261_wp
    real(wp),parameter :: b238 = 0
    real(wp),parameter :: b239 = 0
    real(wp),parameter :: b2310 = 0
    real(wp),parameter :: b2311 = 0
    real(wp),parameter :: b2312 = 0
    real(wp),parameter :: b2313 = 0
    real(wp),parameter :: b2314 = 0
    real(wp),parameter :: b2315 = 0
    real(wp),parameter :: b2316 = 0
    real(wp),parameter :: b2317 = 0
    real(wp),parameter :: b2318 = 0
    real(wp),parameter :: b2319 = 0
    real(wp),parameter :: b2320 = 0
    real(wp),parameter :: b2321 = .3841432262252079340469853296205751191707035527677198880717241749265840367609120672261_wp
    real(wp),parameter :: b2322 = .7003676470588235294117647058823529411764705882352941176470588235294117647058823529412_wp
    real(wp),parameter :: b241 = .4351851851851851851851851851851851851851851851851851851851851851851851851851851851852_wp
    real(wp),parameter :: b242 = 0
    real(wp),parameter :: b243 = -.4244806610219170956648818689597945762515945523075081224120777274585003570971510665439_wp
    real(wp),parameter :: b244 = 0
    real(wp),parameter :: b245 = 0
    real(wp),parameter :: b246 = 0
    real(wp),parameter :: b247 = 0
    real(wp),parameter :: b248 = 0
    real(wp),parameter :: b249 = 0
    real(wp),parameter :: b2410 = 0
    real(wp),parameter :: b2411 = 0
    real(wp),parameter :: b2412 = 0
    real(wp),parameter :: b2413 = 0
    real(wp),parameter :: b2414 = 0
    real(wp),parameter :: b2415 = 0
    real(wp),parameter :: b2416 = 0
    real(wp),parameter :: b2417 = 0
    real(wp),parameter :: b2418 = 0
    real(wp),parameter :: b2419 = 0
    real(wp),parameter :: b2420 = 0
    real(wp),parameter :: b2421 = 0
    real(wp),parameter :: b2422 = 0
    real(wp),parameter :: b2423 = .4244806610219170956648818689597945762515945523075081224120777274585003570971510665439_wp
    real(wp),parameter :: b251 = 14.54990971513478580914495526598245283940024810265087163382498896386800769812364397026_wp
    real(wp),parameter :: b252 = -2.609444444444444444444444444444444444444444444444444444444444444444444444444444444444_wp
    real(wp),parameter :: b253 = -2.016004609236637754870351028563643794559738431497207211298306162299623087053267335725_wp
    real(wp),parameter :: b254 = 0
    real(wp),parameter :: b255 = 0
    real(wp),parameter :: b256 = -1.666875_wp
    real(wp),parameter :: b257 = -1.840010137609049715480551028603992780751854184812582730319893278211480604766702319850_wp
    real(wp),parameter :: b258 = 0
    real(wp),parameter :: b259 = 112.8850026879393650927243496060343039004107661597559344830054269802155904451112935844_wp
    real(wp),parameter :: b2510 = 123.3942086822776167534198661922997568027612066045886506351026838994566451157760161079_wp
    real(wp),parameter :: b2511 = -.7912126656078716671206667350992937018785263673658058027818457299775466324638082312137_wp
    real(wp),parameter :: b2512 = -50.05149873558555185701853277950418028083920384013973228797374663569746882026502756132_wp
    real(wp),parameter :: b2513 = 24.88778291494286023489265916235844707524170235927746532159915897878512972222657954044_wp
    real(wp),parameter :: b2514 = -212.1164197057320847511194729602871445996396409264694299135694995667783941778631290781_wp
    real(wp),parameter :: b2515 = -17.89082255740024165397920079197116480379826062999613139163493416384514698556677051724_wp
    real(wp),parameter :: b2516 = 2.509716434086569985363077956674238286724267256774032757480361264001771307725305753519_wp
    real(wp),parameter :: b2517 = .1162475886937088360535666476065014131326443901611275096836079807311057386735832665899_wp
    real(wp),parameter :: b2518 = .5840328281597421751021178306186155103094570763962374837158281176580253234011237044692_wp
    real(wp),parameter :: b2519 = 1.584417806712708397483834073813399023446091584920600008797782958694409008929749703075_wp
    real(wp),parameter :: b2520 = 1.338635006378392645053446531474068534729248229446179562750186952887872256191439757181_wp
    real(wp),parameter :: b2521 = 1.840010137609049715480551028603992780751854184812582730319893278211480604766702319850_wp
    real(wp),parameter :: b2522 = 1.666875_wp
    real(wp),parameter :: b2523 = 2.016004609236637754870351028563643794559738431497207211298306162299623087053267335725_wp
    real(wp),parameter :: b2524 = 2.609444444444444444444444444444444444444444444444444444444444444444444444444444444444_wp
    real(wp),parameter :: b261 = .4213659219087082450668941175898637044576077571412870755532166081975824639155816025423_wp
    real(wp),parameter :: b262 = 0
    real(wp),parameter :: b263 = 0
    real(wp),parameter :: b264 = 0
    real(wp),parameter :: b265 = 0
    real(wp),parameter :: b266 = 2.360375290413766425107807321597798068298832670876354152068674897991709951684675057597_wp
    real(wp),parameter :: b267 = .7887926811836902144270824477231365437442962856383491851349713361086542561144935964760e-1_wp
    real(wp),parameter :: b268 = -1.881850641776530466652474803895333308262409829541779704822719782535259533895036456964_wp
    real(wp),parameter :: b269 = -1.304700734906095391371228323883517431348033016096248869612640797091741378769699510853_wp
    real(wp),parameter :: b2610 = .1146971532060496506611311517441641422873299688900135018999537656144573134542692746612_wp
    real(wp),parameter :: b2611 = -.5223613182942077907170609676910338480915071906000186351215047373740452623507436746018e-2_wp
    real(wp),parameter :: b2612 = .7134840563194221964556259902880063405282394887795535106616674222878805188799073074497_wp
    real(wp),parameter :: b2613 = 0
    real(wp),parameter :: b2614 = 0
    real(wp),parameter :: b2615 = 0
    real(wp),parameter :: b2616 = 0
    real(wp),parameter :: b2617 = 0
    real(wp),parameter :: b2618 = 0
    real(wp),parameter :: b2619 = 0
    real(wp),parameter :: b2620 = 0
    real(wp),parameter :: b2621 = 0
    real(wp),parameter :: b2622 = 0
    real(wp),parameter :: b2623 = 0
    real(wp),parameter :: b2624 = 0
    real(wp),parameter :: b2625 = 0
    real(wp),parameter :: b271 = -1.016867684065179179311540011641152067739527559831057381007428484343486238202737995538_wp
    real(wp),parameter :: b272 = 0
    real(wp),parameter :: b273 = 0
    real(wp),parameter :: b274 = 0
    real(wp),parameter :: b275 = 0
    real(wp),parameter :: b276 = -7.712044352285817603610736737545203003182107799250377304475316992902646516184627861475_wp
    real(wp),parameter :: b277 = -.4034008409374858753410643280039779311650023076266296210030337872976327092081717454860_wp
    real(wp),parameter :: b278 = 6.739165476490825275476530741799137001781411805688541615342712456896969009075069743117_wp
    real(wp),parameter :: b279 = 6.014994643407224294180918860565568523540180411761624603152605084567461778121724411633_wp
    real(wp),parameter :: b2710 = -1.138427387973993086846707441740657423236451007119331997503516596250482644801941463866_wp
    real(wp),parameter :: b2711 = .5009271973181599563449431362188397685579770251743883251930235916805543539186803311245e-1_wp
    real(wp),parameter :: b2712 = -3.113250932564715585587456369457050044282973515833171733791753110672694235302567710668_wp
    real(wp),parameter :: b2713 = 0
    real(wp),parameter :: b2714 = 0
    real(wp),parameter :: b2715 = 0
    real(wp),parameter :: b2716 = 0
    real(wp),parameter :: b2717 = 0
    real(wp),parameter :: b2718 = 0
    real(wp),parameter :: b2719 = 0
    real(wp),parameter :: b2720 = 0
    real(wp),parameter :: b2721 = 0
    real(wp),parameter :: b2722 = 0
    real(wp),parameter :: b2723 = 0
    real(wp),parameter :: b2724 = 0
    real(wp),parameter :: b2725 = 0
    real(wp),parameter :: b2726 = 1.384086184284282287144691407184059663080846182736441247635994288225760468937471545692_wp
    real(wp),parameter :: b281 = 1.131093475949031458408970675798323789793651098141584053672723684523337489654810813567_wp
    real(wp),parameter :: b282 = 0
    real(wp),parameter :: b283 = 0
    real(wp),parameter :: b284 = 0
    real(wp),parameter :: b285 = 0
    real(wp),parameter :: b286 = -11.30475611955440577592346561419842170756276153921836803995513790246783417109406727002_wp
    real(wp),parameter :: b287 = .8673508908529372037894544277195364499375277491231106383198729858086287127092489733738e-1_wp
    real(wp),parameter :: b288 = 4.971317844154333915807514558966554931901059536397301577948574837963960265600948244094_wp
    real(wp),parameter :: b289 = 14.86493772010299652718002500847699984298963300991290479767277341466986550270499895786_wp
    real(wp),parameter :: b2810 = -5.526130551905351405702373768620234518347747212552226033118137855980983127699224086148_wp
    real(wp),parameter :: b2811 = .1017790491986200061558195486579246543940163889995857946318495764475820398967567810759_wp
    real(wp),parameter :: b2812 = -5.412708567655345677389304794550103449135846140313886894023804835964789906470796078495_wp
    real(wp),parameter :: b2813 = 0
    real(wp),parameter :: b2814 = 0
    real(wp),parameter :: b2815 = 0
    real(wp),parameter :: b2816 = 0
    real(wp),parameter :: b2817 = 0
    real(wp),parameter :: b2818 = 0
    real(wp),parameter :: b2819 = 0
    real(wp),parameter :: b2820 = 0
    real(wp),parameter :: b2821 = 0
    real(wp),parameter :: b2822 = 0
    real(wp),parameter :: b2823 = 0
    real(wp),parameter :: b2824 = 0
    real(wp),parameter :: b2825 = 0
    real(wp),parameter :: b2826 = 2.119905903216124397337756706998226742489167533804374626293399155463227287948398852330_wp
    real(wp),parameter :: b2827 = -.1603789707964253713820928925063521366431305782887091520824325014403564569409562398070_wp
    real(wp),parameter :: b291 = 46.12864603958015905056850990838704062569763496465412763519385202873225169326091596022_wp
    real(wp),parameter :: b292 = 0
    real(wp),parameter :: b293 = 0
    real(wp),parameter :: b294 = 0
    real(wp),parameter :: b295 = 0
    real(wp),parameter :: b296 = 27.91300163119399908845158457840426795358131126287389180101909720743727165524743997096_wp
    real(wp),parameter :: b297 = 16.11362689862451240990975288339484000039234961533373225370151713993919933068137837655_wp
    real(wp),parameter :: b298 = -125.4696763444318726329250646477825268481685879278990572088898486327607587168793188547_wp
    real(wp),parameter :: b299 = 76.57182020120529497684089567511627659347626021577320126137924234653281598557298427556_wp
    real(wp),parameter :: b2910 = -48.97805558723490361747755876313897556229903489002597907918723603724236906416246304792_wp
    real(wp),parameter :: b2911 = -1.242830487244052672528847080627776989497284066852925470069357906835696725179561610962_wp
    real(wp),parameter :: b2912 = 18.85807213383620068645464308722546866214730809725606975205204214283143331463602851729_wp
    real(wp),parameter :: b2913 = 0
    real(wp),parameter :: b2914 = 0
    real(wp),parameter :: b2915 = 0
    real(wp),parameter :: b2916 = 0
    real(wp),parameter :: b2917 = 0
    real(wp),parameter :: b2918 = 0
    real(wp),parameter :: b2919 = 0
    real(wp),parameter :: b2920 = 0
    real(wp),parameter :: b2921 = 0
    real(wp),parameter :: b2922 = 0
    real(wp),parameter :: b2923 = 0
    real(wp),parameter :: b2924 = 0
    real(wp),parameter :: b2925 = 0
    real(wp),parameter :: b2926 = -8.871982194511738170929283011936752038683182063824370821380217475185779921363890310462_wp
    real(wp),parameter :: b2927 = -2.069534982695615656321598541301254038890059866116256078922522295738656144906973495868_wp
    real(wp),parameter :: b2928 = 2.046912691678016537956965912259391642243284658827565955103431482290288593093460219351_wp

    real(wp),parameter :: c1 = .2380952380952380952380952380952380952380952380952380952380952380952380952380952380952e-1_wp
    real(wp),parameter :: c2 = -.11_wp
    real(wp),parameter :: c3 = -.17_wp
    real(wp),parameter :: c4 = 0
    real(wp),parameter :: c5 = 0
    real(wp),parameter :: c6 = -.19_wp
    real(wp),parameter :: c7 = -.21_wp
    real(wp),parameter :: c8 = 0
    real(wp),parameter :: c9 = -.23_wp
    real(wp),parameter :: c10 = -.27_wp
    real(wp),parameter :: c11 = -.29_wp
    real(wp),parameter :: c12 = 0
    real(wp),parameter :: c13 = .1384130236807829740053502031450331467488136400899412345912671194817223119377730668077_wp
    real(wp),parameter :: c14 = .2158726906049313117089355111406811389654720741957730511230185948039919737765126474781_wp
    real(wp),parameter :: c15 = .2438095238095238095238095238095238095238095238095238095238095238095238095238095238095_wp
    real(wp),parameter :: c16 = .2158726906049313117089355111406811389654720741957730511230185948039919737765126474781_wp
    real(wp),parameter :: c17 = .1384130236807829740053502031450331467488136400899412345912671194817223119377730668077_wp
    real(wp),parameter :: c18 = .29_wp
    real(wp),parameter :: c19 = .27_wp
    real(wp),parameter :: c20 = .23_wp
    real(wp),parameter :: c21 = .21_wp
    real(wp),parameter :: c22 = .19_wp
    real(wp),parameter :: c23 = .17_wp
    real(wp),parameter :: c24 = .11_wp
    real(wp),parameter :: c25 = .2380952380952380952380952380952380952380952380952380952380952380952380952380952380952e-1_wp
    real(wp),parameter :: c26 = 0
    real(wp),parameter :: c27 = 0
    real(wp),parameter :: c28 = 0
    real(wp),parameter :: c29 = 0

    real(wp),parameter :: d1 = .1357267366422036691624508570375039921213961405197383540990729271585608575744367221422e-1_wp
    real(wp),parameter :: d2 = 0
    real(wp),parameter :: d3 = 0
    real(wp),parameter :: d4 = 0
    real(wp),parameter :: d5 = 0
    real(wp),parameter :: d6 = 0
    real(wp),parameter :: d7 = 0
    real(wp),parameter :: d8 = 0
    real(wp),parameter :: d9 = .1957242608025905233613155193355413222756712243739956539857744660017249157522413539338_wp
    real(wp),parameter :: d10 = .6188866347435608661616060238380345321740579185698262925072441030810461130118224626650e-1_wp
    real(wp),parameter :: d11 = .2356461254963383884566009189531765712282459575681722654363888116993755654359289982043_wp
    real(wp),parameter :: d12 = .9356981277656948171659125355126822887878370575133833270918675400840988929847712528189e-1_wp
    real(wp),parameter :: d13 = 0
    real(wp),parameter :: d14 = 0
    real(wp),parameter :: d15 = 0
    real(wp),parameter :: d16 = 0
    real(wp),parameter :: d17 = 0
    real(wp),parameter :: d18 = 0
    real(wp),parameter :: d19 = 0
    real(wp),parameter :: d20 = 0
    real(wp),parameter :: d21 = 0
    real(wp),parameter :: d22 = 0
    real(wp),parameter :: d23 = 0
    real(wp),parameter :: d24 = 0
    real(wp),parameter :: d25 = 0
    real(wp),parameter :: d26 = .2788382624223597882496809755901865993371492355322230983832374605891024297949723722505_wp
    real(wp),parameter :: d27 = .4265887719284871852244002213010235951116369449606049493109912099005098392090245967637_wp
    real(wp),parameter :: d28 = -.2878025166474501962999477241233861487151201542084751225511851437798538088192303308430_wp
    real(wp),parameter :: d29 = -.1802605391747162424104685269536402054591231988681564193502526144322952773004003407203e-1_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e2  = c2  - d2
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10 = c10 - d10
    real(wp),parameter :: e11 = c11 - d11
    real(wp),parameter :: e12 = c12 - d12
    real(wp),parameter :: e13 = c13 - d13
    real(wp),parameter :: e14 = c14 - d14
    real(wp),parameter :: e15 = c15 - d15
    real(wp),parameter :: e16 = c16 - d16
    real(wp),parameter :: e17 = c17 - d17
    real(wp),parameter :: e18 = c18 - d18
    real(wp),parameter :: e19 = c19 - d19
    real(wp),parameter :: e20 = c20 - d20
    real(wp),parameter :: e21 = c21 - d21
    real(wp),parameter :: e22 = c22 - d22
    real(wp),parameter :: e23 = c23 - d23
    real(wp),parameter :: e24 = c24 - d24
    real(wp),parameter :: e25 = c25 - d25
    real(wp),parameter :: e26 = c26 - d26
    real(wp),parameter :: e27 = c27 - d27
    real(wp),parameter :: e28 = c28 - d28
    real(wp),parameter :: e29 = c29 - d29

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13), &
               f14 => me%funcs(:,14), &
               f15 => me%funcs(:,15), &
               f16 => me%funcs(:,16), &
               f17 => me%funcs(:,17), &
               f18 => me%funcs(:,18), &
               f19 => me%funcs(:,19), &
               f20 => me%funcs(:,20), &
               f21 => me%funcs(:,21), &
               f22 => me%funcs(:,22), &
               f23 => me%funcs(:,23), &
               f24 => me%funcs(:,24), &
               f25 => me%funcs(:,25), &
               f26 => me%funcs(:,26), &
               f27 => me%funcs(:,27), &
               f28 => me%funcs(:,28), &
               f29 => me%funcs(:,29))

        call me%f(t*h,x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b82*f2+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b92*f2+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b102*f2+b103*f3+b104*f4+b105*f5+b106*f6+b107*f7+&
                            b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b112*f2+b113*f3+b114*f4+b115*f5+b116*f6+b117*f7+&
                            b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b122*f2+b123*f3+b124*f4+b125*f5+b126*f6+b127*f7+&
                            b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b131*f1+b132*f2+b133*f3+b134*f4+b135*f5+b136*f6+b137*f7+&
                            b138*f8+b139*f9+b1310*f10+b1311*f11+b1312*f12),f13)
        call me%f(t+a14*h,x+h*(b141*f1+b142*f2+b143*f3+b144*f4+b145*f5+b146*f6+b147*f7+&
                            b148*f8+b149*f9+b1410*f10+b1411*f11+b1412*f12+b1413*f13),f14)
        call me%f(t+a15*h,x+h*(b151*f1+b152*f2+b153*f3+b154*f4+b155*f5+b156*f6+b157*f7+&
                            b158*f8+b159*f9+b1510*f10+b1511*f11+b1512*f12+b1513*f13+&
                            b1514*f14),f15)
        call me%f(t+a16*h,x+h*(b161*f1+b162*f2+b163*f3+b164*f4+b165*f5+b166*f6+b167*f7+&
                            b168*f8+b169*f9+b1610*f10+b1611*f11+b1612*f12+b1613*f13+&
                            b1614*f14+b1615*f15),f16)
        call me%f(t+a17*h,x+h*(b171*f1+b172*f2+b173*f3+b174*f4+b175*f5+b176*f6+b177*f7+&
                            b178*f8+b179*f9+b1710*f10+b1711*f11+b1712*f12+b1713*f13+&
                            b1714*f14+b1715*f15+b1716*f16),f17)
        call me%f(t+a18*h,x+h*(b181*f1+b182*f2+b183*f3+b184*f4+b185*f5+b186*f6+b187*f7+&
                            b188*f8+b189*f9+b1810*f10+b1811*f11+b1812*f12+b1813*f13+&
                            b1814*f14+b1815*f15+b1816*f16+b1817*f17),f18)
        call me%f(t+a19*h,x+h*(b191*f1+b192*f2+b193*f3+b194*f4+b195*f5+b196*f6+b197*f7+&
                            b198*f8+b199*f9+b1910*f10+b1911*f11+b1912*f12+b1913*f13+&
                            b1914*f14+b1915*f15+b1916*f16+b1917*f17+b1918*f18),f19)
        call me%f(t+a20*h,x+h*(b201*f1+b202*f2+b203*f3+b204*f4+b205*f5+b206*f6+b207*f7+&
                            b208*f8+b209*f9+b2010*f10+b2011*f11+b2012*f12+b2013*f13+&
                            b2014*f14+b2015*f15+b2016*f16+b2017*f17+b2018*f18+b2019*f19),f20)
        call me%f(t+a21*h,x+h*(b211*f1+b212*f2+b213*f3+b214*f4+b215*f5+b216*f6+b217*f7+&
                            b218*f8+b219*f9+b2110*f10+b2111*f11+b2112*f12+b2113*f13+&
                            b2114*f14+b2115*f15+b2116*f16+b2117*f17+b2118*f18+b2119*f19+&
                            b2120*f20),f21)
        call me%f(t+a22*h,x+h*(b221*f1+b222*f2+b223*f3+b224*f4+b225*f5+b226*f6+b227*f7+&
                            b228*f8+b229*f9+b2210*f10+b2211*f11+b2212*f12+b2213*f13+&
                            b2214*f14+b2215*f15+b2216*f16+b2217*f17+b2218*f18+b2219*f19+&
                            b2220*f20+b2221*f21),f22)
        call me%f(t+a23*h,x+h*(b231*f1+b232*f2+b233*f3+b234*f4+b235*f5+b236*f6+b237*f7+&
                            b238*f8+b239*f9+b2310*f10+b2311*f11+b2312*f12+b2313*f13+&
                            b2314*f14+b2315*f15+b2316*f16+b2317*f17+b2318*f18+b2319*f19+&
                            b2320*f20+b2321*f21+b2322*f22),f23)
        call me%f(t+a24*h,x+h*(b241*f1+b242*f2+b243*f3+b244*f4+b245*f5+b246*f6+b247*f7+&
                            b248*f8+b249*f9+b2410*f10+b2411*f11+b2412*f12+b2413*f13+&
                            b2414*f14+b2415*f15+b2416*f16+b2417*f17+b2418*f18+b2419*f19+&
                            b2420*f20+b2421*f21+b2422*f22+b2423*f23),f24)
        call me%f(t+h,    x+h*(b251*f1+b252*f2+b253*f3+b254*f4+b255*f5+b256*f6+b257*f7+&
                            b258*f8+b259*f9+b2510*f10+b2511*f11+b2512*f12+b2513*f13+&
                            b2514*f14+b2515*f15+b2516*f16+b2517*f17+b2518*f18+b2519*f19+&
                            b2520*f20+b2521*f21+b2522*f22+b2523*f23+b2524*f24),f25)
        call me%f(t+a26*h,x+h*(b261*f1+b262*f2+b263*f3+b264*f4+b265*f5+b266*f6+b267*f7+&
                            b268*f8+b269*f9+b2610*f10+b2611*f11+b2612*f12+b2613*f13+&
                            b2614*f14+b2615*f15+b2616*f16+b2617*f17+b2618*f18+b2619*f19+&
                            b2620*f20+b2621*f21+b2622*f22+b2623*f23+b2624*f24+b2625*f25),f26)
        call me%f(t+a27*h,x+h*(b271*f1+b272*f2+b273*f3+b274*f4+b275*f5+b276*f6+b277*f7+&
                            b278*f8+b279*f9+b2710*f10+b2711*f11+b2712*f12+b2713*f13+&
                            b2714*f14+b2715*f15+b2716*f16+b2717*f17+b2718*f18+b2719*f19+&
                            b2720*f20+b2721*f21+b2722*f22+b2723*f23+b2724*f24+b2725*f25+&
                            b2726*f26),f27)
        call me%f(t+a28*h,x+h*(b281*f1+b282*f2+b283*f3+b284*f4+b285*f5+b286*f6+b287*f7+&
                            b288*f8+b289*f9+b2810*f10+b2811*f11+b2812*f12+b2813*f13+&
                            b2814*f14+b2815*f15+b2816*f16+b2817*f17+b2818*f18+b2819*f19+&
                            b2820*f20+b2821*f21+b2822*f22+b2823*f23+b2824*f24+b2825*f25+&
                            b2826*f26+b2827*f27),f28)
        call me%f(t+h,    x+h*(b291*f1+b292*f2+b293*f3+b294*f4+b295*f5+b296*f6+b297*f7+&
                            b298*f8+b299*f9+b2910*f10+b2911*f11+b2912*f12+b2913*f13+&
                            b2914*f14+b2915*f15+b2916*f16+b2917*f17+b2918*f18+b2919*f19+&
                            b2920*f20+b2921*f21+b2922*f22+b2923*f23+b2924*f24+b2925*f25+&
                            b2926*f26+b2927*f27+b2928*f28),f29)

        xf = x+h*(c1*f1+c2*f2+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+&
                c11*f11+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16+c17*f17+c18*f18+&
                c19*f19+c20*f20+c21*f21+c22*f22+c23*f23+c24*f24+c25*f25+c26*f26+&
                c27*f27+c28*f28+c29*f29)

        xerr = h*(e1*f1+e2*f2+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+&
                e11*f11+e12*f12+e13*f13+e14*f14+e15*f15+e16*f16+e17*f17+e18*f18+&
                e19*f19+e20*f20+e21*f21+e22*f22+e23*f23+e24*f24+e25*f25+e26*f26+&
                e27*f27+e28*f28+e29*f29)

    end associate

    end procedure rko129
!*****************************************************************************************

!*****************************************************************************************
!>
!  Stone's 26 stage combined order 10 and 11 Runge-Kutta scheme
!
!### Reference
!  * [Coefficients](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK11/RKcoeff11_a.pdf)

    module procedure rks1110a

    real(wp),parameter :: a2 = .2494510436877389344264053439636521218011835467298426368145366983339993718715373083760_wp
    real(wp),parameter :: a3 = .3741765655316084016396080159454781827017753200947639552218050475009990578073059625640_wp
    real(wp),parameter :: a4 = .5612648482974126024594120239182172740526629801421459328327075712514985867109589438460_wp
    real(wp),parameter :: a5 = .8546995299174592915742289900418092935366740866720183954100852188058051166442961011609e-1_wp
    real(wp),parameter :: a6 = .7310788517640473762691237698106576677338584467168049416410984572334277402413648068397_wp
    real(wp),parameter :: a7 = .1304049416609471516815374056280027453671928620452985586822237474262182566918325326012_wp
    real(wp),parameter :: a8 = .4712636357202847267558933149451482905659499178473869504545717970436656667382069175430_wp
    real(wp),parameter :: a9 = .7730997092193940758126470461455964941402897392485301085000777589597371419303943343132_wp
    real(wp),parameter :: a10 = .2056292271097569160894663338279693801391636310750116019089310570387586267029180604980_wp
    real(wp),parameter :: a11 = .5366419968074299811348135248875344652445218400812654186620229284574082136119576258888_wp
    real(wp),parameter :: a12 = .6950637644720240807995723705690567056506824408065743296638976932495337928615038239116_wp
    real(wp),parameter :: a13 = .8049629952111449717022202873313016978667827601218981279930343926861123204179364388333_wp
    real(wp),parameter :: a14 = .1913760139217099199372138166733346224093015957152232805822661231536077778241788913071_wp
    real(wp),parameter :: a15 = .8894442715562032655697860690613313211578199121477983480401175761840408067430713807992_wp
    real(wp),parameter :: a16 = .4680151547126001101041424878175828807838998612363955308441536427569141462651129332201_wp
    real(wp),parameter :: a17 = .9908386187455954897815362931642001409443269908386187455954897815362931642001409443270_wp
    real(wp),parameter :: a18 = .9000188288457917529655432122010920730559216720015063076633402372434569760873658444737e-1_wp
    real(wp),parameter :: a19 = .4031478884213807074957145083372292348449431198379304971170328814087579865980988000623_wp
    real(wp),parameter :: a20 = .2426620464737056665307786384019567876070118222584590297594781899714635140644109253975_wp
    real(wp),parameter :: a21 = .7208245981830887491264849755415793151642208245981830887491264849755415793151642208246_wp
    real(wp),parameter :: a22 = .1999681579366342939022448654672822798917369845565992676325425887597516319057474924375_wp
    real(wp),parameter :: a23 = .9639977912755383765875207067918277194919933738266151297625621203754831584759801214798_wp
    real(wp),parameter :: a24 = .7436137801734239512538082962268572767752519334426997890789782048277478322006093273963_wp

    real(wp),parameter :: b21 = .2494510436877389344264053439636521218011835467298426368145366983339993718715373083760_wp
    real(wp),parameter :: b31 = .9354414138290210040990200398636954567544383002369098880545126187524976445182649064100e-1_wp
    real(wp),parameter :: b32 = .2806324241487063012297060119591086370263314900710729664163537856257492933554794719230_wp
    real(wp),parameter :: b41 = .1403162120743531506148530059795543185131657450355364832081768928128746466777397359615_wp
    real(wp),parameter :: b42 = 0
    real(wp),parameter :: b43 = .4209486362230594518445590179386629555394972351066094496245306784386239400332192078845_wp
    real(wp),parameter :: b51 = .7019164670454382087443786606405374274908739183454652276598025183935982136630730721870e-1_wp
    real(wp),parameter :: b52 = 0
    real(wp),parameter :: b53 = .2631174841924423444834727510216642181614053183934029855724529099908799066900681707156e-1_wp
    real(wp),parameter :: b54 = -.1103344213204212616536224216203923521156051500668498178221702095786730037088451417417e-1_wp
    real(wp),parameter :: b61 = -.1566211240781183642488518551453027928209011808429926742075994450765608166883593055436_wp
    real(wp),parameter :: b62 = 0
    real(wp),parameter :: b63 = 0
    real(wp),parameter :: b64 = .4022037025415642768212777334976723807953788960750969348079276345083207687523340439985_wp
    real(wp),parameter :: b65 = .4854962733006014636966978914582880797593807314847006810407702678016677881773900683848_wp
    real(wp),parameter :: b71 = .3112249065189454540975938749103688775640860268015765517875901819636788841622883514715e-1_wp
    real(wp),parameter :: b72 = 0
    real(wp),parameter :: b73 = 0
    real(wp),parameter :: b74 = 0
    real(wp),parameter :: b75 = .9925602946470755730976861333746552578010113543371556501156532289522223083686039750464e-1_wp
    real(wp),parameter :: b76 = .2642154434504896200940479950033183068312393142533849189940633462813743874329994944194e-4_wp
    real(wp),parameter :: b81 = .7097221158585919779802904970340267773691990074739280896770286766721241217196318737330_wp
    real(wp),parameter :: b82 = 0
    real(wp),parameter :: b83 = 0
    real(wp),parameter :: b84 = 0
    real(wp),parameter :: b85 = -2.878599621086994245472819276760640913301855954929611359904818127369971660983418731969_wp
    real(wp),parameter :: b86 = .2129446575803202495593011450000675123146819436283352905644086375052608273806710598083e-1_wp
    real(wp),parameter :: b87 = 2.618846675190654969292491980171755675267138670940236691625920383990987123263926669799_wp
    real(wp),parameter :: b91 = .1971162769751622653374279030890901895102271260934752560073298666516820598947694065245e-1_wp
    real(wp),parameter :: b92 = 0
    real(wp),parameter :: b93 = 0
    real(wp),parameter :: b94 = 0
    real(wp),parameter :: b95 = 0
    real(wp),parameter :: b96 = .1314769583961512153957380944910102917010453313237904532896703042317759437588108172284_wp
    real(wp),parameter :: b97 = .2651026426320722003136425762899652531974839391071893983013696279981183473287139122567_wp
    real(wp),parameter :: b98 = .3568084804936544335695235850557119302907377562082027313083048400646746448533926641757_wp
    real(wp),parameter :: b101 = .4520001981236963034666158997158145581787527592168118014437521620068508165483327850010e-1_wp
    real(wp),parameter :: b102 = 0
    real(wp),parameter :: b103 = 0
    real(wp),parameter :: b104 = 0
    real(wp),parameter :: b105 = 0
    real(wp),parameter :: b106 = .4002799368963269424453448726036172129707286035255084033530848251945047638855433381069e-2_wp
    real(wp),parameter :: b107 = .1601500636916154116141648988130496051518236572818790559811874215346710453554656819808_wp
    real(wp),parameter :: b108 = -.6948308961135587660550768108803798325673156163094392313978348939414189769206882520960e-3_wp
    real(wp),parameter :: b109 = -.3028824867077836529758526871817473127675272547494279018764594054601128969315645111849e-2_wp
    real(wp),parameter :: b111 = .6488088874156601542550101009631282425235066317670643554509232342293350370832201888741e-1_wp
    real(wp),parameter :: b112 = 0
    real(wp),parameter :: b113 = 0
    real(wp),parameter :: b114 = 0
    real(wp),parameter :: b115 = 0
    real(wp),parameter :: b116 = 0
    real(wp),parameter :: b117 = 0
    real(wp),parameter :: b118 = .1798223672519538764980931982771653336226267689019203140039759418977874942481337275737_wp
    real(wp),parameter :: b119 = -.1379038853795374658459375939872891452260918957641627984236926962149885290906676389842e-2_wp
    real(wp),parameter :: b1110 = .2933177796677054638696786924539291988218053269602802970971915900988371009464085558176_wp
    real(wp),parameter :: b121 = .6203009498964512152733682602483252800544465468027723696513649660598713938403105468563e-1_wp
    real(wp),parameter :: b122 = 0
    real(wp),parameter :: b123 = 0
    real(wp),parameter :: b124 = 0
    real(wp),parameter :: b125 = 0
    real(wp),parameter :: b126 = 0
    real(wp),parameter :: b127 = 0
    real(wp),parameter :: b128 = 0
    real(wp),parameter :: b129 = .1852494271216683106566732412801252654028709112506185661813162013190347158018146409385e-1_wp
    real(wp),parameter :: b1210 = .3097644836067189024326911104871249008925595428403487952252886015862314098920345337684_wp
    real(wp),parameter :: b1211 = .3047442431634932257738771099290867502123911521608864408553409749254117720052567713637_wp
    real(wp),parameter :: b131 = .6201972260916893080711420739558975045994481991404016270942831042621296675332534575866e-1_wp
    real(wp),parameter :: b132 = 0
    real(wp),parameter :: b133 = 0
    real(wp),parameter :: b134 = 0
    real(wp),parameter :: b135 = 0
    real(wp),parameter :: b136 = 0
    real(wp),parameter :: b137 = 0
    real(wp),parameter :: b138 = 0
    real(wp),parameter :: b139 = .9402397126314593976910409432600376284640745739392856006491693414739377176752496830254e-1_wp
    real(wp),parameter :: b1310 = .3098255331628556417346604747446843311853795538255778076358109181983411673182996479068_wp
    real(wp),parameter :: b1311 = .3036714566645797050698936883147468540428250347499812111630779317183972047785244555068_wp
    real(wp),parameter :: b1312 = .3542231151139475432144782255027699933222589423837038641980029819576720980026202135845e-1_wp
    real(wp),parameter :: b141 = .6969089344825787554650204576987666522661239866526931243865956872566266585875574389090e-1_wp
    real(wp),parameter :: b142 = 0
    real(wp),parameter :: b143 = 0
    real(wp),parameter :: b144 = 0
    real(wp),parameter :: b145 = 0
    real(wp),parameter :: b146 = 0
    real(wp),parameter :: b147 = 0
    real(wp),parameter :: b148 = 0
    real(wp),parameter :: b149 = -1.026068398848221434143517600646828432765332980244057389200313887686341269100904144289_wp
    real(wp),parameter :: b1410 = .1707111685325415936794748807994057363691702390665289122715815578989450621608097471178_wp
    real(wp),parameter :: b1411 = -.2276426317828345244640899307352560724213609963668371269120350778897604268488707131379_wp
    real(wp),parameter :: b1412 = .6469733264162559693641327687577763359541689058709531537823573009011020286605375682754_wp
    real(wp),parameter :: b1413 = .5577116561557104399547116527283603900460440287233664182020166612039997170938506894503_wp
    real(wp),parameter :: b151 = .5627122376293571949132309653445143949301687520262788651750408232740266180616570845510e-1_wp
    real(wp),parameter :: b152 = 0
    real(wp),parameter :: b153 = 0
    real(wp),parameter :: b154 = 0
    real(wp),parameter :: b155 = 0
    real(wp),parameter :: b156 = 0
    real(wp),parameter :: b157 = 0
    real(wp),parameter :: b158 = 0
    real(wp),parameter :: b159 = 0
    real(wp),parameter :: b1510 = 0
    real(wp),parameter :: b1511 = .3678983222850209699517524982831781294240536196092600812178960084767342167600554097921_wp
    real(wp),parameter :: b1512 = -.5195587208803483299356547234278383561187847326601986146678158388835223877380213814230e-1_wp
    real(wp),parameter :: b1513 = .2204300024722190946507889136186101583051939445678464072502258435712922894837870998433_wp
    real(wp),parameter :: b1514 = .2968005951240623144694870329678754295474339460340838345212732256969638774668653008510_wp
    real(wp),parameter :: b161 = .5634846497340329026169333310752166966844307632715132381117641490389593628519826878531e-1_wp
    real(wp),parameter :: b162 = 0
    real(wp),parameter :: b163 = 0
    real(wp),parameter :: b164 = 0
    real(wp),parameter :: b165 = 0
    real(wp),parameter :: b166 = 0
    real(wp),parameter :: b167 = 0
    real(wp),parameter :: b168 = 0
    real(wp),parameter :: b169 = 0
    real(wp),parameter :: b1610 = 0
    real(wp),parameter :: b1611 = .2154493134975348678816921094538633047666156318059356569827382889198219015749935439823_wp
    real(wp),parameter :: b1612 = -.1828038597798339032767249818891038294211236603774749922977117870744959009671532402296_wp
    real(wp),parameter :: b1613 = .1115035565964385314512820745510850004252820806597292903369723111122469246811072381550_wp
    real(wp),parameter :: b1614 = .2961573176911709124764789850722285989373150621850530267865487559044685065325051431373_wp
    real(wp),parameter :: b1615 = -.2863963826611358869027903247801186359263232936399877477557034100902322184153802061016e-1_wp
    real(wp),parameter :: b171 = .5768994198724289036171882589079551200329053317074015736866935607317368058263947892388e-1_wp
    real(wp),parameter :: b172 = 0
    real(wp),parameter :: b173 = 0
    real(wp),parameter :: b174 = 0
    real(wp),parameter :: b175 = 0
    real(wp),parameter :: b176 = 0
    real(wp),parameter :: b177 = 0
    real(wp),parameter :: b178 = 0
    real(wp),parameter :: b179 = 0
    real(wp),parameter :: b1710 = 0
    real(wp),parameter :: b1711 = -.3339594153573746138779657254239687461940518464286089760815893952034552574434211037540e-1_wp
    real(wp),parameter :: b1712 = .4007337036063909756966778946251460958783204419363211999098235920936595305217080757252_wp
    real(wp),parameter :: b1713 = -.2804910522268063380484013221059916164179105523752052409880171087961359378301286147028_wp
    real(wp),parameter :: b1714 = .2860659865450288214940338791467983472419468134474127404785360793939950768656398860364_wp
    real(wp),parameter :: b1715 = .3300830129889967277869093004777042822760017064369869377090469066296160876867312632942_wp
    real(wp),parameter :: b1716 = .2301529673804798738783942876721443945820832328652238487255898956623302521178929654255_wp
    real(wp),parameter :: b181 = .4628349579719024427317175399131733273511297084139209516535128748975595102348466963346e-1_wp
    real(wp),parameter :: b182 = 0
    real(wp),parameter :: b183 = 0
    real(wp),parameter :: b184 = 0
    real(wp),parameter :: b185 = 0
    real(wp),parameter :: b186 = 0
    real(wp),parameter :: b187 = 0
    real(wp),parameter :: b188 = 0
    real(wp),parameter :: b189 = 0
    real(wp),parameter :: b1810 = 0
    real(wp),parameter :: b1811 = 1.630450107397975249685538886816283747362802271995334921665631466207136001357816066251_wp
    real(wp),parameter :: b1812 = -1.519389027110125818833104438792078888366805294218517984315935323689939634346785541266_wp
    real(wp),parameter :: b1813 = 1.451177905612818398158204503782829478002025896171620079182131464615748573245045066933_wp
    real(wp),parameter :: b1814 = .1115858300886827251108540309774994639792980968086331408518485993004009147049002886661_wp
    real(wp),parameter :: b1815 = -.6927910549842217223999151167283689307038087310272710650839151778360002928231938996564_wp
    real(wp),parameter :: b1816 = -1.033967300639154935442351963388333197850853637180263551013370838161871418837432150763_wp
    real(wp),parameter :: b1817 = .9665192672141503474415666456096020214782059380922299431459254579911560328490208464940e-1_wp
    real(wp),parameter :: b191 = -.2030459770939010457179865043457543159890686836908949375313125798640212974640722247139e-1_wp
    real(wp),parameter :: b192 = 0
    real(wp),parameter :: b193 = 0
    real(wp),parameter :: b194 = 0
    real(wp),parameter :: b195 = 0
    real(wp),parameter :: b196 = 0
    real(wp),parameter :: b197 = 0
    real(wp),parameter :: b198 = 0
    real(wp),parameter :: b199 = 0
    real(wp),parameter :: b1910 = 0
    real(wp),parameter :: b1911 = 5.403684273875123798630962329849980016426640735397844799635504992740247972317943415738_wp
    real(wp),parameter :: b1912 = -6.547667878535769638577053800139325829335808492655639368966216657633158529493394018748_wp
    real(wp),parameter :: b1913 = 7.013821350363578373165272471360941682153886129098573011050120208706394847104611280659_wp
    real(wp),parameter :: b1914 = .1870129209851220070917089872488052512865730910524666003812699685036249110639968876355_wp
    real(wp),parameter :: b1915 = -3.586221630318466310345416445340680667399109825679889726092108094877121092019568725906_wp
    real(wp),parameter :: b1916 = -2.796807192734124892412016302005787005622243784243008121795193206999340165569504484433_wp
    real(wp),parameter :: b1917 = .5350966429298067141403865717192268565615462868769074262461851475076297049847405900305_wp
    real(wp),parameter :: b1918 = .2145339995655007603736693460786443623723658483597653704106017814468824679556810775581_wp
    real(wp),parameter :: b201 = -.1742485918290044205614554854196710455479504553190752005798475964042862615568249931979e-1_wp
    real(wp),parameter :: b202 = 0
    real(wp),parameter :: b203 = 0
    real(wp),parameter :: b204 = 0
    real(wp),parameter :: b205 = 0
    real(wp),parameter :: b206 = 0
    real(wp),parameter :: b207 = 0
    real(wp),parameter :: b208 = 0
    real(wp),parameter :: b209 = 0
    real(wp),parameter :: b2010 = 0
    real(wp),parameter :: b2011 = 1.624681387971290566033459737855869526767754661771120924184664522930252560057487715917_wp
    real(wp),parameter :: b2012 = -2.354176134519334554629765394537690359953484690980299665915399924759933249024433201062_wp
    real(wp),parameter :: b2013 = 2.681464673816256355689049206851920401601461184001531947851310974409875440156813229521_wp
    real(wp),parameter :: b2014 = -.3769156802887539884394023979903229235497658274041343397278286356588457999273698837568e-1_wp
    real(wp),parameter :: b2015 = -1.416558131026522083720211494413873707821785422770559705231460034487688457012423059506_wp
    real(wp),parameter :: b2016 = -.6949353599391523142082433283060914200352029386004581317085753581959007030993120101233_wp
    real(wp),parameter :: b2017 = .2176483372495038460388837827561689520877604442565239353657200502175872181787882826377_wp
    real(wp),parameter :: b2018 = .2697602441781337550421531890140860774876933975825845928799740196946913466182988994464_wp
    real(wp),parameter :: b2019 = -.3010654404469406281446127247743328561741318472966391363598843663110743566238944373807e-1_wp
    real(wp),parameter :: b211 = -.5403069633602259108646269644934256957551404599412561453767268074870783669925498016457e-2_wp
    real(wp),parameter :: b212 = 0
    real(wp),parameter :: b213 = 0
    real(wp),parameter :: b214 = 0
    real(wp),parameter :: b215 = 0
    real(wp),parameter :: b216 = 0
    real(wp),parameter :: b217 = 0
    real(wp),parameter :: b218 = 0
    real(wp),parameter :: b219 = 0
    real(wp),parameter :: b2110 = 0
    real(wp),parameter :: b2111 = .4139045320778312241379049484992085439613081939853923740529945110841130228938567470608_wp
    real(wp),parameter :: b2112 = -.6225569512435124987325753043472442306739960604338314680453027342793707050930404307846_wp
    real(wp),parameter :: b2113 = .8734771329093316867790902545034150641322140869164058836843523721471009826005087793734_wp
    real(wp),parameter :: b2114 = -.4460631186034505809671866759278699513326438900110859254673757293843834232949457511945_wp
    real(wp),parameter :: b2115 = -.4770024428429503264261245457799687730095056727976098627656199278228457764648229288301_wp
    real(wp),parameter :: b2116 = .3350194037993700956622540574882038342150917200081189888351569726131444997493018243069_wp
    real(wp),parameter :: b2117 = .7504248468668828335131027827336082864172905052946156155171125652518370192017019285897e-1_wp
    real(wp),parameter :: b2118 = .2868602398412017179778096520767496699075337622961506198525486757858393614463244192166_wp
    real(wp),parameter :: b2119 = -.3105278876109550093409796945938508430591700006793679865817465482241195091701064106306_wp
    real(wp),parameter :: b2120 = .5980742748031364157936282749945094293392110393839614650861749046057502083978432774641_wp
    real(wp),parameter :: b221 = -.8040637783629852670386253408427477129839038039895222958327102384382347854387228140586e-1_wp
    real(wp),parameter :: b222 = 0
    real(wp),parameter :: b223 = 0
    real(wp),parameter :: b224 = 0
    real(wp),parameter :: b225 = 0
    real(wp),parameter :: b226 = 0
    real(wp),parameter :: b227 = 0
    real(wp),parameter :: b228 = 0
    real(wp),parameter :: b229 = 0
    real(wp),parameter :: b2210 = 0
    real(wp),parameter :: b2211 = 14.09003825209192451104055534327205557907144506304998967309786817952326707840468519855_wp
    real(wp),parameter :: b2212 = -16.51095452142836096439418777253125778161824307435047496403861839615632359732680906324_wp
    real(wp),parameter :: b2213 = 17.13478919889388959512706015258229430254357775856472236246578628732774154201548556949_wp
    real(wp),parameter :: b2214 = .5162416754963992955466349330401373481541758099272782749146402071285515617224748797349_wp
    real(wp),parameter :: b2215 = -8.738239054447235555769280431947112005536670926201592356606646679961357074323383907433_wp
    real(wp),parameter :: b2216 = -7.928763291217204762997981382770542063580344552643223907022402362361688687161033776950_wp
    real(wp),parameter :: b2217 = 1.299607095883618341404631332785395816146491714221581325131902479964394904680531866743_wp
    real(wp),parameter :: b2218 = .2985974887158637819626343538772133330364013942516888589723607063393086453938726507764_wp
    real(wp),parameter :: b2219 = .2306899737777788621463004228318904282554044640133218561071121312583287783148540854557_wp
    real(wp),parameter :: b2220 = -.4162475612748372571090619229550372319146077330570330975546626972398912191287912992344_wp
    real(wp),parameter :: b2221 = .3046152792810969736488023713665193266324974471792934717484737567812431778577335699415_wp
    real(wp),parameter :: b231 = .5188944067774664913134965881842338637465167538426812594601866401917453125373780443288e-1_wp
    real(wp),parameter :: b232 = 0
    real(wp),parameter :: b233 = 0
    real(wp),parameter :: b234 = 0
    real(wp),parameter :: b235 = 0
    real(wp),parameter :: b236 = 0
    real(wp),parameter :: b237 = 0
    real(wp),parameter :: b238 = 0
    real(wp),parameter :: b239 = 0
    real(wp),parameter :: b2310 = 0
    real(wp),parameter :: b2311 = -1.283277980375441090269741488559725341590001661469261097548681503281119502510781510671_wp
    real(wp),parameter :: b2312 = 3.197031115178587180279278821590157014124289302868827178494105498343258944734856975984_wp
    real(wp),parameter :: b2313 = -4.925730995528770332329418001427019544036228452951763539493235207541539355193215629707_wp
    real(wp),parameter :: b2314 = .4556842843249686735202095835169824346358162044635023939110786683727545600555997662146_wp
    real(wp),parameter :: b2315 = 2.672937738349530686900845106745665406421870114003531687061025215633257210063413548721_wp
    real(wp),parameter :: b2316 = -1.545486350471516022666399706441233218275810264466541860777249858172892030485892221263_wp
    real(wp),parameter :: b2317 = -.3721987962685219677888284678662933875128991354125828036313277053504464873414697307303_wp
    real(wp),parameter :: b2318 = .1275610106742831413187042273550361954603324559021122028077942983877428065036205767992_wp
    real(wp),parameter :: b2319 = 1.955690675673468524392185835321028843337157689840864927642297223731192314442172030656_wp
    real(wp),parameter :: b2320 = .3191081222387442241092001119049884652557628641382080948098270187047261996513219724153_wp
    real(wp),parameter :: b2321 = 1.121546003224352473639037341460033868464982074905347279402236597305066991912121393629_wp
    real(wp),parameter :: b2322 = -.8107564764218937636489023156262164031679294933798974588613267897756930246095048550015_wp
    real(wp),parameter :: b241 = -.7011708297168838093320073271824231811555516376839339889768780533772204046098056778645e-2_wp
    real(wp),parameter :: b242 = 0
    real(wp),parameter :: b243 = 0
    real(wp),parameter :: b244 = 0
    real(wp),parameter :: b245 = 0
    real(wp),parameter :: b246 = 0
    real(wp),parameter :: b247 = 0
    real(wp),parameter :: b248 = 0
    real(wp),parameter :: b249 = 0
    real(wp),parameter :: b2410 = 0
    real(wp),parameter :: b2411 = 1.191430397515656948185785084406871641958062497543742234300725797836945018037271858084_wp
    real(wp),parameter :: b2412 = -1.828224802245176636695187739870698589078178771180859212263094110856159464053843046402_wp
    real(wp),parameter :: b2413 = 2.587660024998591526533741070898054580235728551249558749988014852858819613806922059203_wp
    real(wp),parameter :: b2414 = -.4609429033777800137012176300671539646289103914701490781289612317160305665226273519358_wp
    real(wp),parameter :: b2415 = -1.325931501693639442980805419646217538577342867896123447497177267595031990967256304102_wp
    real(wp),parameter :: b2416 = .5048216427266366058315103764214691413381451411460686651485736358879672359427752204445_wp
    real(wp),parameter :: b2417 = .2194255008759965600775266235750828151587685421380286198819468613425068758315070811985_wp
    real(wp),parameter :: b2418 = .2580087037350744276193737624019399564007443016182003098695556923721330317034140914297_wp
    real(wp),parameter :: b2419 = -.7820644366919761953621998768725631028114098091524728093576852041863328545044120664888_wp
    real(wp),parameter :: b2420 = .4981236112353889237817573857026156554544669576903254922783923277260844155420240361688_wp
    real(wp),parameter :: b2421 = -.2818659230370184110487250730841700637198882434735413160224991411608463306147813033701_wp
    real(wp),parameter :: b2422 = .1975630649870201406787729762339550948339749510061280462627577046361617264269398213974_wp
    real(wp),parameter :: b2423 = -.2737789055818164357320317060050411797735340939936712549180293178469667438122671145178e-1_wp
    real(wp),parameter :: b251 = .6971223157985491441889543463058968335832800922315696226973214263599946413187133757441e-1_wp
    real(wp),parameter :: b252 = 0
    real(wp),parameter :: b253 = 0
    real(wp),parameter :: b254 = 0
    real(wp),parameter :: b255 = 0
    real(wp),parameter :: b256 = 0
    real(wp),parameter :: b257 = 0
    real(wp),parameter :: b258 = 0
    real(wp),parameter :: b259 = 0
    real(wp),parameter :: b2510 = 0
    real(wp),parameter :: b2511 = .8485072933915981598683071881067011723176540155630964375777672731824955304237165117739_wp
    real(wp),parameter :: b2512 = -1.504790599380062029550479292004278750516583328735111485609300227641118522265117453228_wp
    real(wp),parameter :: b2513 = 2.971661962491066226601618850699945230551598000336030549022007199066569229936095601265_wp
    real(wp),parameter :: b2514 = .2189212281135321618639152301919425667361998064642432899287016486436734890226912997813_wp
    real(wp),parameter :: b2515 = -1.363504916136495083863504916136495083863504916136495083863504916136495083863504916136_wp
    real(wp),parameter :: b2516 = 1.535256651208129648912110087693939955244136470648889182343875125782209378315005545886_wp
    real(wp),parameter :: b2517 = .2255130922481182477299374148203865138344993252214597927624562120740486384043245934395_wp
    real(wp),parameter :: b2518 = -.1252487725356194666272520523138563196763786809348751670013078003838735310461344402659_wp
    real(wp),parameter :: b2519 = -1.433947471579772638181105448843590748725989807918463347706781654253234025872206977656_wp
    real(wp),parameter :: b2520 = -.2653692486145121794045624436138677664647506121922928212398504961979636551101946127078_wp
    real(wp),parameter :: b2521 = -1.057083353895191276322695714066286715782319071838176809244435034528746563196511305239_wp
    real(wp),parameter :: b2522 = .5958992342288874679818860848371024166473432615899828093938621236679595663138268737944_wp
    real(wp),parameter :: b2523 = .5640763974002360165840684759509622393779174110889642877127279949714649497561235682228e-1_wp
    real(wp),parameter :: b2524 = .2280650291404422449145227284026716224019757875996592625955056045913295898305255848960_wp
    real(wp),parameter :: b261 = .5847724407501707800509964027663995519804001870472850379216053933603321448937320382873e-1_wp
    real(wp),parameter :: b262 = 0
    real(wp),parameter :: b263 = 0
    real(wp),parameter :: b264 = 0
    real(wp),parameter :: b265 = 0
    real(wp),parameter :: b266 = 0
    real(wp),parameter :: b267 = 0
    real(wp),parameter :: b268 = 0
    real(wp),parameter :: b269 = 0
    real(wp),parameter :: b2610 = 0
    real(wp),parameter :: b2611 = .8482433234853130399472226781516394610296267520150369993929425712656780053034142790857_wp
    real(wp),parameter :: b2612 = -1.818778337229596341437999685014432821665031881560984595319320029619105771387827231293_wp
    real(wp),parameter :: b2613 = 3.360950949051261986329911895518526801628857587757616068272635548399517354842337530945_wp
    real(wp),parameter :: b2614 = .1508779656009112841743785043192199515750215073024373320382690391467484988490972888214_wp
    real(wp),parameter :: b2615 = -1.535625970104296058211732653128688607443092018291498032244748403611478465692632705036_wp
    real(wp),parameter :: b2616 = 1.689710648596447426456373308903510510459502815029443327222583037297243278873766013804_wp
    real(wp),parameter :: b2617 = .2753763916359371389555315070580325732352897337628204267122697551960198539312003364345_wp
    real(wp),parameter :: b2618 = -.8663993203603434555880694559625317699151588603507943515655135799742692170332074641456e-1_wp
    real(wp),parameter :: b2619 = -1.514109179619982965788582168222695819982704233094238060182049846398573819503074487947_wp
    real(wp),parameter :: b2620 = -.2731051919432418524342072455220752640028637692479133311886846927223922900302148881433_wp
    real(wp),parameter :: b2621 = -.7386383863141516263346011834004149434481207761507984915319537912098408949133636097415_wp
    real(wp),parameter :: b2622 = .6343578675072522696688679414505709253269883958978404947025376176286935002340724649840_wp
    real(wp),parameter :: b2623 = .2219915268009860510634252213682523269888255016876797090075141605397073674221546364928e-1_wp
    real(wp),parameter :: b2624 = -.7329654538493563887779811693040477761888079625817917741084140276508628003504291297654e-1_wp
    real(wp),parameter :: b2625 = 0

    real(wp),parameter :: c1 = .2560300422448818640056973733027381517733824076443136953634798402199579680290993667785e-1_wp
    real(wp),parameter :: c2 = 0
    real(wp),parameter :: c3 = 0
    real(wp),parameter :: c4 = 0
    real(wp),parameter :: c5 = 0
    real(wp),parameter :: c6 = 0
    real(wp),parameter :: c7 = 0
    real(wp),parameter :: c8 = 0
    real(wp),parameter :: c9 = 0
    real(wp),parameter :: c10 = 0
    real(wp),parameter :: c11 = 0
    real(wp),parameter :: c12 = 0
    real(wp),parameter :: c13 = 0
    real(wp),parameter :: c14 = 0
    real(wp),parameter :: c15 = .4518777011767389121126592910724741440985372571238852320538195476794258345326944681173_wp
    real(wp),parameter :: c16 = .1571644835050410342525892218785914023945906781166129113400861908694541485546171465980_wp
    real(wp),parameter :: c17 = 1.587056574024945279877130621604901273153398820663935100541595105228850005414302839405_wp
    real(wp),parameter :: c18 = .1457164101207840798598697293357824068221719394674815441268910066025070033320249231167_wp
    real(wp),parameter :: c19 = .6756690992874993365704777906957166121534838869284524431811474400358351375791575269125e-1_wp
    real(wp),parameter :: c20 = .2263702829370636470261171189724674060331565741885612750461048803721449531302274282834_wp
    real(wp),parameter :: c21 = .9585657493693419066100150421626762898713313017234629340723230525017707321578125838875_wp
    real(wp),parameter :: c22 = -.5737510952253720251243392241397701824182886145229967660377157325379672132470673999182e-1_wp
    real(wp),parameter :: c23 = -.8207294363127299380342111455087248111504033000422844250308430626629630466499721876523_wp
    real(wp),parameter :: c24 = -.8199946795810496733757980690172928576129340372072025495094374198905045106663577622820_wp
    real(wp),parameter :: c25 = -.9218218898708361658735554044867437117607070020394289598912304554724677090414683888511_wp
    real(wp),parameter :: c26 = 0

    real(wp),parameter :: d1 = .2611453567558047174327245053028336466626806615883961914681281336363189972915245544873e-1_wp
    real(wp),parameter :: d2 = 0
    real(wp),parameter :: d3 = 0
    real(wp),parameter :: d4 = 0
    real(wp),parameter :: d5 = 0
    real(wp),parameter :: d6 = 0
    real(wp),parameter :: d7 = 0
    real(wp),parameter :: d8 = 0
    real(wp),parameter :: d9 = 0
    real(wp),parameter :: d10 = 0
    real(wp),parameter :: d11 = 0
    real(wp),parameter :: d12 = 0
    real(wp),parameter :: d13 = 0
    real(wp),parameter :: d14 = 0
    real(wp),parameter :: d15 = .6130789590044024420375347235457514503322168389154611016924109584567258667669428084328e-1_wp
    real(wp),parameter :: d16 = .2225176387207912085954871898906978473895949131807976650805573845853338932568562016822_wp
    real(wp),parameter :: d17 = -.6452541374509687528609607022885324791758461431623241957335732784511248267767614258399_wp
    real(wp),parameter :: d18 = .1428813432326775237316522259538641672308222297639395170211160113665720641033224975230_wp
    real(wp),parameter :: d19 = .1998503703511764116053495106126866164319741023221332167221899831502033833897226507607e-1_wp
    real(wp),parameter :: d20 = .2271215662186804051109071916317932366756181626009529050272447503184463361820614170681_wp
    real(wp),parameter :: d21 = .3555653412963184679461951194576569724124866534036012519458466746801575903667969215385_wp
    real(wp),parameter :: d22 = -.5084144204984584161908312330177564954877067154001079429559777704302946848644928736984e-1_wp
    real(wp),parameter :: d23 = .3451056668568212037936618506732414889276840710489545476476791628380007277211978725437_wp
    real(wp),parameter :: d24 = -.1182965488838884338743861432044520655991039617854064994056837668531638997325328536863_wp
    real(wp),parameter :: d25 = 0
    real(wp),parameter :: d26 = .4137931034482758620689655172413793103448275862068965517241379310344827586206896551724_wp

    real(wp),parameter :: e1  = c1  - d1
    real(wp),parameter :: e2  = c2  - d2
    real(wp),parameter :: e3  = c3  - d3
    real(wp),parameter :: e4  = c4  - d4
    real(wp),parameter :: e5  = c5  - d5
    real(wp),parameter :: e6  = c6  - d6
    real(wp),parameter :: e7  = c7  - d7
    real(wp),parameter :: e8  = c8  - d8
    real(wp),parameter :: e9  = c9  - d9
    real(wp),parameter :: e10 = c10 - d10
    real(wp),parameter :: e11 = c11 - d11
    real(wp),parameter :: e12 = c12 - d12
    real(wp),parameter :: e13 = c13 - d13
    real(wp),parameter :: e14 = c14 - d14
    real(wp),parameter :: e15 = c15 - d15
    real(wp),parameter :: e16 = c16 - d16
    real(wp),parameter :: e17 = c17 - d17
    real(wp),parameter :: e18 = c18 - d18
    real(wp),parameter :: e19 = c19 - d19
    real(wp),parameter :: e20 = c20 - d20
    real(wp),parameter :: e21 = c21 - d21
    real(wp),parameter :: e22 = c22 - d22
    real(wp),parameter :: e23 = c23 - d23
    real(wp),parameter :: e24 = c24 - d24
    real(wp),parameter :: e25 = c25 - d25
    real(wp),parameter :: e26 = c26 - d26

    associate (f1 => me%funcs(:,1), &
               f2 => me%funcs(:,2), &
               f3 => me%funcs(:,3), &
               f4 => me%funcs(:,4), &
               f5 => me%funcs(:,5), &
               f6 => me%funcs(:,6), &
               f7 => me%funcs(:,7), &
               f8 => me%funcs(:,8), &
               f9 => me%funcs(:,9), &
               f10 => me%funcs(:,10), &
               f11 => me%funcs(:,11), &
               f12 => me%funcs(:,12), &
               f13 => me%funcs(:,13), &
               f14 => me%funcs(:,14), &
               f15 => me%funcs(:,15), &
               f16 => me%funcs(:,16), &
               f17 => me%funcs(:,17), &
               f18 => me%funcs(:,18), &
               f19 => me%funcs(:,19), &
               f20 => me%funcs(:,20), &
               f21 => me%funcs(:,21), &
               f22 => me%funcs(:,22), &
               f23 => me%funcs(:,23), &
               f24 => me%funcs(:,24), &
               f25 => me%funcs(:,25), &
               f26 => me%funcs(:,26))

        call me%f(t+h,   x,f1)
        call me%f(t+a2*h,x+h*(b21*f1),f2)
        call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
        call me%f(t+a4*h,x+h*(b41*f1+b42*f2+b43*f3),f4)
        call me%f(t+a5*h,x+h*(b51*f1+b52*f2+b53*f3+b54*f4),f5)
        call me%f(t+a6*h,x+h*(b61*f1+b62*f2+b63*f3+b64*f4+b65*f5),f6)
        call me%f(t+a7*h,x+h*(b71*f1+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
        call me%f(t+a8*h,x+h*(b81*f1+b82*f2+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
        call me%f(t+a9*h,x+h*(b91*f1+b92*f2+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
        call me%f(t+a10*h,x+h*(b101*f1+b102*f2+b103*f3+b104*f4+b105*f5+b106*f6+b107*f7+&
                            b108*f8+b109*f9),f10)
        call me%f(t+a11*h,x+h*(b111*f1+b112*f2+b113*f3+b114*f4+b115*f5+b116*f6+b117*f7+&
                            b118*f8+b119*f9+b1110*f10),f11)
        call me%f(t+a12*h,x+h*(b121*f1+b122*f2+b123*f3+b124*f4+b125*f5+b126*f6+b127*f7+&
                            b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
        call me%f(t+a13*h,x+h*(b131*f1+b132*f2+b133*f3+b134*f4+b135*f5+b136*f6+b137*f7+&
                            b138*f8+b139*f9+b1310*f10+b1311*f11+b1312*f12),f13)
        call me%f(t+a14*h,x+h*(b141*f1+b142*f2+b143*f3+b144*f4+b145*f5+b146*f6+b147*f7+&
                            b148*f8+b149*f9+b1410*f10+b1411*f11+b1412*f12+b1413*f13),f14)
        call me%f(t+a15*h,x+h*(b151*f1+b152*f2+b153*f3+b154*f4+b155*f5+b156*f6+b157*f7+&
                            b158*f8+b159*f9+b1510*f10+b1511*f11+b1512*f12+b1513*f13+&
                            b1514*f14),f15)
        call me%f(t+a16*h,x+h*(b161*f1+b162*f2+b163*f3+b164*f4+b165*f5+b166*f6+b167*f7+&
                            b168*f8+b169*f9+b1610*f10+b1611*f11+b1612*f12+b1613*f13+&
                            b1614*f14+b1615*f15),f16)
        call me%f(t+a17*h,x+h*(b171*f1+b172*f2+b173*f3+b174*f4+b175*f5+b176*f6+b177*f7+&
                            b178*f8+b179*f9+b1710*f10+b1711*f11+b1712*f12+b1713*f13+&
                            b1714*f14+b1715*f15+b1716*f16),f17)
        call me%f(t+a18*h,x+h*(b181*f1+b182*f2+b183*f3+b184*f4+b185*f5+b186*f6+b187*f7+&
                            b188*f8+b189*f9+b1810*f10+b1811*f11+b1812*f12+b1813*f13+&
                            b1814*f14+b1815*f15+b1816*f16+b1817*f17),f18)
        call me%f(t+a19*h,x+h*(b191*f1+b192*f2+b193*f3+b194*f4+b195*f5+b196*f6+b197*f7+&
                            b198*f8+b199*f9+b1910*f10+b1911*f11+b1912*f12+b1913*f13+&
                            b1914*f14+b1915*f15+b1916*f16+b1917*f17+b1918*f18),f19)
        call me%f(t+a20*h,x+h*(b201*f1+b202*f2+b203*f3+b204*f4+b205*f5+b206*f6+b207*f7+&
                            b208*f8+b209*f9+b2010*f10+b2011*f11+b2012*f12+b2013*f13+&
                            b2014*f14+b2015*f15+b2016*f16+b2017*f17+b2018*f18+b2019*f19),f20)
        call me%f(t+a21*h,x+h*(b211*f1+b212*f2+b213*f3+b214*f4+b215*f5+b216*f6+b217*f7+&
                            b218*f8+b219*f9+b2110*f10+b2111*f11+b2112*f12+b2113*f13+&
                            b2114*f14+b2115*f15+b2116*f16+b2117*f17+b2118*f18+b2119*f19+&
                            b2120*f20),f21)
        call me%f(t+a22*h,x+h*(b221*f1+b222*f2+b223*f3+b224*f4+b225*f5+b226*f6+b227*f7+&
                            b228*f8+b229*f9+b2210*f10+b2211*f11+b2212*f12+b2213*f13+&
                            b2214*f14+b2215*f15+b2216*f16+b2217*f17+b2218*f18+b2219*f19+&
                            b2220*f20+b2221*f21),f22)
        call me%f(t+a23*h,x+h*(b231*f1+b232*f2+b233*f3+b234*f4+b235*f5+b236*f6+b237*f7+&
                            b238*f8+b239*f9+b2310*f10+b2311*f11+b2312*f12+b2313*f13+&
                            b2314*f14+b2315*f15+b2316*f16+b2317*f17+b2318*f18+b2319*f19+&
                            b2320*f20+b2321*f21+b2322*f22),f23)
        call me%f(t+a24*h,x+h*(b241*f1+b242*f2+b243*f3+b244*f4+b245*f5+b246*f6+b247*f7+&
                            b248*f8+b249*f9+b2410*f10+b2411*f11+b2412*f12+b2413*f13+&
                            b2414*f14+b2415*f15+b2416*f16+b2417*f17+b2418*f18+b2419*f19+&
                            b2420*f20+b2421*f21+b2422*f22+b2423*f23),f24)
        call me%f(t+h,    x+h*(b251*f1+b252*f2+b253*f3+b254*f4+b255*f5+b256*f6+b257*f7+&
                            b258*f8+b259*f9+b2510*f10+b2511*f11+b2512*f12+b2513*f13+&
                            b2514*f14+b2515*f15+b2516*f16+b2517*f17+b2518*f18+b2519*f19+&
                            b2520*f20+b2521*f21+b2522*f22+b2523*f23+b2524*f24),f25)
        call me%f(t+h,    x+h*(b261*f1+b262*f2+b263*f3+b264*f4+b265*f5+b266*f6+b267*f7+&
                            b268*f8+b269*f9+b2610*f10+b2611*f11+b2612*f12+b2613*f13+&
                            b2614*f14+b2615*f15+b2616*f16+b2617*f17+b2618*f18+b2619*f19+&
                            b2620*f20+b2621*f21+b2622*f22+b2623*f23+b2624*f24+b2625*f25),f26)

        xf = x+h*(c1*f1+c2*f2+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11+c12*f12+&
                c13*f13+c14*f14+c15*f15+c16*f16+c17*f17+c18*f18+c19*f19+c20*f20+c21*f21+c22*f22+&
                c23*f23+c24*f24+c25*f25+c26*f26)

        xerr = h*(e1*f1+e2*f2+e3*f3+e4*f4+e5*f5+e6*f6+e7*f7+e8*f8+e9*f9+e10*f10+e11*f11+e12*f12+&
                e13*f13+e14*f14+e15*f15+e16*f16+e17*f17+e18*f18+e19*f19+e20*f20+e21*f21+e22*f22+&
                e23*f23+e24*f24+e25*f25+e26*f26)

    end associate

    end procedure rks1110a
!*****************************************************************************************

!*****************************************************************************************
    end submodule rklib_variable_steps
!*****************************************************************************************