!*****************************************************************************************
!>
!  Fixed-step RK formulas.

    submodule(rklib_module) rklib_fixed_steps

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Euler (1st order) integration method.

    module procedure euler

    real(wp),dimension(me%n) :: f1

    call me%f(t,x,f1)

    xf = x + h*f1

    end procedure euler
!*****************************************************************************************

!*****************************************************************************************
!>
!  Midpoint (2nd order) integration method.

    module procedure midpoint

    real(wp),dimension(me%n) :: f1,f2

    call me%f(t,x,f1)
    call me%f(t+0.5_wp*h,x+0.5_wp*h*f1,f2)

    xf = x + h*f2

    end procedure midpoint
!*****************************************************************************************

!*****************************************************************************************
!>
!  Heun's (2nd order) integration method

    module procedure heun

    real(wp),dimension(me%n) :: f1,f2

    call me%f(t,x,f1)
    call me%f(t+h,x+h*f1,f2)

    xf = x + 0.5_wp*h*(f1+f2)

    end procedure heun
!*****************************************************************************************

!*****************************************************************************************
!>
!  2-stage, 2nd order TVD Runge-Kutta method of Shu and Osher (1988). CFL=1.0.
!
!### Reference
!  * C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory
!    shock-capturing schemes", Journal of Computational Physics, 77, 1988, 439-471.
!    https://doi.org/10.1016/0021-9991(88)90177-5.

    module procedure rkssp22

    real(wp),dimension(me%n) :: fs

    call me%f(t, x, fs)
    xf = x + h*fs
    call me%f(t + h, xf, fs)

    xf = (x + xf + h*fs) / 2.0_wp

    end procedure rkssp22
!*****************************************************************************************

!*****************************************************************************************
!>
!  3rd order, 3 steps RK integration method

    module procedure rk3

    real(wp),dimension(me%n) :: f1,f2,f3

    real(wp),parameter :: a1  =  1.0_wp/6.0_wp
    real(wp),parameter :: a2  =  2.0_wp/3.0_wp
    real(wp),parameter :: a3  =  1.0_wp/6.0_wp
    real(wp),parameter :: b2  =  1.0_wp/2.0_wp
    !real(wp),parameter :: b3  =  1.0_wp
    real(wp),parameter :: c21 =  1.0_wp/2.0_wp
    !real(wp),parameter :: c31 = -1.0_wp
    real(wp),parameter :: c32 =  2.0_wp

    call me%f(t,      x,                f1)
    call me%f(t+b2*h, x+h*c21*f1,       f2)
    call me%f(t+h,    x+h*(-f1+c32*f2), f3)

    xf = x + h*( a1*f1 + a2*f2 + a3*f3 )

    end procedure rk3
!*****************************************************************************************

!*****************************************************************************************
!>
!  3-stage, 3rd order TVD Runge-Kutta method of Shu and Osher (1988). CFL=1.0.
!
!### Reference
!  * C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory
!    shock-capturing schemes", Journal of Computational Physics, 77, 1988, 439-471.
!    https://doi.org/10.1016/0021-9991(88)90177-5.

    module procedure rkssp33

    real(wp),dimension(me%n) :: fs

    call me%f(t, x, fs)
    xf = x + h*fs
    call me%f(t + h, xf, fs)
    xf = (3.0_wp*x + xf + h*fs)/4.0_wp
    call me%f(t + h/2.0_wp, xf, fs)
    xf = (x + 2.0_wp*xf + 2.0_wp*h*fs)/3.0_wp

     end procedure rkssp33
!*****************************************************************************************

!*****************************************************************************************
!>
!   5-stage, 3rd order SSP Runge-Kutta method of Spiteri and Ruuth (2005). CFL=2.65.
!
!### Reference
!   * Ruuth, Steven. "Global optimization of explicit strong-stability-preserving Runge-Kutta
!   methods." Mathematics of Computation 75.253 (2006): 183-207.
!   https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf
!
!@note the coefficients here are only 15 digits of precision.

    module procedure rkssp53

    real(wp), parameter :: a30 = 0.355909775063327_wp
    real(wp), parameter :: a32 = 0.644090224936674_wp
    real(wp), parameter :: a40 = 0.367933791638137_wp
    real(wp), parameter :: a43 = 0.632066208361863_wp
    real(wp), parameter :: a52 = 0.237593836598569_wp
    real(wp), parameter :: a54 = 0.762406163401431_wp
    real(wp), parameter :: b10 = 0.377268915331368_wp
    real(wp), parameter :: b21 = 0.377268915331368_wp
    real(wp), parameter :: b32 = 0.242995220537396_wp
    real(wp), parameter :: b43 = 0.238458932846290_wp
    real(wp), parameter :: b54 = 0.287632146308408_wp
    real(wp), parameter :: c1  = 0.377268915331368_wp
    real(wp), parameter :: c2  = 0.754537830662736_wp
    real(wp), parameter :: c3  = 0.728985661612188_wp
    real(wp), parameter :: c4  = 0.699226135931670_wp

    real(wp), dimension(me%n) :: xs, fs

    call me%f(t, x, fs)
    ! x1 as xs
    xs = x + b10*h*fs
    call me%f(t + c1*h, xs, fs)
    ! x2 as xf
    xf = xs + b21*h*fs
    call me%f(t + c2*h, xf, fs)
    ! x3 as xs
    xs = a30*x + a32*xf + b32*h*fs
    call me%f(t + c3*h, xs, fs)
    ! x4 as xs
    xs = a40*x + a43*xs + b43*h*fs
    call me%f(t + c4*h, xs, fs)

    xf = a52*xf + a54*xs + b54*h*fs

    end procedure rkssp53
!*****************************************************************************************

!*****************************************************************************************
!>
!  Take one Runge Kutta 4 integration step: `t -> t+h (x -> xf)`

    module procedure rk4

    real(wp),dimension(me%n) :: f1,f2,f3,f4

    associate (h2 => 0.5_wp*h)
        call me%f(t,x,f1)
        call me%f(t+h2,x+h2*f1,f2)
        call me%f(t+h2,x+h2*f2,f3)
        call me%f(t+h,x+h*f3,f4)
    end associate

    xf = x + h*(f1+f2+f2+f3+f3+f4)/6.0_wp

    end procedure rk4
!*****************************************************************************************

!*****************************************************************************************
!>
!  4th order Runge Kutta Shanks (4 points)
!
!### Reference
!  * E. B. Shanks, "Solutions of Differential Equations by Evaluations of Functions"
!    Math. Comp. 20 (1966).

    module procedure rks4

    real(wp),dimension(me%n) :: f0,f1,f2,f3

    real(wp),parameter :: a1  =  1.0_wp / 100.0_wp
    real(wp),parameter :: a2  =  3.0_wp / 5.0_wp
    !real(wp),parameter :: a3  =  1.0_wp
    real(wp),parameter :: c   =  1.0_wp / 70092.0_wp
    real(wp),parameter :: c0  = -179124.0_wp
    real(wp),parameter :: c1  =  200000.0_wp
    real(wp),parameter :: c2  =  40425.0_wp
    real(wp),parameter :: c3  =  8791.0_wp
    real(wp),parameter :: aa1 =  1.0_wp / 100.0_wp
    real(wp),parameter :: aa2 =  1.0_wp / 245.0_wp
    real(wp),parameter :: aa3 =  1.0_wp / 8791.0_wp
    real(wp),parameter :: b20 = -4278.0_wp
    real(wp),parameter :: b21 =  4425.0_wp
    real(wp),parameter :: b30 =  524746.0_wp
    real(wp),parameter :: b31 = -532125.0_wp
    real(wp),parameter :: b32 =  16170.0_wp

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+aa1*h*f0,f1)
    call me%f(t+a2*h,x+aa2*h*(b20*f0+b21*f1),f2)
    call me%f(t+h,   x+aa3*h*(b30*f0+b31*f1+b32*f2),f3)

    xf = x + h*c*(c0*f0+c1*f1+c2*f2+c3*f3)

    end procedure rks4
!*****************************************************************************************

!*****************************************************************************************
!>
!  4-stage, 4th order low storage non-TVD Runge-Kutta method of Jiang and Shu (1988).
!
!### Reference
!  * Method: Jiang, Guang-Shan, and Chi-Wang Shu. "Efficient implementation of weighted ENO
!    schemes." Journal of computational physics 126.1 (1996): 202-228.
!    https://ntrs.nasa.gov/api/citations/19960007052/downloads/19960007052.pdf
!  * Implementation: J. M. F. Donnert et al 2019 ApJS 241 23.
!    https://iopscience.iop.org/article/10.3847/1538-4365/ab09fb

    module procedure rkls44

    real(wp), dimension(me%n) :: xs, fs

    xf = x
    xs = -4.0_wp*x/3.0_wp
    call me%f(t, xf, fs)
    xf = x - h*fs/2.0_wp
    xs = xs + xf/3.0_wp
    call me%f(t + h/2.0_wp, xf, fs)
    xf = x - h*fs/2.0_wp
    xs = xs + 2.0_wp*xf/3.0_wp
    call me%f(t + h/2.0_wp, xf, fs)
    xf = x - h*fs
    xs = xs + xf/3.0_wp
    call me%f(t + h, xf, fs)
    xf = x - h*fs/6.0_wp
    xf = xf + xs

    end procedure rkls44
!*****************************************************************************************

!*****************************************************************************************
!>
!   5-stage, 4th order SSP Runge-Kutta method of Spiteri and Ruuth (2005). CFL=1.508.
!
!### Reference
!   * Ruuth, Steven. "Global optimization of explicit strong-stability-preserving Runge-Kutta
!   methods." Mathematics of Computation 75.253 (2006): 183-207.
!   https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf

    module procedure rkssp54

    real(wp), parameter :: b10 = 0.391752226571890_wp
    real(wp), parameter :: a20 = 0.444370493651235_wp
    real(wp), parameter :: a21 = 0.555629506348765_wp
    real(wp), parameter :: b21 = 0.368410593050371_wp
    real(wp), parameter :: a30 = 0.620101851488403_wp
    real(wp), parameter :: a32 = 0.379898148511597_wp
    real(wp), parameter :: b32 = 0.251891774271694_wp
    real(wp), parameter :: a40 = 0.178079954393132_wp
    real(wp), parameter :: a43 = 0.821920045606868_wp
    real(wp), parameter :: b43 = 0.544974750228521_wp
    real(wp), parameter :: a52 = 0.517231671970585_wp
    real(wp), parameter :: a53 = 0.096059710526147_wp
    real(wp), parameter :: b53 = 0.063692468666290_wp
    real(wp), parameter :: a54 = 0.386708617503269_wp
    real(wp), parameter :: b54 = 0.226007483236906_wp
    real(wp), parameter :: c1  = 0.391752226571890_wp
    real(wp), parameter :: c2  = 0.586079689311540_wp
    real(wp), parameter :: c3  = 0.474542363121400_wp
    real(wp), parameter :: c4  = 0.935010630967653_wp

    real(wp), dimension(me%n) :: x2, x3, f3, fs

    call me%f(t, x, fs)

    ! x2 as x1
    x2 = x + b10*h*fs
    call me%f(t + c1*h, x2, fs)

    x2 = a20*x + a21*x2 + b21*h*fs
    call me%f(t + c2*h, x2, fs)

    x3 = a30*x + a32*x2 + b32*h*fs
    call me%f(t + c3*h, x3, f3)

    ! xf as x4
    xf = a40*x + a43*x3 + b43*h*f3
    call me%f(t + c4*h, xf, fs)

    xf = a52*x2 + a53*x3 + b53*h*f3 + a54*xf + b54*h*fs

    end procedure rkssp54
!*****************************************************************************************

!*****************************************************************************************
!>
!  Runge Kutta Shanks (5th order)
!
!### Reference
!  * E. B. Shanks, "Solutions of Differential Equations by Evaluations of Functions"
!    Math. Comp. 20 (1966).

    module procedure rks5

    real(wp),dimension(me%n) :: f0,f1,f2,f3,f4

    real(wp),parameter :: a1  =  1.0_wp / 9000.0_wp
    real(wp),parameter :: a2  =  3.0_wp / 10.0_wp
    real(wp),parameter :: a3  =  3.0_wp / 4.0_wp
    !real(wp),parameter :: a4  =  1.0_wp
    real(wp),parameter :: c   =  1.0_wp / 1134.0_wp
    real(wp),parameter :: c0  =  105.0_wp
    real(wp),parameter :: c2  =  500.0_wp
    real(wp),parameter :: c3  =  448.0_wp
    real(wp),parameter :: c4  =  81.0_wp
    real(wp),parameter :: aa1 =  1.0_wp / 9000.0_wp
    real(wp),parameter :: aa2 =  1.0_wp / 10.0_wp
    real(wp),parameter :: aa3 =  1.0_wp / 8.0_wp
    real(wp),parameter :: aa4 =  1.0_wp / 81.0_wp
    real(wp),parameter :: b20 = -4047.0_wp
    real(wp),parameter :: b21 =  4050.0_wp
    real(wp),parameter :: b30 =  20241.0_wp
    real(wp),parameter :: b31 = -20250.0_wp
    real(wp),parameter :: b32 =  15.0_wp
    real(wp),parameter :: b40 = -931041.0_wp
    real(wp),parameter :: b41 =  931500.0_wp
    real(wp),parameter :: b42 = -490.0_wp
    real(wp),parameter :: b43 =  112.0_wp

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+aa1*h*f0,f1)
    call me%f(t+a2*h,x+aa2*h*(b20*f0+b21*f1),f2)
    call me%f(t+a3*h,x+aa3*h*(b30*f0+b31*f1+b32*f2),f3)
    call me%f(t+h,   x+aa4*h*(b40*f0+b41*f1+b42*f2+b43*f3),f4)

    xf = x + h * c * (c0*f0 + c2*f2 + c3*f3 + c4*f4)

    end procedure rks5
!*****************************************************************************************

!*****************************************************************************************
!>
!  Butcher's 6th order method. 7 function evaluations.
!
!### References
!  * Butcher, J. (1964). On Runge-Kutta processes of high order.
!    Journal of the Australian Mathematical Society, 4(2), 179-194.

    module procedure rkb6

    real(wp),dimension(me%n) :: f1,f2,f3,f4,f5,f6,f7

    real(wp),parameter :: a2 = 1.0_wp / 3.0_wp
    real(wp),parameter :: a3 = 2.0_wp / 3.0_wp
    real(wp),parameter :: a4 = 1.0_wp / 3.0_wp
    real(wp),parameter :: a5 = 1.0_wp / 2.0_wp
    real(wp),parameter :: a6 = 1.0_wp / 2.0_wp

    real(wp),parameter :: b21 =  1.0_wp  / 3.0_wp
    real(wp),parameter :: b32 =  2.0_wp  / 3.0_wp
    real(wp),parameter :: b41 =  1.0_wp  / 12.0_wp
    real(wp),parameter :: b42 =  1.0_wp  / 3.0_wp
    real(wp),parameter :: b43 = -1.0_wp  / 12.0_wp
    real(wp),parameter :: b51 = -1.0_wp  / 16.0_wp
    real(wp),parameter :: b52 =  9.0_wp  / 8.0_wp
    real(wp),parameter :: b53 = -3.0_wp  / 16.0_wp
    real(wp),parameter :: b54 = -3.0_wp  / 8.0_wp
    real(wp),parameter :: b62 =  9.0_wp  / 8.0_wp
    real(wp),parameter :: b63 = -3.0_wp  / 8.0_wp
    real(wp),parameter :: b64 = -3.0_wp  / 4.0_wp
    real(wp),parameter :: b65 =  1.0_wp  / 2.0_wp
    real(wp),parameter :: b71 =  9.0_wp  / 44.0_wp
    real(wp),parameter :: b72 = -9.0_wp  / 11.0_wp
    real(wp),parameter :: b73 =  63.0_wp / 44.0_wp
    real(wp),parameter :: b74 =  18.0_wp / 11.0_wp
    real(wp),parameter :: b76 = -16.0_wp / 11.0_wp

    real(wp),parameter :: c1 =  11.0_wp / 120.0_wp
    real(wp),parameter :: c3 =  27.0_wp / 40.0_wp
    real(wp),parameter :: c4 =  27.0_wp / 40.0_wp
    real(wp),parameter :: c5 = -4.0_wp  / 15.0_wp
    real(wp),parameter :: c6 = -4.0_wp  / 15.0_wp
    real(wp),parameter :: c7 =  11.0_wp / 120.0_wp

    call me%f(t,       x,f1)
    call me%f(t+a2*h,  x+h*(b21*f1),f2)
    call me%f(t+a3*h,  x+h*(         b32*f2),f3)
    call me%f(t+a4*h,  x+h*(b41*f1 + b42*f2 + b43*f3),f4)
    call me%f(t+a5*h,  x+h*(b51*f1 + b52*f2 + b53*f3 + b54*f4),f5)
    call me%f(t+a6*h,  x+h*(         b62*f2 + b63*f3 + b64*f4 + b65*f5),f6)
    call me%f(t+h,     x+h*(b71*f1 + b72*f2 + b73*f3 + b74*f4 + b76*f6),f7)

    xf = x + h*(c1*f1 + c3*f3 + c4*f4 + c5*f5 + c6*f6 + c7*f7)

    end procedure rkb6
!*****************************************************************************************

!*****************************************************************************************
!>
!  Take one Runge Kutta 7 integration step: `t -> t+h (x -> xf)`
!
!### Reference
!  * E. B. Shanks, "Solutions of Differential Equations by Evaluations of Functions"
!    Math. Comp. 20 (1966).

    module procedure rk7

    real(wp),dimension(me%n) :: f0,f1,f2,f3,f4,f5,f6,f7,f8

    real(wp),parameter :: a1 = 2.0_wp / 9.0_wp
    real(wp),parameter :: a2 = 1.0_wp / 3.0_wp
    real(wp),parameter :: a3 = 1.0_wp / 2.0_wp
    real(wp),parameter :: a4 = 1.0_wp / 6.0_wp
    real(wp),parameter :: a5 = 8.0_wp / 9.0_wp
    real(wp),parameter :: a6 = 1.0_wp / 9.0_wp
    real(wp),parameter :: a7 = 5.0_wp / 6.0_wp
    !real(wp),parameter :: a8 = 1.0_wp

    real(wp),parameter :: c = 1.0_wp / 2140320.0_wp

    real(wp),parameter :: c0 = 110201.0_wp
    real(wp),parameter :: c3 = 767936.0_wp
    real(wp),parameter :: c4 = 635040.0_wp
    real(wp),parameter :: c5 = -59049.0_wp
    real(wp),parameter :: c6 = -59049.0_wp
    real(wp),parameter :: c7 = 635040.0_wp
    real(wp),parameter :: c8 = 110201.0_wp

    real(wp),parameter :: aa1 = 2.0_wp / 9.0_wp
    real(wp),parameter :: aa2 = 1.0_wp / 12.0_wp
    real(wp),parameter :: aa3 = 1.0_wp / 8.0_wp
    real(wp),parameter :: aa4 = 1.0_wp / 216.0_wp
    real(wp),parameter :: aa5 = 1.0_wp / 729.0_wp
    real(wp),parameter :: aa6 = 1.0_wp / 151632.0_wp
    real(wp),parameter :: aa7 = 1.0_wp / 1375920.0_wp
    real(wp),parameter :: aa8 = 1.0_wp / 251888.0_wp

    !real(wp),parameter :: b20 = 1.0_wp
    real(wp),parameter :: b21 = 3.0_wp

    !real(wp),parameter :: b30 = 1.0_wp
    real(wp),parameter :: b32 = 3.0_wp

    real(wp),parameter :: b40 = 23.0_wp
    real(wp),parameter :: b42 = 21.0_wp
    real(wp),parameter :: b43 = -8.0_wp

    real(wp),parameter :: b50 = -4136.0_wp
    real(wp),parameter :: b52 = -13584.0_wp
    real(wp),parameter :: b53 = 5264.0_wp
    real(wp),parameter :: b54 = 13104.0_wp

    real(wp),parameter :: b60 = 105131.0_wp
    real(wp),parameter :: b62 = 302016.0_wp
    real(wp),parameter :: b63 = -107744.0_wp
    real(wp),parameter :: b64 = -284256.0_wp
    real(wp),parameter :: b65 = 1701.0_wp

    real(wp),parameter :: b70 = -775229.0_wp
    real(wp),parameter :: b72 = -2770950.0_wp
    real(wp),parameter :: b73 = 1735136.0_wp
    real(wp),parameter :: b74 = 2547216.0_wp
    real(wp),parameter :: b75 = 81891.0_wp
    real(wp),parameter :: b76 = 328536.0_wp

    real(wp),parameter :: b80 = 23569.0_wp
    real(wp),parameter :: b82 = -122304.0_wp
    real(wp),parameter :: b83 = -20384.0_wp
    real(wp),parameter :: b84 = 695520.0_wp
    real(wp),parameter :: b85 = -99873.0_wp
    real(wp),parameter :: b86 = -466560.0_wp
    real(wp),parameter :: b87 = 241920.0_wp

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+aa1*h*(f0),f1)
    call me%f(t+a2*h,x+aa2*h*(    f0+b21*f1),f2)
    call me%f(t+a3*h,x+aa3*h*(    f0+b32*f2),f3)
    call me%f(t+a4*h,x+aa4*h*(b40*f0+b42*f2+b43*f3),f4)
    call me%f(t+a5*h,x+aa5*h*(b50*f0+b52*f2+b53*f3+b54*f4),f5)
    call me%f(t+a6*h,x+aa6*h*(b60*f0+b62*f2+b63*f3+b64*f4+b65*f5),f6)
    call me%f(t+a7*h,x+aa7*h*(b70*f0+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
    call me%f(t+h,   x+aa8*h*(b80*f0+b82*f2+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)

    xf = x + h * c * (c0*f0 + c3*f3 + c4*f4 + c5*f5 + c6*f6 + c7*f7 + c8*f8)

    end procedure rk7
!*****************************************************************************************

!*****************************************************************************************
!>
!  Take one Runge Kutta 8 integration step: `t -> t+h (x -> xf)`
!  This is Formula (8-10) from Reference [1].
!
!# Reference
!  1. E. B. Shanks, "[Higher Order Approximations of Runge-Kutta Type](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)",
!     NASA Technical Note, NASA TN D-2920, Sept. 1965.

    module procedure rk8_10

    real(wp),dimension(me%n) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9

    !parameters:
    real(wp),parameter :: a1  = 4.0_wp/27.0_wp
    real(wp),parameter :: a2  = 2.0_wp/9.0_wp
    real(wp),parameter :: a3  = 1.0_wp/3.0_wp
    real(wp),parameter :: a4  = 1.0_wp/2.0_wp
    real(wp),parameter :: a5  = 2.0_wp/3.0_wp
    real(wp),parameter :: a6  = 1.0_wp/6.0_wp
    real(wp),parameter :: a8  = 5.0_wp/6.0_wp
    real(wp),parameter :: c   = 1.0_wp/840.0_wp
    real(wp),parameter :: c0  = 41.0_wp
    real(wp),parameter :: c3  = 27.0_wp
    real(wp),parameter :: c4  = 272.0_wp
    real(wp),parameter :: c5  = 27.0_wp
    real(wp),parameter :: c6  = 216.0_wp
    real(wp),parameter :: c8  = 216.0_wp
    real(wp),parameter :: c9  = 41.0_wp
    real(wp),parameter :: aa1 = 4.0_wp/27.0_wp
    real(wp),parameter :: aa2 = 1.0_wp/18.0_wp
    real(wp),parameter :: aa3 = 1.0_wp/12.0_wp
    real(wp),parameter :: aa4 = 1.0_wp/8.0_wp
    real(wp),parameter :: aa5 = 1.0_wp/54.0_wp
    real(wp),parameter :: aa6 = 1.0_wp/4320.0_wp
    real(wp),parameter :: aa7 = 1.0_wp/20.0_wp
    real(wp),parameter :: aa8 = 1.0_wp/288.0_wp
    real(wp),parameter :: aa9 = 1.0_wp/820.0_wp
    real(wp),parameter :: b21 = 3.0_wp
    real(wp),parameter :: b32 = 3.0_wp
    real(wp),parameter :: b43 = 3.0_wp
    real(wp),parameter :: b50 = 13.0_wp
    real(wp),parameter :: b52 = -27.0_wp
    real(wp),parameter :: b53 = 42.0_wp
    real(wp),parameter :: b54 = 8.0_wp
    real(wp),parameter :: b60 = 389.0_wp
    real(wp),parameter :: b62 = -54.0_wp
    real(wp),parameter :: b63 = 966.0_wp
    real(wp),parameter :: b64 = -824.0_wp
    real(wp),parameter :: b65 = 243.0_wp
    real(wp),parameter :: b70 = -231.0_wp
    real(wp),parameter :: b72 = 81.0_wp
    real(wp),parameter :: b73 = -1164.0_wp
    real(wp),parameter :: b74 = 656.0_wp
    real(wp),parameter :: b75 = -122.0_wp
    real(wp),parameter :: b76 = 800.0_wp
    real(wp),parameter :: b80 = -127.0_wp
    real(wp),parameter :: b82 = 18.0_wp
    real(wp),parameter :: b83 = -678.0_wp
    real(wp),parameter :: b84 = 456.0_wp
    real(wp),parameter :: b85 = -9.0_wp
    real(wp),parameter :: b86 = 576.0_wp
    real(wp),parameter :: b87 = 4.0_wp
    real(wp),parameter :: b90 = 1481.0_wp
    real(wp),parameter :: b92 = -81.0_wp
    real(wp),parameter :: b93 = 7104.0_wp
    real(wp),parameter :: b94 = -3376.0_wp
    real(wp),parameter :: b95 = 72.0_wp
    real(wp),parameter :: b96 = -5040.0_wp
    real(wp),parameter :: b97 = -60.0_wp
    real(wp),parameter :: b98 = 720.0_wp

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+aa1*h*f0,f1)
    call me%f(t+a2*h,x+aa2*h*(f0+b21*f1),f2)
    call me%f(t+a3*h,x+aa3*h*(f0+b32*f2),f3)
    call me%f(t+a4*h,x+aa4*h*(f0+b43*f3),f4)
    call me%f(t+a5*h,x+aa5*h*(b50*f0+b52*f2+b53*f3+b54*f4),f5)
    call me%f(t+a6*h,x+aa6*h*(b60*f0+b62*f2+b63*f3+b64*f4+b65*f5),f6)
    call me%f(t+h,x+aa7*h*(b70*f0+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
    call me%f(t+a8*h,x+aa8*h*(b80*f0+b82*f2+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
    call me%f(t+h,x+aa9*h*(b90*f0+b92*f2+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)

    xf = x + h*c*(c0*f0+c3*f3+c4*f4+c5*f5+c6*f6+c8*f8+c9*f9)

    end procedure rk8_10
!*****************************************************************************************

!*****************************************************************************************
!>
!  8th order Shanks, 12 function evaluations.
!
!### Reference
!  * E. B. Shanks, "Solutions of Differential Equations by Evaluations of Functions"
!    Math. Comp. 20 (1966).

    module procedure rk8_12

    real(wp),dimension(me%n) :: f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11

    real(wp),parameter :: a1    =  1.0_wp / 9.0_wp
    real(wp),parameter :: a2    =  1.0_wp / 6.0_wp
    real(wp),parameter :: a3    =  1.0_wp / 4.0_wp
    real(wp),parameter :: a4    =  1.0_wp / 10.0_wp
    real(wp),parameter :: a5    =  1.0_wp / 6.0_wp
    real(wp),parameter :: a6    =  1.0_wp / 2.0_wp
    real(wp),parameter :: a7    =  2.0_wp / 3.0_wp
    real(wp),parameter :: a8    =  1.0_wp / 3.0_wp
    real(wp),parameter :: a9    =  5.0_wp / 6.0_wp
    real(wp),parameter :: a10   =  5.0_wp / 6.0_wp
    ! real(wp),parameter :: a11   =  1.0_wp
    real(wp),parameter :: c     =  1.0_wp / 840.0_wp
    real(wp),parameter :: c0    =  41.0_wp
    real(wp),parameter :: c5    =  216.0_wp
    real(wp),parameter :: c6    =  272.0_wp
    real(wp),parameter :: c7    =  27.0_wp
    real(wp),parameter :: c8    =  27.0_wp
    real(wp),parameter :: c9    =  36.0_wp
    real(wp),parameter :: c10   =  180.0_wp
    real(wp),parameter :: c11   =  41.0_wp
    real(wp),parameter :: aa1   =  1.0_wp / 9.0_wp
    real(wp),parameter :: aa2   =  1.0_wp / 24.0_wp
    real(wp),parameter :: aa3   =  1.0_wp / 16.0_wp
    real(wp),parameter :: aa4   =  1.0_wp / 500.0_wp
    real(wp),parameter :: aa5   =  1.0_wp / 972.0_wp
    real(wp),parameter :: aa6   =  1.0_wp / 36.0_wp
    real(wp),parameter :: aa7   =  1.0_wp / 243.0_wp
    real(wp),parameter :: aa8   =  1.0_wp / 324.0_wp
    real(wp),parameter :: aa9   =  1.0_wp / 324.0_wp
    real(wp),parameter :: aa10  =  1.0_wp / 1620.0_wp
    real(wp),parameter :: aa11  =  1.0_wp / 4428.0_wp
    !real(wp),parameter :: b20   =  1.0_wp
    real(wp),parameter :: b21   =  3.0_wp
    !real(wp),parameter :: b30   =  1.0_wp
    real(wp),parameter :: b32   =  3.0_wp
    real(wp),parameter :: b40   =  29.0_wp
    real(wp),parameter :: b42   =  33.0_wp
    real(wp),parameter :: b43   = -12.0_wp
    real(wp),parameter :: b50   =  33.0_wp
    real(wp),parameter :: b53   =  4.0_wp
    real(wp),parameter :: b54   =  125.0_wp
    real(wp),parameter :: b60   = -21.0_wp
    real(wp),parameter :: b63   =  76.0_wp
    real(wp),parameter :: b64   =  125.0_wp
    real(wp),parameter :: b65   = -162.0_wp
    real(wp),parameter :: b70   = -30.0_wp
    real(wp),parameter :: b73   = -32.0_wp
    real(wp),parameter :: b74   =  125.0_wp
    real(wp),parameter :: b76   =  99.0_wp
    real(wp),parameter :: b80   =  1175.0_wp
    real(wp),parameter :: b83   = -3456.0_wp
    real(wp),parameter :: b84   = -6250.0_wp
    real(wp),parameter :: b85   =  8424.0_wp
    real(wp),parameter :: b86   =  242.0_wp
    real(wp),parameter :: b87   = -27.0_wp
    real(wp),parameter :: b90   =  293.0_wp
    real(wp),parameter :: b93   = -852.0_wp
    real(wp),parameter :: b94   = -1375.0_wp
    real(wp),parameter :: b95   =  1836.0_wp
    real(wp),parameter :: b96   = -118.0_wp
    real(wp),parameter :: b97   =  162.0_wp
    real(wp),parameter :: b98   =  324.0_wp
    real(wp),parameter :: b100  =  1303.0_wp
    real(wp),parameter :: b103  = -4260.0_wp
    real(wp),parameter :: b104  = -6875.0_wp
    real(wp),parameter :: b105  =  9990.0_wp
    real(wp),parameter :: b106  =  1030.0_wp
    real(wp),parameter :: b109  =  162.0_wp
    real(wp),parameter :: b110  = -8595.0_wp
    real(wp),parameter :: b113  =  30720.0_wp
    real(wp),parameter :: b114  =  48750.0_wp
    real(wp),parameter :: b115  = -66096.0_wp
    real(wp),parameter :: b116  =  378.0_wp
    real(wp),parameter :: b117  = -729.0_wp
    real(wp),parameter :: b118  = -1944.0_wp
    real(wp),parameter :: b119  = -1296.0_wp
    real(wp),parameter :: b1110 =  3240.0_wp

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+aa1*h*(f0),f1)
    call me%f(t+a2*h,x+aa2*h*(    f0+b21*f1),f2)
    call me%f(t+a3*h,x+aa3*h*(    f0+b32*f2),f3)
    call me%f(t+a4*h,x+aa4*h*(b40*f0+b42*f2+b43*f3),f4)
    call me%f(t+a5*h,x+aa5*h*(b50*f0+b53*f3+b54*f4),f5)
    call me%f(t+a6*h,x+aa6*h*(b60*f0+b63*f3+b64*f4+b65*f5),f6)
    call me%f(t+a7*h,x+aa7*h*(b70*f0+b73*f3+b74*f4+b76*f6),f7)
    call me%f(t+a8*h,x+aa8*h*(b80*f0+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
    call me%f(t+a9*h,x+aa9*h*(b90*f0+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
    call me%f(t+a10*h,x+aa10*h*(b100*f0+b103*f3+b104*f4+b105*f5+b106*f6+b109*f9),f10)
    call me%f(t+h,    x+aa11*h*(b110*f0+b113*f3+b114*f4+b115*f5+b116*f6+b117*f7+&
                                b118*f8+b119*f9+b1110*f10),f11)

    xf = x + h*c*(c0*f0+c5*f5+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+c11*f11)

    end procedure rk8_12
!*****************************************************************************************

!*****************************************************************************************
!>
!  Cooper-Verner 11 stage, 8th order Runge-Kutta method.
!
!### Reference
!  * Some Explicit Runge-Kutta Methods of High Order, by G. J. Cooper and J. H. Verner,
!    SIAM Journal on Numerical Analysis, Vol. 9, No. 3, (September 1972), pages 389 to 405
!  * http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK8/RKcoeff8b_1.pdf

    module procedure rkcv8

    real(wp),dimension(me%n) :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11

    real(wp),parameter :: s = sqrt(21.0_wp)

    real(wp),parameter :: a2 = 1.0_wp / 2.0_wp
    real(wp),parameter :: a3 = 1.0_wp / 2.0_wp
    real(wp),parameter :: a4 = 1.0_wp / 2.0_wp - 1.0_wp / 14.0_wp*s
    real(wp),parameter :: a5 = 1.0_wp / 2.0_wp - 1.0_wp / 14.0_wp*s
    real(wp),parameter :: a6 = 1.0_wp / 2.0_wp
    real(wp),parameter :: a7 = 1.0_wp / 2.0_wp + 1.0_wp / 14.0_wp*s
    real(wp),parameter :: a8 = 1.0_wp / 2.0_wp + 1.0_wp / 14.0_wp*s
    real(wp),parameter :: a9 = 1.0_wp / 2.0_wp
    real(wp),parameter :: a10 = 1.0_wp / 2.0_wp - 1.0_wp / 14.0_wp*s

    real(wp),parameter :: b21 = 1.0_wp / 2.0_wp
    real(wp),parameter :: b31 = 1.0_wp / 4.0_wp
    real(wp),parameter :: b32 = 1.0_wp / 4.0_wp
    real(wp),parameter :: b41 = 1.0_wp / 7.0_wp
    real(wp),parameter :: b42 = -1.0_wp / 14.0_wp + 3.0_wp / 98.0_wp*s
    real(wp),parameter :: b43 = 3.0_wp / 7.0_wp - 5.0_wp / 49.0_wp*s
    real(wp),parameter :: b51 = 11.0_wp / 84.0_wp - 1.0_wp / 84.0_wp*s
    real(wp),parameter :: b53 = 2.0_wp / 7.0_wp - 4.0_wp / 63.0_wp*s
    real(wp),parameter :: b54 = 1.0_wp / 12.0_wp + 1.0_wp / 252.0_wp*s
    real(wp),parameter :: b61 = 5.0_wp / 48.0_wp - 1.0_wp / 48.0_wp*s
    real(wp),parameter :: b63 = 1.0_wp / 4.0_wp - 1.0_wp / 36.0_wp*s
    real(wp),parameter :: b64 = -77.0_wp / 120.0_wp-7.0_wp / 180.0_wp*s
    real(wp),parameter :: b65 = 63.0_wp / 80.0_wp + 7.0_wp / 80.0_wp*s
    real(wp),parameter :: b71 = 5.0_wp / 21.0_wp + 1.0_wp / 42.0_wp*s
    real(wp),parameter :: b73 = -48.0_wp / 35.0_wp - 92.0_wp / 315.0_wp*s
    real(wp),parameter :: b74 = 211.0_wp / 30.0_wp + 29.0_wp / 18.0_wp*s
    real(wp),parameter :: b75 = -36.0_wp / 5.0_wp - 23.0_wp / 14.0_wp*s
    real(wp),parameter :: b76 = 9.0_wp / 5.0_wp + 13.0_wp / 35.0_wp*s
    real(wp),parameter :: b81 = 1.0_wp / 14.0_wp
    real(wp),parameter :: b85 = 1.0_wp / 9.0_wp + 1.0_wp / 42.0_wp*s
    real(wp),parameter :: b86 = 13.0_wp / 63.0_wp + 1.0_wp / 21.0_wp*s
    real(wp),parameter :: b87 = 1.0_wp / 9.0_wp
    real(wp),parameter :: b91 = 1.0_wp / 32.0_wp
    real(wp),parameter :: b95 = 91.0_wp / 576.0_wp + 7.0_wp / 192.0_wp*s
    real(wp),parameter :: b96 = 11.0_wp / 72.0_wp
    real(wp),parameter :: b97 = -385.0_wp / 1152.0_wp + 25.0_wp / 384.0_wp*s
    real(wp),parameter :: b98 = 63.0_wp / 128.0_wp - 13.0_wp / 128.0_wp*s
    real(wp),parameter :: b101 = 1.0_wp / 14.0_wp
    real(wp),parameter :: b105 = 1.0_wp / 9.0_wp
    real(wp),parameter :: b106 = -733.0_wp / 2205.0_wp + 1.0_wp / 15.0_wp*s
    real(wp),parameter :: b107 = 515.0_wp / 504.0_wp - 37.0_wp / 168.0_wp*s
    real(wp),parameter :: b108 = -51.0_wp / 56.0_wp + 11.0_wp / 56.0_wp*s
    real(wp),parameter :: b109 = 132.0_wp / 245.0_wp - 4.0_wp / 35.0_wp*s
    real(wp),parameter :: b115 = -7.0_wp / 3.0_wp - 7.0_wp / 18.0_wp*s
    real(wp),parameter :: b116 = -2.0_wp / 5.0_wp - 28.0_wp / 45.0_wp*s
    real(wp),parameter :: b117 = -91.0_wp / 24.0_wp + 53.0_wp / 72.0_wp*s
    real(wp),parameter :: b118 = 301.0_wp / 72.0_wp - 53.0_wp / 72.0_wp*s
    real(wp),parameter :: b119 = 28.0_wp / 45.0_wp + 28.0_wp / 45.0_wp*s
    real(wp),parameter :: b1110 = 49.0_wp / 18 + 7.0_wp / 18.0_wp*s

    real(wp),parameter :: c1 = 1.0_wp / 20.0_wp
    real(wp),parameter :: c8 = 49.0_wp / 180.0_wp
    real(wp),parameter :: c9 = 16.0_wp / 45.0_wp
    real(wp),parameter :: c10 = 49.0_wp / 180.0_wp
    real(wp),parameter :: c11 = 1.0_wp / 20.0_wp

    call me%f(t,      x,f1)
    call me%f(t+a2*h, x+h*(b21*f1),f2)
    call me%f(t+a3*h, x+h*(b31*f1  + b32*f2),f3)
    call me%f(t+a4*h, x+h*(b41*f1  + b42*f2  + b43*f3),f4)
    call me%f(t+a5*h, x+h*(b51*f1            + b53*f3  + b54*f4),f5)
    call me%f(t+a6*h, x+h*(b61*f1            + b63*f3  + b64*f4  + b65*f5),f6)
    call me%f(t+a7*h, x+h*(b71*f1            + b73*f3  + b74*f4  + b75*f5 + b76*f6),f7)
    call me%f(t+a8*h, x+h*(b81*f1                                + b85*f5  + b86*f6  + b87*f7),f8)
    call me%f(t+a9*h, x+h*(b91*f1                                + b95*f5  + b96*f6  + b97*f7 + b98*f8),f9)
    call me%f(t+a10*h,x+h*(b101*f1                               + b105*f5 + b106*f6 + b107*f7 + b108*f8 + b109*f9),f10)
    call me%f(t+h,    x+h*(                                        b115*f5 + b116*f6 + b117*f7 + b118*f8 + b119*f9 + b1110*f10),f11)

    xf = x + h*(c1*f1+c8*f8+c9*f9+c10*f10+c11*f11)

    end procedure rkcv8
!*****************************************************************************************

!*****************************************************************************************
!>
!  Zhang 10th order
!
!### Reference
!  * David K. Zhang, "Discovering New Runge-Kutta Methods Using Unstructured Numerical Search",
!    Thesis, April 16, 2019 [1911.00318.pdf](https://arxiv.org/pdf/1911.00318.pdf)
!  * [Coefficients](https://github.com/dzhang314/rktk/blob/master/methods/Zhang10.txt)

    module procedure rkz10

    real(wp),parameter :: b21   = +0.06888096612188652230677098661632935381315159322698285980682460033484176156636_wp
    real(wp),parameter :: b31   = -0.83810520353364237535186366202804838694106968892100493081079119861616316483888_wp
    real(wp),parameter :: b32   = +1.25369004871695465344923436259290210937215845259642921936068038992469719487344_wp
    real(wp),parameter :: b41   = -0.00490948675932680175369634467679043271190185159760893197940412354439119879747_wp
    real(wp),parameter :: b42   = +0.08232173030768021571147603044143500716498050758632346872918243517328221749373_wp
    real(wp),parameter :: b43   = -0.00853127742646689165100869914831522063992706276173167694295371129404925712991_wp
    real(wp),parameter :: b51   = +1.04893195376435953811751460952804681519865641127045728585681123905425821150864_wp
    real(wp),parameter :: b52   = -0.75382817731759581241904280322911590179250916990960557195131331600659134298994_wp
    real(wp),parameter :: b53   = +0.80522815974064911780476652756884437914316713728791543021500322792680878839117_wp
    real(wp),parameter :: b54   = -0.38446322074238561063992159166768203951623890872336238198515458977800812023865_wp
    real(wp),parameter :: b61   = -0.23992383433329995270018501071487605020336614587127799309879281304981537313900_wp
    real(wp),parameter :: b62   = -0.06364261163929229571107164643776483490956349179570579724023318438638668084701_wp
    real(wp),parameter :: b63   = +0.19967543135197895864283978317448663244372194705250026204952108330132418721869_wp
    real(wp),parameter :: b64   = +0.61035379509547002158981202467825146906358429601618359949777530100370262806814_wp
    real(wp),parameter :: b65   = +0.37859947224999048426403753064329493462141566215033797838496558400426835557965_wp
    real(wp),parameter :: b71   = +0.01778833946316172550067667113299992644020408256145083580468064240809624146756_wp
    real(wp),parameter :: b72   = -0.01105021905501825472152093430479940021746699824126398377648788328289961044531_wp
    real(wp),parameter :: b73   = -0.00439342855052892929633342674635339618911779936473248023015437836908839819799_wp
    real(wp),parameter :: b74   = +0.10597527290509018731505658395393495306363828909501538659771103682791206047696_wp
    real(wp),parameter :: b75   = +0.00405089069638330732472929903714130352675532845602736909932582009874330236167_wp
    real(wp),parameter :: b76   = -0.00167929709134223461022001825411827908084682305838572957913124250043819202374_wp
    real(wp),parameter :: b81   = +0.23566046225418537925692690356824973135127585829120749799176100081366132717900_wp
    real(wp),parameter :: b82   = +0.08893362689701559656112930868794335701480724068673657440558996917556014206620_wp
    real(wp),parameter :: b83   = +0.04138834709876858545169518341465722051417932518424807026039922164100253591914_wp
    real(wp),parameter :: b84   = -0.85290303603263069912309614093599359240031750375681439893926580845622126359984_wp
    real(wp),parameter :: b85   = -0.02308087548113625563008379864870135277272977633339373075859414529485111031087_wp
    real(wp),parameter :: b86   = +0.00968592188591852514590807745954989979539211465004178878459403530669594442891_wp
    real(wp),parameter :: b87   = +0.80144977934196893864833794701047679699141990680205169668196460489197395944287_wp
    real(wp),parameter :: b91   = +0.09048693760705330681795027653195090842362561521235737211709949607524998063814_wp
    real(wp),parameter :: b92   = +0.03204427202479226242378001797390171191984804295279358947221505285617826137779_wp
    real(wp),parameter :: b93   = +0.12433768215891056280095053314827014992247446913700152962717366852615035330919_wp
    real(wp),parameter :: b94   = -0.30731521755038158104946232011583084977900951861660916480042444618524837737411_wp
    real(wp),parameter :: b95   = +0.16447843681160268810163161551627427448581444120771873237010859382161886567365_wp
    real(wp),parameter :: b96   = -0.04100681167344768439646557295106150922485658224571650562701715147123044716022_wp
    real(wp),parameter :: b97   = +0.40661989930773198456977153545461924102851176247811339291440854490078881632775_wp
    real(wp),parameter :: b98   = +0.18669987124631939180018777607871869034096434382893787663176153320303681978038_wp
    real(wp),parameter :: b101  = -0.12815849137721058028379122684114729957595452282110310092889040855742427673083_wp
    real(wp),parameter :: b102  = -0.09494292242532245543748881606192770346771760827967185000597507680447193729080_wp
    real(wp),parameter :: b103  = -0.14450344511826259894917045729900439186304081642314445476990786282061867102157_wp
    real(wp),parameter :: b104  = +0.92134940704406912363606156819205798706915260363096464763106069459083167107327_wp
    real(wp),parameter :: b105  = -0.13265301054949046320484467982507968172620145761107284691851092173475314294014_wp
    real(wp),parameter :: b106  = +0.01661485872631464348319012586897860991564771318706900087871886151957318825574_wp
    real(wp),parameter :: b107  = -0.64461773243991825568826538659502156175792556790721387489238296660149971176174_wp
    real(wp),parameter :: b108  = +0.47148463683422024098060807066053294761028203774917809473770514989820710158152_wp
    real(wp),parameter :: b109  = +0.15101154448891262356107150246546481622684638215041867281807172181868980886911_wp
    real(wp),parameter :: b111  = -0.35394262889926481874318272132638227200275604000107564690424492655250744148526_wp
    real(wp),parameter :: b112  = -0.12999250061102426001102206164524015480707041252098090235235458563228391897841_wp
    real(wp),parameter :: b113  = +0.07296717474130959358123820183761890201581557542115580687780801792533382974931_wp
    real(wp),parameter :: b114  = +1.23409288203437673148811178545199717478084359792767549031487434785241641765067_wp
    real(wp),parameter :: b115  = +0.24325479314414790350916184694938182359002523242070889454561579643424180130675_wp
    real(wp),parameter :: b116  = -0.04886415340347670065787549059652108014832442895629947389530738182689787996728_wp
    real(wp),parameter :: b117  = -0.51442371407978914516038594645191171303742863244306956194414424856852330710239_wp
    real(wp),parameter :: b118  = +0.07088500449199304304537010242765865457286791828303557337358385517356109210099_wp
    real(wp),parameter :: b119  = -0.20555333994396552677698205301698855743095009068049025897065388713681681672772_wp
    real(wp),parameter :: b1110 = +0.04716132770900545782293703693524094489806604422476436750471220364001025348790_wp
    real(wp),parameter :: b121  = -1.28883596401187309154174375754588391722337985378291577953064993833085325077145_wp
    real(wp),parameter :: b122  = -0.28020869761589692854933312826845810393637896090137224607343344880154309805585_wp
    real(wp),parameter :: b123  = +2.95502531060581466025398140893693575414694538631086569687435523799636000133369_wp
    real(wp),parameter :: b124  = +2.68729452804277601048934070442766694554427130448420986662047690713132366345266_wp
    real(wp),parameter :: b125  = +4.94062342991935573838476417426809163760371481941041460313006887242753416771859_wp
    real(wp),parameter :: b126  = -1.10669748003212502578007105340045335120159611347391447923604138591653946190627_wp
    real(wp),parameter :: b127  = -1.01499184284449721513261150169983810954496707583685170059698090784083909089787_wp
    real(wp),parameter :: b128  = +2.53750982609373480524394613209930567128221625331336116097749979112883387569908_wp
    real(wp),parameter :: b129  = -3.05085435525303312587778339284072228780737127406865428592625671405448243077510_wp
    real(wp),parameter :: b1210 = -4.48248369012964780236263511347668602256466957962418723733767273268597925464354_wp
    real(wp),parameter :: b1211 = -1.21257843188510861623323738621796488386816954981559478788233918659317859600661_wp
    real(wp),parameter :: b131  = -0.26440344286774078027313975052115341819194404740705505988980388884472351746613_wp
    real(wp),parameter :: b132  = -0.06681118811399606825394685412417487195704333031683770523718698572711060614338_wp
    real(wp),parameter :: b133  = +0.22983178806827308655694675660046084355033982772309698316158066683955449428729_wp
    real(wp),parameter :: b134  = +0.64074149645736181719485655167408177784648088775776468738158085006283347730738_wp
    real(wp),parameter :: b135  = +0.42942911081576828545053652299580217770451868606900245607955202491949872419495_wp
    real(wp),parameter :: b136  = -0.00967078876662303036856558223558993132212111489612543886857914521712390529304_wp
    real(wp),parameter :: b137  = +0.01613666672999681392501812041282804272240296713357286534780333216134004281864_wp
    real(wp),parameter :: b138  = -0.03596314318487727365094853269870982256844289210631058204493077377818835633949_wp
    real(wp),parameter :: b139  = -0.04902670091896965612083428201101205091198148478269054303875029017005640576667_wp
    real(wp),parameter :: b1310 = +0.01388537673498896790224810676002758506392838309962714417674884945696557987538_wp
    real(wp),parameter :: b1311 = -0.02329475109524907760636494276820464664967810607545806961371703044079964798003_wp
    real(wp),parameter :: b1312 = +0.00411071472206196739547606137318514638902786727009559929394539277815884070579_wp
    real(wp),parameter :: b141  = +1.28995927881082034516023610876593719136951515200841564716720963602983766539309_wp
    real(wp),parameter :: b142  = +0.13503265202676397192358271319486347590315093450300003507414087911915346743493_wp
    real(wp),parameter :: b143  = -1.57533735338626298015011591165511659112640458688783385930436332029909777715227_wp
    real(wp),parameter :: b144  = -1.13614144125767380945133294734170682018114470434355261134291967031467884444496_wp
    real(wp),parameter :: b145  = -2.60629976327340284860734294234508518773200013764798796731373995419345584010989_wp
    real(wp),parameter :: b146  = -0.14225514857221282929776869245957057127882456911938845124276329363250774479820_wp
    real(wp),parameter :: b147  = -1.30050279397530175137190231972314623324238865480815606289913301764254981931424_wp
    real(wp),parameter :: b148  = +2.70797163114645615543070920761013389431867941711011998785845670786144007733081_wp
    real(wp),parameter :: b149  = +2.43554296946750176014394493560037156798756615586273835416220476865949138763199_wp
    real(wp),parameter :: b1410 = -0.26458182516446999987007272710827096690076304619434648405874330108389046658365_wp
    real(wp),parameter :: b1411 = +0.28431747275350668428305940521170347967298804085514279367418757904848425107603_wp
    real(wp),parameter :: b1412 = -0.01567336501605444936187349334932741547111318613954509638698908027760861917463_wp
    real(wp),parameter :: b1413 = +0.60355253162364202926624736416406789911182794847681800316234125803391629274555_wp
    real(wp),parameter :: b151  = -0.82170145239485721248182614307298014363501501247665665859542963907661741852252_wp
    real(wp),parameter :: b152  = +1.25591221619821840505806060191670198289896021578755350597122093077684012752648_wp
    real(wp),parameter :: b153  = -0.01615284116525821921698589364201132298146948312403205059324494512949466638820_wp
    real(wp),parameter :: b154  = -0.02297038987201048533939731478067813642773684981394747854111883904855646952919_wp
    real(wp),parameter :: b155  = -0.02796034376707111482574969609615809805957457730932703561247060660045527963337_wp
    real(wp),parameter :: b156  = -0.75094298727080692288567852962945843784728392190426468013732662611274129611030_wp
    real(wp),parameter :: b157  = -0.00831029089278008603544533798969572897337060350821384521268253868513966100453_wp
    real(wp),parameter :: b158  = +0.02785523862077057819395231768788127899618354611871665730619924470585495216854_wp
    real(wp),parameter :: b159  = +0.04021421539333651503730432703988977024358721548426886753518601253708363296062_wp
    real(wp),parameter :: b1510 = +1.61713646360095425893804116572314412328254887955289755326352595104935489454125_wp
    real(wp),parameter :: b1511 = -1.39628377741243214317623603663032161218084680028714877542197062209713837494256_wp
    real(wp),parameter :: b1512 = -0.01424439159061064208304848726946155039179802539115130281665541763757472383980_wp
    real(wp),parameter :: b1513 = +0.75709054336643265169970374917936196757291644158102449813751290304179069353124_wp
    real(wp),parameter :: b1514 = -0.22405735763057330478532402187136037006601226103429496673285661641467238072308_wp
    real(wp),parameter :: b161  = +0.27726401853185531071095250844208662954769209843860654782270439586107555912462_wp
    real(wp),parameter :: b162  = +0.09453818370728229619312145290868518019954577268055790438905422269619151226603_wp
    real(wp),parameter :: b163  = +0.84172754295454286541980788975065897179825972249436132684274103676922879908590_wp
    real(wp),parameter :: b164  = -0.90665259832844475943298457713954683980649043330939335614048943465366310987303_wp
    real(wp),parameter :: b165  = -0.09334858033923263720339222938292231543826754868386677793091178055855225779084_wp
    real(wp),parameter :: b166  = +4.08887758914114010965574538791054120551785765841976165061624483740399180627400_wp
    real(wp),parameter :: b167  = +0.79539986499428990687079492545982314739495812782573436535546220822826278953650_wp
    real(wp),parameter :: b168  = -0.04859175145875636251218702854395529418032114353274620641007859711365783956136_wp
    real(wp),parameter :: b169  = +0.14819851457498430887143064575970769402240035004018486609649448109702272035773_wp
    real(wp),parameter :: b1610 = -1.77021148062189597056755123808303838576982978034562033343027660487016998355672_wp
    real(wp),parameter :: b1611 = +1.83821084492008174732030666893446362369060575712213328902105414298283276556307_wp
    real(wp),parameter :: b1612 = +0.04909344463430916722179909282844531508454281205967923461622573420987424382584_wp
    real(wp),parameter :: b1613 = -3.80378823067007538071243467410680278604792596795503781745130930797860484235966_wp
    real(wp),parameter :: b1614 = +0.33157239872339854149205035234217673029656296754907585248267975563293121409049_wp
    real(wp),parameter :: b1615 = -0.84228976076347914332745917708032287630959039280343054587959508970676337698255_wp

    real(wp),parameter :: c1 = +0.03181927458023409759419419944088926839595967458804889313841865851348230504492_wp
    real(wp),parameter :: c3 = +0.04681369289018421954398607025172424191865836856120072277976273594277845007872_wp
    real(wp),parameter :: c6 = +1.37553536170545749013126339747598151094852489823041688723220548486290269295630_wp
    real(wp),parameter :: c7 = +0.17506567143964248943590740397548086595186840317795589108615338820052040464001_wp
    real(wp),parameter :: c8 = +0.14924798653008455234489054168544610975862154045444166485218499565185924479229_wp
    real(wp),parameter :: c9 = +0.27180992312662372142355988162582830282467844949919584507610421296773792690656_wp
    real(wp),parameter :: c10 = -0.17077940511361075387191294440611902706762757201779959741914461935949309724569_wp
    real(wp),parameter :: c11 = +0.30035519768848521819453255694344499143765502380409939822945376435190034136303_wp
    real(wp),parameter :: c12 = -0.01014254403202798636405231311418830930667910487160410011307114999518762983823_wp
    real(wp),parameter :: c13 = -1.19177815782093190956613115811061881440838744821819017576359936663422556606341_wp
    real(wp),parameter :: c14 = +0.03639139354170713653095246097302347593652682408895579802073390992985525665017_wp
    real(wp),parameter :: c15 = -0.04683316086510645556720041443823709643011325850704077997225219431585033242635_wp
    real(wp),parameter :: c16 = +0.03249476632925818017001031769734448004031420121031955285305017988372000314167_wp

    ! see Equation 1.20 in reference:
    real(wp),parameter :: a1  = 0
    real(wp),parameter :: a2  = b21
    real(wp),parameter :: a3  = b31 +b32
    real(wp),parameter :: a4  = b41 +b42 +b43
    real(wp),parameter :: a5  = b51 +b52 +b53 +b54
    real(wp),parameter :: a6  = b61 +b62 +b63 +b64 +b65
    real(wp),parameter :: a7  = b71 +b72 +b73 +b74 +b75 +b76
    real(wp),parameter :: a8  = b81 +b82 +b83 +b84 +b85 +b86 +b87
    real(wp),parameter :: a9  = b91 +b92 +b93 +b94 +b95 +b96 +b97 +b98
    real(wp),parameter :: a10 = b101+b102+b103+b104+b105+b106+b107+b108+b109
    real(wp),parameter :: a11 = b111+b112+b113+b114+b115+b116+b117+b118+b119+b1110
    real(wp),parameter :: a12 = b121+b122+b123+b124+b125+b126+b127+b128+b129+b1210+b1211
    real(wp),parameter :: a13 = b131+b132+b133+b134+b135+b136+b137+b138+b139+b1310+b1311+b1312
    real(wp),parameter :: a14 = b141+b142+b143+b144+b145+b146+b147+b148+b149+b1410+b1411+b1412+b1413
    real(wp),parameter :: a15 = b151+b152+b153+b154+b155+b156+b157+b158+b159+b1510+b1511+b1512+b1513+b1514
    real(wp),parameter :: a16 = b161+b162+b163+b164+b165+b166+b167+b168+b169+b1610+b1611+b1612+b1613+b1614+b1615

    real(wp),dimension(me%n) :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16

    call me%f(t+a1*h,  x,f1)
    call me%f(t+a2*h,  x+h*(b21*f1),f2)
    call me%f(t+a3*h,  x+h*(b31*f1  + b32*f2),f3)
    call me%f(t+a4*h,  x+h*(b41*f1  + b42*f2  + b43*f3),f4)
    call me%f(t+a5*h,  x+h*(b51*f1  + b52*f2  + b53*f3  + b54*f4),f5)
    call me%f(t+a6*h,  x+h*(b61*f1  + b62*f2  + b63*f3  + b64*f4  + b65*f5),f6)
    call me%f(t+a7*h,  x+h*(b71*f1  + b72*f2  + b73*f3  + b74*f4  + b75*f5 + b76*f6),f7)
    call me%f(t+a8*h,  x+h*(b81*f1  + b82*f2  + b83*f3  + b84*f4  + b85*f5  + b86*f6  + &
                            b87*f7),f8)
    call me%f(t+a9*h,  x+h*(b91*f1  + b92*f2  + b93*f3  + b94*f4  + b95*f5  + b96*f6  + &
                            b97*f7 + b98*f8),f9)
    call me%f(t+a10*h, x+h*(b101*f1 + b102*f2 + b103*f3 + b104*f4 + b105*f5 + b106*f6 + &
                            b107*f7 + b108*f8 + b109*f9),f10)
    call me%f(t+a11*h, x+h*(b111*f1 + b112*f2 + b113*f3 + b114*f4 + b115*f5 + b116*f6 + &
                            b117*f7 + b118*f8 + b119*f9 + b1110*f10),f11)
    call me%f(t+a12*h, x+h*(b121*f1 + b122*f2 + b123*f3 + b124*f4 + b125*f5 + b126*f6 + &
                            b127*f7 + b128*f8 + b129*f9 + b1210*f10 + b1211*f11),f12)
    call me%f(t+a13*h, x+h*(b131*f1 + b132*f2 + b133*f3 + b134*f4 + b135*f5 + b136*f6 + &
                            b137*f7 + b138*f8 + b139*f9 + b1310*f10 + b1311*f11 + &
                            b1312*f12),f13)
    call me%f(t+a14*h, x+h*(b141*f1 + b142*f2 + b143*f3 + b144*f4 + b145*f5 + b146*f6 + &
                            b147*f7 + b148*f8 + b149*f9 + b1410*f10 + b1411*f11 + &
                            b1412*f12 + b1413*f13),f14)
    call me%f(t+a15*h, x+h*(b151*f1 + b152*f2 + b153*f3 + b154*f4 + b155*f5 + b156*f6 + &
                            b157*f7 + b158*f8 + b159*f9 + b1510*f10 + b1511*f11 + &
                            b1512*f12 + b1513*f13 + b1514*f14),f15)
    call me%f(t+a16*h, x+h*(b161*f1 + b162*f2 + b163*f3 + b164*f4 + b165*f5 + b166*f6 + &
                            b167*f7 + b168*f8 + b169*f9 + b1610*f10 + b1611*f11 + &
                            b1612*f12 + b1613*f13 + b1614*f14 + b1615*f15),f16)

    xf = x+h*(c1*f1+c3*f3+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+&
              c11*f11+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16)

    end procedure rkz10
!*****************************************************************************************

!*****************************************************************************************
!>
!  Ono's 10th order method
!
!### References
!  * Hiroshi Ono, "A Runge-Kutta method of order 10 which minimizes truncation error",
!    The Japan Society for Industrial and Applied Mathematics,
!    Vol. 13, No. 1, 2003, pp 35 - 44.

    module procedure rko10

    real(wp),parameter :: a2 = .3357505083417036184129939488963472892525539577157696170900805997207182929518116563365_wp
    real(wp),parameter :: a3 = .5263563553500217821118595221339767434067222948292617539626688113669112965383197614021_wp
    real(wp),parameter :: a4 = .7895345330250326731677892832009651151100834422438926309440032170503669448074796421031_wp
    real(wp),parameter :: a5 = .1852155685265047076970656199437875208987680209449849941608024849050604421517806064228_wp
    real(wp),parameter :: a6 = .2895345330250326731677892832009651151100834422438926309440032170503669448074796421031_wp
    real(wp),parameter :: a7 = .7659879027055932401898717206953675811356898851157703878511505293452324618973555436361_wp
    real(wp),parameter :: a8 = .1080739009578824490100240661758266719014599665988034811870817465871751958605805741080_wp
    real(wp),parameter :: a9 = .3573842417596774518429245029795604640404982636367873040901247917361510345429002009092_wp
    real(wp),parameter :: a10 = .8825276619647323464255014869796690751828678442680521196637911779185276585194132570617_wp
    real(wp),parameter :: a11 = .6426157582403225481570754970204395359595017363632126959098752082638489654570997990908_wp
    real(wp),parameter :: a12 = .1174723380352676535744985130203309248171321557319478803362088220814723414805867429383_wp
    real(wp),parameter :: a13 = .7659879027055932401898717206953675811356898851157703878511505293452324618973555436361_wp
    real(wp),parameter :: a14 = .2895345330250326731677892832009651151100834422438926309440032170503669448074796421031_wp
    real(wp),parameter :: a15 = .5263563553500217821118595221339767434067222948292617539626688113669112965383197614021_wp
    real(wp),parameter :: a16 = .3357505083417036184129939488963472892525539577157696170900805997207182929518116563365_wp

    real(wp),parameter :: b21 = .3357505083417036184129939488963472892525539577157696170900805997207182929518116563365_wp
    real(wp),parameter :: b31 = .1137717040478783056449396184425018992536986324028992632761012508174661576282091865815_wp
    real(wp),parameter :: b32 = .4125846513021434764669199036914748441530236624263624906865675605494451389101105748206_wp
    real(wp),parameter :: b41 = .1973836332562581682919473208002412787775208605609731577360008042625917362018699105258_wp
    real(wp),parameter :: b43 = .5921508997687745048758419624007238363325625816829194732080024127877752086056097315774_wp
    real(wp),parameter :: b51 = .1360001717992528326665261557210531883399922632350163155393510462488093550801937826726_wp
    real(wp),parameter :: b53 = .8247208052028112223790799775745048180379960003630910417779053602715946236310693542837e-1_wp
    real(wp),parameter :: b54 = -.3325668379302924720736853353471614924502384232634042555633909737090837529152011167815e-1_wp
    real(wp),parameter :: b61 = .6546785199481110579756630335868870724165088583024509548938434492848601278421108871818e-1_wp
    real(wp),parameter :: b64 = .6858715620423033993966826640550281985093888681236135859756994188871385539662876478859e-3_wp
    real(wp),parameter :: b65 = .2233808094681792639708262971782213796699231675455239218686431727029937934693022657371_wp
    real(wp),parameter :: b71 = .2379439999135994223360182301765245366877534054337310606780367192080205084533052094979_wp
    real(wp),parameter :: b74 = .1285793626493546102687479560301629068933431333159680850577212924587597635005521433930_wp
    real(wp),parameter :: b75 = -.7303763776380165799305971651206183012312272113455920027358381669467410091843571268521_wp
    real(wp),parameter :: b76 = 1.129840917780655787515702699609298438785820557711663244851230684625193199127855317597_wp
    real(wp),parameter :: b81 = .6062780896441720418288517632894493561013839292138719356621888281245740094338544165775e-1_wp
    real(wp),parameter :: b85 = .7888222248692022386398029561573557081956734138939980123787302990185683517283562747179e-1_wp
    real(wp),parameter :: b86 = -.3213213504125939005101626782074483911778772517654658420946839341553511609992151857815e-1_wp
    real(wp),parameter :: b87 = .6960045478044110141748620518910045895419574645630705924582272883960758442810235566394e-3_wp
    real(wp),parameter :: b91 = .3104415865438721889084097538819995930358479868352200893591644921124828667005686536920e-1_wp
    real(wp),parameter :: b96 = .1571656120644442171909172634188670081467914525944196307657031602142991040495820182979_wp
    real(wp),parameter :: b97 = .1117638541384084302989687338271505642661663918417037071041752801553969715970722155609e-3_wp
    real(wp),parameter :: b98 = .1690627071867076073308672954386663460258558459670039606814010070304482468516642450265_wp
    real(wp),parameter :: b101 = .1430614201422677732943411456202705981582929884075538519224332811759391256045628709601e-1_wp
    real(wp),parameter :: b106 = -.3914725335338579386439799990026912235591541327277623987271673700062264924050571982953_wp
    real(wp),parameter :: b107 = .2848400967322984698670378405298733357388045123736554424201702352859211869408079104816_wp
    real(wp),parameter :: b108 = .2559426777170804941167590501634194225203620020234561863154326261989141389268124787652_wp
    real(wp),parameter :: b109 = .7189112790349845437562504807270404806670261637579475044631123583223249124963937790142_wp
    real(wp),parameter :: b111 = .1006668457426631788172270822620227260906121023931371723741841060852730989522549641648e-1_wp
    real(wp),parameter :: b116 = -.3683271809355001194608455001594963737510628435069995678571139196388186857259494651043_wp
    real(wp),parameter :: b117 = .8575502624022299903546939609254513705068159712672738891650679949767497156363396523266e-1_wp
    real(wp),parameter :: b118 = .2633792404483482304677774849110341109669789276529274236259927084605103685123897043694_wp
    real(wp),parameter :: b119 = .6783120744401823352568978305425238981627472606453781308422134412090277276134222894570_wp
    real(wp),parameter :: b1110 = -.2657008652719721502394642259236950907890441579413439685514223187307272640162219128045e-1_wp
    real(wp),parameter :: b121 = .4212573015701458022875431871453412827491044618882839077863171149554249469159410384126e-1_wp
    real(wp),parameter :: b126 = -.1140730257394986628323587415013320453327540775054758894945046140858856585220720793160e-1_wp
    real(wp),parameter :: b127 = -.3577323583881041642639854910390843838477094752248779125994771974507179107770816673622e-2_wp
    real(wp),parameter :: b128 = .8478085069047181595883795669333443274198965196004767054783136203447118275436721594692e-1_wp
    real(wp),parameter :: b129 = -.8236897403845190534516723638464155343698796681468739910260266595384301707513187151823e-4_wp
    real(wp),parameter :: b1210 = .7920174936649160445251535005941104644271462692404428580786438040911208546142377592512e-3_wp
    real(wp),parameter :: b1211 = .4840734825985701173601980408776943260994401783442431626214940796417131157064341867563e-2_wp
    real(wp),parameter :: b131 = .1866168584194642618748985099397434698267694490527249444362505121514680479505442836574_wp
    real(wp),parameter :: b134 = .1285793626493546102687479560301629068933431333159680850577212924587597635005521433930_wp
    real(wp),parameter :: b135 = -.7303763776380165799305971651206183012312272113455920027358381669467410091843571268521_wp
    real(wp),parameter :: b136 = .3365620746651354218129539461265307620933082614096471567813938850533189571406833598679_wp
    real(wp),parameter :: b137 = .3439152195696137654283681672732838006806635738118480999635325742827372086353048385257_wp
    real(wp),parameter :: b138 = .6116809614283795957124860804122425718138169583157463948158033894717678649195939132136_wp
    real(wp),parameter :: b139 = .8555642651045021202243366091904917631763787495189342130253621820535903092162656518481_wp
    real(wp),parameter :: b1310 = -.8913928524866960048456401557668818573165120479499861285407465844438950550960191679472e-1_wp
    real(wp),parameter :: b1311 = -.4263317531098201582091498458499961586420389075833366648914840378325891404712881232497_wp
    real(wp),parameter :: b1312 = -.4510834231343501965076085217297850477436729165851712257475164429026900343003414799732_wp
    real(wp),parameter :: b141 = .8132722174581177948072809791683141502457873136125978954741528825822682868334020231374e-1_wp
    real(wp),parameter :: b144 = .6858715620423033993966826640550281985093888681236135859756994188871385539662876478859e-3_wp
    real(wp),parameter :: b145 = .2233808094681792639708262971782213796699231675455239218686431727029937934693022657371_wp
    real(wp),parameter :: b146 = -.3543152669595614715438074770044728248424634682417120916726014150358797064660842419218_wp
    real(wp),parameter :: b147 = -.1614228337471357340842753214971350545154909424968569413717974406833347832344548170007_wp
    real(wp),parameter :: b148 = -1.267989518319081979614362221100086229540210555088918164010191288219362963604567702803_wp
    real(wp),parameter :: b149 = .2482505660965057784745955362600844051763121832196360142449926255422280563592831868088_wp
    real(wp),parameter :: b1410 = .1894015072204329862304305088267495664543558350848310970103331708369628629196601961145e-2_wp
    real(wp),parameter :: b1411 = -.2367108504250404697299275875477514623653536202773466252372015412173547424356865161882e-1_wp
    real(wp),parameter :: b1412 = 1.376373544047006195976245547649688122109906103559164922810806152088122419417100442009_wp
    real(wp),parameter :: b1413 = .1650212091015662542191305948002865244010106371945579174943772453918520072439660689694_wp
    real(wp),parameter :: b151 = .1137717040478783056449396184425018992536986324028992632761012508174661576282091865815_wp
    real(wp),parameter :: b152 = .4125846513021434764669199036914748441530236624263624906865675605494451389101105748206_wp
    real(wp),parameter :: b156 = -.5266141894222268880524004361760397955773267936583064281992495328485812564705031059898_wp
    real(wp),parameter :: b157 = .4270299995350464166415042303891280056028480962227556089080302280248036279775533169729_wp
    real(wp),parameter :: b1513 = -.4270299995350464166415042303891280056028480962227556089080302280248036279775533169729_wp
    real(wp),parameter :: b1514 = .5266141894222268880524004361760397955773267936583064281992495328485812564705031059898_wp
    real(wp),parameter :: b161 = .3357505083417036184129939488963472892525539577157696170900805997207182929518116563365_wp
    real(wp),parameter :: b163 = -.4368367083891407610190158064846822677554138099678168504498524386288907219928146135373_wp
    real(wp),parameter :: b1615 = .4368367083891407610190158064846822677554138099678168504498524386288907219928146135373_wp
    real(wp),parameter :: b171 = .3528400914539065757695374250020526983792042094778706633574048575138598395445375459964e-1_wp
    real(wp),parameter :: b172 = -.4368588562339554617519938616027395568153523718226810275461527253101515442644071737416_wp
    real(wp),parameter :: b173 = -.5185253025751911700815354370658368214653404020923461413204927007907407344767758710502_wp
    real(wp),parameter :: b176 = .8353882146318347708823979907871444293693024897533008980286851134950209556868285481998e-1_wp
    real(wp),parameter :: b177 = .3357324883823615793514438806416832609180686549969009829861766772233990660034579152038_wp
    real(wp),parameter :: b178 = -.1180346853290197678148440843338719378456952414935754110833043790640005049977678056530_wp
    real(wp),parameter :: b179 = -.2028121524999718341813208304021045923675061150408674489110338814798518734278169865166_wp
    real(wp),parameter :: b1710 = .4018734442526045674154192467995252408036818247069631361197899796078897150680466096292_wp
    real(wp),parameter :: b1711 = .6982808681238487516481594542254805552510541794271997632640921307356093400075120505480_wp
    real(wp),parameter :: b1712 = .2745953340127460413630142920026546525788337913278252157685398122819348957290706395825_wp
    real(wp),parameter :: b1713 = -.8071051158813288531951186779278104056732084407804403857002131168492985419909905100904_wp
    real(wp),parameter :: b1714 = .2986469883301853807480531774155235135599206769328769914173437804434298240853514778770_wp
    real(wp),parameter :: b1715 = .5185253025751911700815354370658368214653404020923461413204927007907407344767758710502_wp
    real(wp),parameter :: b1716 = .4368588562339554617519938616027395568153523718226810275461527253101515442644071737416_wp

    real(wp),parameter :: c1 = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1_wp
    real(wp),parameter :: c2 = -.2192242833052276559865092748735244519392917369308600337268128161888701517706576728499e-1_wp
    real(wp),parameter :: c3 = -.5671077504725897920604914933837429111531190926275992438563327032136105860113421550095e-1_wp
    real(wp),parameter :: c6 = -.5604719764011799410029498525073746312684365781710914454277286135693215339233038348083e-1_wp
    real(wp),parameter :: c7 = .1789297658862876254180602006688963210702341137123745819397993311036789297658862876254_wp
    real(wp),parameter :: c9 = .2774291885177431765083602625606543404285043197180408363394722409866844803871713937960_wp
    real(wp),parameter :: c10 = .1892374781489234901583064041060123262381623469486258303271944256799821862794952728707_wp
    real(wp),parameter :: c11 = .2774291885177431765083602625606543404285043197180408363394722409866844803871713937960_wp
    real(wp),parameter :: c12 = .1892374781489234901583064041060123262381623469486258303271944256799821862794952728707_wp
    real(wp),parameter :: c13 = -.1789297658862876254180602006688963210702341137123745819397993311036789297658862876254_wp
    real(wp),parameter :: c14 = .5604719764011799410029498525073746312684365781710914454277286135693215339233038348083e-1_wp
    real(wp),parameter :: c15 = .5671077504725897920604914933837429111531190926275992438563327032136105860113421550095e-1_wp
    real(wp),parameter :: c16 = .2192242833052276559865092748735244519392917369308600337268128161888701517706576728499e-1_wp
    real(wp),parameter :: c17 = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1_wp

    real(wp),dimension(me%n) :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17

    call me%f(t+h,   x,f1)
    call me%f(t+a2*h,x+h*(b21*f1),f2)
    call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
    call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
    call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
    call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
    call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
    call me%f(t+a8*h,x+h*(b81*f1+b85*f5+b86*f6+b87*f7),f8)
    call me%f(t+a9*h,x+h*(b91*f1+b96*f6+b97*f7+b98*f8),f9)
    call me%f(t+a10*h,x+h*(b101*f1+b106*f6+b107*f7+&
                           b108*f8+b109*f9),f10)
    call me%f(t+a11*h,x+h*(b111*f1+b116*f6+b117*f7+&
                           b118*f8+b119*f9+b1110*f10),f11)
    call me%f(t+a12*h,x+h*(b121*f1+b126*f6+b127*f7+&
                           b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
    call me%f(t+a13*h,x+h*(b131*f1+b134*f4+b135*f5+b136*f6+b137*f7+&
                           b138*f8+b139*f9+b1310*f10+b1311*f11+b1312*f12),f13)
    call me%f(t+a14*h,x+h*(b141*f1+b144*f4+b145*f5+b146*f6+b147*f7+&
                           b148*f8+b149*f9+b1410*f10+b1411*f11+b1412*f12+b1413*f13),f14)
    call me%f(t+a15*h,x+h*(b151*f1+b152*f2+b156*f6+b157*f7+&
                           b1513*f13+b1514*f14),f15)
    call me%f(t+a16*h,x+h*(b161*f1+b163*f3+b1615*f15),f16)
    call me%f(t+h,    x+h*(b171*f1+b172*f2+b173*f3+b176*f6+b177*f7+&
                           b178*f8+b179*f9+b1710*f10+b1711*f11+b1712*f12+b1713*f13+&
                           b1714*f14+b1715*f15+b1716*f16),f17)

    xf = x+h*(c1*f1+c2*f2+c3*f3+c6*f6+c7*f7+c9*f9+c10*f10+&
              c11*f11+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16+c17*f17)

    end procedure rko10
!*****************************************************************************************

!*****************************************************************************************
!>
!  Hairer 10th order method.
!
!### References
!  * Ernst Hairer, "A Runge-Kutta Method of Order 10"
!    January 1978, IMA Journal of Applied Mathematics 21(1)

    module procedure rkh10

    real(wp),parameter :: a2 = .5233584004620047139632937023215170497515953383496781610502942213792083195573343614542_wp
    real(wp),parameter :: a3 = .5265091001416125727329516734775408259523008085442588621240322448552569517961160776999_wp
    real(wp),parameter :: a4 = .7897636502124188590994275102163112389284512128163882931860483672828854276941741165498_wp
    real(wp),parameter :: a5 = .3939235701256720143227738119996521367488350094732736567444747389961242969755195414769_wp
    real(wp),parameter :: a6 = .7666539862535505911932668693560686601417881683220066628197494388337786524595845606945_wp
    real(wp),parameter :: a7 = .2897636502124188590994275102163112389284512128163882931860483672828854276941741165498_wp
    real(wp),parameter :: a8 = .1084776892195672933536461100396272126536495082872530661076180009074872219675135232227_wp
    real(wp),parameter :: a9 = .3573842417596774518429245029795604640404982636367873040901247917361510345429002009092_wp
    real(wp),parameter :: a10 = .8825276619647323464255014869796690751828678442680521196637911779185276585194132570617_wp
    real(wp),parameter :: a11 = .6426157582403225481570754970204395359595017363632126959098752082638489654570997990908_wp
    real(wp),parameter :: a12 = .1174723380352676535744985130203309248171321557319478803362088220814723414805867429383_wp
    real(wp),parameter :: a13 = .7666539862535505911932668693560686601417881683220066628197494388337786524595845606945_wp
    real(wp),parameter :: a14 = .2897636502124188590994275102163112389284512128163882931860483672828854276941741165498_wp
    real(wp),parameter :: a15 = .5265091001416125727329516734775408259523008085442588621240322448552569517961160776999_wp
    real(wp),parameter :: a16 = .5233584004620047139632937023215170497515953383496781610502942213792083195573343614542_wp

    real(wp),parameter :: b21 = .5233584004620047139632937023215170497515953383496781610502942213792083195573343614542_wp
    real(wp),parameter :: b31 = .2616697163778127283312402097548997641973627614039327706102102491953717704187699219879_wp
    real(wp),parameter :: b32 = .2648393837637998444017114637226410617549380471403260915138219956598851813773461557120_wp
    real(wp),parameter :: b41 = .1974409125531047147748568775540778097321128032040970732965120918207213569235435291375_wp
    real(wp),parameter :: b43 = .5923227376593141443245706326622334291963384096122912198895362754621640707706305874124_wp
    real(wp),parameter :: b51 = .1973205486287023067036649485978952117578825584337199991656132317825743833021058860612_wp
    real(wp),parameter :: b53 = .2950833340926721918228255598274359230145509784596048092593866299668586683096812321792_wp
    real(wp),parameter :: b54 = -.9848031259570248420371669642567899802359852742005115168052512275330875463626757676351e-1_wp
    real(wp),parameter :: b61 = .1313134173444616536130177999345909470542768366990469474228776563495603985387834598888_wp
    real(wp),parameter :: b64 = .1101544395386396206773696967716892932905881833462590372574439152868807925960672597269_wp
    real(wp),parameter :: b65 = .5251861293704493169028793726497884197969231482767006781394278671973374613247338410788_wp
    real(wp),parameter :: b71 = .1342003418463226002727476951680931444018781918996634475042938727232827990267510411664_wp
    real(wp),parameter :: b74 = .6960887032881160802299824047678314828416281106463280160258778625630871045892118300249_wp
    real(wp),parameter :: b75 = .2504977215703398097125518092509800218006823479445679226772611578145189247821334145350_wp
    real(wp),parameter :: b76 = -.7910231164923596311158543989705934101157374376741710930213845258180034007039221691766_wp
    real(wp),parameter :: b81 = .7221827418966261942005081845517049770474117717860435538167259617193885317787239789517e-1_wp
    real(wp),parameter :: b85 = -.5833632293645610716380606138930651874161115931850840194371530896953184153028603826158e-1_wp
    real(wp),parameter :: b86 = .3047557668574525220174950070294036092519082552287015783049982467857131629284936232933e-2_wp
    real(wp),parameter :: b87 = .9154818029778625587722640290346919759800040787487009688661073123722307869064222735619e-1_wp
    real(wp),parameter :: b91 = .3125500813516617947050120528263476612150928811800844295622424453129979711857118942140e-1_wp
    real(wp),parameter :: b96 = .1091238215424128929294834955207716494619546474783251185719596191189550956073995218558e-3_wp
    real(wp),parameter :: b97 = .1567257586309938356246107465648794708917441337700138495832023584660674808897051732313_wp
    real(wp),parameter :: b98 = .1692943511719750238548830676365254553777828871012866864321262291196648014390164387346_wp
    real(wp),parameter :: b101 = .1190660441466861924216884258080925099509546061156163781761478570338134316043617322173e-1_wp
    real(wp),parameter :: b106 = .2834370820246027860255992266982188665428702252125274101591962479489267620092644600260_wp
    real(wp),parameter :: b107 = -.4163121675706282353724276181356300212164192944619403444080980777942870933794472685216_wp
    real(wp),parameter :: b108 = .2646463339497663668210902091085361027053564926326924812994390366777834273576276339310_wp
    real(wp),parameter :: b109 = .7388498091463228097090708267277348761559649602732109347956391853827232193715322584046_wp
    real(wp),parameter :: b111 = .2340657369133197891470838377984007842503946857756845416362339865999662377057916960290e-1_wp
    real(wp),parameter :: b116 = .9449313018949365401300253095605614324982516623777346096448800625130954018295939047740e-1_wp
    real(wp),parameter :: b117 = -.2728720559019952606363092580665963250433705067252372208829562550632597611313198545757_wp
    real(wp),parameter :: b118 = .2240220461156057997944315522518131846261124699330640294272458923970695162204120486767_wp
    real(wp),parameter :: b119 = .6043814410751657569719347222576085340011863610739072982875214128849988434755301302413_wp
    real(wp),parameter :: b1110 = -.3081537692927938090069243415828207929929122273386332605004724686626579706106108533173e-1_wp
    real(wp),parameter :: b121 = .4544377531017616315765389908153096498645890941991961177836633137902353666760633687088e-1_wp
    real(wp),parameter :: b126 = -.1187996671864028586765254219285356343376285990176386474891150929245660355742130860183e-2_wp
    real(wp),parameter :: b127 = .1203565499092261097966188217234362058515446695047694116231895919296441480090125575596e-1_wp
    real(wp),parameter :: b128 = .7512690298764966821627521371565572140275315500656240513504528112598150181393445048666e-1_wp
    real(wp),parameter :: b129 = -.1822092409888012403141186105974838892761579859192182074696828369088959767419077567178e-1_wp
    real(wp),parameter :: b1210 = -.2571528540841043468806376221771396205460354181513400383106090430345956162981031278207e-3_wp
    real(wp),parameter :: b1211 = .4532078371347468185965270952011502734303744355238469520648294046672741844375709484540e-2_wp
    real(wp),parameter :: b131 = .1767137782592772030958798765711993346076326211800572275450227165783753236705910865492_wp
    real(wp),parameter :: b134 = .1101544395386396206773696967716892932905881833462590372574439152868807925960672597269_wp
    real(wp),parameter :: b135 = .5251861293704493169028793726497884197969231482767006781394278671973374613247338410788_wp
    real(wp),parameter :: b136 = -.4716207672801957948798217912152359376250630852495511063738116933651587031904328351457_wp
    real(wp),parameter :: b137 = .8990310498491875266368990071875152922763468480002185650326986125011485318362907529907_wp
    real(wp),parameter :: b138 = -.7467230306916289638599602008088168117750310724922743198498253813592425510843163068237_wp
    real(wp),parameter :: b139 = -1.017101516756146040853186972006065972987027196800421553809421717321497529906933631477_wp
    real(wp),parameter :: b1310 = .1263508715195988962951307827687648346421985369266969430473204298972536422365713122404_wp
    real(wp),parameter :: b1311 = .5660138272355064270682732249907470012763799581315503842554078250210353407723389384909_wp
    real(wp),parameter :: b1312 = .5986492052088624001098038724464832066388402270027708075754868643976463442046741430643_wp
    real(wp),parameter :: b141 = .1277534947480869822694777006880571541639616513225826576695303067404023367054772185702_wp
    real(wp),parameter :: b144 = .6960887032881160802299824047678314828416281106463280160258778625630871045892118300249_wp
    real(wp),parameter :: b145 = .2504977215703398097125518092509800218006823479445679226772611578145189247821334145350_wp
    real(wp),parameter :: b146 = -.7368246436028416867609246757454535374296880219263938462439002090823944915566264811824_wp
    real(wp),parameter :: b147 = -.2778578777108241826773273374900723250222301109862216853553157201018147214465526588169_wp
    real(wp),parameter :: b148 = -.5997526313598403501296884799197753021563938240370770948150479630779446286262003432092_wp
    real(wp),parameter :: b149 = .2024692338910704693500237585621903123505161701229471467587157451308903694383321235511_wp
    real(wp),parameter :: b1410 = .5432036982363849780600684652634443601468189969678666775046718813224989883416104871445e-2_wp
    real(wp),parameter :: b1411 = -.1074472474155047920101206919894381337125444664272205024314798936418733258920769563370e-1_wp
    real(wp),parameter :: b1412 = .6951688484570234004700591858164146072357628221597117426434839740273190245052113679250_wp
    real(wp),parameter :: b1413 = -.6246651130952503394431547116755180508600167575701318270645551618021614799102076408562e-1_wp
    real(wp),parameter :: b151 = .2616697163778127283312402097548997641973627614039327706102102491953717704187699219879_wp
    real(wp),parameter :: b152 = .2648393837637998444017114637226410617549380471403260915138219956598851813773461557120_wp
    real(wp),parameter :: b156 = -.1998011270205324791079663580830885049848273745422651189682301346802905866051733476638_wp
    real(wp),parameter :: b157 = -.6510499873052827124921914489683813643155863882516440645794556633240216912803403931627_wp
    real(wp),parameter :: b1513 = .1998011270205324791079663580830885049848273745422651189682301346802905866051733476638_wp
    real(wp),parameter :: b1514 = .6510499873052827124921914489683813643155863882516440645794556633240216912803403931627_wp
    real(wp),parameter :: b161 = .5233584004620047139632937023215170497515953383496781610502942213792083195573343614542_wp
    real(wp),parameter :: b163 = -.5558812136754302060726143105309293455559184141943321053532734480099926250948077261183_wp
    real(wp),parameter :: b1615 = .5558812136754302060726143105309293455559184141943321053532734480099926250948077261183_wp
    real(wp),parameter :: b171 = .5732079543206559103114261705103983656495216504867462310285994428078568043160654439795e-1_wp
    real(wp),parameter :: b172 = -.5499710763899945608115841896290187887481592249811405834035066676393750158953834290913_wp
    real(wp),parameter :: b173 = -.6499374174008749135116607420010890619711618624173024222960650740195874521599402439688_wp
    real(wp),parameter :: b176 = -1.061667370401756207240019539023157074172524666307437022389776456477183230723296269940_wp
    real(wp),parameter :: b177 = -.4040156689806358294269682234212183308262562023912486365220642577870402491555711062480e-1_wp
    real(wp),parameter :: b178 = -.1828302366407607254710272774065261039379052622607190097473388370699414811305446343873_wp
    real(wp),parameter :: b179 = -.3336592706492786845666575661828162687906558601961826440714525336287466822150370633233_wp
    real(wp),parameter :: b1710 = .3956485423760567568801345107166015519577734440834727480004748180136901286634710478955_wp
    real(wp),parameter :: b1711 = .6950570494599735891002099282005158129027126868215679095299345058137097320818106877162_wp
    real(wp),parameter :: b1712 = .2714873764573748588377263058539220945263829691804714618529052530298982146739754552950_wp
    real(wp),parameter :: b1713 = .6071810560414041202873774349794680164722661545496003750296400378855628528787164400954_wp
    real(wp),parameter :: b1714 = .5918636248229842840838104081530739675596239893196764223449596939309288102548549028752_wp
    real(wp),parameter :: b1715 = .6499374174008749135116607420010890619711618624173024222960650740195874521599402439688_wp
    real(wp),parameter :: b1716 = .5499710763899945608115841896290187887481592249811405834035066676393750158953834290913_wp

    real(wp),parameter :: c1 = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1_wp
    real(wp),parameter :: c2 = -.3846153846153846153846153846153846153846153846153846153846153846153846153846153846154e-1_wp
    real(wp),parameter :: c3 = -.9090909090909090909090909090909090909090909090909090909090909090909090909090909090909e-1_wp
    real(wp),parameter :: c6 = -.1348314606741573033707865168539325842696629213483146067415730337078651685393258426966_wp
    real(wp),parameter :: c7 = -.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111_wp
    real(wp),parameter :: c9 = .2774291885177431765083602625606543404285043197180408363394722409866844803871713937960_wp
    real(wp),parameter :: c10 = .1892374781489234901583064041060123262381623469486258303271944256799821862794952728707_wp
    real(wp),parameter :: c11 = .2774291885177431765083602625606543404285043197180408363394722409866844803871713937960_wp
    real(wp),parameter :: c12 = .1892374781489234901583064041060123262381623469486258303271944256799821862794952728707_wp
    real(wp),parameter :: c13 = .1348314606741573033707865168539325842696629213483146067415730337078651685393258426966_wp
    real(wp),parameter :: c14 = .1111111111111111111111111111111111111111111111111111111111111111111111111111111111111_wp
    real(wp),parameter :: c15 = .9090909090909090909090909090909090909090909090909090909090909090909090909090909090909e-1_wp
    real(wp),parameter :: c16 = .3846153846153846153846153846153846153846153846153846153846153846153846153846153846154e-1_wp
    real(wp),parameter :: c17 = .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1_wp

    real(wp),dimension(me%n) :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17

    call me%f(t+h,   x,f1)
    call me%f(t+a2*h,x+h*(b21*f1),f2)
    call me%f(t+a3*h,x+h*(b31*f1+b32*f2),f3)
    call me%f(t+a4*h,x+h*(b41*f1+b43*f3),f4)
    call me%f(t+a5*h,x+h*(b51*f1+b53*f3+b54*f4),f5)
    call me%f(t+a6*h,x+h*(b61*f1+b64*f4+b65*f5),f6)
    call me%f(t+a7*h,x+h*(b71*f1+b74*f4+b75*f5+b76*f6),f7)
    call me%f(t+a8*h,x+h*(b81*f1+b85*f5+b86*f6+b87*f7),f8)
    call me%f(t+a9*h,x+h*(b91*f1+b96*f6+b97*f7+b98*f8),f9)
    call me%f(t+a10*h,x+h*(b101*f1+b106*f6+b107*f7+b108*f8+b109*f9),f10)
    call me%f(t+a11*h,x+h*(b111*f1+b116*f6+b117*f7+b118*f8+b119*f9+b1110*f10),f11)
    call me%f(t+a12*h,x+h*(b121*f1+b126*f6+b127*f7+b128*f8+b129*f9+b1210*f10+b1211*f11),f12)
    call me%f(t+a13*h,x+h*(b131*f1+b134*f4+b135*f5+b136*f6+b137*f7+b138*f8+b139*f9+&
                           b1310*f10+b1311*f11+b1312*f12),f13)
    call me%f(t+a14*h,x+h*(b141*f1+b144*f4+b145*f5+b146*f6+b147*f7+b148*f8+b149*f9+&
                           b1410*f10+b1411*f11+b1412*f12+b1413*f13),f14)
    call me%f(t+a15*h,x+h*(b151*f1+b152*f2+b156*f6+b157*f7+b1513*f13+b1514*f14),f15)
    call me%f(t+a16*h,x+h*(b161*f1+b163*f3+b1615*f15),f16)
    call me%f(t+h,    x+h*(b171*f1+b172*f2+b173*f3+b176*f6+b177*f7+b178*f8+b179*f9+&
                           b1710*f10+b1711*f11+b1712*f12+b1713*f13+b1714*f14+b1715*f15+b1716*f16),f17)

    xf = x+h*(c1*f1+c2*f2+c3*f3+c6*f6+c7*f7+c9*f9+c10*f10+&
              c11*f11+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16+c17*f17)

    end procedure rkh10
!*****************************************************************************************

!*****************************************************************************************
    end submodule rklib_fixed_steps
!*****************************************************************************************