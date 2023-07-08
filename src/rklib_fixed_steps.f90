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

    if (h==zero) then
        xf = x
        return
    end if

    call me%f(t,x,f1)

    xf = x + h*f1

    end procedure euler
!*****************************************************************************************

!*****************************************************************************************
!>
!  Midpoint (2nd order) integration method.

    module procedure midpoint

    real(wp),dimension(me%n) :: f1,f2

    if (h==zero) then
        xf = x
        return
    end if

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

    if (h==zero) then
        xf = x
        return
    end if

    call me%f(t,x,f1)
    call me%f(t+h,x+h*f1,f2)

    xf = x + 0.5_wp*h*(f1+f2)

    end procedure heun
!*****************************************************************************************

!*****************************************************************************************
!>
!  Shu and Osher's (2nd order) TVD integration method
!
!### Reference
!   * C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory
!   shock-capturing schemes", Journal of Computational Physics, 77, 1988, 439-471.
!   https://doi.org/10.1016/0021-9991(88)90177-5.
    module procedure rktvd2

    real(wp),dimension(me%n) :: fs

    if (h==zero) then
        xf = x
        return
    end if

    call me%f(t, x, fs)
    xf = x + h*fs
    call me%f(t + h, xf, fs)

    xf = (x + xf + h*fs)/2

    end procedure rktvd2
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
    real(wp),parameter :: b3  =  1.0_wp
    real(wp),parameter :: c21 =  1.0_wp/2.0_wp
    real(wp),parameter :: c31 = -1.0_wp
    real(wp),parameter :: c32 =  2.0_wp

    if (h==zero) then
        xf = x
        return
    end if

    call me%f(t,      x,                   f1)
    call me%f(t+b2*h, x+h*c21*f1,          f2)
    call me%f(t+b3*h, x+h*(c31*f1+c32*f2), f3)

    xf = x + h*( a1*f1 + a2*f2 + a3*f3 )

    end procedure rk3
!*****************************************************************************************

!*****************************************************************************************
!>
!  Shu and Osher's (3rd order) TVD integration method
!
!### Reference
!   * C.-W. Shu, S. Osher, "Efficient implementation of essentially non-oscillatory
!   shock-capturing schemes", Journal of Computational Physics, 77, 1988, 439-471.
!   https://doi.org/10.1016/0021-9991(88)90177-5.
    module procedure rktvd3

        real(wp),dimension(me%n) :: xs, fs

        if (h==zero) then
            xf = x
            return
        end if

        call me%f(t, x, fs)
        xs = x + h*fs
        call me%f(t + h, xs, fs)
        xs = (3*x + xs + h*fs)/4
        call me%f(t + h, xs, fs)

        xf = (x + 2*xs + 2*h*fs)/3

     end procedure rktvd3
!*****************************************************************************************

!*****************************************************************************************
!>
!  Take one Runge Kutta 4 integration step: `t -> t+h (x -> xf)`

    module procedure rk4

    real(wp),dimension(me%n) :: f1,f2,f3,f4

    if (h==zero) then
        xf = x
        return
    end if

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
    real(wp),parameter :: a3  =  1.0_wp
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

    if (h==zero) then
        xf = x
        return
    end if

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+aa1*h*f0,f1)
    call me%f(t+a2*h,x+aa2*h*(b20*f0+b21*f1),f2)
    call me%f(t+a3*h,x+aa3*h*(b30*f0+b31*f1+b32*f2),f3)

    xf = x + h*c*(c0*f0+c1*f1+c2*f2+c3*f3)

    end procedure rks4
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
    real(wp),parameter :: a4  =  1.0_wp
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

    if (h==zero) then
        xf = x
        return
    end if

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+aa1*h*f0,f1)
    call me%f(t+a2*h,x+aa2*h*(b20*f0+b21*f1),f2)
    call me%f(t+a3*h,x+aa3*h*(b30*f0+b31*f1+b32*f2),f3)
    call me%f(t+a4*h,x+aa4*h*(b40*f0+b41*f1+b42*f2+b43*f3),f4)

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
    real(wp),parameter :: a8 = 1.0_wp

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

    real(wp),parameter :: b20 = 1.0_wp
    real(wp),parameter :: b21 = 3.0_wp

    real(wp),parameter :: b30 = 1.0_wp
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

    if (h==zero) then
        xf = x
        return
    end if

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+aa1*h*(f0),f1)
    call me%f(t+a2*h,x+aa2*h*(b20*f0+b21*f1),f2)
    call me%f(t+a3*h,x+aa3*h*(b30*f0+b32*f2),f3)
    call me%f(t+a4*h,x+aa4*h*(b40*f0+b42*f2+b43*f3),f4)
    call me%f(t+a5*h,x+aa5*h*(b50*f0+b52*f2+b53*f3+b54*f4),f5)
    call me%f(t+a6*h,x+aa6*h*(b60*f0+b62*f2+b63*f3+b64*f4+b65*f5),f6)
    call me%f(t+a7*h,x+aa7*h*(b70*f0+b72*f2+b73*f3+b74*f4+b75*f5+b76*f6),f7)
    call me%f(t+a8*h,x+aa8*h*(b80*f0+b82*f2+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)

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

    if (h==zero) then
        xf = x
        return
    end if

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
    real(wp),parameter :: a11   =  1.0_wp
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
    real(wp),parameter :: b20   =  1.0_wp
    real(wp),parameter :: b21   =  3.0_wp
    real(wp),parameter :: b30   =  1.0_wp
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
    call me%f(t+a2*h,x+aa2*h*(b20*f0+b21*f1),f2)
    call me%f(t+a3*h,x+aa3*h*(b30*f0+b32*f2),f3)
    call me%f(t+a4*h,x+aa4*h*(b40*f0+b42*f2+b43*f3),f4)
    call me%f(t+a5*h,x+aa5*h*(b50*f0+b53*f3+b54*f4),f5)
    call me%f(t+a6*h,x+aa6*h*(b60*f0+b63*f3+b64*f4+b65*f5),f6)
    call me%f(t+a7*h,x+aa7*h*(b70*f0+b73*f3+b74*f4+b76*f6),f7)
    call me%f(t+a8*h,x+aa8*h*(b80*f0+b83*f3+b84*f4+b85*f5+b86*f6+b87*f7),f8)
    call me%f(t+a9*h,x+aa9*h*(b90*f0+b93*f3+b94*f4+b95*f5+b96*f6+b97*f7+b98*f8),f9)
    call me%f(t+a10*h,x+aa10*h*(b100*f0+b103*f3+b104*f4+b105*f5+b106*f6+b109*f9),f10)
    call me%f(t+a11*h,x+aa11*h*(b110*f0+b113*f3+b114*f4+b115*f5+b116*f6+b117*f7+b118*f8+b119*f9+b1110*f10),f11)

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
    real(wp),parameter :: a11 = 1

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

    if (h==zero) then
        xf = x
        return
    end if

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
    real(wp),parameter :: c2 = +0.00000000000000000000000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c3 = +0.04681369289018421954398607025172424191865836856120072277976273594277845007872_wp
    real(wp),parameter :: c4 = +0.00000000000000000000000000000000000000000000000000000000000000000000000000000_wp
    real(wp),parameter :: c5 = +0.00000000000000000000000000000000000000000000000000000000000000000000000000000_wp
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

    if (h==zero) then
        xf = x
        return
    end if

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

    xf = x+h*(c1*f1+c2*f2+c3*f3+c4*f4+c5*f5+c6*f6+c7*f7+c8*f8+c9*f9+c10*f10+&
              c11*f11+c12*f12+c13*f13+c14*f14+c15*f15+c16*f16)

    end procedure rkz10
!*****************************************************************************************

!*****************************************************************************************
    end submodule rklib_fixed_steps
!*****************************************************************************************