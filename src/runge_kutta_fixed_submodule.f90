!*****************************************************************************************
!>
!  Fixed-step RK formulas

    submodule(runge_kutta_module) runge_kutta_fixed_submodule

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

!**********************************************************************************
!>
! Runge Kutta Shanks (5th order)
!
!### Reference
!  * E. B. Shanks, "Solutions of Differential Equations by Evaluations of Functions"
!     Math. Comp. 20 (1966).

    module procedure rks5

    real(wp),dimension(me%n) :: f0,f1,f2,f3,f4

    real(wp),parameter :: a1 = 1.0_wp / 9000.0_wp
    real(wp),parameter :: a2 = 3.0_wp / 10.0_wp
    real(wp),parameter :: a3 = 3.0_wp / 4.0_wp
    real(wp),parameter :: a4 = 1.0_wp

    real(wp),parameter :: c = 1.0_wp / 1134.0_wp

    real(wp),parameter :: c0 = 105.0_wp
    real(wp),parameter :: c2 = 500.0_wp
    real(wp),parameter :: c3 = 448.0_wp
    real(wp),parameter :: c4 = 81.0_wp

    real(wp),parameter :: aa1 = 1.0_wp / 9000.0_wp
    real(wp),parameter :: aa2 = 1.0_wp / 10.0_wp
    real(wp),parameter :: aa3 = 1.0_wp / 8.0_wp
    real(wp),parameter :: aa4 = 1.0_wp / 81.0_wp

    real(wp),parameter :: b20 = -4047.0_wp
    real(wp),parameter :: b21 = 4050.0_wp

    real(wp),parameter :: b30 = 20241.0_wp
    real(wp),parameter :: b31 = -20250.0_wp
    real(wp),parameter :: b32 = 15.0_wp

    real(wp),parameter :: b40 = -931041.0_wp
    real(wp),parameter :: b41 = 931500.0_wp
    real(wp),parameter :: b42 = -490.0_wp
    real(wp),parameter :: b43 = 112.0_wp

    call me%f(t,x,f0)
    call me%f(t+a1*h,x+aa1*h*f0,f1)
    call me%f(t+a2*h,x+aa2*h*(b20*f0+b21*f1),f2)
    call me%f(t+a3*h,x+aa3*h*(b30*f0+b31*f1+b32*f2),f3)
    call me%f(t+a4*h,x+aa4*h*(b40*f0+b41*f1+b42*f2+b43*f3),f4)

    xf = x + h * c * (c0*f0 + c2*f2 + c3*f3 + c4*f4)

    end procedure rks5
!**********************************************************************************

!*****************************************************************************************
!>
!  Take one Runge Kutta 7 integration step: `t -> t+h (x -> xf)`
!
!### Reference
!  * E. B. Shanks, "Solutions of Differential Equations by Evaluations of Functions"
!     Math. Comp. 20 (1966).

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
    end submodule runge_kutta_fixed_submodule
!*****************************************************************************************