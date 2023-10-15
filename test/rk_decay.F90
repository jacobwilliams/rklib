program rk_decay
!! Unit test of the [[rk_module]].
!! This test complements the test in file [[rk_test_order]]; that test
!! checks the correctness of the coefficients other than b-ij in formulae such as
!!
!!     call me%f(t+a6*h, x+h*(b61*f1 + b62*f2 + b63*f3 + b64*f4 + b65*f5), f6)
!!
!! Errors in the second argument are not caught in that example since
!! the R.H.S. of the ODE used in it depends only on the independent variable, i.e., (\dot x = t^n).
!!
!! The test in this file checks if the b-ij coefficients and the expressions in the second arguments
!! to the calls to me%f are likewise free from typographical errors.
!! The test problem is the autonomous equation \dot x = -x, with the initial condition
!! x(0) = 1, integrated until t = 1.
!!
!!### History
!! * Original Author: `@mecej4` on https://fortran-lang.discourse.group/t/a-new-runge-kutta-fortran-library/6055/30

    use rklib_module, wp => rk_module_rk
    implicit none

    real(wp), parameter :: zero = 0.0_wp, one = 1.0_wp
    real(wp), parameter :: eps = epsilon(one)

    class(rk_class), allocatable :: s

    ! Test all methods
    write(*, '(a10, a10, a24, a9)') 'Method', "Variable", "error", "Success"
    write(*, *) "----------------------------------------------------"

#include "rklib_allocate_and_test.inc"

contains

    subroutine run_test()
        !! Order test
        use face, only: colorize

        integer, parameter :: n = 1
        real(wp), parameter :: t0 = zero, tf = one, dt = 1.0e-4_wp, x0(n) = one, tol=1.0e-12_wp
        real(wp) :: xf(n)
        type(rklib_properties) :: p
        type(stepsize_class) :: sz
        character(:), allocatable :: method
        logical :: check, isvariable
        character(len=100) :: tmp

        p = s%properties()
        method = p%short_name

        ! step size constructor
        call sz%initialize(fixed_step_mode = .false.)

        select type (s)
        class is (rk_fixed_step_class)
            isvariable = .false.
            call s%initialize(n=n, f=derivative_test_order)
            call s%integrate(t0, x0, dt, tf, xf)

        class is (rk_variable_step_class)
            isvariable = .true.
            call s%initialize(n=n, f=derivative_test_order, rtol=[tol], atol=[tol], &
                              stepsize_method=sz , hinit_method=2)
            call s%integrate(t0, x0, dt, tf, xf)

        end select

        check = abs(xf(1) - exp(-one)) <= 1.0e-4_wp
        write(tmp, '(a10, l10, es24.15e2, l9)') method, isvariable, abs(xf(1) - exp(-one)), check
        if (check) then
            write(*,'(a)') tmp
        else
            write(*,'(a)') colorize(trim(tmp), color_fg='red')
        end if

        deallocate(s)

    end subroutine run_test

    subroutine derivative_test_order(me, t, x, xdot)
        class(rk_class), intent(inout)      :: me
        real(wp), intent(in)                :: t
        real(wp), dimension(:), intent(in)  :: x
        real(wp), dimension(:), intent(out) :: xdot
        type(rklib_properties) :: p
        integer :: order

        xdot(1) =  -x(1) ! Analytical solution is x = exp(-t)

    end subroutine

end program rk_decay
