program rk_test_order
!! Unit test of the [[rk_module]].
!!
!! An ODE integration method of order \(n\) has a local truncation error of order \(n+1\).
!! Thus, a method of order \(n\) should be able to integrate exactly (within machine precision)
!! a polynomial function of the same order, i.e., \(x(t)=t^n\), over a given interval with a
!! single step.
!!
!! This test checks if this behavior is fullfilled by integrating \(\dot{x}(t)=n t^{n-1}\) over
!! \([0,1]\), which should give \(x(1)=1\) for all \(n\).
    use rklib_module, wp => rk_module_rk
    implicit none

    real(wp), parameter :: zero = 0.0_wp, one = 1.0_wp
    real(wp), parameter :: eps = epsilon(one)

    class(rk_class), allocatable :: s

    ! Test all methods
    write(*, '(a10, a10, a24, a9)') 'Method', "Variable", "xf", "Success"
    write(*, *) "----------------------------------------------------"

#include "rklib_allocate_and_test.inc"

contains

    subroutine run_test()
        !! Order test
        use face, only: colorize

        integer, parameter :: n = 1
        real(wp), parameter :: t0 = zero, tf = one, dt = (tf - t0), x0(n) = zero
        real(wp) :: xf(n)
        type(rklib_properties) :: p
        type(stepsize_class) :: sz
        character(:), allocatable :: method
        logical :: check, isvariable
        character(len=100) :: tmp

        p = s%properties()
        method = p%short_name

        ! step size constructor
        call sz%initialize(fixed_step_mode = .true.)

        select type (s)
        class is (rk_fixed_step_class)
            isvariable = .false.
            call s%initialize(n=n, f=derivative_test_order)
            call s%integrate(t0, x0, dt, tf, xf)

        class is (rk_variable_step_class)
            isvariable = .true.
            call s%initialize(n=n, f=derivative_test_order, rtol=[one], atol=[one], &
                              stepsize_method=sz , hinit_method=2)
            call s%integrate(t0, x0, dt, tf, xf)

        end select

        ! 20*eps is ok for real32 and real64, but not always for real128, likely because some
        ! methods only have constants with double precision
        check = abs(xf(1) - one) <= 10000*eps
        write(tmp, '(a10, l10, es24.15e2, l9)') method, isvariable, xf, check
        if (check) then
            write(*,'(a)') tmp
        else
            write(*,'(a)') colorize(trim(tmp), color_fg='red')
        end if

        deallocate(s)

    end subroutine run_test

    subroutine derivative_test_order(me, t, x, xdot)
    !! Derivative for order test.
    !! $$
    !! \dot{x}(t)=n t^{n-1}
    !! $$
    !! where \(n \ge 1\) is the global order of the method.
        class(rk_class), intent(inout)      :: me
        real(wp), intent(in)                :: t
        real(wp), dimension(:), intent(in)  :: x
        real(wp), dimension(:), intent(out) :: xdot
        type(rklib_properties) :: p
        integer :: order

        p = me%properties()
        order = p%order
        xdot(1) = order*t**(order - 1)

    end subroutine

end program rk_test_order
