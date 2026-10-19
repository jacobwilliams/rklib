program rk_test_nonautonomous
!! Unit test of the [[rk_module]].
!! This test complements [[rk_decay]] (an autonomous equation) and [[rk_test_order]]:
!! it checks that the times of the stages are correct for all the methods,
!! using the non-autonomous equation \dot x = cos(t), with x(1) = sin(1),
!! integrated from t = 1 to 3 (starting at t = 0 would hide errors such as
!! using `t*h` instead of `t` for a stage).

    use rklib_module, wp => rk_module_rk
    implicit none

    real(wp), parameter :: t0 = 1.0_wp, tf = 3.0_wp

    class(rk_class), allocatable :: s
    logical :: all_ok

    all_ok = .true.
    write(*, '(a10, a10, a24, a9)') 'Method', "Variable", "error", "Success"
    write(*, *) "----------------------------------------------------"

#include "rklib_allocate_and_test.inc"

    if (.not. all_ok) error stop 'some methods failed'
    write(*,*) ''
    write(*,*) 'PASSED'

contains

    subroutine run_test()

        real(wp), parameter :: tol = 1.0e-10_wp, dt = 1.0e-3_wp
        real(wp) :: xf(1), err
        type(rklib_properties) :: p
        logical :: isvariable, check

        p = s%properties()

        select type (s)
        class is (rk_fixed_step_class)
            isvariable = .false.
            call s%initialize(n=1, f=deriv)
            call s%integrate(t0, [sin(t0)], dt, tf, xf)
        class is (rk_variable_step_class)
            isvariable = .true.
            call s%initialize(n=1, f=deriv, rtol=[tol], atol=[tol])
            call s%integrate(t0, [sin(t0)], dt, tf, xf)
        end select

        err = abs(xf(1) - sin(tf))
        ! (the euler method is first order, so its error is ~dt)
        check = err <= merge(1.0e-2_wp, 1.0e-6_wp, p%order == 1)
        write(*, '(a10, l10, es24.15e2, l9)') p%short_name, isvariable, err, check
        if (.not. check) all_ok = .false.

        deallocate(s)

    end subroutine run_test

    subroutine deriv(me, t, x, xdot)
        class(rk_class), intent(inout)      :: me
        real(wp), intent(in)                :: t
        real(wp), dimension(:), intent(in)  :: x
        real(wp), dimension(:), intent(out) :: xdot
        xdot(1) = cos(t)
    end subroutine deriv

end program rk_test_nonautonomous
