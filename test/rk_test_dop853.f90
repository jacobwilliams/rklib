!*****************************************************************************************
!>
!  Comparison of the [[rkdp85]] method (with and without dense output) with the
!  reference DOP853 implementation (https://github.com/jacobwilliams/dop853).
!
!  1. With identical steps (constant step size in both): the solution at each step,
!     and the dense output within each step, should agree to roundoff.
!  2. With the usual adaptive step size control (which is different in the two codes):
!     both should be accurate to about the requested tolerance.
!
!  The test problem is a circular two-body orbit (`mu=1`), which has the exact
!  solution `[cos(t), sin(t), -sin(t), cos(t)]`.

    module rk_test_dop853_types
    !! Integrator classes for [[rk_test_dop853]], extended to save the data for the test.

    use rklib_module, wp => rk_module_rk
    use dop853_module, only: dop853_class

    implicit none

    integer,parameter :: n = 4 !! number of state variables
    integer,parameter :: max_out = 5000 !! max number of saved points

    type :: recorder
        !! saved points from an integration
        integer :: n_out = 0
        real(wp),dimension(max_out) :: t = 0.0_wp
        real(wp),dimension(n,max_out) :: x = 0.0_wp
    contains
        procedure :: add => recorder_add
        procedure :: reset => recorder_reset
    end type recorder

    type,extends(rkdp85_class) :: rkdp85_test_class
        !! rklib integrator, with the data for the test
        integer :: fevals = 0
        logical :: interior = .false. !! save the dense output at interior points of each step
        type(recorder) :: steps !! points passed to `report`
        type(recorder) :: dense !! dense output at interior points
    end type rkdp85_test_class

    type,extends(dop853_class) :: dop853_test_class
        !! DOP853 integrator, with the data for the test
        logical :: interior = .false. !! save the dense output at interior points of each step
        real(wp) :: dt_out = 0.0_wp !! if > 0, save the dense output at this fixed interval instead
        integer :: k_out = 0 !! index of the last fixed interval point
        type(recorder) :: steps !! points passed to `solout` (or the fixed interval points)
        type(recorder) :: dense !! dense output at interior points
    end type dop853_test_class

    contains

!*****************************************************************************************
    subroutine recorder_add(me, t, x)
        class(recorder),intent(inout) :: me
        real(wp),intent(in) :: t
        real(wp),dimension(:),intent(in) :: x
        if (me%n_out >= max_out) error stop 'too many points'
        me%n_out = me%n_out + 1
        me%t(me%n_out) = t
        me%x(:,me%n_out) = x
    end subroutine recorder_add
!*****************************************************************************************

!*****************************************************************************************
    subroutine recorder_reset(me)
        class(recorder),intent(inout) :: me
        me%n_out = 0
    end subroutine recorder_reset
!*****************************************************************************************

    end module rk_test_dop853_types
!*****************************************************************************************

    program rk_test_dop853

    use rklib_module, wp => rk_module_rk
    use dop853_module, only: dop853_class, dop853_wp
    use rk_test_dop853_types

    implicit none

    integer,parameter :: n_interior = 9 !! interior points per step for the dense output comparison
    real(wp),parameter :: t0 = 0.0_wp
    real(wp),parameter :: tf = 10.0_wp
    real(wp),dimension(n),parameter :: x0 = [1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]


    if (dop853_wp /= wp) error stop 'rklib and dop853 must be compiled with the same real kind'

    call test_identical_steps()
    call test_adaptive()

    write(*,*) ''
    write(*,*) 'PASSED'

    contains

!*****************************************************************************************
    subroutine test_identical_steps()
        !! constant step size in both codes, so that the steps are identical

        real(wp),parameter :: h = 0.125_wp !! step size (exactly representable)
        real(wp),parameter :: tol_match = 1.0e-12_wp !! for agreement between the two codes

        type(rkdp85_test_class) :: rk
        type(dop853_test_class) :: dop
        type(stepsize_class) :: sz
        logical :: dense
        integer :: i, nfcn
        real(wp) :: diff_steps, diff_dense
        character(len=*),parameter :: fmt = '(1X,A,L2,A,I5,A,ES10.3,A,ES10.3,A,I6,A,I6)'

        write(*,*) ''
        write(*,*) 'Identical steps (h = 0.125):'

        do i = 1, 2
            dense = i==2

            ! rklib:
            call sz%initialize(fixed_step_mode=.true.)
            call rk%initialize(n=n, f=rk_f, report=rk_report, stepsize_method=sz, dense_output=dense)
            rk%interior = dense
            call run_rklib(rk, h)

            ! dop853: never reject or change the step:
            call run_dop853(dop, h, rtol=1.0e10_wp, hmax=h, interior=dense, nfcn=nfcn)

            if (rk%steps%n_out /= dop%steps%n_out) error stop 'different number of steps'
            if (any(rk%steps%t(1:rk%steps%n_out) /= dop%steps%t(1:dop%steps%n_out))) &
                error stop 'the steps are not identical'
            diff_steps = maxval(abs(rk%steps%x(:,1:rk%steps%n_out) - dop%steps%x(:,1:dop%steps%n_out)))
            diff_dense = 0.0_wp
            if (dense) then
                if (rk%dense%n_out /= dop%dense%n_out .or. rk%dense%n_out /= n_interior*(rk%steps%n_out-1)) &
                    error stop 'different number of dense output points'
                diff_dense = maxval(abs(rk%dense%x(:,1:rk%dense%n_out) - dop%dense%x(:,1:dop%dense%n_out)))
            end if

            write(*,fmt) 'dense:', dense, '  steps:', rk%steps%n_out-1, &
                         '  max diff (steps):', diff_steps, '  max diff (dense):', diff_dense, &
                         '  fevals rklib:', rk%fevals, '  dop853:', nfcn
            if (diff_steps > tol_match) error stop 'rkdp85 does not match dop853 at the steps'
            if (diff_dense > tol_match) error stop 'rkdp85 dense output does not match dop853'
        end do

    end subroutine test_identical_steps
!*****************************************************************************************

!*****************************************************************************************
    subroutine test_adaptive()
        !! the usual step size control in both codes, with output at a fixed interval

        real(wp),parameter :: tol = 1.0e-10_wp !! integration tolerance
        real(wp),parameter :: dt_out = 0.1_wp !! output interval
        real(wp),parameter :: tol_err = 1.0e-8_wp !! max allowed error vs. the exact solution

        type(rkdp85_test_class) :: rk
        type(dop853_test_class) :: dop
        integer :: nfcn, num_steps, naccpt
        real(wp) :: err_rk, err_dop
        character(len=*),parameter :: fmt = '(1X,A,A,I5,A,I5,A,I6,A,ES10.3)'

        write(*,*) ''
        write(*,*) 'Adaptive steps (tol = 1e-10), output every 0.1:'

        ! rklib, without dense output (points at the steps):
        call rk%initialize(n=n, f=rk_f, rtol=[tol], atol=[tol], report=rk_report)
        call run_rklib(rk, 0.0_wp)
        call rk%info(num_steps=num_steps)
        write(*,fmt) 'rklib  (steps) ', ' points:', rk%steps%n_out, '  steps:', num_steps, &
                     '  fevals:', rk%fevals, '  max error:', max_error(rk%steps)
        if (max_error(rk%steps) > tol_err) error stop 'rklib error too large'

        ! rklib, with dense output:
        call rk%initialize(n=n, f=rk_f, rtol=[tol], atol=[tol], report=rk_report, dense_output_step=dt_out)
        call run_rklib(rk, 0.0_wp)
        call rk%info(num_steps=num_steps)
        err_rk = max_error(rk%steps)
        write(*,fmt) 'rklib  (dense) ', ' points:', rk%steps%n_out, '  steps:', num_steps, &
                     '  fevals:', rk%fevals, '  max error:', err_rk

        ! dop853, with dense output:
        dop%dt_out = dt_out
        call run_dop853(dop, 0.0_wp, rtol=tol, hmax=0.0_wp, interior=.false., nfcn=nfcn, naccpt=naccpt)
        err_dop = max_error(dop%steps)
        write(*,fmt) 'dop853 (dense) ', ' points:', dop%steps%n_out, '  steps:', naccpt, &
                     '  fevals:', nfcn, '  max error:', err_dop

        if (rk%steps%n_out /= dop%steps%n_out) error stop 'different number of output points'
        if (any(abs(rk%steps%t(1:rk%steps%n_out) - dop%steps%t(1:dop%steps%n_out)) > 0.0_wp)) &
            error stop 'different output times'
        if (err_rk > tol_err) error stop 'rklib dense output error too large'
        if (err_dop > tol_err) error stop 'dop853 dense output error too large'

    end subroutine test_adaptive
!*****************************************************************************************

!*****************************************************************************************
    subroutine run_rklib(rk, h)
        !! integrate with rklib from `t0` to `tf`
        type(rkdp85_test_class),intent(inout) :: rk
        real(wp),intent(in) :: h !! initial step size (0 to compute it)
        real(wp),dimension(n) :: xf
        rk%fevals = 0
        call rk%steps%reset()
        call rk%dense%reset()
        call rk%integrate(t0, x0, h, tf, xf)
        if (rk%failed()) error stop 'rklib integration failed'
    end subroutine run_rklib
!*****************************************************************************************

!*****************************************************************************************
    subroutine run_dop853(dop, h, rtol, hmax, interior, nfcn, naccpt)
        !! integrate with dop853 from `t0` to `tf` (with dense output)
        type(dop853_test_class),intent(inout) :: dop
        real(wp),intent(in) :: h !! initial step size (0 to compute it)
        real(wp),intent(in) :: rtol !! relative and absolute tolerance
        real(wp),intent(in) :: hmax !! max step size (0 for no limit)
        logical,intent(in) :: interior !! save the dense output at interior points of each step
        integer,intent(out) :: nfcn !! number of function evaluations
        integer,intent(out),optional :: naccpt !! number of accepted steps

        real(wp) :: x
        real(wp),dimension(n) :: y
        real(wp) :: dt_out
        integer :: idid
        logical :: status_ok

        dt_out = dop%dt_out ! (reset by initialize)
        call dop%initialize(n=n, fcn=dop_f, solout=dop_solout, nstiff=-1, &
                            hinitial=h, hmax=hmax, icomp=[1,2,3,4], status_ok=status_ok)
        if (.not. status_ok) error stop 'dop853 initialization failed'
        dop%interior = interior
        dop%dt_out = dt_out
        dop%k_out = 0
        call dop%steps%reset()
        call dop%dense%reset()

        x = t0
        y = x0
        call dop%integrate(x, y, tf, [rtol], [rtol], iout=2, idid=idid)
        if (idid /= 1) error stop 'dop853 integration failed'
        call dop%info(nfcn=nfcn, naccpt=naccpt)

    end subroutine run_dop853
!*****************************************************************************************

!*****************************************************************************************
    function max_error(rec) result(err)
        !! max error of the saved points vs. the exact solution
        type(recorder),intent(in) :: rec
        real(wp) :: err
        integer :: i
        err = 0.0_wp
        do i = 1, rec%n_out
            err = max(err, maxval(abs(rec%x(:,i) - exact(rec%t(i)))))
        end do
    end function max_error
!*****************************************************************************************

!*****************************************************************************************
    pure function exact(t) result(x)
        !! exact solution of the circular orbit
        real(wp),intent(in) :: t
        real(wp),dimension(n) :: x
        x = [cos(t), sin(t), -sin(t), cos(t)]
    end function exact
!*****************************************************************************************

!*****************************************************************************************
    pure function twobody(x) result(xdot)
        !! two-body equations of motion (`mu=1`)
        real(wp),dimension(:),intent(in) :: x
        real(wp),dimension(n) :: xdot
        real(wp) :: r3
        r3 = norm2(x(1:2))**3
        xdot = [x(3), x(4), -x(1)/r3, -x(2)/r3]
    end function twobody
!*****************************************************************************************

!*****************************************************************************************
    subroutine rk_f(me,t,x,xdot)
        !! rklib derivative function
        class(rk_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot
        xdot = twobody(x)
        select type (me)
        class is (rkdp85_test_class)
            me%fevals = me%fevals + 1
        end select
    end subroutine rk_f
!*****************************************************************************************

!*****************************************************************************************
    subroutine rk_report(me,t,x)
        !! rklib report function: save the point, and the dense output within the step
        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x

        logical :: available
        real(wp) :: t_prev, ti
        real(wp),dimension(n) :: xi
        integer :: k

        select type (me)
        class is (rkdp85_test_class)
            if (me%interior) then
                call me%dense_interval(available, t0=t_prev)
                if (available) then
                    do k = 1, n_interior
                        ti = t_prev + (t-t_prev)*real(k,wp)/real(n_interior+1,wp)
                        call me%interpolate(ti, xi)
                        call me%dense%add(ti, xi)
                    end do
                end if
            end if
            call me%steps%add(t, x)
        end select

    end subroutine rk_report
!*****************************************************************************************

!*****************************************************************************************
    subroutine dop_f(me,x,y,f)
        !! dop853 derivative function
        class(dop853_class),intent(inout) :: me
        real(wp),intent(in)               :: x
        real(wp),dimension(:),intent(in)  :: y
        real(wp),dimension(:),intent(out) :: f
        f = twobody(y)
    end subroutine dop_f
!*****************************************************************************************

!*****************************************************************************************
    subroutine dop_solout(me,nr,xold,x,y,irtrn,xout)
        !! dop853 output function: save the point (or the fixed interval
        !! points), and the dense output within the step
        class(dop853_class),intent(inout) :: me
        integer,intent(in)                :: nr
        real(wp),intent(in)               :: xold
        real(wp),intent(in)               :: x
        real(wp),dimension(:),intent(in)  :: y
        integer,intent(inout)             :: irtrn
        real(wp),intent(out)              :: xout

        real(wp) :: ti
        real(wp),dimension(n) :: yi
        integer :: i, k

        xout = 0.0_wp ! not used (iout=2)

        select type (me)
        class is (dop853_test_class)

            if (me%dt_out > 0.0_wp) then
                ! points at the fixed interval (like rklib's `dense_output_step`):
                if (nr==1) then
                    call me%steps%add(x, y)
                else
                    do
                        ti = t0 + real(me%k_out+1,wp)*me%dt_out
                        if (x - ti <= 2*spacing(abs(x))) exit ! (points at the end of the step are done next time)
                        do i = 1, n
                            yi(i) = me%contd8(i, ti)
                        end do
                        call me%steps%add(ti, yi)
                        me%k_out = me%k_out + 1
                    end do
                    if (x == tf) call me%steps%add(x, y) ! final point
                end if
            else
                if (me%interior .and. nr > 1) then
                    do k = 1, n_interior
                        ti = xold + (x-xold)*real(k,wp)/real(n_interior+1,wp)
                        do i = 1, n
                            yi(i) = me%contd8(i, ti)
                        end do
                        call me%dense%add(ti, yi)
                    end do
                end if
                call me%steps%add(x, y)
            end if

        end select

    end subroutine dop_solout
!*****************************************************************************************


    end program rk_test_dop853
!*****************************************************************************************
