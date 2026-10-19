!*****************************************************************************************
!>
!  Comparison of rklib with the FLINT library (https://github.com/princemahajan/FLINT),
!  which has four of the same methods:
!
!   FLINT      | rklib  | interpolant order (FLINT / rklib)
!   ---------- | ------ | ---------------------------------
!   DOP54      | rkdp54 | 4 / 4 (the same)
!   DOP853     | rkdp85 | 7 / 7 (the same)
!   Verner65E  | rkv65e | 5 / 6
!   Verner98R  | rkv98r | 8 / 9
!
!  1. With identical (constant) steps: the steps of the two libraries should agree to roundoff,
!     and so should the dense output for the methods with the same interpolant. For the
!     Verner methods (where rklib uses the higher order interpolant), the error of rklib's
!     dense output should not be larger than FLINT's.
!  2. With each library's own step size control (which differ), and output at a fixed
!     interval: the results are printed for comparison (with only a loose accuracy check).
!  3. With rklib using the step size control options of Hairer's codes (which FLINT uses,
!     see [[stepsize_class]]): rklib should take (almost) the same steps as FLINT.
!
!  The test problem is a circular two-body orbit (`mu=1`), which has the exact
!  solution `[cos(t), sin(t), -sin(t), cos(t)]`.

    module rk_test_flint_module

    !! FLINT problem definition for [[rk_test_flint]].

    use FLINT, FWP => FLINT_WP

    implicit none

    type,extends(DiffEqSys) :: orbit_sys
    contains
        procedure :: F => orbit_f
    end type orbit_sys

    contains

    function orbit_f(me, X, Y, Params)
        class(orbit_sys),intent(inout) :: me
        real(FWP),intent(in) :: X
        real(FWP),dimension(me%n),intent(in) :: Y
        real(FWP),dimension(:),intent(in),optional :: Params
        real(FWP),dimension(size(Y)) :: orbit_f
        real(FWP) :: r3
        r3 = hypot(Y(1),Y(2))**3
        orbit_f = [Y(3), Y(4), -Y(1)/r3, -Y(2)/r3]
    end function orbit_f

    end module rk_test_flint_module
!*****************************************************************************************

!*****************************************************************************************
    program rk_test_flint

    use rklib_module, wp => rk_module_rk
    use FLINT, FWP => FLINT_WP
    use rk_test_flint_module
    use iso_fortran_env, only: int64

    implicit none

    integer,parameter :: n = 4 !! number of state variables
    integer,parameter :: max_out = 2000 !! max number of output points
    integer,parameter :: flint_max_steps = 5000 !! max number of steps for FLINT
    integer,parameter :: reps = 20 !! number of runs for the timing
    real(wp),parameter :: t0 = 0.0_wp, tf = 10.0_wp
    real(wp),parameter :: h = 0.125_wp !! step size for the identical steps test (exactly representable)
    real(wp),parameter :: dt_out = 0.01_wp !! output interval for the adaptive test
    real(wp),dimension(n),parameter :: x0 = [1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]
    real(wp),dimension(3),parameter :: tols = [1.0e-8_wp, 1.0e-10_wp, 1.0e-12_wp]
    integer,parameter :: n_steps = nint((tf-t0)/h) !! number of steps for the identical steps test
    integer,parameter :: n_grid = nint((tf-t0)/dt_out) + 1 !! number of output points for the adaptive test
    character(len=9),dimension(4),parameter :: flint_names = ['DOP54    ','DOP853   ','Verner65E','Verner98R']

    class(rk_variable_step_class),allocatable :: s
    type(ERK_class) :: erk
    type(orbit_sys),target :: diffeq
    integer :: im, i, k, fevals, n_out
    real(wp),dimension(9*n_steps) :: t_int !! interior points of the steps
    real(wp),dimension(n,9*n_steps) :: x_int_rk, x_int_fl
    real(wp),dimension(n_grid) :: t_grid
    real(wp),dimension(n,n_grid) :: x_grid_fl
    real(wp),dimension(max_out) :: t_out
    real(wp),dimension(n,max_out) :: x_out

    if (FWP /= wp) then
        ! [rklib uses the REAL32/REAL128 preprocessor flags and FLINT uses _REAL32/_REAL128]
        write(*,*) 'Skipping the FLINT comparison: rklib and FLINT are compiled with different real kinds'
        stop
    end if

    diffeq%n = n
    diffeq%m = 0
    do k = 0, n_steps-1
        do i = 1, 9
            t_int(9*k+i) = (t0 + k*h) + h*real(i,wp)/10.0_wp
        end do
    end do
    t_grid = [(t0 + i*dt_out, i = 0, n_grid-1)]

    do im = 1, 4
        select case (im)
        case(1); allocate(rkdp54_class :: s)
        case(2); allocate(rkdp85_class :: s)
        case(3); allocate(rkv65e_class :: s)
        case(4); allocate(rkv98r_class :: s)
        end select
        call identical_steps()
        call adaptive()
        deallocate(s)
    end do

    call hairer_settings()

    write(*,*) ''
    write(*,*) 'PASSED'

    contains

!*****************************************************************************************
    function method_id(i) result(id)
        !! FLINT method for case `i`
        integer,intent(in) :: i
        integer :: id
        select case (i)
        case(1); id = ERK_DOP54
        case(2); id = ERK_DOP853
        case(3); id = ERK_VERNER65E
        case(4); id = ERK_VERNER98R
        end select
    end function method_id
!*****************************************************************************************

!*****************************************************************************************
    subroutine identical_steps()
        !! constant step size in both libraries

        type(stepsize_class) :: sz
        type(rklib_properties) :: p
        real(wp),dimension(n) :: xf, yf
        real(wp) :: xf_flint, step, diff_step, diff_dense, err_rk, err_fl
        integer(kind(FLINT_SUCCESS)) :: status
        logical :: same_interpolant

        p = s%properties()
        same_interpolant = im <= 2

        ! rklib:
        call sz%initialize(fixed_step_mode=.true.)
        call s%initialize(n=n, f=f, stepsize_method=sz, report=report_interior, dense_output=.true.)
        n_out = 0
        call s%integrate(t0, x0, h, tf, xf)
        if (s%failed() .or. n_out /= size(t_int)) error stop 'rklib integration failed'

        ! FLINT (the tolerances are not used with a constant step size):
        call erk%Init(diffeq, flint_max_steps, Method=method_id(im), ATol=[1.0_wp], RTol=[1.0_wp], InterpOn=.true.)
        xf_flint = tf
        step = h
        call erk%Integrate(t0, x0, xf_flint, yf, StepSz=step, UseConstStepSz=.true.)
        call erk%Interpolate(t_int, x_int_fl)
        call erk%Info(status)
        if (status /= FLINT_SUCCESS) error stop 'FLINT integration failed'

        diff_step = maxval(abs(xf - yf))
        diff_dense = maxval(abs(x_int_rk - x_int_fl))
        err_rk = max_error(t_int, x_int_rk)
        err_fl = max_error(t_int, x_int_fl)
        write(*,*) ''
        write(*,'(1X,A,A,A)') trim(flint_names(im)), ' vs ', p%short_name
        write(*,'(3X,A,ES10.3)') 'identical steps: max diff (final state):      ', diff_step
        write(*,'(3X,A,ES10.3)') 'identical steps: max diff (dense output):     ', diff_dense
        write(*,'(3X,A,2ES11.3)') 'identical steps: dense output error rklib/FLINT:', err_rk, err_fl

        ! the steps should agree to roundoff (rkv65e has large coefficients, and so a larger roundoff):
        if (diff_step > merge(1.0e-10_wp, 1.0e-12_wp, im==3)) error stop 'the steps do not agree with FLINT'
        if (same_interpolant) then
            if (diff_dense > 1.0e-12_wp) error stop 'the dense output does not agree with FLINT'
        else
            if (err_rk > err_fl) error stop 'the dense output is less accurate than FLINT''s'
        end if

    end subroutine identical_steps
!*****************************************************************************************

!*****************************************************************************************
    subroutine adaptive()
        !! each library's own step size control, with output at a fixed interval

        real(wp),dimension(n) :: xf, yf
        real(wp) :: xf_flint, step, err_rk, err_fl, time_rk, time_fl
        integer :: i, r, num_steps, flint_steps, flint_fevals
        integer(kind(FLINT_SUCCESS)) :: status
        integer(int64) :: c1, c2, rate

        call system_clock(count_rate=rate)
        write(*,'(3X,A)') 'adaptive (rklib/FLINT):   tol     steps         fevals        output error            final error             time (us)'
        do i = 1, size(tols)

            ! rklib:
            call s%initialize(n=n, f=f, rtol=[tols(i)], atol=[tols(i)], report=report_grid, dense_output_step=dt_out)
            time_rk = huge(1.0_wp)
            do r = 1, reps
                fevals = 0
                n_out = 0
                call system_clock(c1)
                call s%integrate(t0, x0, 0.0_wp, tf, xf)
                call system_clock(c2)
                time_rk = min(time_rk, real(c2-c1,wp)/real(rate,wp))
            end do
            if (s%failed() .or. n_out /= n_grid) error stop 'rklib integration failed'
            call s%info(num_steps=num_steps)
            err_rk = max_error(t_out(1:n_out), x_out(:,1:n_out))

            ! FLINT (interpolation after the integration):
            ! [the interpolation storage is kept for the next run: if it is deallocated,
            !  FLINT must be initialized again before the next integration]
            call erk%Init(diffeq, flint_max_steps, Method=method_id(im), ATol=[tols(i)], RTol=[tols(i)], InterpOn=.true.)
            time_fl = huge(1.0_wp)
            do r = 1, reps
                xf_flint = tf
                step = 0.0_wp
                call system_clock(c1)
                call erk%Integrate(t0, x0, xf_flint, yf, StepSz=step)
                call erk%Interpolate(t_grid, x_grid_fl, .true.)
                call system_clock(c2)
                time_fl = min(time_fl, real(c2-c1,wp)/real(rate,wp))
            end do
            call erk%Info(status, nAccept=flint_steps, nFCalls=flint_fevals)
            if (status /= FLINT_SUCCESS) error stop 'FLINT integration failed'
            err_fl = max_error(t_grid, x_grid_fl)

            write(*,'(3X,A,ES8.1,2(I7,A,I6),4(ES11.3,A),2(F8.1,A))') '                     ', tols(i), &
                num_steps, '/', flint_steps, fevals, '/', flint_fevals, &
                err_rk, ' /', err_fl, ' ', maxval(abs(xf-exact(tf))), ' /', maxval(abs(yf-exact(tf))), ' ', &
                1.0e6_wp*time_rk, ' /', 1.0e6_wp*time_fl, ''

            ! loose accuracy check (the two libraries use different error norms and step size control):
            if (err_rk > max(1.0e3_wp*tols(i), 1.0e-9_wp)) error stop 'rklib output error too large'

        end do

    end subroutine adaptive
!*****************************************************************************************

!*****************************************************************************************
    subroutine hairer_settings()
        !! rklib with the step size control options of Hairer's `DOPRI5` and `DOP853` codes
        !! (which FLINT uses) should take (almost) the same steps as FLINT's DOP54 and DOP853.
        !! [the treatment of the last step is slightly different, so the number of
        !!  steps may differ by one, which also changes the final error slightly]

        real(wp),dimension(5),parameter :: tol = [1.0e-6_wp, 1.0e-8_wp, 1.0e-10_wp, 1.0e-12_wp, 1.0e-14_wp]

        type(rkdp54_class) :: s54
        type(rkdp85_class) :: s85
        type(stepsize_class) :: sz
        real(wp),dimension(n) :: xf, yf
        real(wp) :: xf_flint, step, err_rk, err_fl
        integer :: i, j, num_steps, num_rejected, flint_steps, flint_rejected
        integer(kind(FLINT_SUCCESS)) :: status

        write(*,*) ''
        write(*,*) 'rklib with the step size control options of DOPRI5/DOP853 (no output):'
        do j = 1, 2
            do i = 1, size(tol)
                if (j==1) then
                    ! DOPRI5: with the PI (Lund) stabilization
                    call sz%initialize(norm=rms_func, rtol_scaling=2, p_exponent_offset=0, &
                                       hfactor_accept=6.0_wp, hfactor_reject=0.333_wp, safety_factor=0.9_wp, &
                                       no_increase_after_reject=.true., pi_beta=0.04_wp, pi_beta_mult=0.75_wp, &
                                       hmin=0.0_wp)
                    call s54%initialize(n=n, f=f, rtol=[tol(i)], atol=[tol(i)], stepsize_method=sz, hinit_method=2)
                    call s54%integrate(t0, x0, 0.0_wp, tf, xf)
                    call s54%info(num_steps=num_steps, num_rejected_steps=num_rejected)
                    if (s54%failed()) error stop 'rklib integration failed'
                    call erk%Init(diffeq, flint_max_steps, Method=ERK_DOP54, ATol=[tol(i)], RTol=[tol(i)])
                else
                    ! DOP853: with the combined error estimate using the norms
                    call sz%initialize(norm=rms_func, rtol_scaling=2, p_exponent_offset=0, &
                                       hfactor_accept=6.0_wp, hfactor_reject=0.333_wp, safety_factor=0.9_wp, &
                                       no_increase_after_reject=.true., error_combination=2, hmin=0.0_wp)
                    call s85%initialize(n=n, f=f, rtol=[tol(i)], atol=[tol(i)], stepsize_method=sz, hinit_method=2)
                    call s85%integrate(t0, x0, 0.0_wp, tf, xf)
                    call s85%info(num_steps=num_steps, num_rejected_steps=num_rejected)
                    if (s85%failed()) error stop 'rklib integration failed'
                    call erk%Init(diffeq, flint_max_steps, Method=ERK_DOP853, ATol=[tol(i)], RTol=[tol(i)])
                end if
                xf_flint = tf
                step = 0.0_wp
                call erk%Integrate(t0, x0, xf_flint, yf, StepSz=step)
                call erk%Info(status, nAccept=flint_steps, nReject=flint_rejected)
                if (status /= FLINT_SUCCESS) error stop 'FLINT integration failed'
                err_rk = maxval(abs(xf - exact(tf)))
                err_fl = maxval(abs(yf - exact(tf)))
                write(*,'(3X,A,A,ES8.1,A,2I6,A,2I6,A,2ES11.3)') flint_names(j), ' tol', tol(i), &
                    '  accepted/rejected rklib:', num_steps, num_rejected, '  FLINT:', flint_steps, flint_rejected, &
                    '  error rklib/FLINT:', err_rk, err_fl
                if (abs(num_steps - flint_steps) > 1 .or. num_rejected /= flint_rejected) &
                    error stop 'not the same steps as FLINT'
                if (abs(err_rk - err_fl) > max(0.25_wp*err_fl, 100*epsilon(1.0_wp))) &
                    error stop 'not the same error as FLINT'
            end do
        end do

    end subroutine hairer_settings
!*****************************************************************************************

!*****************************************************************************************
    subroutine f(me,t,x,xdot)
        !! rklib derivative function
        class(rk_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot
        real(wp) :: r3
        r3 = hypot(x(1),x(2))**3
        xdot = [x(3), x(4), -x(1)/r3, -x(2)/r3]
        fevals = fevals + 1
    end subroutine f
!*****************************************************************************************

!*****************************************************************************************
    subroutine report_interior(me,t,x)
        !! rklib report function: dense output at the 9 interior points of each step
        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x
        logical :: available
        integer :: j
        call me%dense_interval(available)
        if (.not. available) return
        do j = 1, 9
            n_out = n_out + 1
            call me%interpolate(t_int(n_out), x_int_rk(:,n_out))
        end do
    end subroutine report_interior
!*****************************************************************************************

!*****************************************************************************************
    subroutine report_grid(me,t,x)
        !! rklib report function: save the points
        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x
        if (n_out >= max_out) error stop 'too many output points'
        n_out = n_out + 1
        t_out(n_out) = t
        x_out(:,n_out) = x
    end subroutine report_grid
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
    function max_error(ts, xs) result(err)
        !! max error of the points vs. the exact solution
        real(wp),dimension(:),intent(in) :: ts
        real(wp),dimension(:,:),intent(in) :: xs
        real(wp) :: err
        integer :: j
        err = 0.0_wp
        do j = 1, size(ts)
            err = max(err, maxval(abs(xs(:,j) - exact(ts(j)))))
        end do
    end function max_error
!*****************************************************************************************

    end program rk_test_flint
!*****************************************************************************************
