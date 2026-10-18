!*****************************************************************************************
!>
!  Unit test of the dense output of the variable-step methods.

    program rk_test_dense

    use rklib_module, wp => rk_module_rk
    use iso_fortran_env

    implicit none

    real(wp),parameter :: pi = acos(-1.0_wp)
    real(wp),parameter :: t0_step = 0.5_wp !! start time for the order test
                                           !! (not 0, where the solution is symmetric)
    integer,parameter :: n_interp = 10 !! number of interpolated points per step

    class(rk_variable_step_class),allocatable :: s
    real(wp) :: max_err, max_endpoint_err
    real(wp),dimension(1) :: y0, yf, y
    integer :: fevals
    integer :: f_interpolate_status !! status from calling `interpolate` in [[nonlinear_interpolate]]

    ! points reported by [[report_grid]]:
    integer,parameter :: max_rep = 2000
    integer :: n_rep
    real(wp),dimension(max_rep) :: t_rep
    real(wp),dimension(2,max_rep) :: x_rep

    !                      method           | order | extra fevals | h for order test
    allocate(rkdp54_class :: s); call run_tests(4,      0,             0.0125_wp)
    allocate(rkdp85_class :: s); call run_tests(7,      4,             0.1_wp)
    allocate(rkv65e_class :: s); call run_tests(6,      3,        0.1_wp)
    allocate(rkv65r_class :: s); call run_tests(6,      3,        0.1_wp)
    allocate(rkv76e_class :: s); call run_tests(7,      6,        0.1_wp)
    allocate(rkv76r_class :: s); call run_tests(7,      6,        0.1_wp)
    allocate(rkv87e_class :: s); call run_tests(8,      8,        0.2_wp)
    allocate(rkv87r_class :: s); call run_tests(8,      8,        0.1_wp)
    allocate(rkv98e_class :: s); call run_tests(9,      10,        0.2_wp)
    allocate(rkv98r_class :: s); call run_tests(9,      10,        0.2_wp)

    call test_no_dense()
    call test_interpolate_in_f()

    write(*,*) 'PASSED'

    contains

    !*********************************************************
        subroutine run_tests(expected_order, extra_fevals, h)
        !! run all the tests for the method `s`

        integer,intent(in) :: expected_order !! expected order of the interpolant
        integer,intent(in) :: extra_fevals !! number of extra function evaluations
                                           !! per step required for the dense output
        real(wp),intent(in) :: h !! step size for the order test (should be
                                 !! in the asymptotic range for the method)

        type(rklib_properties) :: p
        type(stepsize_class) :: sz
        real(wp) :: err1, err2, order, tf_actual, gf
        real(wp),dimension(2) :: xf
        integer :: istatus, fevals_no_dense, num_steps

        p = s%properties()
        write(*,*) ''
        write(*,*) 'Method: ', p%long_name, ' dense: ', p%dense
        if (.not. p%dense) error stop 'method should have dense output'

        !-----------------------------------------------------
        ! function evaluations without using the dense output:
        fevals = 0
        call s%initialize(n=2, f=oscillator, rtol=[1.0e-10_wp], atol=[1.0e-10_wp])
        call s%integrate(0.0_wp, [1.0_wp, 0.0_wp], 0.0_wp, 10.0_wp, xf)
        if (s%failed()) error stop 'integration failed'
        fevals_no_dense = fevals
        call s%info(num_steps=num_steps)

        !-----------------------------------------------------
        ! interpolate within every step, from the report function:
        fevals = 0
        call s%initialize(n=2, f=oscillator, rtol=[1.0e-10_wp], atol=[1.0e-10_wp], report=report, dense_output=.true.)
        max_err = 0.0_wp
        max_endpoint_err = 0.0_wp
        call s%integrate(0.0_wp, [1.0_wp, 0.0_wp], 0.0_wp, 10.0_wp, xf)
        if (s%failed()) error stop 'integration failed'
        write(*,*) 'max dense output error:       ', max_err
        write(*,*) 'max endpoint mismatch:        ', max_endpoint_err
        write(*,*) 'error at final point:         ', maxval(abs(xf - exact_oscillator(10.0_wp)))
        write(*,*) 'function evaluations:         ', fevals_no_dense, '(without dense output)', &
                                                     fevals, '(with dense output)'
        if (max_err > 1.0e-8_wp) error stop 'dense output error too large'
        if (max_endpoint_err > 100*epsilon(1.0_wp)) error stop 'dense output does not match step endpoints'
        if (fevals - fevals_no_dense /= extra_fevals*num_steps) &
            error stop 'unexpected number of extra function evaluations'

        !-----------------------------------------------------
        ! same thing, integrating backwards:
        max_err = 0.0_wp
        max_endpoint_err = 0.0_wp
        call s%integrate(10.0_wp, exact_oscillator(10.0_wp), 0.0_wp, 0.0_wp, xf)
        if (s%failed()) error stop 'backward integration failed'
        write(*,*) 'backward: max dense output error: ', max_err
        write(*,*) 'backward: max endpoint mismatch:  ', max_endpoint_err
        if (max_err > 1.0e-8_wp) error stop 'backward dense output error too large'
        if (max_endpoint_err > 100*epsilon(1.0_wp)) &
            error stop 'backward dense output does not match step endpoints'

        !-----------------------------------------------------
        ! reporting at a fixed step size using the dense output:
        call test_dense_output_step(fevals_no_dense)

        !-----------------------------------------------------
        ! order of the interpolant: single step of a nonlinear, non-autonomous problem
        ! (local error should be O(h^(p+1)) for a p-th order interpolant)
        call sz%initialize(fixed_step_mode=.true.)
        call s%initialize(n=1, f=nonlinear, stepsize_method=sz, dense_output=.true.)
        y0 = [exact_nonlinear(t0_step)]
        err1 = single_step_error(h)
        err2 = single_step_error(h/2.0_wp)
        order = log(err1/err2)/log(2.0_wp) - 1.0_wp
        write(*,*) 'interpolant errors (h, h/2):  ', err1, err2
        write(*,*) 'observed interpolant order:   ', order
        ! (the tolerance allows for the step size not being exactly in the asymptotic range,
        !  but would detect an interpolant that is an order lower than expected)
        if (abs(order-expected_order) > 0.5_wp) error stop 'unexpected order of the interpolant'

        !-----------------------------------------------------
        ! out of range of the last step:
        call s%interpolate(t0_step + h, y)
        call s%status(istatus)
        if (istatus /= RKLIB_ERROR_DENSE_OUT_OF_RANGE) error stop 'expected out of range error'

        !-----------------------------------------------------
        ! event location using the dense output:
        ! x(1) = cos(t) = 0.5  -->  t = pi/3
        ! (without and with the dense output enabled)
        fevals = 0
        call s%initialize(n=2, f=oscillator, g=event, rtol=[1.0e-12_wp], atol=[1.0e-12_wp])
        call s%integrate_to_event(0.0_wp, [1.0_wp, 0.0_wp], 0.0_wp, 10.0_wp, 1.0e-12_wp, tf_actual, xf, gf)
        if (s%failed()) error stop 'event integration failed'
        write(*,*) 'event time error:             ', tf_actual - pi/3.0_wp, '(without dense output)'
        write(*,*) 'event function evaluations:   ', fevals, '(without dense output)'
        if (abs(tf_actual - pi/3.0_wp) > 1.0e-10_wp) error stop 'event time error too large'
        fevals = 0
        call s%initialize(n=2, f=oscillator, g=event, rtol=[1.0e-12_wp], atol=[1.0e-12_wp], dense_output=.true.)
        call s%integrate_to_event(0.0_wp, [1.0_wp, 0.0_wp], 0.0_wp, 10.0_wp, 1.0e-12_wp, tf_actual, xf, gf)
        if (s%failed()) error stop 'event integration failed'
        write(*,*) 'event time error:             ', tf_actual - pi/3.0_wp, '(with dense output)'
        write(*,*) 'event function evaluations:   ', fevals, '(with dense output)'
        if (abs(tf_actual - pi/3.0_wp) > 1.0e-10_wp) error stop 'event time error too large'

        !-----------------------------------------------------
        call test_dense_options()

        deallocate(s)

        end subroutine run_tests
    !*********************************************************

    !*********************************************************
        subroutine test_dense_options()
        !! test the `dense_output` and `dense_output_step` options for the method `s`

        integer :: istatus
        logical :: available
        real(wp),dimension(2) :: x, xf
        real(wp) :: t0_last

        ! not enabled (the default):
        call s%initialize(n=2, f=oscillator)
        call s%integrate(0.0_wp, [1.0_wp, 0.0_wp], 0.1_wp, 1.0_wp, xf)
        call s%interpolate(0.95_wp, x)
        call s%status(istatus)
        if (istatus /= RKLIB_ERROR_DENSE_NOT_ENABLED) error stop 'expected dense not enabled error from interpolate'
        call s%initialize(n=2, f=oscillator)
        call s%integrate(0.0_wp, [1.0_wp, 0.0_wp], 0.1_wp, 1.0_wp, xf)
        call s%dense_interval(available)
        call s%status(istatus)
        if (available .or. istatus /= RKLIB_ERROR_DENSE_NOT_ENABLED) &
            error stop 'expected dense not enabled error from dense_interval'

        ! explicitly disabled, but with dense_output_step:
        call s%initialize(n=2, f=oscillator, dense_output=.false., dense_output_step=0.1_wp)
        call s%status(istatus)
        if (istatus /= RKLIB_ERROR_DENSE_STEP_DISABLED) error stop 'expected dense step disabled error'

        ! explicitly enabled, with dense_output_step:
        call s%initialize(n=2, f=oscillator, dense_output=.true., dense_output_step=0.1_wp)
        if (s%failed()) error stop 'dense_output=.true. with dense_output_step should be ok'

        ! dense_output_step enables it automatically:
        call s%initialize(n=2, f=oscillator, dense_output_step=0.1_wp)
        call s%integrate(0.0_wp, [1.0_wp, 0.0_wp], 0.1_wp, 1.0_wp, xf)
        call s%dense_interval(available, t0=t0_last)
        if (.not. available) error stop 'dense_output_step should enable the dense output'
        call s%interpolate((t0_last+1.0_wp)/2.0_wp, x)
        if (s%failed()) error stop 'dense_output_step should enable the dense output'

        write(*,*) 'dense output options: ok'

        end subroutine test_dense_options
    !*********************************************************

    !*********************************************************
        subroutine test_dense_output_step(fevals_no_dense)
        !! test reporting at a fixed step size (`dense_output_step`)

        integer,intent(in) :: fevals_no_dense !! function evaluations for the
                                              !! integration from 0 to 10 without dense output

        real(wp),parameter :: dt_out = 0.01_wp !! output step size
        real(wp),parameter :: tol = 1.0e-10_wp

        real(wp),dimension(2) :: xf
        real(wp) :: tf_actual, gf
        integer :: num_steps

        ! without a report function, it should integrate normally:
        fevals = 0
        call s%initialize(n=2, f=oscillator, rtol=[tol], atol=[tol], dense_output_step=dt_out)
        call s%integrate(0.0_wp, [1.0_wp, 0.0_wp], 0.0_wp, 10.0_wp, xf)
        if (s%failed()) error stop 'integration failed'
        if (fevals /= fevals_no_dense) error stop 'dense_output_step without report should not change the integration'

        ! forward (the final time is on the output grid, so it shouldn't be duplicated):
        call s%initialize(n=2, f=oscillator, rtol=[tol], atol=[tol], report=report_grid, dense_output_step=dt_out)
        n_rep = 0
        call s%integrate(0.0_wp, [1.0_wp, 0.0_wp], 0.0_wp, 10.0_wp, xf)
        call s%info(num_steps=num_steps)
        call check_grid('forward', 0.0_wp, 10.0_wp, dt_out, 1001, num_steps)

        ! backward (the final time is not on the output grid):
        n_rep = 0
        call s%integrate(10.0_wp, exact_oscillator(10.0_wp), 0.0_wp, 0.255_wp, xf)
        call s%info(num_steps=num_steps)
        call check_grid('backward', 10.0_wp, 0.255_wp, dt_out, 976, num_steps)

        ! event finding: x(1) = cos(t) = 0.5  -->  t = pi/3
        call s%initialize(n=2, f=oscillator, g=event, rtol=[tol], atol=[tol], &
                          report=report_grid, dense_output_step=dt_out)
        n_rep = 0
        call s%integrate_to_event(0.0_wp, [1.0_wp, 0.0_wp], 0.0_wp, 10.0_wp, 1.0e-12_wp, tf_actual, xf, gf)
        call s%info(num_steps=num_steps)
        if (abs(tf_actual - pi/3.0_wp) > 1.0e-10_wp) error stop 'event time error too large'
        call check_grid('event', 0.0_wp, tf_actual, dt_out, 106, num_steps)

        end subroutine test_dense_output_step
    !*********************************************************

    !*********************************************************
        subroutine check_grid(name, t0, tf, dt_out, n_expected, num_steps)
        !! check the points reported by [[report_grid]]

        character(len=*),intent(in) :: name
        real(wp),intent(in) :: t0, tf, dt_out
        integer,intent(in) :: n_expected !! expected number of reported points
        integer,intent(in) :: num_steps !! number of steps taken by the integrator

        integer :: i
        real(wp) :: direction, max_err, max_grid_err

        if (s%failed()) error stop 'integration failed'
        direction = sign(1.0_wp, tf-t0)
        max_err = 0.0_wp
        max_grid_err = 0.0_wp
        do i = 1, n_rep
            max_err = max(max_err, maxval(abs(x_rep(:,i) - exact_oscillator(t_rep(i)))))
            if (i < n_rep) max_grid_err = max(max_grid_err, abs(t_rep(i) - (t0 + direction*(i-1)*dt_out)))
        end do
        write(*,'(1X,A,A,I5,A,I5,A,ES10.3,A,ES10.3)') name, ': points: ', n_rep, &
            '  steps: ', num_steps, '  max error: ', max_err, '  max grid error: ', max_grid_err

        if (n_rep /= n_expected) error stop 'unexpected number of dense output points'
        if (t_rep(1) /= t0 .or. t_rep(n_rep) /= tf) error stop 'first and last points should be t0 and tf'
        if (max_grid_err > 0.0_wp) error stop 'dense output points are not on the grid'
        if (direction*(t_rep(n_rep) - t_rep(n_rep-1)) <= 0.0_wp .or. &
            direction*(t_rep(n_rep) - t_rep(n_rep-1)) > dt_out*(1.0_wp+1.0e-12_wp)) &
            error stop 'invalid last interval'
        if (max_err > 1.0e-8_wp) error stop 'dense output error too large'
        if (num_steps >= n_rep) error stop 'expected fewer steps than reported points'

        end subroutine check_grid
    !*********************************************************

    !*********************************************************
        subroutine report_grid(me,t,x)
        !! save the reported points
        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x
        if (n_rep >= max_rep) error stop 'too many reported points'
        n_rep = n_rep + 1
        t_rep(n_rep) = t
        x_rep(:,n_rep) = x
        end subroutine report_grid
    !*********************************************************

    !*********************************************************
        subroutine test_no_dense()
        !! a method without a dense output should raise an error

        type(rkf45_class) :: s_nodense
        integer :: istatus

        call s_nodense%initialize(n=1, f=nonlinear)
        call s_nodense%integrate(0.0_wp, [1.0_wp], 0.1_wp, 0.1_wp, yf)
        call s_nodense%interpolate(0.05_wp, y)
        call s_nodense%status(istatus)
        if (istatus /= RKLIB_ERROR_NO_DENSE_OUTPUT) error stop 'expected no dense output error'

        call s_nodense%initialize(n=1, f=nonlinear, dense_output_step=0.1_wp)
        call s_nodense%status(istatus)
        if (istatus /= RKLIB_ERROR_NO_DENSE_OUTPUT) error stop 'expected no dense output error for dense_output_step'

        call s_nodense%initialize(n=1, f=nonlinear, dense_output=.true.)
        call s_nodense%status(istatus)
        if (istatus /= RKLIB_ERROR_NO_DENSE_OUTPUT) error stop 'expected no dense output error for dense_output'
        write(*,*) ''
        write(*,*) 'no dense output error: ok'

        end subroutine test_no_dense
    !*********************************************************

    !*********************************************************
        subroutine test_interpolate_in_f()
        !! calling `interpolate` from the derivative function should raise an error

        type(rkdp54_class) :: s54
        type(rkdp85_class) :: s85
        type(stepsize_class) :: sz
        integer :: istatus

        call sz%initialize(fixed_step_mode=.true.)

        ! during a step (there is no dense output for the step that is in progress):
        call s54%initialize(n=1, f=nonlinear_interpolate, stepsize_method=sz, dense_output=.true.)
        f_interpolate_status = 0
        call s54%integrate(0.0_wp, [1.0_wp], 0.1_wp, 0.1_wp, yf)
        if (f_interpolate_status /= RKLIB_ERROR_NO_DENSE_OUTPUT) &
            error stop 'expected no dense output error from interpolate in f'

        ! during the extra function evaluations for the dense output
        ! (this would be a recursive call to interpolate):
        call s85%initialize(n=1, f=nonlinear, stepsize_method=sz, dense_output=.true.)
        call s85%integrate(0.0_wp, [1.0_wp], 0.1_wp, 0.1_wp, yf)
        if (s85%failed()) error stop 'integration failed'
        call s85%initialize(n=1, f=nonlinear_interpolate, stepsize_method=sz, dense_output=.true.)
        call s85%integrate(0.0_wp, [1.0_wp], 0.1_wp, 0.1_wp, yf)
        f_interpolate_status = 0
        call s85%interpolate(0.05_wp, y)
        call s85%status(istatus)
        if (f_interpolate_status /= RKLIB_ERROR_DENSE_RECURSIVE .or. &
            istatus /= RKLIB_ERROR_DENSE_RECURSIVE) &
            error stop 'expected recursive interpolate error'

        write(*,*) 'interpolate called from f error: ok'

        end subroutine test_interpolate_in_f
    !*********************************************************

    !*********************************************************
        subroutine oscillator(me,t,x,xdot)
        !! x'' = -x
        class(rk_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot
        xdot = [x(2), -x(1)]
        fevals = fevals + 1
        end subroutine oscillator

        pure function exact_oscillator(t) result(x)
        real(wp),intent(in) :: t
        real(wp),dimension(2) :: x
        x = [cos(t), -sin(t)]
        end function exact_oscillator

        subroutine event(me,t,x,g)
        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x
        real(wp),intent(out)             :: g
        g = x(1) - 0.5_wp
        end subroutine event

        subroutine nonlinear(me,t,x,xdot)
        !! y' = -2 t y^2  -->  y = 1/(1+t^2)
        class(rk_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot
        xdot = -2.0_wp*t*x**2
        end subroutine nonlinear

        subroutine nonlinear_interpolate(me,t,x,xdot)
        !! same as [[nonlinear]], but also calls `interpolate` (which isn't allowed)
        class(rk_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot
        real(wp),dimension(1) :: xi
        call nonlinear(me,t,x,xdot)
        if (f_interpolate_status == 0) then
            call me%interpolate(t, xi)
            call me%status(f_interpolate_status)
        end if
        end subroutine nonlinear_interpolate

        pure real(wp) function exact_nonlinear(t)
        real(wp),intent(in) :: t
        exact_nonlinear = 1.0_wp/(1.0_wp+t**2)
        end function exact_nonlinear
    !*********************************************************

    !*********************************************************
        subroutine report(me,t,x)
        !! check the dense output at points within the step that was just taken
        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x

        integer :: i
        real(wp) :: ti, t_prev
        real(wp),dimension(2) :: xi, x_prev
        logical :: available

        ! the step that was just taken (not available for the first point):
        call me%dense_interval(available, t0=t_prev, x0=x_prev)

        if (available) then
            do i = 0, n_interp
                ti = t_prev + (t-t_prev)*real(i,wp)/real(n_interp,wp)
                call me%interpolate(ti, xi)
                if (me%failed()) error stop 'interpolate failed'
                max_err = max(max_err, maxval(abs(xi - exact_oscillator(ti))))
                if (i==0)        max_endpoint_err = max(max_endpoint_err, maxval(abs(xi - x_prev)))
                if (i==n_interp) max_endpoint_err = max(max_endpoint_err, maxval(abs(xi - x)))
            end do
        end if

        end subroutine report
    !*********************************************************

    !*********************************************************
        function single_step_error(h) result(err)
        !! take one step of size `h`, return the max interpolation error inside it
        real(wp),intent(in) :: h
        real(wp) :: err
        integer :: i
        real(wp) :: ti
        call s%integrate(t0_step, y0, h, t0_step+h, yf)
        err = 0.0_wp
        do i = 1, n_interp-1
            ti = t0_step + h*real(i,wp)/real(n_interp,wp)
            call s%interpolate(ti, y)
            err = max(err, abs(y(1) - exact_nonlinear(ti)))
        end do
        end function single_step_error
    !*********************************************************

    end program rk_test_dense
!*****************************************************************************************
