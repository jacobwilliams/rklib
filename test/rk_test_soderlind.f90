!*****************************************************************************************
!>
!  Comparison of the step size controllers of [[stepsize_class]]: the standard
!  (elementary) controller, the PI (Lund) controller, and Söderlind's digital filters
!  (with the hard limits and with the smooth arctangent limiter).
!
!  The test problems are:
!
!  1. `kepler`: an eccentric Kepler orbit (`e=0.9`, two revolutions), which requires large
!     changes in the step size. The exact solution is computed by solving Kepler's equation.
!  2. `arenstorf`: the Arenstorf orbit (one period), a restricted three-body problem
!     (Hairer, Norsett & Wanner, Solving ODEs I, section II.0). The exact solution at the
!     final time is the initial state.
!  3. `vanderpol`: the Van der Pol oscillator with `mu=10` (`t` from 0 to 20), which has fast
!     transitions between slow phases. The reference solution is computed with [[rkv98e]] at a
!     tight tolerance.
!  4. `prothero`: the Prothero-Robinson problem `x' = lambda*(x-cos(t)) - sin(t)`, with
!     `lambda = -1000` and `x(0)=1` (exact solution `cos(t)`), for `t` from 0 to 10. This is a
!     stiff problem, so the step size of an explicit method is limited by its stability
!     rather than by the accuracy, which typically gives an oscillating step size sequence
!     with many rejected steps for the standard controller (the case that the PI controllers
!     are designed for, see Hairer & Wanner, Solving ODEs II, section IV.2, and Söderlind (2003)).
!
!  Each problem is integrated with three methods ([[rkdp54]], [[rkv65e]], [[rkdp85]]) for
!  four tolerances. The other step size options are the ones used in `DOPRI5` and `DOP853`
!  (except for `pi_beta`): the `rms` norm, `rtol_scaling=2`, `p_exponent_offset=0`,
!  `hfactor_accept=6`, `hfactor_reject=0.333`, and `no_increase_after_reject`.
!
!  For each problem, the results are summarized in a markdown table (`rk_test_soderlind.md`),
!  relative to the standard controller (geometric means over the methods and the tolerances):
!
!  * `rel. cost`: the number of function evaluations,
!  * `rel. error`: the error at the final time,
!  * `equiv. cost`: the relative cost for the same accuracy as the standard controller,
!    assuming that the error is proportional to `fevals**(-p)` (`p` is the order of the method):
!    `rel_cost * rel_error**(1/p)`. A value less than 1 is more efficient than the standard
!    controller,
!  * `rejected`: the total number of rejected steps,
!  * `roughness`: the mean of `|log(h(i+1)/h(i))|` over the accepted steps
!    (a measure of the smoothness of the step size sequence).
!
!  The test checks that:
!
!  * all the integrations succeed,
!  * all the controllers have a similar efficiency (`equiv. cost` within a factor of 1.5),
!  * for the stiff problem, the PI controllers (Lund, `PI3040`, `PI3333`, `PI4020`) have much
!    fewer rejected steps and a smoother step size sequence than the standard controller
!    (for this problem, the low-pass filters `H211b` and `H312b` do not reduce the number of
!    rejected steps much), and
!  * for the nonstiff problems, the low-pass filters (`H211b`, `H211PI`, `H312b`, `H312PID`)
!    give a smoother step size sequence than the standard controller.

    module rk_test_soderlind_module

    !! Test problems for [[rk_test_soderlind]].

    use rklib_module, wp => rk_module_rk

    implicit none

    real(wp),parameter :: pi = acos(-1.0_wp)

    real(wp),parameter :: ecc = 0.9_wp !! eccentricity of the Kepler orbit (`a = mu = 1`)

    ! Arenstorf orbit:
    real(wp),parameter :: amu = 0.012277471_wp
    real(wp),parameter :: amup = 1.0_wp - amu
    real(wp),parameter :: arenstorf_period = 17.0652165601579625588917206249_wp
    real(wp),dimension(4),parameter :: arenstorf_x0 = [0.994_wp, 0.0_wp, 0.0_wp, &
                                                       -2.00158510637908252240537862224_wp]

    real(wp),parameter :: vdp_mu = 10.0_wp !! Van der Pol parameter
    real(wp),parameter :: pr_lambda = -1000.0_wp !! Prothero-Robinson parameter

    integer :: problem = 1 !! the current problem
    integer :: fevals = 0  !! number of function evaluations

    ! for the step sizes:
    real(wp) :: t_last = 0.0_wp   !! the last reported time
    real(wp) :: h_last = 0.0_wp   !! the last step size
    real(wp) :: rough_sum = 0.0_wp !! sum of `|log(h(i+1)/h(i))|`
    real(wp) :: rough_last = 0.0_wp !! the last term of `rough_sum`
    integer  :: n_rough = 0        !! number of terms in `rough_sum`
    integer  :: n_reports = 0      !! number of reported points

    contains

    subroutine deriv(me,t,x,xdot)
        !! the derivative function for all the problems
        class(rk_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot
        real(wp) :: r1, r2
        fevals = fevals + 1
        select case (problem)
        case (1) ! kepler
            r1 = hypot(x(1),x(2))**3
            xdot = [x(3), x(4), -x(1)/r1, -x(2)/r1]
        case (2) ! arenstorf
            r1 = ((x(1)+amu)**2 + x(2)**2)**1.5_wp
            r2 = ((x(1)-amup)**2 + x(2)**2)**1.5_wp
            xdot(1) = x(3)
            xdot(2) = x(4)
            xdot(3) = x(1) + 2.0_wp*x(4) - amup*(x(1)+amu)/r1 - amu*(x(1)-amup)/r2
            xdot(4) = x(2) - 2.0_wp*x(3) - amup*x(2)/r1 - amu*x(2)/r2
        case (3) ! van der pol
            xdot = [x(2), vdp_mu*(1.0_wp-x(1)**2)*x(2) - x(1)]
        case (4) ! prothero-robinson
            xdot = [pr_lambda*(x(1)-cos(t)) - sin(t)]
        end select
    end subroutine deriv

    subroutine report(me,t,x)
        !! accumulate the roughness of the step size sequence
        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x
        real(wp) :: h
        n_reports = n_reports + 1
        if (n_reports > 1) then
            h = abs(t - t_last)
            if (h > 0.0_wp) then
                if (h_last > 0.0_wp) then
                    rough_last = abs(log(h/h_last))
                    rough_sum = rough_sum + rough_last
                    n_rough = n_rough + 1
                end if
                h_last = h
            end if
        end if
        t_last = t
    end subroutine report

    pure function kepler_exact(t) result(x)
        !! exact solution of the Kepler problem (starting at periapsis at `t=0`)
        real(wp),intent(in) :: t
        real(wp),dimension(4) :: x
        real(wp) :: ea, dea, r
        integer :: i
        ea = t + ecc*sin(t)
        do i = 1, 100
            dea = (ea - ecc*sin(ea) - t) / (1.0_wp - ecc*cos(ea))
            ea = ea - dea
            if (abs(dea) <= 4*epsilon(1.0_wp)*max(1.0_wp,abs(ea))) exit
        end do
        r = 1.0_wp - ecc*cos(ea)
        x = [cos(ea) - ecc, sqrt(1.0_wp-ecc**2)*sin(ea), -sin(ea)/r, sqrt(1.0_wp-ecc**2)*cos(ea)/r]
    end function kepler_exact

    end module rk_test_soderlind_module
!*****************************************************************************************

!*****************************************************************************************
    program rk_test_soderlind

    use rklib_module, wp => rk_module_rk
    use rk_test_soderlind_module

    implicit none

    integer,parameter :: n_problems = 4
    character(len=*),dimension(n_problems),parameter :: problem_names = &
        [character(len=10) :: 'kepler', 'arenstorf', 'vanderpol', 'prothero']
    character(len=*),dimension(n_problems),parameter :: problem_desc = &
        [character(len=80) :: 'Kepler orbit (e=0.9, two revolutions)', &
                              'Arenstorf orbit (one period)', &
                              'Van der Pol oscillator (mu=10, t=0..20)', &
                              'Prothero-Robinson (lambda=-1000, t=0..10, stiff)']
    integer,dimension(n_problems),parameter :: problem_n = [4, 4, 2, 1]

    integer,parameter :: n_methods = 3
    character(len=*),dimension(n_methods),parameter :: method_names = &
        [character(len=6) :: 'rkdp54', 'rkv65e', 'rkdp85']

    real(wp),dimension(4),parameter :: tols = [1.0e-4_wp, 1.0e-6_wp, 1.0e-8_wp, 1.0e-10_wp]

    ! the controllers: (name, filter, limiter, pi_beta)
    integer,parameter :: n_ctrl = 17
    character(len=*),dimension(n_ctrl),parameter :: ctrl_filter = [character(len=7) :: &
        '', '', 'PI3040', 'PI3333', 'PI4020', 'PC11', 'H0211', 'H211b', 'H211PI', &
        'H0312', 'H312b', 'H312PID', 'H0321', 'H321', 'I', 'H211b', 'H312b']
    integer,dimension(n_ctrl),parameter :: ctrl_limiter = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2]
    real(wp),dimension(n_ctrl),parameter :: ctrl_pi_beta = [0.0_wp, 0.04_wp, spread(0.0_wp,1,n_ctrl-2)]

    character(len=*),parameter :: report_file = 'rk_test_soderlind.md'

    ! results [controller, method, tolerance, problem]:
    integer,dimension(n_ctrl,n_methods,size(tols),n_problems) :: res_fevals = 0, res_rejected = 0, res_steps = 0
    real(wp),dimension(n_ctrl,n_methods,size(tols),n_problems) :: res_error = 0.0_wp, res_rough = 0.0_wp

    ! summary [controller, problem]:
    real(wp),dimension(n_ctrl,n_problems) :: rel_cost, rel_error, equiv_cost, mean_rough
    integer,dimension(n_ctrl,n_problems) :: total_rejected

    real(wp),dimension(:),allocatable :: x_ref
    integer :: iprob, im, itol, ic, iunit, ip
    logical :: all_ok

    all_ok = .true.

    do iprob = 1, n_problems
        problem = iprob
        x_ref = reference_solution(iprob)
        do im = 1, n_methods
            do itol = 1, size(tols)
                do ic = 1, n_ctrl
                    call run(iprob, im, tols(itol), ic, &
                             res_fevals(ic,im,itol,iprob), res_steps(ic,im,itol,iprob), &
                             res_rejected(ic,im,itol,iprob), res_error(ic,im,itol,iprob), &
                             res_rough(ic,im,itol,iprob))
                end do
            end do
        end do
    end do

    ! summary relative to the standard controller (geometric means):
    do iprob = 1, n_problems
        do ic = 1, n_ctrl
            rel_cost(ic,iprob)  = 0.0_wp
            rel_error(ic,iprob) = 0.0_wp
            equiv_cost(ic,iprob) = 0.0_wp
            do im = 1, n_methods
                ip = method_order(im)
                do itol = 1, size(tols)
                    associate (c  => log(real(res_fevals(ic,im,itol,iprob),wp) / real(res_fevals(1,im,itol,iprob),wp)), &
                               er => log(err_floor(res_error(ic,im,itol,iprob)) / err_floor(res_error(1,im,itol,iprob))))
                        rel_cost(ic,iprob)   = rel_cost(ic,iprob)   + c
                        rel_error(ic,iprob)  = rel_error(ic,iprob)  + er
                        equiv_cost(ic,iprob) = equiv_cost(ic,iprob) + c + er/ip
                    end associate
                end do
            end do
            rel_cost(ic,iprob)   = exp(rel_cost(ic,iprob)   / (n_methods*size(tols)))
            rel_error(ic,iprob)  = exp(rel_error(ic,iprob)  / (n_methods*size(tols)))
            equiv_cost(ic,iprob) = exp(equiv_cost(ic,iprob) / (n_methods*size(tols)))
            total_rejected(ic,iprob) = sum(res_rejected(ic,:,:,iprob))
            mean_rough(ic,iprob) = sum(res_rough(ic,:,:,iprob)) / (n_methods*size(tols))
        end do
    end do

    call write_report()

    ! print the summary:
    do iprob = 1, n_problems
        write(*,'(/,a)') trim(problem_desc(iprob))
        write(*,'(a24,5a12)') 'controller', 'rel cost', 'rel error', 'equiv cost', 'rejected', 'roughness'
        do ic = 1, n_ctrl
            write(*,'(a24,3f12.3,i12,f12.4)') controller_name(ic), rel_cost(ic,iprob), rel_error(ic,iprob), &
                                              equiv_cost(ic,iprob), total_rejected(ic,iprob), mean_rough(ic,iprob)
        end do
    end do

    ! checks:
    !  * all the controllers have a similar efficiency as the standard one (within a factor of 1.5)
    do iprob = 1, n_problems
        do ic = 2, n_ctrl
            if (equiv_cost(ic,iprob) > 1.5_wp .or. equiv_cost(ic,iprob) < 1.0_wp/1.5_wp) then
                write(*,*) 'unexpected efficiency: ', trim(problem_names(iprob)), ' ', controller_name(ic)
                all_ok = .false.
            end if
        end do
    end do
    !  * for the stiff problem, the PI controllers (Lund, PI3040, PI3333, PI4020) have much fewer
    !    rejected steps (less than 1/3), and a smoother step size sequence, than the standard controller:
    do ic = 2, n_ctrl
        if ((ctrl_pi_beta(ic) > 0.0_wp .or. any(ctrl_filter(ic) == ['PI3040', 'PI3333', 'PI4020'])) &
            .and. ctrl_limiter(ic) == 1) then
            if (3*total_rejected(ic,4) >= total_rejected(1,4)) then
                write(*,*) 'the PI controller should reduce the number of rejected steps: ', controller_name(ic)
                all_ok = .false.
            end if
            if (mean_rough(ic,4) >= mean_rough(1,4)) then
                write(*,*) 'the PI controller should reduce the roughness of the step sizes: ', controller_name(ic)
                all_ok = .false.
            end if
        end if
    end do
    !  * for the nonstiff problems, the low-pass filters (H211b, H211PI, H312b, H312PID) give a smoother
    !    step size sequence than the standard controller:
    do ic = 2, n_ctrl
        if (any(ctrl_filter(ic) == ['H211b  ', 'H211PI ', 'H312b  ', 'H312PID']) .and. ctrl_limiter(ic) == 1) then
            do iprob = 1, 3
                if (mean_rough(ic,iprob) >= mean_rough(1,iprob)) then
                    write(*,*) 'the filter should reduce the roughness of the step sizes: ', &
                               controller_name(ic), ' ', trim(problem_names(iprob))
                    all_ok = .false.
                end if
            end do
        end if
    end do

    if (.not. all_ok) error stop 'FAILED'
    write(*,'(/,a)') 'results written to '//report_file
    write(*,*) 'PASSED'

    contains

    function controller_name(ic) result(name)
        !! description of a controller
        integer,intent(in) :: ic
        character(len=:),allocatable :: name
        if (ctrl_filter(ic) == '') then
            if (ctrl_pi_beta(ic) > 0.0_wp) then
                name = 'PI (Lund, beta=0.04)'
            else
                name = 'standard'
            end if
        else
            name = trim(ctrl_filter(ic))
        end if
        if (ctrl_limiter(ic) == 2) name = name//' + smooth limiter'
    end function controller_name

    pure real(wp) function err_floor(e)
        !! the error, limited to a small value (for the ratios)
        real(wp),intent(in) :: e
        err_floor = max(e, 1.0e-15_wp)
    end function err_floor

    pure integer function method_order(im)
        integer,intent(in) :: im
        integer,dimension(n_methods),parameter :: orders = [5, 6, 8]
        method_order = orders(im)
    end function method_order

    function reference_solution(iprob) result(x)
        !! the solution at the final time
        integer,intent(in) :: iprob
        real(wp),dimension(:),allocatable :: x
        type(rkv98e_class) :: s
        type(stepsize_class) :: sz
        real(wp),dimension(2) :: xf
        select case (iprob)
        case (1); x = kepler_exact(final_time(iprob))
        case (2); x = arenstorf_x0
        case (3)
            ! reference solution with a high order method at a tight tolerance
            call sz%initialize(hmin=1.0e-12_wp, hmax=1.0_wp, norm=rms_func, rtol_scaling=2, p_exponent_offset=0)
            call s%initialize(n=2, f=deriv, rtol=[1.0e-15_wp,1.0e-15_wp], atol=[1.0e-15_wp,1.0e-15_wp], &
                              stepsize_method=sz, hinit_method=2)
            call s%integrate(0.0_wp, initial_state(iprob), 0.0_wp, final_time(iprob), xf)
            x = xf
        case (4); x = [cos(final_time(iprob))]
        end select
    end function reference_solution

    pure function initial_state(iprob) result(x)
        integer,intent(in) :: iprob
        real(wp),dimension(:),allocatable :: x
        select case (iprob)
        case (1); x = [1.0_wp-ecc, 0.0_wp, 0.0_wp, sqrt((1.0_wp+ecc)/(1.0_wp-ecc))]
        case (2); x = arenstorf_x0
        case (3); x = [2.0_wp, 0.0_wp]
        case (4); x = [1.0_wp]
        end select
    end function initial_state

    pure real(wp) function final_time(iprob)
        integer,intent(in) :: iprob
        select case (iprob)
        case (1); final_time = 4.0_wp*pi
        case (2); final_time = arenstorf_period
        case (3); final_time = 20.0_wp
        case (4); final_time = 10.0_wp
        end select
    end function final_time

    subroutine run(iprob, im, tol, ic, nfevals, nsteps, nrejected, error, roughness)
        !! integrate a problem with a method, tolerance, and controller
        integer,intent(in) :: iprob, im, ic
        real(wp),intent(in) :: tol
        integer,intent(out) :: nfevals, nsteps, nrejected
        real(wp),intent(out) :: error, roughness

        class(rk_variable_step_class),allocatable :: s
        type(stepsize_class) :: sz
        real(wp),dimension(problem_n(iprob)) :: xf
        integer :: n, istatus

        n = problem_n(iprob)
        select case (im)
        case (1); allocate(rkdp54_class :: s)
        case (2); allocate(rkv65e_class :: s)
        case (3); allocate(rkdp85_class :: s)
        end select

        ! the DOPRI5/DOP853 options, with the controller:
        if (ctrl_filter(ic) == '') then
            call sz%initialize(hmin=1.0e-12_wp, hmax=10.0_wp, norm=rms_func, rtol_scaling=2, &
                               p_exponent_offset=0, hfactor_accept=6.0_wp, hfactor_reject=0.333_wp, &
                               no_increase_after_reject=.true., pi_beta=ctrl_pi_beta(ic), &
                               limiter=ctrl_limiter(ic), max_attempts=100)
        else
            call sz%initialize(hmin=1.0e-12_wp, hmax=10.0_wp, norm=rms_func, rtol_scaling=2, &
                               p_exponent_offset=0, hfactor_accept=6.0_wp, hfactor_reject=0.333_wp, &
                               no_increase_after_reject=.true., filter=trim(ctrl_filter(ic)), &
                               limiter=ctrl_limiter(ic), max_attempts=100)
        end if
        call s%initialize(n=n, f=deriv, rtol=spread(tol,1,n), atol=spread(tol,1,n), &
                          stepsize_method=sz, hinit_method=2, report=report, &
                          max_number_of_steps=1000000)

        fevals = 0
        n_reports = 0
        h_last = 0.0_wp
        rough_sum = 0.0_wp
        rough_last = 0.0_wp
        n_rough = 0
        call s%integrate(0.0_wp, initial_state(iprob), 0.0_wp, final_time(iprob), xf)
        call s%status(istatus)
        if (istatus /= 0) then
            write(*,*) 'integration failed: ', trim(problem_names(iprob)), ' ', method_names(im), &
                       ' tol=', tol, ' ', controller_name(ic), ' status=', istatus
            all_ok = .false.
        end if
        call s%info(num_steps=nsteps, num_rejected_steps=nrejected)

        nfevals = fevals
        error = maxval(abs(xf - x_ref))
        ! (the last step is not included, since it is limited by the final time)
        roughness = 0.0_wp
        if (n_rough > 1) roughness = (rough_sum - rough_last) / (n_rough-1)

    end subroutine run

    subroutine write_report()
        !! write the markdown tables
        integer :: iprob, ic, im, itol
        character(len=32) :: c1, c2, c3, c4, c5

        open(newunit=iunit, file=report_file, status='replace')
        write(iunit,'(a)') '# Comparison of the step size controllers'
        write(iunit,'(a)') ''
        write(iunit,'(a)') 'Generated by `test/rk_test_soderlind.f90` (see that file for the details).'
        write(iunit,'(a)') 'Methods: `rkdp54`, `rkv65e`, `rkdp85`. Tolerances (`rtol=atol`): 1e-4, 1e-6, 1e-8, 1e-10.'
        write(iunit,'(a)') 'Other options: the `DOPRI5`/`DOP853` settings (`rms` norm, `rtol_scaling=2`, `p_exponent_offset=0`, '//&
                           '`hfactor_accept=6`, `hfactor_reject=0.333`, `no_increase_after_reject`).'
        write(iunit,'(a)') ''
        write(iunit,'(a)') 'The values are relative to the standard controller (geometric means over the methods and '//&
                           'tolerances). `equiv. cost` is the relative cost for the same accuracy (less than 1 is '//&
                           'more efficient). `roughness` is the mean `|log(h(i+1)/h(i))|`.'
        do iprob = 1, n_problems
            write(iunit,'(a)') ''
            write(iunit,'(a)') '## '//trim(problem_desc(iprob))
            write(iunit,'(a)') ''
            write(iunit,'(a)') 'Controller | rel. cost | rel. error | equiv. cost | rejected | roughness'
            write(iunit,'(a)') '--- | --- | --- | --- | --- | ---'
            do ic = 1, n_ctrl
                write(c1,'(f0.3)') rel_cost(ic,iprob)
                write(c2,'(f0.3)') rel_error(ic,iprob)
                write(c3,'(f0.3)') equiv_cost(ic,iprob)
                write(c4,'(i0)')   total_rejected(ic,iprob)
                write(c5,'(f0.4)') mean_rough(ic,iprob)
                write(iunit,'(a)') controller_name(ic)//' | '//trim(c1)//' | '//trim(c2)//' | '//trim(c3)//&
                                   ' | '//trim(c4)//' | '//trim(c5)
            end do
        end do

        ! details for the stiff problem:
        write(iunit,'(a)') ''
        write(iunit,'(a)') '## Details: '//trim(problem_desc(4))
        write(iunit,'(a)') ''
        write(iunit,'(a)') 'Function evaluations / rejected steps for each method and tolerance.'
        write(iunit,'(a)') ''
        write(iunit,'(a)', advance='no') 'Controller'
        do im = 1, n_methods
            do itol = 1, size(tols)
                write(c1,'(es7.0)') tols(itol)
                write(iunit,'(a)', advance='no') ' | `'//trim(method_names(im))//'` '//trim(adjustl(c1))
            end do
        end do
        write(iunit,'(a)') ''
        write(iunit,'(a)', advance='no') '---'
        do im = 1, n_methods*size(tols)
            write(iunit,'(a)', advance='no') ' | ---'
        end do
        write(iunit,'(a)') ''
        do ic = 1, n_ctrl
            write(iunit,'(a)', advance='no') controller_name(ic)
            do im = 1, n_methods
                do itol = 1, size(tols)
                    write(c1,'(i0,a,i0)') res_fevals(ic,im,itol,4), ' / ', res_rejected(ic,im,itol,4)
                    write(iunit,'(a)', advance='no') ' | '//trim(c1)
                end do
            end do
            write(iunit,'(a)') ''
        end do
        close(iunit)

    end subroutine write_report

    end program rk_test_soderlind
!*****************************************************************************************
