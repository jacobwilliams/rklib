!*****************************************************************************************
!>
!  Comparison of many combinations of the step size control options of rklib
!  (see [[stepsize_class]]) for the Dormand-Prince 8(5) method ([[rkdp85]]),
!  with the reference DOP853 implementation (https://github.com/jacobwilliams/dop853)
!  and FLINT's DOP853 (https://github.com/princemahajan/FLINT).
!
!  The test problem is an eccentric Kepler orbit (`e=0.6`, three revolutions),
!  which requires large changes in the step size. The exact solution is computed by
!  solving Kepler's equation.
!
!  Two tests are done for each tolerance:
!
!  1. Integration to the final time (no output): the error at the final time.
!  2. Dense output at a fixed interval: the RMS (over all the output points) of the
!     max abs error of each point.
!
!  The results (accuracy and the total number of function evaluations) are written to
!  a markdown file (`rk_test_dop853_options.md`). The relative cost and relative error
!  columns are the geometric means (over the tolerances) of the number of function
!  evaluations and of the error, divided by those of the reference DOP853. The
!  equivalent cost combines them: it is the relative cost for the same accuracy as DOP853,
!  assuming that the error is proportional to `fevals**(-8)` (for an 8th order method),
!  i.e., `rel_cost * rel_error**(1/8)`. A value less than 1 is more efficient than DOP853.
!
!  The error criterion is one of four combinations of `relative_err` and `p_exponent_offset`:
!  the error per step with `p_exponent_offset = 1` (the rklib default) or `0` (the exponent
!  `1/(q+1)` for this `8(7)` error estimate, as in DOP853), and the error per unit step with
!  `p_exponent_offset = 0` or `-1` (`-1` is the exponent consistent with `0` for the error
!  per step, since the error per unit step is proportional to one less power of `h`).
!  Note that for the error per unit step, `tol` is the allowed error per unit time,
!  so the same `tol` gives a different accuracy (the equivalent cost accounts for this).
!
!  The test checks that all the integrations succeed, and that the combination of the
!  options that corresponds to DOP853's step size control takes the same steps as DOP853.

    module rk_test_dop853_options_module

    !! The test problem (Kepler orbit) for [[rk_test_dop853_options]].

    use FLINT, FWP => FLINT_WP

    implicit none

    integer,parameter :: kp = FWP !! real kind (checked to be the same as rklib's and dop853's)
    real(kp),parameter :: ecc = 0.6_kp !! eccentricity of the orbit (`a = mu = 1`)

    type,extends(DiffEqSys) :: kepler_sys
        !! the problem for FLINT
    contains
        procedure :: F => kepler_flint
    end type kepler_sys

    contains

    pure function kepler_deriv(x) result(xdot)
        !! equations of motion (`mu=1`)
        real(kp),dimension(:),intent(in) :: x
        real(kp),dimension(4) :: xdot
        real(kp) :: r3
        r3 = hypot(x(1),x(2))**3
        xdot = [x(3), x(4), -x(1)/r3, -x(2)/r3]
    end function kepler_deriv

    pure function kepler_exact(t) result(x)
        !! exact solution (starting at periapsis at `t=0`), from Kepler's equation
        real(kp),intent(in) :: t
        real(kp),dimension(4) :: x
        real(kp) :: ea, dea, r
        integer :: i
        ea = t + ecc*sin(t) ! initial guess for the eccentric anomaly (mean motion = 1)
        do i = 1, 100
            dea = (ea - ecc*sin(ea) - t) / (1.0_kp - ecc*cos(ea))
            ea = ea - dea
            if (abs(dea) <= 4*epsilon(1.0_kp)*max(1.0_kp,abs(ea))) exit
        end do
        r = 1.0_kp - ecc*cos(ea)
        x = [cos(ea) - ecc, sqrt(1.0_kp-ecc**2)*sin(ea), -sin(ea)/r, sqrt(1.0_kp-ecc**2)*cos(ea)/r]
    end function kepler_exact

    function kepler_flint(me, X, Y, Params)
        class(kepler_sys),intent(inout) :: me
        real(kp),intent(in) :: X
        real(kp),dimension(me%n),intent(in) :: Y
        real(kp),dimension(:),intent(in),optional :: Params
        real(kp),dimension(size(Y)) :: kepler_flint
        kepler_flint = kepler_deriv(Y)
    end function kepler_flint

    end module rk_test_dop853_options_module
!*****************************************************************************************

!*****************************************************************************************
    program rk_test_dop853_options

    use rklib_module, wp => rk_module_rk
    use dop853_module, only: dop853_class, dop853_wp
    use FLINT, FWP => FLINT_WP
    use rk_test_dop853_options_module

    implicit none

    integer,parameter :: n = 4 !! number of state variables
    real(wp),parameter :: pi = acos(-1.0_wp)
    real(wp),parameter :: t0 = 0.0_wp
    real(wp),parameter :: tf = 6.0_wp*pi !! three revolutions
    integer,parameter :: n_dense = 1000 !! number of intervals for the dense output
    real(wp),parameter :: dt_out = (tf-t0)/n_dense !! dense output interval
    real(wp),dimension(n),parameter :: x0 = [1.0_wp-ecc, 0.0_wp, 0.0_wp, sqrt((1.0_wp+ecc)/(1.0_wp-ecc))]
    real(wp),dimension(4),parameter :: tols = [1.0e-6_wp, 1.0e-8_wp, 1.0e-10_wp, 1.0e-12_wp]
    character(len=*),parameter :: report_file = 'rk_test_dop853_options.md'

    ! the options that are varied (see [[stepsize_class]]):
    character(len=4),dimension(3),parameter :: norm_names = ['max ','2   ','rms ']
    ! the error criterion: (relative_err, p_exponent_offset)
    logical,dimension(4),parameter :: crit_relerr = [.false., .false., .true., .true.]
    integer,dimension(4),parameter :: crit_offset = [1, 0, 0, -1]
    integer,parameter :: n_options = 3*2*4*2*2*2*2 !! number of option combinations

    integer,parameter :: n_ref = 3 !! number of reference rows (DOP853, FLINT, rklib defaults)
    integer,parameter :: n_rows = n_ref + n_options

    ! results [row, tolerance, test (1=final time, 2=dense)]:
    integer,dimension(n_rows,size(tols),2) :: res_fevals = 0
    real(wp),dimension(n_rows,size(tols),2) :: res_error = 0.0_wp
    character(len=:),allocatable :: row_label(:) !! description of each row (markdown table cells)
    integer :: hairer_row = 0 !! the row with the options of DOP853's step size control
    integer :: hairer_steps, hairer_rejected !! accepted/rejected steps for `hairer_row`
    integer,dimension(2) :: dop_steps, dop_rejected !! accepted/rejected steps for DOP853 [test]

    ! output points (the same for all the methods): t0, t0+dt_out, ..., and tf
    real(wp),dimension(n_dense+1) :: t_grid
    real(wp),dimension(n,n_dense+1) :: x_grid
    integer :: n_grid_out !! number of points saved

    integer :: fevals !! function evaluation counter for rklib
    integer :: i, j, itol, itest, row
    integer :: inorm, iscale, icrit, icomb, ilim, inoinc, ipi

    if (wp /= kp .or. wp /= dop853_wp) then
        write(*,*) 'Skipping the DOP853 options comparison: the libraries are compiled with different real kinds'
        stop
    end if

    t_grid(1:n_dense) = [(t0 + real(i,wp)*dt_out, i = 0, n_dense-1)]
    t_grid(n_dense+1) = tf

    allocate(character(len=200) :: row_label(n_rows))
    row_label(1) = '**DOP853** (reference) | | | | | | | |'
    row_label(2) = '**FLINT** DOP853 | | | | | | | |'
    row_label(3) = '**rklib defaults** (`hinit_method=1`) | max | 1 | F | 1 | 1 | 0.5, 2 | F | 0'
    row = n_ref
    do inorm = 1, 3
    do iscale = 1, 2
    do icrit = 1, 4
    do icomb = 1, 2
    do ilim = 1, 2
    do inoinc = 1, 2
    do ipi = 1, 2
        row = row + 1
        write(row_label(row),'(*(G0))') &
            'rklib ', row - n_ref, ' | ', trim(norm_names(inorm)), ' | ', iscale, ' | ', &
            merge('T', 'F', crit_relerr(icrit)), ' | ', crit_offset(icrit), ' | ', icomb, ' | ', &
            trim(merge('0.5, 2    ', '0.333, 6  ', ilim==1)), ' | ', merge('F', 'T', inoinc==1), ' | ', &
            trim(merge('0   ', '0.04', ipi==1))
        if (inorm==3 .and. iscale==2 .and. icrit==2 .and. icomb==2 .and. ilim==2 .and. inoinc==2 .and. ipi==1) &
            hairer_row = row
    end do; end do; end do; end do; end do; end do; end do
    row_label(hairer_row) = trim(row_label(hairer_row)) // ' (= DOP853)'

    do itest = 1, 2
        do itol = 1, size(tols)
            call run_dop853(1, tols(itol), itest==2)
            call run_flint(2, tols(itol), itest==2)
            call run_rklib_defaults(3, tols(itol), itest==2)
            row = n_ref
            do inorm = 1, 3
            do iscale = 1, 2
            do icrit = 1, 4
            do icomb = 1, 2
            do ilim = 1, 2
            do inoinc = 1, 2
            do ipi = 1, 2
                row = row + 1
                call run_rklib(row, tols(itol), itest==2)
            end do; end do; end do; end do; end do; end do; end do
        end do
    end do

    call write_report()

    write(*,*) 'Results written to: '//report_file
    write(*,*) ''
    write(*,*) 'PASSED'

    contains

!*****************************************************************************************
    subroutine run_rklib(row, tol, dense)
        !! rkdp85 with the combination of options for this row
        integer,intent(in) :: row
        real(wp),intent(in) :: tol
        logical,intent(in) :: dense
        type(stepsize_class) :: sz
        select case (inorm)
        case(1); call sz%initialize(norm=maxval_func)
        case(2); call sz%initialize(norm=norm2_func)
        case(3); call sz%initialize(norm=rms_func)
        end select
        call set_options(sz)
        call integrate_rklib(row, tol, dense, sz, hinit_method=2)
        if (row == hairer_row .and. .not. dense) call check_hairer(row, tol)
    end subroutine run_rklib
!*****************************************************************************************

!*****************************************************************************************
    subroutine set_options(sz)
        !! set the other options of the combination (keeping the norm, which is set by the caller)
        !! [stepsize_class%initialize resets the options that are not specified, so the norm
        !!  has to be specified again]
        type(stepsize_class),intent(inout) :: sz
        select case (inorm)
        case(1); call sz%initialize(norm=maxval_func, rtol_scaling=iscale, relative_err=crit_relerr(icrit), p_exponent_offset=crit_offset(icrit), &
                                    error_combination=icomb, hfactor_reject=merge(0.5_wp, 0.333_wp, ilim==1), &
                                    hfactor_accept=merge(2.0_wp, 6.0_wp, ilim==1), no_increase_after_reject=inoinc==2, &
                                    pi_beta=merge(0.0_wp, 0.04_wp, ipi==1), pi_beta_mult=0.2_wp, hmin=0.0_wp)
        case(2); call sz%initialize(norm=norm2_func, rtol_scaling=iscale, relative_err=crit_relerr(icrit), p_exponent_offset=crit_offset(icrit), &
                                    error_combination=icomb, hfactor_reject=merge(0.5_wp, 0.333_wp, ilim==1), &
                                    hfactor_accept=merge(2.0_wp, 6.0_wp, ilim==1), no_increase_after_reject=inoinc==2, &
                                    pi_beta=merge(0.0_wp, 0.04_wp, ipi==1), pi_beta_mult=0.2_wp, hmin=0.0_wp)
        case(3); call sz%initialize(norm=rms_func, rtol_scaling=iscale, relative_err=crit_relerr(icrit), p_exponent_offset=crit_offset(icrit), &
                                    error_combination=icomb, hfactor_reject=merge(0.5_wp, 0.333_wp, ilim==1), &
                                    hfactor_accept=merge(2.0_wp, 6.0_wp, ilim==1), no_increase_after_reject=inoinc==2, &
                                    pi_beta=merge(0.0_wp, 0.04_wp, ipi==1), pi_beta_mult=0.2_wp, hmin=0.0_wp)
        end select
    end subroutine set_options
!*****************************************************************************************

!*****************************************************************************************
    subroutine run_rklib_defaults(row, tol, dense)
        !! rkdp85 with all the default options
        integer,intent(in) :: row
        real(wp),intent(in) :: tol
        logical,intent(in) :: dense
        type(stepsize_class) :: sz
        call sz%initialize()
        call integrate_rklib(row, tol, dense, sz, hinit_method=1)
    end subroutine run_rklib_defaults
!*****************************************************************************************

!*****************************************************************************************
    subroutine integrate_rklib(row, tol, dense, sz, hinit_method)
        integer,intent(in) :: row
        real(wp),intent(in) :: tol
        logical,intent(in) :: dense
        type(stepsize_class),intent(in) :: sz
        integer,intent(in) :: hinit_method
        type(rkdp85_class) :: s
        real(wp),dimension(n) :: xf
        integer :: itest
        itest = merge(2, 1, dense)
        fevals = 0
        n_grid_out = 0
        if (dense) then
            call s%initialize(n=n, f=f, rtol=[tol], atol=[tol], stepsize_method=sz, hinit_method=hinit_method, &
                              report=report_grid, dense_output_step=dt_out)
        else
            call s%initialize(n=n, f=f, rtol=[tol], atol=[tol], stepsize_method=sz, hinit_method=hinit_method)
        end if
        call s%integrate(t0, x0, 0.0_wp, tf, xf)
        if (s%failed()) error stop 'rklib integration failed'
        call save_result(row, tol, itest, fevals, xf)
        if (row == hairer_row .and. .not. dense) then
            call s%info(num_steps=hairer_steps, num_rejected_steps=hairer_rejected)
        end if
    end subroutine integrate_rklib
!*****************************************************************************************

!*****************************************************************************************
    subroutine run_dop853(row, tol, dense)
        !! the reference DOP853 implementation
        integer,intent(in) :: row
        real(wp),intent(in) :: tol
        logical,intent(in) :: dense
        type(dop853_class) :: dop
        real(wp) :: x
        real(wp),dimension(n) :: y
        integer :: idid, nfcn
        logical :: status_ok
        n_grid_out = 0
        if (dense) then
            call dop%initialize(n=n, fcn=dop_f, solout=dop_solout, nstiff=-1, icomp=[1,2,3,4], status_ok=status_ok)
        else
            call dop%initialize(n=n, fcn=dop_f, nstiff=-1, status_ok=status_ok)
        end if
        if (.not. status_ok) error stop 'dop853 initialization failed'
        x = t0
        y = x0
        call dop%integrate(x, y, tf, [tol], [tol], iout=merge(2, 0, dense), idid=idid)
        if (idid /= 1) error stop 'dop853 integration failed'
        call dop%info(nfcn=nfcn, naccpt=dop_steps(merge(2,1,dense)), nrejct=dop_rejected(merge(2,1,dense)))
        call save_result(row, tol, merge(2, 1, dense), nfcn, y)
    end subroutine run_dop853
!*****************************************************************************************

!*****************************************************************************************
    subroutine run_flint(row, tol, dense)
        !! FLINT's DOP853
        integer,intent(in) :: row
        real(wp),intent(in) :: tol
        logical,intent(in) :: dense
        type(ERK_class) :: erk
        type(kepler_sys),target :: diffeq
        real(wp) :: xf, step
        real(wp),dimension(n) :: yf
        integer :: nfcn
        integer(kind(FLINT_SUCCESS)) :: status
        diffeq%n = n
        diffeq%m = 0
        if (.not. dense) then
            ! [workaround for a FLINT issue: Integrate uses the size of an array that is only
            !  allocated if the interpolation is on, so initialize it once with the interpolation on]
            call erk%Init(diffeq, 20000, Method=ERK_DOP853, ATol=[tol], RTol=[tol], InterpOn=.true.)
        end if
        call erk%Init(diffeq, 20000, Method=ERK_DOP853, ATol=[tol], RTol=[tol], InterpOn=dense)
        xf = tf
        step = 0.0_wp
        call erk%Integrate(t0, x0, xf, yf, StepSz=step)
        if (dense) then
            call erk%Interpolate(t_grid, x_grid)
            n_grid_out = size(t_grid)
        end if
        call erk%Info(status, nFCalls=nfcn)
        if (status /= FLINT_SUCCESS) error stop 'FLINT integration failed'
        call save_result(row, tol, merge(2, 1, dense), nfcn, yf)
    end subroutine run_flint
!*****************************************************************************************

!*****************************************************************************************
    subroutine save_result(row, tol, itest, nfcn, xf)
        !! save the number of function evaluations and the error
        integer,intent(in) :: row, itest, nfcn
        real(wp),intent(in) :: tol
        real(wp),dimension(n),intent(in) :: xf !! final state
        integer :: it, k
        real(wp) :: err
        it = findloc(tols, tol, dim=1)
        res_fevals(row, it, itest) = nfcn
        if (itest == 1) then
            ! error at the final time:
            res_error(row, it, itest) = maxval(abs(xf - kepler_exact(tf)))
        else
            ! RMS of the errors at the output points:
            if (n_grid_out /= size(t_grid)) error stop 'unexpected number of dense output points'
            err = 0.0_wp
            do k = 1, n_grid_out
                err = err + maxval(abs(x_grid(:,k) - kepler_exact(t_grid(k))))**2
            end do
            res_error(row, it, itest) = sqrt(err/n_grid_out)
        end if
        if (.not. (res_error(row, it, itest) < 1.0e-2_wp)) error stop 'error too large'  ! (also catches NaN)
    end subroutine save_result
!*****************************************************************************************

!*****************************************************************************************
    subroutine check_hairer(row, tol)
        !! the options of DOP853's step size control should give the same steps as DOP853
        integer,intent(in) :: row
        real(wp),intent(in) :: tol
        if (hairer_steps /= dop_steps(1) .or. hairer_rejected /= dop_rejected(1)) then
            write(*,*) 'tol: ', tol, ' rklib accepted/rejected: ', hairer_steps, hairer_rejected, &
                       ' DOP853: ', dop_steps(1), dop_rejected(1)
            error stop 'the DOP853 options do not give the same steps as DOP853'
        end if
    end subroutine check_hairer
!*****************************************************************************************

!*****************************************************************************************
    subroutine write_report()
        !! write the results as markdown tables

        integer,parameter :: n_top = 10 !! number of rows in the ranking table
        integer :: iunit, itest, it, row, best_cost, best_err, k, irank
        real(wp),dimension(n_rows) :: rel_cost, rel_err, equiv_cost
        logical,dimension(n_rows) :: ranked
        character(len=*),parameter :: hdr_options = &
            'Method | Norm | `rtol_scaling` | `relative_err` | `p_exponent_offset` | `error_combination` | '//&
            '`hfactor_reject`, `hfactor_accept` | `no_increase_after_reject` | `pi_beta`'

        open(newunit=iunit, file=report_file, status='replace')
        write(iunit,'(A)') '# rkdp85 step size control options vs. DOP853'
        write(iunit,'(A)') ''
        write(iunit,'(A)') 'Generated by `test/rk_test_dop853_options.f90`.'
        write(iunit,'(A)') ''
        write(iunit,'(A)') 'Problem: Kepler orbit with eccentricity 0.6 (`mu = a = 1`), three revolutions '//&
                           '(`t` from 0 to 6&pi;), starting at periapsis. `rtol = atol = tol`.'
        write(iunit,'(A)') ''
        write(iunit,'(A)') 'For all the rklib combinations, `hinit_method=2` (the initial step size of DOP853), '//&
                           '`safety_factor=0.9`, and `pi_beta_mult=0.2` (as in DOP853) are used. '//&
                           'See `stepsize_class` for a description of the options. With `relative_err=T` '//&
                           '(the error per unit step), `tol` is the allowed error per unit time; '//&
                           '`p_exponent_offset=-1` is then the exponent consistent with `0` for the error per step.'
        write(iunit,'(A)') ''
        write(iunit,'(A)') '* **fevals**: total number of function evaluations.'
        write(iunit,'(A)') '* **rel. cost**, **rel. error**: geometric means (over the tolerances) of fevals '//&
                           'and of the error, relative to the reference DOP853 (lower is better).'
        write(iunit,'(A)') '* **equiv. cost**: the relative cost for the same accuracy as DOP853, '//&
                           'assuming that the error is proportional to fevals<sup>-8</sup> (8th order method): '//&
                           'rel. cost &times; rel. error<sup>1/8</sup>. Less than 1 is more efficient than DOP853.'

        do itest = 1, 2
            do row = 1, n_rows
                rel_cost(row) = exp(sum(log(real(res_fevals(row,:,itest),wp)/real(res_fevals(1,:,itest),wp)))/size(tols))
                rel_err(row)  = exp(sum(log(res_error(row,:,itest)/res_error(1,:,itest)))/size(tols))
                equiv_cost(row) = rel_cost(row) * rel_err(row)**(1.0_wp/8.0_wp)
            end do
            write(iunit,'(A)') ''
            if (itest == 1) then
                write(iunit,'(A)') '## Integration to the final time'
                write(iunit,'(A)') ''
                write(iunit,'(A)') 'Error: max abs error of the state at the final time.'
            else
                write(iunit,'(A)') '## Dense output'
                write(iunit,'(A)') ''
                write(iunit,'(A,I0,A)') 'Output at ', n_dense+1, ' evenly spaced points (using the dense output '//&
                    'of each method; `dense_output_step` for rklib). Error: RMS over all the output points '//&
                    'of the max abs error of each point.'
            end if
            write(iunit,'(A)') ''
            best_cost = minloc(rel_cost(n_ref+1:), dim=1) + n_ref
            best_err  = minloc(rel_err(n_ref+1:), dim=1) + n_ref
            write(iunit,'(A)') '* Lowest rel. cost: '//trim(plain(row_label(best_cost)))//' ('//&
                               trim(fmt_real(rel_cost(best_cost),'F6.3'))//', rel. error '//&
                               trim(fmt_real(rel_err(best_cost),'F6.3'))//')'
            write(iunit,'(A)') '* Lowest rel. error: '//trim(plain(row_label(best_err)))//' ('//&
                               trim(fmt_real(rel_err(best_err),'F6.3'))//', rel. cost '//&
                               trim(fmt_real(rel_cost(best_err),'F6.3'))//')'
            ! ranking by the equivalent cost:
            write(iunit,'(A)') ''
            write(iunit,'(A,I0,A)') 'Top ', n_top, ' by equivalent cost (all the rows):'
            write(iunit,'(A)') ''
            write(iunit,'(A)') '| Rank | '//hdr_options//' | rel. cost | rel. error | equiv. cost |'
            write(iunit,'(A)') '|--:|---|---|---|---|---|---|---|---|---|--:|--:|--:|'
            ranked = .false.
            do irank = 1, n_top
                k = minloc(equiv_cost, mask=.not. ranked, dim=1)
                ranked(k) = .true.
                write(iunit,'(A,I0,A)') '| ', irank, ' | '//trim(row_label(k))//' | '//&
                    trim(fmt_real(rel_cost(k),'F6.3'))//' | '//trim(fmt_real(rel_err(k),'F6.3'))//' | '//&
                    trim(fmt_real(equiv_cost(k),'F6.3'))//' |'
            end do
            write(iunit,'(A)') ''
            write(iunit,'(A)') 'All the results:'
            write(iunit,'(A)') ''
            write(iunit,'(A)', advance='no') '| '//hdr_options
            do it = 1, size(tols)
                write(iunit,'(A)', advance='no') ' | fevals (tol='//tol_str(tols(it))//') | error (tol='//&
                                                 tol_str(tols(it))//')'
            end do
            write(iunit,'(A)') ' | rel. cost | rel. error | equiv. cost |'
            write(iunit,'(A)', advance='no') '|---|---|---|---|---|---|---|---|---'
            do it = 1, size(tols)
                write(iunit,'(A)', advance='no') '|--:|--:'
            end do
            write(iunit,'(A)') '|--:|--:|--:|'
            do row = 1, n_rows
                write(iunit,'(A)', advance='no') '| '//trim(row_label(row))
                do it = 1, size(tols)
                    write(iunit,'(A,I0,A,A)', advance='no') ' | ', res_fevals(row,it,itest), ' | ', &
                        trim(fmt_real(res_error(row,it,itest),'ES9.2'))
                end do
                write(iunit,'(A)') ' | '//trim(fmt_real(rel_cost(row),'F6.3'))//' | '//&
                                   trim(fmt_real(rel_err(row),'F6.3'))//' | '//&
                                   trim(fmt_real(equiv_cost(row),'F6.3'))//' |'
            end do
        end do
        close(iunit)

    end subroutine write_report
!*****************************************************************************************

!*****************************************************************************************
    function tol_str(tol) result(str)
        !! a tolerance as a short string, e.g., `1e-6`
        real(wp),intent(in) :: tol
        character(len=:),allocatable :: str
        character(len=10) :: tmp
        write(tmp,'(A,I0)') '1e', nint(log10(tol))
        str = trim(tmp)
    end function tol_str
!*****************************************************************************************

!*****************************************************************************************
    function fmt_real(x, fmt) result(str)
        real(wp),intent(in) :: x
        character(len=*),intent(in) :: fmt
        character(len=:),allocatable :: str
        character(len=30) :: tmp
        write(tmp,'('//fmt//')') x
        str = trim(adjustl(tmp))
    end function fmt_real
!*****************************************************************************************

!*****************************************************************************************
    function plain(label) result(str)
        !! a row label as text: "rklib 12 | rms | 2 | ..." --> "rklib 12 (norm=rms, rtol_scaling=2, ...)"
        character(len=*),intent(in) :: label
        character(len=:),allocatable :: str
        character(len=*),dimension(8),parameter :: names = [character(len=25) :: 'norm', 'rtol_scaling', &
            'relative_err', 'p_exponent_offset', 'error_combination', 'hfactor_reject/accept', 'no_increase_after_reject', 'pi_beta']
        integer :: i, k
        character(len=:),allocatable :: rest, val
        k = index(label, '|')
        str = trim(label(1:k-1)) // ' ('
        rest = label(k+1:)
        do i = 1, 8
            k = index(rest, '|')
            if (k == 0) then
                val = trim(adjustl(rest))
            else
                val = trim(adjustl(rest(1:k-1)))
                rest = rest(k+1:)
            end if
            str = str // trim(names(i)) // '=' // val
            if (i < 8) str = str // ', '
        end do
        str = str // ')'
    end function plain
!*****************************************************************************************

!*****************************************************************************************
    subroutine f(me,t,x,xdot)
        !! rklib derivative function
        class(rk_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot
        xdot = kepler_deriv(x)
        fevals = fevals + 1
    end subroutine f
!*****************************************************************************************

!*****************************************************************************************
    subroutine report_grid(me,t,x)
        !! rklib report function (dense output points)
        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x
        n_grid_out = n_grid_out + 1
        if (n_grid_out > size(t_grid)) error stop 'too many output points'
        if (abs(t - t_grid(n_grid_out)) > 10*spacing(tf)) error stop 'unexpected output time'
        x_grid(:,n_grid_out) = x
    end subroutine report_grid
!*****************************************************************************************

!*****************************************************************************************
    subroutine dop_f(me,x,y,f)
        !! dop853 derivative function
        class(dop853_class),intent(inout) :: me
        real(wp),intent(in)               :: x
        real(wp),dimension(:),intent(in)  :: y
        real(wp),dimension(:),intent(out) :: f
        f = kepler_deriv(y)
    end subroutine dop_f
!*****************************************************************************************

!*****************************************************************************************
    subroutine dop_solout(me,nr,xold,x,y,irtrn,xout)
        !! dop853 output function: the dense output points within the step
        class(dop853_class),intent(inout) :: me
        integer,intent(in)                :: nr
        real(wp),intent(in)               :: xold
        real(wp),intent(in)               :: x
        real(wp),dimension(:),intent(in)  :: y
        integer,intent(inout)             :: irtrn
        real(wp),intent(out)              :: xout
        integer :: i, k
        logical :: last
        xout = 0.0_wp
        if (nr == 1) then
            n_grid_out = 1
            x_grid(:,1) = y
            return
        end if
        ! [dop853 computes the final time as x+(xend-x), which may differ from tf by roundoff]
        last = abs(x - tf) <= 10*spacing(tf)
        ! the points within this step (points at the end of the step are done for the next step,
        ! except for the final time):
        do while (n_grid_out < size(t_grid))
            k = n_grid_out + 1
            if (k == size(t_grid)) then
                if (.not. last) exit
                x_grid(:,k) = y  ! final point
            else
                if (t_grid(k) >= x) exit
                do i = 1, n
                    x_grid(i,k) = me%contd8(i, t_grid(k))
                end do
            end if
            n_grid_out = k
        end do
    end subroutine dop_solout
!*****************************************************************************************

    end program rk_test_dop853_options
!*****************************************************************************************
