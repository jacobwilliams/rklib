
!*****************************************************************************************
!>
!  Unit tests for step size adjustment routines.

    program step_size_test

    use rklib_module, wp => rk_module_rk

    implicit none

    type(stepsize_class)   :: s1     !! for testing the different methods
    type(stepsize_class)   :: s2     !! for testing the different methods
    real(wp)               :: h      !! current step size
    real(wp),dimension(1)  :: tol    !! abs error tolerance
    real(wp),dimension(1)  :: err    !! truncation error estimate
    integer                :: p      !! order of the method
    real(wp)               :: hnew   !! new step size
    logical                :: accept !! if the step is accepted
    real(wp)               :: expected !! expected new step size
    real(wp),dimension(1)  :: tol2   !! tolerance from the [[stepsize_class]]

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' step_size_test'
    write(*,*) '---------------'
    write(*,*) ''

    h   = 10.0_wp
    tol = 1.0e-9_wp
    err = 1.0e-7_wp
    p   = 4

    call s1%initialize(accept_mode=1)
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=1,relative_err=F : hnew = ', hnew

    call s2%initialize(accept_mode=2)
    call s2%compute_stepsize(1,h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=2,relative_err=F : hnew = ', hnew

    call s1%initialize(accept_mode=1, relative_err=.true.)
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=1,relative_err=T : hnew = ', hnew

    call s2%initialize(accept_mode=2, relative_err=.true.)
    call s2%compute_stepsize(1,h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=2,relative_err=T : hnew = ', hnew

    ! a case where the error per step and the error per unit step criteria differ:
    !  error per step:      err/tol     = 5   --> reject
    !  error per unit step: err/(h*tol) = 0.5 --> accept
    write(*,*) ''
    h   = 10.0_wp
    tol = 1.0e-9_wp
    err = 5.0e-9_wp
    p   = 4
    call s1%initialize(accept_mode=1, relative_err=.false.)
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=1,relative_err=F : accept = ', accept, ' hnew = ', hnew
    if (accept .or. hnew >= h) error stop 'error per step: step should be rejected'
    call s2%initialize(accept_mode=2, relative_err=.false.)
    call s2%compute_stepsize(1,h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=2,relative_err=F : accept = ', accept, ' hnew = ', hnew
    if (accept .or. hnew >= h) error stop 'error per step: step should be rejected'
    call s1%initialize(accept_mode=1, relative_err=.true.)
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=1,relative_err=T : accept = ', accept, ' hnew = ', hnew
    if (.not. accept .or. hnew <= h) error stop 'error per unit step: step should be accepted'
    call s2%initialize(accept_mode=2, relative_err=.true.)
    call s2%compute_stepsize(1,h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=2,relative_err=T : accept = ', accept, ' hnew = ', hnew
    if (.not. accept .or. hnew <= h) error stop 'error per unit step: step should be accepted'

    ! initialize resets the options that are not specified (relative_err is false again):
    call s1%initialize(accept_mode=1, relative_err=.true.)
    call s1%initialize(accept_mode=1)
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    if (accept) error stop 'initialize should reset the options that are not specified'

    ! rms norm:
    if (abs(rms_func([3.0_wp, 4.0_wp]) - sqrt(12.5_wp)) > 10*epsilon(1.0_wp)) error stop 'rms_func'

    ! scaling of the relative tolerance:
    call s1%initialize()
    tol2 = s1%tolerance([1.0e-3_wp], [1.0e-6_wp], [2.0_wp], [-1.0_wp])
    if (abs(tol2(1) - (1.0e-3_wp*1.0_wp + 1.0e-6_wp)) > 1.0e-18_wp) error stop 'rtol_scaling=1'
    call s1%initialize(rtol_scaling=2)
    tol2 = s1%tolerance([1.0e-3_wp], [1.0e-6_wp], [2.0_wp], [-1.0_wp])
    if (abs(tol2(1) - (1.0e-3_wp*2.0_wp + 1.0e-6_wp)) > 1.0e-18_wp) error stop 'rtol_scaling=2'

    ! PI control: two accepted steps, with errors 0.5 and 0.25 (relative to the tolerance)
    ! [exponent: 1/(4+1) - 0.75*0.04 = 0.17, and the error of the previous step starts at 1e-4]
    h   = 1.0_wp
    tol = 1.0_wp
    p   = 4
    call s1%initialize(pi_beta=0.04_wp)
    err = 0.5_wp
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    expected = h * 0.9_wp * (1.0_wp/0.5_wp)**0.17_wp * 1.0e-4_wp**0.04_wp
    write(*,*) 'PI control: hnew = ', hnew, ' expected = ', expected
    if (.not. accept .or. abs(hnew-expected) > 1.0e-14_wp) error stop 'PI control (first step)'
    err = 0.25_wp
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    expected = h * 0.9_wp * (1.0_wp/0.25_wp)**0.17_wp * 0.5_wp**0.04_wp
    write(*,*) 'PI control: hnew = ', hnew, ' expected = ', expected
    if (.not. accept .or. abs(hnew-expected) > 1.0e-14_wp) error stop 'PI control (second step)'

    ! no increase of the step size after a rejected step:
    call s1%initialize(no_increase_after_reject=.true.)
    call s2%initialize()
    err = 5.0_wp
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    if (accept) error stop 'the step should be rejected'
    call s2%compute_stepsize(1,h,tol,err,p,hnew,accept)
    err = 0.01_wp
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    if (.not. accept .or. hnew > h) error stop 'the step size should not increase after a rejected step'
    call s2%compute_stepsize(1,h,tol,err,p,hnew,accept)
    if (.not. accept .or. hnew <= h) error stop 'the step size should increase (default)'

    ! a negative p_exponent_offset (e.g., for the error per unit step with a p(p-1) pair):
    ! exponent 1/(4-1) = 1/3
    h   = 1.0_wp
    tol = 1.0_wp
    p   = 4
    err = 0.5_wp
    call s1%initialize(p_exponent_offset=-1)
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    expected = h * 0.9_wp * (1.0_wp/0.5_wp)**(1.0_wp/3.0_wp)
    write(*,*) 'p_exponent_offset=-1: hnew = ', hnew, ' expected = ', expected
    if (.not. accept .or. abs(hnew-expected) > 1.0e-14_wp) error stop 'p_exponent_offset=-1'
    ! (p + p_exponent_offset is limited to at least 1: exponent 1)
    p = 1
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    expected = min(2.0_wp, h * 0.9_wp * (1.0_wp/0.5_wp))
    if (.not. accept .or. abs(hnew-expected) > 1.0e-14_wp) error stop 'p_exponent_offset=-1 with p=1'

    ! the "extremely small error" check is relative to the tolerance:
    ! (an error of 1e-15 is not small if the tolerance is 1e-16)
    call s1%initialize()
    tol = 1.0e-16_wp
    err = 1.0e-15_wp
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    if (accept) error stop 'a small error that is larger than the tolerance should be rejected'

    call test_filters()

    write(*,*) ''
    write(*,*) 'PASSED'

    contains

    subroutine test_filters()

    !! Tests of the digital filters (see [[stepsize_class]]).

    character(len=*),dimension(*),parameter :: filters = [character(len=7) :: 'I', 'PI3040', 'PI3333', 'PI4020', &
                                                          'PC11', 'H0211', 'H211b', 'H211PI', 'H0312', 'H312b', &
                                                          'H312PID', 'H0321', 'H321']
    real(wp),dimension(*),parameter :: errs = [0.5_wp, 2.0_wp, 0.3_wp, 0.01_wp, 1.0e-20_wp, 0.8_wp, 0.05_wp]
    real(wp),parameter :: sf = 0.9_wp   !! safety factor (default)
    real(wp),parameter :: k = 5.0_wp    !! p + p_exponent_offset (for p=4)
    real(wp),parameter :: theta = sf**k !! target error of the filters
    real(wp),parameter :: c = 1.0e3_wp  !! for the error model: err = c*h**k

    real(wp) :: ha, hb, hnew_a, hnew_b, h1, h2, h3, hstar
    logical :: accept_a, accept_b, found
    real(wp),dimension(3) :: beta
    real(wp),dimension(2) :: alpha
    integer :: i, j

    write(*,*) ''
    write(*,*) 'digital filters:'

    tol = 1.0_wp
    p   = 4

    ! the 'I' filter is the standard controller (also for rejected and very small errors):
    call s1%initialize()
    call s2%initialize(filter='I')
    ha = 1.0_wp
    hb = 1.0_wp
    do i = 1, size(errs)
        call s1%compute_stepsize(1,ha,tol,[errs(i)],p,hnew_a,accept_a)
        call s2%compute_stepsize(1,hb,tol,[errs(i)],p,hnew_b,accept_b)
        if (accept_a .neqv. accept_b .or. abs(hnew_a-hnew_b) > 1.0e-14_wp*hnew_a) &
            error stop 'the I filter should be the same as the standard controller'
        ha = hnew_a
        hb = hnew_b
    end do

    ! H211b: (theta/e_n)^(1/(4k)) * (theta/e_{n-1})^(1/(4k)) * (h_n/h_{n-1})^(-1/4)
    ! [the previous step is not available for the first step, so it's the same as the current one]
    call s1%initialize(filter='H211b')
    h1 = 1.0_wp
    call s1%compute_stepsize(1,h1,tol,[0.5_wp],p,h2,accept_a)
    expected = h1 * (theta/0.5_wp)**(0.5_wp/k)
    write(*,*) '  H211b (step 1): hnew = ', h2, ' expected = ', expected
    if (.not. accept_a .or. abs(h2-expected) > 1.0e-14_wp) error stop 'H211b (step 1)'
    call s1%compute_stepsize(1,h2,tol,[0.25_wp],p,h3,accept_a)
    expected = h2 * (theta/0.25_wp)**(0.25_wp/k) * (theta/0.5_wp)**(0.25_wp/k) * (h2/h1)**(-0.25_wp)
    write(*,*) '  H211b (step 2): hnew = ', h3, ' expected = ', expected
    if (.not. accept_a .or. abs(h3-expected) > 1.0e-14_wp) error stop 'H211b (step 2)'

    ! H312b: three accepted steps; the third one uses the full history
    call s1%initialize(filter='H312b')
    call s2%initialize(filter_beta=[1.0_wp, 2.0_wp, 1.0_wp]/8.0_wp, filter_alpha=[3.0_wp, 1.0_wp]/8.0_wp)
    ha = 1.0_wp
    hb = 1.0_wp
    call s1%compute_stepsize(1,ha,tol,[0.5_wp],p,h1,accept_a)
    call s1%compute_stepsize(1,h1,tol,[0.25_wp],p,h2,accept_a)
    call s1%compute_stepsize(1,h2,tol,[0.1_wp],p,h3,accept_a)
    expected = h2 * (theta/0.1_wp)**(1.0_wp/(8*k)) * (theta/0.25_wp)**(2.0_wp/(8*k)) * &
               (theta/0.5_wp)**(1.0_wp/(8*k)) * (h2/h1)**(-3.0_wp/8) * (h1/ha)**(-1.0_wp/8)
    write(*,*) '  H312b (step 3): hnew = ', h3, ' expected = ', expected
    if (.not. accept_a .or. abs(h3-expected) > 1.0e-14_wp) error stop 'H312b (step 3)'
    ! the same filter, with user-specified coefficients:
    call s2%compute_stepsize(1,hb,tol,[0.5_wp],p,hnew_b,accept_b)
    call s2%compute_stepsize(1,hnew_b,tol,[0.25_wp],p,hnew_b,accept_b)
    call s2%compute_stepsize(1,hnew_b,tol,[0.1_wp],p,hnew_b,accept_b)
    if (abs(hnew_b-h3) > 1.0e-14_wp) error stop 'user-specified filter coefficients'

    ! initialize resets the history (the first step is the same as before):
    call s1%initialize(filter='H211b')
    call s1%compute_stepsize(1,1.0_wp,tol,[0.5_wp],p,hnew_a,accept_a)
    if (abs(hnew_a - (theta/0.5_wp)**(0.5_wp/k)) > 1.0e-14_wp) error stop 'initialize should reset the filter history'

    ! a rejected step uses the standard controller:
    call s1%initialize(filter='H312b')
    call s2%initialize()
    call s1%compute_stepsize(1,1.0_wp,tol,[0.5_wp],p,hnew_a,accept_a)
    call s1%compute_stepsize(1,1.0_wp,tol,[3.0_wp],p,hnew_a,accept_a)
    call s2%compute_stepsize(1,1.0_wp,tol,[3.0_wp],p,hnew_b,accept_b)
    if (accept_a .or. abs(hnew_a-hnew_b) > 1.0e-14_wp) error stop 'filter: rejected step'

    ! smooth limiter: 1 + kappa*atan((f-1)/kappa), then the hard limits
    call s1%initialize(limiter=2)
    call s1%compute_stepsize(1,1.0_wp,tol,[0.01_wp],p,hnew_a,accept_a)
    expected = 1.0_wp + atan(sf*(1.0_wp/0.01_wp)**(1.0_wp/k) - 1.0_wp)
    write(*,*) '  smooth limiter: hnew = ', hnew_a, ' expected = ', expected
    if (abs(hnew_a-expected) > 1.0e-14_wp) error stop 'smooth limiter'
    call s1%initialize(limiter=2, limiter_kappa=0.5_wp)
    call s1%compute_stepsize(1,1.0_wp,tol,[2.0_wp],p,hnew_a,accept_a)
    expected = 1.0_wp + 0.5_wp*atan((sf*(1.0_wp/2.0_wp)**(1.0_wp/k) - 1.0_wp)/0.5_wp)
    if (accept_a .or. abs(hnew_a-expected) > 1.0e-14_wp) error stop 'smooth limiter (kappa=0.5, rejected)'

    ! the names are not case sensitive, and unknown names are not found:
    beta = -1.0_wp; alpha = -1.0_wp
    call stepsize_filter_coefficients('h211B', beta, alpha, found)
    if (.not. found .or. any(beta /= [0.25_wp, 0.25_wp, 0.0_wp]) .or. any(alpha /= [0.25_wp, 0.0_wp])) &
        error stop 'stepsize_filter_coefficients(h211B)'
    beta = -1.0_wp; alpha = -1.0_wp
    call stepsize_filter_coefficients('H999', beta, alpha, found)
    if (found .or. any(beta /= -1.0_wp) .or. any(alpha /= -1.0_wp)) error stop 'stepsize_filter_coefficients(H999)'

    ! for the asymptotic error model err = c*h^k, the step size of all the filters converges to the
    ! step size for which the error is the target error theta: h* = (theta/c)^(1/k)
    ! [with wide hard limits, so that they are not used]
    hstar = (theta/c)**(1.0_wp/k)
    do j = 1, size(filters)
        call s1%initialize(filter=trim(filters(j)), hfactor_accept=100.0_wp, hfactor_reject=0.01_wp, hmin=1.0e-12_wp)
        ha = 0.01_wp
        do i = 1, 300
            call s1%compute_stepsize(1,ha,tol,[c*ha**k],p,hnew_a,accept_a)
            ha = hnew_a
        end do
        write(*,'(a,a8,a,es10.2)') '   convergence: ', trim(filters(j)), ': h/h* - 1 = ', ha/hstar - 1.0_wp
        if (abs(ha/hstar - 1.0_wp) > 1.0e-10_wp) error stop 'the step size of the filter does not converge'
    end do

    end subroutine test_filters

    end program step_size_test
!*****************************************************************************************
