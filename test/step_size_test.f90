
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

    ! the "extremely small error" check is relative to the tolerance:
    ! (an error of 1e-15 is not small if the tolerance is 1e-16)
    call s1%initialize()
    tol = 1.0e-16_wp
    err = 1.0e-15_wp
    call s1%compute_stepsize(1,h,tol,err,p,hnew,accept)
    if (accept) error stop 'a small error that is larger than the tolerance should be rejected'

    write(*,*) ''
    write(*,*) 'PASSED'

    end program step_size_test
!*****************************************************************************************
