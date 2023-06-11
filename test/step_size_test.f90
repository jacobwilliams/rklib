
!*****************************************************************************************
!>
!  Unit tests for step size adjustment routines.

    program step_size_test

    use runge_kutta_module, wp => rk_module_rk

    implicit none

    type(stepsize_class)   :: s1     !! for testing the different methods
    type(stepsize_class)   :: s2     !! for testing the different methods
    type(stepsize_class)   :: s3     !! for testing the different methods
    real(wp)               :: h      !! current step size
    real(wp)               :: tol    !! abs error tolerance
    real(wp)               :: err    !! truncation error estimate
    integer                :: p      !! order of the method
    real(wp)               :: hnew   !! new step size
    logical                :: accept !! if the step is accepted

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
    call s1%compute_stepsize(h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=1,relative_err=F : hnew = ', hnew

    call s2%initialize(accept_mode=2)
    call s2%compute_stepsize(h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=2,relative_err=F : hnew = ', hnew

    call s1%initialize(accept_mode=1, relative_err=.true.)
    call s1%compute_stepsize(h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=1,relative_err=T : hnew = ', hnew

    call s2%initialize(accept_mode=2, relative_err=.true.)
    call s2%compute_stepsize(h,tol,err,p,hnew,accept)
    write(*,*) 'accept_mode=2,relative_err=T : hnew = ', hnew

    end program step_size_test
!*****************************************************************************************
