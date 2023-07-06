!*****************************************************************************************
!>
!  Compare rk78 to another implementation.

    program rk78_validation

    use rklib_module, wp => rk_module_rk
    use rk78_module
    use test_support

    implicit none

    integer,parameter  :: n = 6            !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp !! event location tolerance
    real(wp),parameter :: mu = 3.9860043543609593E+05_wp  !! central body gravitational parameter (km^3/s^2) - Earth

    type(rkf78_class) :: s
    type(rkv78_class) :: s2
    type(stepsize_class) :: sz
    logical :: first
    real(wp) :: t
    real(wp) :: x(n), x0(n), xf(n), x02(n)
    real(wp),dimension(3) :: r,v
    integer :: n_func_evals

    real(wp),parameter :: rtol = 1.0e-12_wp
    real(wp),parameter :: atol = 1.0e-12_wp
    real(wp),parameter :: t0  = 0.0_wp       !! initial time (sec)
    real(wp),parameter :: dt  = 10.0_wp      !! initial time step (sec)
    real(wp),parameter :: tf  = 10000.0_wp   !! final time (sec)

    !initial conditions:
    real(wp),parameter :: a    = 8000.0_wp ! km
    real(wp),parameter :: ecc  = 0.4_wp
    real(wp),parameter :: inc  = 45.0_wp * deg2rad
    real(wp),parameter :: raan = 45.0_wp * deg2rad
    real(wp),parameter :: aop  = 45.0_wp * deg2rad
    real(wp),parameter :: tru  = 45.0_wp * deg2rad
    real(wp),parameter :: p = a * (1.0_wp - ecc**2)

    call orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)
    x0 = [r,v]

    write(*,*) ''
    write(*,*) 'independent rk78'

    first = .true.
    n_func_evals = 0
    t = t0
    x = x0
    call rk78(der,t,tf,n,x,dt,rtol,atol) ! forward
    call rk78(der,t,t0,n,x,dt,rtol,atol) ! reverse
    write(*,'(i5,1x,*(d15.6,1X))') n_func_evals,x-x0

    !-----------------------------------

    ! step size constructor:
    call sz%initialize( hmin              = 1.0e-6_wp,    &
                        hmax              = 1.0e+6_wp,    &
                        hfactor_reject    = 0.5_wp,       &
                        hfactor_accept    = 2.0_wp,       &
                        max_attempts      = 10000,        &
                        accept_mode       = 2,            &
                        norm              = maxval_func,  &
                        relative_err      = .false.,      &
                        safety_factor     = 0.7_wp,       &
                        p_exponent_offset = 1             )

    write(*,*) ''
    write(*,*) 'rkf78 from this library'

    ! integrator constructor:
    call s%initialize(n=n,f=twobody,rtol=[rtol],atol=[atol],stepsize_method=sz)

    ! integrate:
    first = .true.
    n_func_evals = 0
    call s%integrate(t0,x0,dt,tf,xf)     !forward
    call s%integrate(tf,xf,dt,t0,x02)    !reverse
    write(*,'(i5,1x,*(d15.6,1X))') n_func_evals,x02-x0

    write(*,*) ''
    write(*,*) 'rkv78 from this library'

    ! integrator constructor:
    call s2%initialize(n=n,f=twobody,rtol=[rtol],atol=[atol],stepsize_method=sz)

    ! integrate:
    first = .true.
    n_func_evals = 0
    call s2%integrate(t0,x0,dt,tf,xf)     !forward
    call s2%integrate(tf,xf,dt,t0,x02)    !reverse
    write(*,'(i5,1x,*(d15.6,1X))') n_func_evals,x02-x0

    write(*,*) ''

    contains

    !*********************************************************
    subroutine twobody(me,t,x,xdot)

        !! derivative routine for [[rkf78_class]]

        class(rk_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot

        call der(t,x,xdot)

        end subroutine twobody
    !*********************************************************

    !*********************************************************
    subroutine der(t,x,xdot)

        !! derivative routine for two-body orbit propagation

        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(*),intent(out) :: xdot

        real(wp),dimension(3) :: r,v,a_grav
        real(wp) :: rmag

        r = x(1:3)
        v = x(4:6)
        rmag = norm2(r)
        a_grav = -mu/rmag**3 * r !acceleration due to gravity

        xdot(1:3) = v
        xdot(4:6) = a_grav

        n_func_evals = n_func_evals + 1

        end subroutine der
    !*********************************************************

    end program rk78_validation
!*****************************************************************************************