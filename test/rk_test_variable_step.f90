
!*****************************************************************************************
!>
!  Unit test of the [[rk_module]].
!  Integrate a two-body orbit around the Earth.

    program rk_test_variable_step

    use runge_kutta_module, wp => rk_module_rk
    use test_support

    implicit none

    real(wp),parameter :: mu = 3.9860043543609593E+05_wp  !! central body gravitational parameter (km^3/s^2) - Earth
    integer,parameter  :: n = 6            !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp !! event location tolerance

    class(rk_class),allocatable :: s, s2
    type(stepsize_class) :: sz
    integer :: fevals
    integer :: ierr !! error flag
    integer :: icase
    integer :: p_exponent_offset
    logical :: relative_err
    real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n),gf,tf_actual,rtol,atol
    real(wp) :: safety_factor, hfactor_accept
    real(wp) :: a,p,ecc,inc,raan,aop,tru
    real(wp),dimension(3) :: r,v
    logical :: first

    ! test all the methods:
    allocate(rkf78_class   :: s);  allocate(rkf78_class   :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(rkf89_class   :: s);  allocate(rkf89_class   :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(rkv89_class   :: s);  allocate(rkv89_class   :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(rkf108_class  :: s);  allocate(rkf108_class  :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(rkf1210_class :: s);  allocate(rkf1210_class :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(rkf1412_class :: s);  allocate(rkf1412_class :: s2); call run_test(); deallocate(s); deallocate(s2)

    contains
!*****************************************************************************************

    subroutine run_test()

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' rk_variable_step_test'
    write(*,*) '---------------'
    write(*,*) ''

    do icase = 1, 4

        write(*,*) ''
        write(*,*) '***************'
        write(*,*) ' case' , icase
        write(*,*) '***************'
        write(*,*) ''

        ! ... relative_err,safety_factor,p_exponent_offset don't
        !     seem to change the results at all ....

        select case (icase)
        case(1)
            ! defaults
            relative_err = .false.
            safety_factor = 0.8_wp
            p_exponent_offset = 1
            hfactor_accept = 2.0_wp    ! changing this does change result
        case(2)
            relative_err = .false.
            safety_factor = 0.9_wp
            p_exponent_offset = 1
            hfactor_accept = 5.0_wp
        case(3)
            relative_err = .false.
            safety_factor = 0.95_wp
            p_exponent_offset = 1
            hfactor_accept = 10.0_wp
        case(4)
            relative_err = .false.
            safety_factor = 0.9_wp
            p_exponent_offset = 0
            hfactor_accept = 2.0_wp
        end select

        !step size constructor:
        !call sz%initialize(hmin=1.0e-6_wp,hmax=1.0e6_wp)
        call sz%initialize( hmin              = 1.0e-6_wp,        &
                            hmax              = 1.0e+6_wp,        &
                            hfactor_reject    = 0.5_wp,           &
                            hfactor_accept    = hfactor_accept,   &
                            max_attempts      = 1000,              &
                            accept_mode       = 2,                &
                            norm              = maxval_func,      &
                            relative_err      = relative_err,     &
                            safety_factor     = safety_factor,    &
                            p_exponent_offset = p_exponent_offset )

        select type (s)
        class is (rk_variable_step_class)

            !integrator constructor:
            call s%initialize(n=n,f=twobody,rtol=[1.0e-15_wp],atol=[1.0e-12_wp],&
                              stepsize_method=sz,report=twobody_report)

            !initial conditions:
            !write(*,*) 'general elliptical:'
            a    = 8000.0_wp ! km
            ecc  = 0.1_wp
            inc  = 45.0_wp * deg2rad
            raan = 45.0_wp * deg2rad
            aop  = 45.0_wp * deg2rad
            tru  = 45.0_wp * deg2rad
            p = a * (1.0_wp - ecc**2)

            call orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)
            x0 = [r,v]

            !x0   = [10000.0_wp,10000.0_wp,10000.0_wp,&   ! initial state [r,v] (km,km/s)
            !        1.0_wp,2.0_wp,3.0_wp]
            t0   = 0.0_wp              ! initial time (sec)
            !dt   = 0.0_wp           ! automatically compute an initial time step (sec)
            dt   = 10.0_wp           ! initial time step (sec)
            tf   = 10000.0_wp         ! final time (sec)

            !s%num_rejected_steps = 0
            fevals = 0
            first = .true.
            call s%integrate(t0,x0,dt,tf,xf,ierr)    !forward
            write(*,*) ''
            write(*,*) 'ierr = ', ierr
            write(*,'(A/,*(F15.6/))') 'Final state:',xf
            write(*,'(A,I5)') 'Function evaluations:', fevals
            !write(*,'(A,I5)') 'Number of rejected steps:',s%num_rejected_steps   ! why is this 0 when ierr = -3 ???

            !s%num_rejected_steps = 0
            fevals = 0
            !s%report => null()    !disable reporting
            call s%integrate(tf,xf,-dt,t0,x02,ierr)  !backwards

            write(*,*) 'ierr = ', ierr
            write(*,'(A/,*(E20.12/))') 'Error:',x02-x0
            write(*,'(A,I5)') 'Function evaluations:', fevals
            !write(*,'(A,I5)') 'Number of rejected steps:',s%num_rejected_steps
            write(*,*) ''

        end select

    end do

    !***************************************************************************
    !event finding test:

    write(*,*) ' Event test - integrate until z = 12,000'

    select type (s2)
    class is (rk_variable_step_class)

        call s2%initialize(n=n,f=twobody,g=twobody_event,&
                           rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
                           stepsize_method=sz,&
                           report=twobody_report)

        fevals = 0
        first = .true.
        x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
              1.0_wp,2.0_wp,3.0_wp]
        t0 = 0.0_wp       !initial time (sec)
        dt = 10.0_wp    !time step (sec)
        tf = 1000.0_wp  !final time (sec)

        call s2%integrate_to_event(t0,x0,dt,tf,tol,tf_actual,xf,gf,ierr)
        write(*,*) ''
        write(*,'(A,I5)')         'ierr:       ',ierr
        write(*,'(A/,*(F15.6/))') 'Final time: ',tf_actual
        write(*,'(A/,*(F15.6/))') 'Final state:',xf
        write(*,'(A/,*(F15.6/))') 'Event func :',gf
        write(*,'(A,I5)') 'Function evaluations:', fevals
        !write(*,'(A,I5)') 'Number of rejected steps:',s2%num_rejected_steps

    end select

    end subroutine run_test
    !*********************************************************

    !*********************************************************
        subroutine twobody(me,t,x,xdot)

        !! derivative routine for two-body orbit propagation

        implicit none

        class(rk_class),intent(inout)     :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot

        real(wp),dimension(3) :: r,v,a_grav
        real(wp) :: rmag

        r = x(1:3)
        v = x(4:6)
        rmag = norm2(r)
        a_grav = -mu/rmag**3 * r !acceleration due to gravity

        xdot(1:3) = v
        xdot(4:6) = a_grav

        fevals = fevals + 1

        end subroutine twobody
    !*********************************************************

    !*********************************************************
        subroutine twobody_report(me,t,x)

        !! report function - write time,state to console

        implicit none

        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x

        if (first) then  !print header
            write(*,*) ''
            write(*,'(*(A15,1X))')  'time (sec)','x (km)','y (km)','z (km)',&
                                    'vx (km/s)','vy (km/s)','vz (km/s)'
            first = .false.
        end if

        write(*,'(*(F15.6,1X))') t,x

        end subroutine twobody_report
    !*********************************************************

    !*********************************************************
        subroutine twobody_event(me,t,x,g)

        !! event function (z = 12,000)

        implicit none

        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x
        real(wp),intent(out)             :: g

        g = x(3) - 12000.0_wp

        end subroutine twobody_event
    !*********************************************************

    end program rk_test_variable_step
!*****************************************************************************************
