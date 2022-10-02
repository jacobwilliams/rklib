
!*****************************************************************************************
!>
!  Unit test of the [[rk_module]].
!  Integrate a two-body orbit around the Earth.

    program rk_test

    use runge_kutta_module, wp => rk_module_rk

    real(wp),parameter :: mu = 398600.436233_wp  !! central body gravitational parameter (km3/s2) - Earth
    integer,parameter :: n = 6  !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp  !! event location tolerance

    class(rk_class),allocatable :: s, s2
    integer :: fevals   !! number of function evaluations
    logical :: first    !! first point is being exported
    real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n),gf,tf_actual
    character(len=:),allocatable :: method !! method name

    ! test all the methods:
    allocate(euler_class :: s);    method='euler';    allocate(euler_class :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(midpoint_class :: s); method='midpoint'; allocate(midpoint_class :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(heun_class :: s);     method='heun';     allocate(heun_class :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(rk4_class :: s);      method='rk4';      allocate(rk4_class :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(rk7_class :: s);      method='rk7';      allocate(rk7_class :: s2); call run_test(); deallocate(s); deallocate(s2)
    allocate(rk8_10_class :: s);   method='rk8_10';   allocate(rk8_10_class :: s2); call run_test(); deallocate(s); deallocate(s2)

    contains
!*****************************************************************************************

    subroutine run_test()

        write(*,*) ''
        write(*,*) '---------------'
        write(*,*) ' rk_test: '//method
        write(*,*) '---------------'
        write(*,*) ''

        !***************************************************************************

        select type (s)
        class is (rk_fixed_step_class)

            !constructor (main body is Earth):
            call s%initialize(n=n,f=twobody,report=twobody_report)

            !initial conditions:
            x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
                    1.0_wp,2.0_wp,3.0_wp]
            t0 = 0.0_wp     !initial time (sec)
            dt = 10.0_wp    !time step (sec)
            tf = 1000.0_wp  !final time (sec)

            fevals = 0
            first = .true.
            call s%integrate(t0,x0,dt,tf,xf)    !forward
            write(*,*) ''
            write(*,'(A/,*(F15.6/))') 'Final state:',xf

            fevals = 0
            !s%report => null()    !disable reporting
            call s%integrate(tf,xf,-dt,t0,x02)  !backwards

        end select

        write(*,'(A/,*(E20.12/))') 'Error:',x02-x0
        write(*,'(A,I5)') 'Function evaluations:', fevals
        write(*,*) ''

        !***************************************************************************
        !event finding test:
        select type (s2)
        class is (rk_fixed_step_class)
            write(*,*) ' Event test - integrate until z = 12,000'
            call s2%initialize(n=n,f=twobody,g=twobody_event,report=twobody_report)
            x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
                    1.0_wp,2.0_wp,3.0_wp]
            t0 = 0.0_wp       !initial time (sec)
            dt = 10.0_wp    !time step (sec)
            tf = 1000.0_wp  !final time (sec)
            call s2%integrate_to_event(t0,x0,dt,tf,tol,tf_actual,xf,gf)
        end select
        write(*,*) ''
        write(*,'(A/,*(F15.6/))') 'Final time: ',tf_actual
        write(*,'(A/,*(F15.6/))') 'Final state:',xf
        write(*,'(A/,*(F15.6/))') 'Event func :',gf

    end subroutine run_test

    !*********************************************************
        subroutine twobody(me,t,x,xdot)

        !! derivative routine for two-body orbit propagation

        implicit none

        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(:),intent(in)     :: x
        real(wp),dimension(:),intent(out)    :: xdot

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

        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(:),intent(in)     :: x

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

        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(:),intent(in)     :: x
        real(wp),intent(out)                 :: g

        g = x(3) - 12000.0_wp

        end subroutine twobody_event
    !*********************************************************

    end program rk_test
!*****************************************************************************************
