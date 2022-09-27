
!*****************************************************************************************
!>
!  Unit test of the [[rk_module]].
!  Integrate a two-body orbit around the Earth.

    program rk_test

    use runge_kutta_module
    use rk_kind_module
    use rk_numbers_module

    type,extends(rk4_class) :: spacecraft
        !! spacecraft propagation type.
        !! extends the [[rk4_class]] to include data used in the deriv routine
        real(wp) :: mu = zero      !! central body gravitational parameter (km3/s2)
        integer :: fevals = 0      !! number of function evaluations
        logical :: first = .true.  !! first point is being exported
    end type spacecraft

    integer,parameter :: n=6    !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp !! event location tolerance

    type(spacecraft) :: s, s2
    real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n),gf,tf_actual

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' rk_test'
    write(*,*) '---------------'
    write(*,*) ''

    !***************************************************************************

    !constructor (main body is Earth):
    s = spacecraft(n=n,f=twobody,mu=398600.436233_wp,report=twobody_report)

    !initial conditions:
    x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
            1.0_wp,2.0_wp,3.0_wp]
    t0 = zero       !initial time (sec)
    dt = 10.0_wp    !time step (sec)
    tf = 1000.0_wp  !final time (sec)

    s%fevals = 0
    s%first = .true.
    call s%integrate(t0,x0,dt,tf,xf)    !forward
    write(*,*) ''
    write(*,'(A/,*(F15.6/))') 'Final state:',xf

    s%fevals = 0
    !s%report => null()    !disable reporting
    call s%integrate(tf,xf,-dt,t0,x02)  !backwards

    write(*,'(A/,*(E20.12/))') 'Error:',x02-x0
    write(*,'(A,I5)') 'Function evaluations:', s%fevals
    write(*,*) ''

    !***************************************************************************
    !event finding test:

    write(*,*) ' Event test - integrate until z = 12,000'
    s2 = spacecraft(n=n,f=twobody,g=twobody_event,mu=398600.436233_wp,report=twobody_report)
    x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
            1.0_wp,2.0_wp,3.0_wp]
    t0 = zero       !initial time (sec)
    dt = 10.0_wp    !time step (sec)
    tf = 1000.0_wp  !final time (sec)
    call s2%integrate_to_event(t0,x0,dt,tf,tol,tf_actual,xf,gf)
    write(*,*) ''
    write(*,'(A/,*(F15.6/))') 'Final time: ',tf_actual
    write(*,'(A/,*(F15.6/))') 'Final state:',xf
    write(*,'(A/,*(F15.6/))') 'Event func :',gf

    contains
!*****************************************************************************************

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

        select type (me)
        class is (spacecraft)

            r = x(1:3)
            v = x(4:6)
            rmag = norm2(r)
            a_grav = -me%mu/rmag**3 * r !acceleration due to gravity

            xdot(1:3) = v
            xdot(4:6) = a_grav

            me%fevals = me%fevals + 1

        end select

        end subroutine twobody
    !*********************************************************

    !*********************************************************
        subroutine twobody_report(me,t,x)

        !! report function - write time,state to console

        implicit none

        class(rk_class),intent(inout)        :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(:),intent(in)     :: x

        select type (me)
        class is (spacecraft)
            if (me%first) then  !print header
                write(*,*) ''
                write(*,'(*(A15,1X))')  'time (sec)','x (km)','y (km)','z (km)',&
                                        'vx (km/s)','vy (km/s)','vz (km/s)'
                me%first = .false.
            end if
        end select

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
