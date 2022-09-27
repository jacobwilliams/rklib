
!*****************************************************************************************
!>
!  Unit test of the [[rk_module]].
!  Integrate a two-body orbit around the Earth.

    program rk_test_variable_step

    use runge_kutta_module
    use rk_kind_module
    use rk_numbers_module

    implicit none

    real(wp),parameter :: deg2rad = acos(-1.0_wp) / 180.0_wp

    !type,extends(rkf78_class) :: spacecraft
    !type,extends(rkf89_class) :: spacecraft
    !type,extends(rkv89_class) :: spacecraft
    !type,extends(rkf108_class) :: spacecraft
    !type,extends(rkf1210_class) :: spacecraft
    type,extends(rkf1412_class) :: spacecraft
        !! spacecraft propagation type.
        real(wp) :: mu     = zero     !! central body gravitational parameter (km3/s2)
        integer  :: fevals = 0        !! number of function evaluations
        logical  :: first  = .true.   !! first point is being exported
    end type spacecraft

    integer,parameter  :: n = 6            !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp !! event location tolerance

    type(spacecraft) :: s, s2
    real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n),gf,tf_actual,rtol,atol
    integer :: ierr !! error flag
    type(stepsize_class) :: sz
    integer :: icase
    logical :: relative_err
    real(wp) :: safety_factor, hfactor_accept
    integer :: p_exponent_offset
    real(wp) :: mu
    real(wp) :: a,p,ecc,inc,raan,aop,tru
    real(wp),dimension(3) :: r,v

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' rk_variable_step_test'
    write(*,*) '---------------'
    write(*,*) ''

    !***************************************************************************

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

        !integrator constructor:
        call s%initialize(n=n,f=twobody,rtol=[1.0e-15_wp],atol=[1.0e-12_wp],&
                            stepsize_method=sz,report=twobody_report)

        !initial conditions:
        !write(*,*) 'general elliptical:'
        mu = 3.9860043543609593E+05_wp ! for earth
        a    = 8000.0_wp ! km
        ecc  = 0.1_wp
        inc  = 45.0_wp * deg2rad
        raan = 45.0_wp * deg2rad
        aop  = 45.0_wp * deg2rad
        tru  = 45.0_wp * deg2rad
        p = a * (one - ecc**2)

        call orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)
        x0 = [r,v]

        !x0   = [10000.0_wp,10000.0_wp,10000.0_wp,&   ! initial state [r,v] (km,km/s)
        !        1.0_wp,2.0_wp,3.0_wp]
        t0   = zero              ! initial time (sec)
        !dt   = 0.0_wp           ! automatically compute an initial time step (sec)
        dt   = 10.0_wp           ! initial time step (sec)
        tf   = 10000.0_wp         ! final time (sec)
        s%mu = 398600.436233_wp  ! main body is Earth

        !s%num_rejected_steps = 0
        s%fevals = 0
        s%first = .true.
        call s%integrate(t0,x0,dt,tf,xf,ierr)    !forward
        write(*,*) ''
        write(*,*) 'ierr = ', ierr
        write(*,'(A/,*(F15.6/))') 'Final state:',xf
        write(*,'(A,I5)') 'Function evaluations:', s%fevals
        !write(*,'(A,I5)') 'Number of rejected steps:',s%num_rejected_steps   ! why is this 0 when ierr = -3 ???

        !s%num_rejected_steps = 0
        s%fevals = 0
        !s%report => null()    !disable reporting
        call s%integrate(tf,xf,-dt,t0,x02,ierr)  !backwards

        write(*,*) 'ierr = ', ierr
        write(*,'(A/,*(E20.12/))') 'Error:',x02-x0
        write(*,'(A,I5)') 'Function evaluations:', s%fevals
        !write(*,'(A,I5)') 'Number of rejected steps:',s%num_rejected_steps
        write(*,*) ''

    end do

    !***************************************************************************
    !event finding test:

    write(*,*) ' Event test - integrate until z = 12,000'

   ! NOTE: the following causes an ICE in gfortran 7.1, but works with ifort:
   ! s2 = spacecraft(n=n,f=twobody,g=twobody_event,mu=398600.436233_wp,&
   !                  rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
   !                  stepsize_method=sz,report=twobody_report)
   ! do it this way instead:
   call s2%initialize(n=n,f=twobody,g=twobody_event,&
                      rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
                      stepsize_method=sz,&
                      report=twobody_report)
   s2%mu = 398600.436233_wp

   s2%fevals = 0
   s2%first = .true.
   x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
           1.0_wp,2.0_wp,3.0_wp]
   t0 = zero       !initial time (sec)
   dt = 10.0_wp    !time step (sec)
   tf = 1000.0_wp  !final time (sec)

   call s2%integrate_to_event(t0,x0,dt,tf,tol,tf_actual,xf,gf,ierr)
   write(*,*) ''
   write(*,'(A,I5)')         'ierr:       ',ierr
   write(*,'(A/,*(F15.6/))') 'Final time: ',tf_actual
   write(*,'(A/,*(F15.6/))') 'Final state:',xf
   write(*,'(A/,*(F15.6/))') 'Event func :',gf
   write(*,'(A,I5)') 'Function evaluations:', s2%fevals
   !write(*,'(A,I5)') 'Number of rejected steps:',s2%num_rejected_steps

    contains
!*****************************************************************************************

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

        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x

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

        class(rk_class),intent(inout)    :: me
        real(wp),intent(in)              :: t
        real(wp),dimension(:),intent(in) :: x
        real(wp),intent(out)             :: g

        g = x(3) - 12000.0_wp

        end subroutine twobody_event
    !*********************************************************


!*****************************************************************************************
!>
!  Convert orbital elements to position and velocity vectors.

        pure subroutine orbital_elements_to_rv(mu, p, ecc, inc, raan, aop, tru, r, v)

        implicit none
    
        real(wp),intent(in)               :: mu   !! gravitational parameter [\(km^{3}/s^{2}\)]
        real(wp),intent(in)               :: p    !! semiparameter \(a(1-e^{2})\) [km]
        real(wp),intent(in)               :: ecc  !! eccentricity [--]
        real(wp),intent(in)               :: inc  !! inclination [rad]
        real(wp),intent(in)               :: raan !! raan [rad]
        real(wp),intent(in)               :: aop  !! argument of peripsis [rad]
        real(wp),intent(in)               :: tru  !! true anomaly [rad]
        real(wp),dimension(3),intent(out) :: r    !! position vector [km]
        real(wp),dimension(3),intent(out) :: v    !! velocity vector [km/s]
    
        real(wp),dimension(3,2) :: rotmat
        real(wp),dimension(2)   :: r_pqw,v_pqw
        logical                 :: circular,equatorial
        real(wp)                :: ctru,stru,sr,cr,si,ci,sa,ca,raan_tmp,aop_tmp
    
        call orbit_check(ecc,inc,circular,equatorial)
    
        if (circular) then   ! periapsis undefined
            aop_tmp = zero
        else
            aop_tmp = aop
        end if
    
        if (equatorial) then   ! node undefined
            raan_tmp = zero
        else
            raan_tmp = raan
        end if
    
        ! perifocal position and velocity:
        ctru   = cos(tru)
        stru   = sin(tru)
        r_pqw  = [ctru, stru] * p/(one+ecc*ctru)
        v_pqw  = [-stru, (ecc+ctru)] * sqrt(mu/p)
    
        ! perifocal to cartesian:
        sr          = sin(raan_tmp)
        cr          = cos(raan_tmp)
        si          = sin(inc)
        ci          = cos(inc)
        sa          = sin(aop_tmp)
        ca          = cos(aop_tmp)
        rotmat(1,:) = [cr*ca-sr*sa*ci, -cr*sa-sr*ca*ci]
        rotmat(2,:) = [sr*ca+cr*sa*ci, -sr*sa+cr*ca*ci]
        rotmat(3,:) = [sa*si, ca*si]
    
        ! transform:
        r = matmul(rotmat,r_pqw)
        v = matmul(rotmat,v_pqw)
    
        end subroutine orbital_elements_to_rv
    !*****************************************************************************************
    
!*****************************************************************************************
!>
!  Check the orbit for singularities.

        pure subroutine orbit_check(ecc,inc,circular,equatorial)

        implicit none
    
        real(wp),intent(in) :: ecc        !! eccentricity
        real(wp),intent(in) :: inc        !! inclination [rad]
        logical,intent(out) :: circular   !! is the orbit circular?
        logical,intent(out) :: equatorial !! is the orbit equatorial?
    
        real(wp),parameter :: tol = 1.0e-10_wp !! tolerance for circular & equatorial checks
    
        circular   = ecc < tol
        equatorial = (one - abs(cos(inc))) < tol  ! 0 or 180 deg
    
        end subroutine orbit_check
    !*****************************************************************************************
    
    end program rk_test_variable_step
!*****************************************************************************************
