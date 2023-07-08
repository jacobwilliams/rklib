
!*****************************************************************************************
!>
!  Unit test of the [[rk_module]].
!  Integrate a two-body orbit around the Earth.

    program rk_test

    use rklib_module, wp => rk_module_rk
    use pyplot_module

    real(wp),parameter :: mu = 398600.436233_wp  !! central body gravitational parameter (km3/s2) - Earth
    integer,parameter :: n = 6  !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp  !! event location tolerance

    class(rk_class),allocatable :: s, s2
    integer :: fevals   !! number of function evaluations
    logical :: first    !! first point is being exported
    real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n),gf,tf_actual
    type(pyplot) :: plt
    integer :: istat
    character(len=3) :: rstr

    ! initialize plot
    call plt%initialize(grid=.true.,xlabel='Relative Error',&
                        ylabel='Number of Function Evaluations',&
                        figsize=[20,10],font_size=20,axes_labelsize=20,&
                        xtick_labelsize=20, ytick_labelsize=20,&
                        legend_fontsize=20,&
                        title='Fixed-Step Runge Kutta Methods',legend=.true.)

    ! test all the methods:
    allocate(euler_class :: s);    allocate(s2, source=s); call run_all_tests('euler'   ,[255,0,0]);      call finish()
    allocate(midpoint_class :: s); allocate(s2, source=s); call run_all_tests('midpoint',[235, 110, 52]); call finish()
    allocate(heun_class :: s);     allocate(s2, source=s); call run_all_tests('heun'    ,[235, 165, 52]); call finish()
    allocate(rkssp22_class :: s);  allocate(s2, source=s); call run_all_tests('rkssp22' ,[235, 195, 52]); call finish()
    allocate(rk3_class :: s);      allocate(s2, source=s); call run_all_tests('rk3'     ,[220, 235, 52]); call finish()
    allocate(rkssp33_class :: s);  allocate(s2, source=s); call run_all_tests('rkssp33' ,[220, 255, 52]); call finish()
    allocate(rk4_class :: s);      allocate(s2, source=s); call run_all_tests('rk4'     ,[0,255,0]);      call finish()
    allocate(rks4_class :: s);     allocate(s2, source=s); call run_all_tests('rks4'    ,[52, 235, 186]); call finish()
    allocate(rkssp54_class :: s);  allocate(s2, source=s); call run_all_tests('rkssp54' ,[52, 220, 210]); call finish()
    allocate(rks5_class :: s);     allocate(s2, source=s); call run_all_tests('rks5'    ,[52, 198, 235]); call finish()
    allocate(rkb6_class :: s);     allocate(s2, source=s); call run_all_tests('rkb6'    ,[0, 0, 0]);      call finish()
    allocate(rk7_class :: s);      allocate(s2, source=s); call run_all_tests('rk7'     ,[52, 64, 235]);  call finish()
    allocate(rk8_10_class :: s);   allocate(s2, source=s); call run_all_tests('rk8_10'  ,[122, 52, 235]); call finish()
    allocate(rk8_12_class :: s);   allocate(s2, source=s); call run_all_tests('rk8_12'  ,[229, 52, 235]); call finish()
    allocate(rkcv8_class :: s);    allocate(s2, source=s); call run_all_tests('rkcv8'   ,[217, 163, 163]);call finish()
    allocate(rkz10_class :: s);    allocate(s2, source=s); call run_all_tests('rkz10'   ,[222, 115, 73]); call finish()
    allocate(rko10_class :: s);    allocate(s2, source=s); call run_all_tests('rko10'   ,[200, 200, 200]); call finish()

    ! save plot:
    write(rstr,'(I3)') wp
    call plt%savefig(figfile='rk_test_R'//trim(adjustl(rstr))//'.png',istat=istat)

    contains
!*****************************************************************************************

    subroutine finish()
        deallocate(s); deallocate(s2)
    end subroutine finish

    subroutine run_all_tests(method,color)
        !! run all the tests
        character(len=*),intent(in) :: method !! name of the RK method to use
        integer,dimension(3),intent(in) :: color !! color for the plot
        call performance_test(method,color)
        call run_test(method)
    end subroutine run_all_tests

    subroutine performance_test(method,color)
        !! generate a performance plot for all the methods
        character(len=*),intent(in) :: method !! name of the RK method to use
        integer,dimension(3),intent(in) :: color !! color for the plot

        integer,parameter :: factor = 1
        integer,parameter :: n_cases = factor*1000  !! used for `dt`

        integer :: i !! counter
        real(wp),dimension(n_cases) :: r_error, v_error
        integer,dimension(n_cases) :: feval
        real(wp) :: xerror(n)

        do i = 1, n_cases

            !initial conditions:
            x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
                  1.0_wp,2.0_wp,3.0_wp]
            t0 = 0.0_wp            ! initial time (sec)
            dt = real(i,wp)/factor ! time step (sec)
            tf = real(n_cases,wp)/factor  ! final time (sec)
            fevals = 0
            first = .true.

            select type (s)
            class is (rk_fixed_step_class)

                !constructor (main body is Earth):
                call s%initialize(n=n,f=twobody)
                call s%integrate(t0,x0,dt,tf,xf)   ! forward
                fevals = 0
                call s%integrate(tf,xf,-dt,t0,x02) ! backwards

            end select

            ! compute relative error:
            where (x0 /= 0.0_wp)
                xerror = (x02-x0)/(x0)
            else where
                xerror = (x02-x0)
            end where
            r_error(i) = norm2(xerror(1:3))
            v_error(i) = norm2(xerror(4:6))
            feval(i) = fevals

        end do

        ! add to the plot:
        call plt%add_plot(r_error,real(feval,wp),&
                            label=method,&
                            linestyle='-',color=real(color/255.0_wp,wp),&
                            markersize=5,linewidth=4,istat=istat,&
                            xscale='log',yscale='log')

    end subroutine performance_test

    subroutine run_test(method)
        !! basic test

        character(len=*),intent(in) :: method !! name of the RK method to use

        write(*,*) ''
        write(*,*) '---------------------------------------------'
        write(*,*) ' rk_test: '//method
        write(*,*) '---------------------------------------------'
        write(*,*) ''

        !***************************************************************************

        select type (s)
        class is (rk_fixed_step_class)

            !constructor (main body is Earth):
            call s%initialize(n=n,f=twobody,report=twobody_report,report_rate=10)

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
            call s2%initialize(n=n,f=twobody,g=twobody_event,report=twobody_report,report_rate=10)
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
