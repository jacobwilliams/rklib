
!*****************************************************************************************
!>
!  Unit test of the [[rk_module]].
!  Integrate a two-body orbit around the Earth.

    program rk_test_variable_step

    use rklib_module, wp => rk_module_rk
    use test_support
    use pyplot_module

    implicit none

    real(wp),parameter :: mu = 3.9860043543609593E+05_wp  !! central body gravitational parameter (km^3/s^2) - Earth
    integer,parameter  :: n = 6            !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp !! event location tolerance

    class(rk_class),allocatable :: s, s2
    integer :: fevals
    logical :: first

    type(pyplot) :: plt, plt2
    integer :: istat
    character(len=3) :: rstr

    integer,parameter :: font_size = 30
    integer,parameter :: legend_fontsize = 20

    ! initialize plot
    call plt%initialize(grid=.true.,xlabel='Relative Error',&
                        ylabel='Number of Function Evaluations',&
                        figsize=[30,15],font_size=font_size,axes_labelsize=font_size,&
                        xtick_labelsize=font_size, ytick_labelsize=font_size,&
                        legend_fontsize=legend_fontsize,&
                        title='Variable-Step Runge Kutta Methods',legend=.true.)
    call plt2%initialize(grid=.true.,xlabel='Relative Error',&
                        ylabel='Number of Function Evaluations',&
                        figsize=[30,15],font_size=font_size,axes_labelsize=font_size,&
                        xtick_labelsize=font_size, ytick_labelsize=font_size,&
                        legend_fontsize=legend_fontsize,&
                        title='Variable-Step Runge Kutta Methods [Fixed-Step Mode]',legend=.true.)

    ! test all the methods:
    allocate(rkbs32_class  :: s);  allocate(s2, source=s); call run_all_tests('rkbs32',  [255,0,0]);  call finish()

    allocate(rkf45_class   :: s);  allocate(s2, source=s); call run_all_tests('rkf45',   [245, 152, 152]);  call finish()

    allocate(rkck54_class  :: s);  allocate(s2, source=s); call run_all_tests('rkck54',  [255, 102, 0]); call finish()
    allocate(rkdp54_class  :: s);  allocate(s2, source=s); call run_all_tests('rkdp54',  [189, 90, 25]); call finish()
    allocate(rkt54_class   :: s);  allocate(s2, source=s); call run_all_tests('rkt54',   [143, 78, 36]); call finish()
    allocate(rks54_class   :: s);  allocate(s2, source=s); call run_all_tests('rks54',   [243, 78, 36]); call finish()

    allocate(rkdp65_class  :: s);  allocate(s2, source=s); call run_all_tests('rkdp65',  [251, 255, 0]);  call finish()
    allocate(rkc65_class   :: s);  allocate(s2, source=s); call run_all_tests('rkc65',   [207, 194, 145]);call finish()
    allocate(rktp64_class  :: s);  allocate(s2, source=s); call run_all_tests('rktp64',  [187, 189, 49]); call finish()
    allocate(rkv65e_class  :: s);  allocate(s2, source=s); call run_all_tests('rkv65e',  [149, 150, 63]); call finish()
    allocate(rkv65r_class  :: s);  allocate(s2, source=s); call run_all_tests('rkv65r',  [65, 71, 41]); call finish()
    allocate(rkv65_class   :: s);  allocate(s2, source=s); call run_all_tests('rkv65',   [75, 81, 51]); call finish()
    allocate(rktf65_class  :: s);  allocate(s2, source=s); call run_all_tests('rktf65',  [2, 2, 2]); call finish()

    allocate(rktp75_class  :: s);  allocate(s2, source=s); call run_all_tests('rktp75',  [0, 255, 38]);    call finish()
    allocate(rktmy7_class  :: s);  allocate(s2, source=s); call run_all_tests('rktmy7',  [102, 247, 255]); call finish()
    allocate(rkv76e_class  :: s);  allocate(s2, source=s); call run_all_tests('rkv76e',  [38, 189, 60]);   call finish()
    allocate(rkv76r_class  :: s);  allocate(s2, source=s); call run_all_tests('rkv76r',  [149, 163, 93]);   call finish()
    allocate(rkf78_class   :: s);  allocate(s2, source=s); call run_all_tests('rkf78',   [66, 143, 77]);   call finish()
    allocate(rkv78_class   :: s);  allocate(s2, source=s); call run_all_tests('rkv78',   [77, 105, 81]);   call finish()

    allocate(rktp86_class  :: s);  allocate(s2, source=s); call run_all_tests('rktp86',  [0, 47, 255]);    call finish()
    allocate(rkdp87_class  :: s);  allocate(s2, source=s); call run_all_tests('rkdp87',  [51, 83, 222]);   call finish()
    allocate(rkv87e_class  :: s);  allocate(s2, source=s); call run_all_tests('rkv87e',  [90, 116, 230]);  call finish()
    allocate(rkv87r_class  :: s);  allocate(s2, source=s); call run_all_tests('rkv87r',  [0,0,0]);  call finish()
    allocate(rkk87_class  :: s);  allocate(s2, source=s); call run_all_tests('rkk87',  [10,10,10]);  call finish()
    allocate(rkf89_class   :: s);  allocate(s2, source=s); call run_all_tests('rkf89',   [116, 133, 207]); call finish()
    allocate(rkv89_class   :: s);  allocate(s2, source=s); call run_all_tests('rkv89',   [169, 176, 219]); call finish()

    allocate(rkt98a_class  :: s);  allocate(s2, source=s); call run_all_tests('rkt98a',  [195, 0, 255]);   call finish()
    allocate(rkv98e_class  :: s);  allocate(s2, source=s); call run_all_tests('rkv98e',  [192, 52, 235]);  call finish()
    allocate(rkv98r_class  :: s);  allocate(s2, source=s); call run_all_tests('rkv98r',  [79, 5, 153]);    call finish()
    allocate(rkf108_class  :: s);  allocate(s2, source=s); call run_all_tests('rkf108',  [198, 149, 245]);  call finish()
    allocate(rkc108_class  :: s);  allocate(s2, source=s); call run_all_tests('rkc108',  [232, 207, 255]); call finish()
    allocate(rkb109_class  :: s);  allocate(s2, source=s); call run_all_tests('rkb109',  [200, 180, 230]); call finish()

    allocate(rks1110a_class :: s); allocate(s2, source=s); call run_all_tests('rks1110a',[0,0,0]);          call finish()
    allocate(rkf1210_class :: s);  allocate(s2, source=s); call run_all_tests('rkf1210', [94,94,94]);       call finish()
    allocate(rko129_class :: s);   allocate(s2, source=s); call run_all_tests('rko129',  [145, 145, 145]);  call finish()
    allocate(rkf1412_class :: s);  allocate(s2, source=s); call run_all_tests('rkf1412', [225, 230, 230]);  call finish()

    ! save plot:
    write(rstr,'(I3)') wp
    call plt%savefig(figfile='rk_test_variable_step_R'//trim(adjustl(rstr))//'.png',istat=istat)
    call plt2%savefig(figfile='rk_test_variable_step_FIXED_MODE_R'//trim(adjustl(rstr))//'.png',istat=istat)

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
        call performance_test_fixed(method,color)
        call run_test(method)
    end subroutine run_all_tests

    subroutine performance_test(method,color)
        !! generate a performance plot for all the methods
        character(len=*),intent(in) :: method !! name of the RK method to use
        integer,dimension(3),intent(in) :: color !! color for the plot

        integer,parameter :: exp_start = 8
        integer,parameter :: exp_stop = 15
        integer,parameter :: factor = 3

        integer :: i !! counter
        real(wp),dimension(factor*exp_start:factor*exp_stop) :: r_error, v_error
        integer,dimension(factor*exp_start:factor*exp_stop) :: feval
        real(wp) :: xerror(n)
        real(wp) :: rtol, atol

        integer :: icase
        integer :: p_exponent_offset
        logical :: relative_err
        real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n),gf,tf_actual
        real(wp) :: safety_factor, hfactor_accept
        real(wp) :: a,p,ecc,inc,raan,aop,tru
        real(wp),dimension(3) :: r,v
        type(stepsize_class) :: sz

        !initial conditions:
        x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
              1.0_wp,2.0_wp,3.0_wp]
        t0 = 0.0_wp      ! initial time (sec)
        dt = 10.0_wp     ! initial time step (sec)
        tf = 10000.0_wp  ! final time (sec)

        do i = exp_start*factor, exp_stop*factor

            rtol = 10.0_wp ** (-i/real(factor,wp))
            atol = 10.0_wp ** (-i/real(factor,wp))

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

            select type (s)
            class is (rk_variable_step_class)

                ! integrator constructor:
                call s%initialize(n=n,f=twobody,rtol=[rtol],atol=[atol],stepsize_method=sz)

                ! integrate:
                first = .true.
                fevals = 0
                call s%integrate(t0,x0,dt,tf,xf)     !forward
                if (s%failed()) error stop 'integration failed'
                feval(i) = fevals
                call s%integrate(tf,xf,dt,t0,x02)    !reverse
                !write(*,'(i5,1x,*(d15.6,1X))') n_func_evals,x02-x0

            end select

            ! compute relative error:
            where (x0 /= 0.0_wp)
                xerror = (x02-x0)/(x0)
            else where
                xerror = (x02-x0)
            end where
            r_error(i) = norm2(xerror(1:3))
            v_error(i) = norm2(xerror(4:6))

        end do

        ! add to the plot:
        call plt%add_plot(r_error,real(feval,wp),&
                            label=method,&
                            linestyle='.-',color=real(color/255.0_wp,wp),&
                            markersize=5,linewidth=4,istat=istat,&
                            xscale='log',yscale='log')

    end subroutine performance_test

    subroutine performance_test_fixed(method,color)
        !! generate a performance plot for all the methods
        character(len=*),intent(in) :: method !! name of the RK method to use
        integer,dimension(3),intent(in) :: color !! color for the plot

        integer,parameter :: factor = 1
        integer,parameter :: n_cases = factor*10000  !! used for `dt`

        integer :: i !! counter
        real(wp),dimension(n_cases) :: r_error, v_error
        integer,dimension(n_cases) :: feval
        real(wp) :: xerror(n)

        integer :: icase
        logical :: relative_err
        real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n),gf,tf_actual
        real(wp) :: safety_factor, hfactor_accept
        real(wp) :: a,p,ecc,inc,raan,aop,tru
        real(wp),dimension(3) :: r,v
        type(stepsize_class) :: sz

        !initial conditions:
        x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
              1.0_wp,2.0_wp,3.0_wp]
        t0 = 0.0_wp           ! initial time (sec)
        tf = real(n_cases,wp) ! final time (sec)

        do i = 1, n_cases

            dt = real(i,wp)/factor ! time step (sec)

            ! step size constructor:
            call sz%initialize(fixed_step_mode = .true.)

            select type (s)
            class is (rk_variable_step_class)

                ! integrator constructor:
                ! [rtol,atol not used in fixed step mode]
                call s%initialize(n=n,f=twobody,rtol=[1.0_wp],atol=[1.0_wp],stepsize_method=sz)

                ! integrate:
                first = .true.
                fevals = 0
                call s%integrate(t0,x0,dt,tf,xf)     !forward
                feval(i) = fevals
                call s%integrate(tf,xf,dt,t0,x02)    !reverse
                !write(*,'(i5,1x,*(d15.6,1X))') n_func_evals,x02-x0

            end select

            ! compute relative error:
            where (x0 /= 0.0_wp)
                xerror = (x02-x0)/(x0)
            else where
                xerror = (x02-x0)
            end where
            r_error(i) = norm2(xerror(1:3))
            v_error(i) = norm2(xerror(4:6))

        end do

        ! add to the plot:
        call plt2%add_plot(r_error,real(feval,wp),&
                            label=method,&
                            linestyle='.-',color=real(color/255.0_wp,wp),&
                            markersize=5,linewidth=4,istat=istat,&
                            xscale='log',yscale='log')

    end subroutine performance_test_fixed

    subroutine run_test(name)

    character(len=*),intent(in) :: name !! name of the method

    integer :: ierr !! status code
    integer :: icase
    integer :: p_exponent_offset
    logical :: relative_err
    real(wp) :: t0,tf,x0(n),dt,xf(n),x02(n),gf,tf_actual
    real(wp) :: safety_factor, hfactor_accept
    real(wp) :: a,p,ecc,inc,raan,aop,tru
    real(wp),dimension(3) :: r,v
    type(stepsize_class) :: sz

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) name//' : rk_variable_step_test'
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
                              stepsize_method=sz,report=twobody_report,report_rate=10)

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
            t0   = 0.0_wp            ! initial time (sec)
            !dt   = 0.0_wp           ! automatically compute an initial time step (sec)
            dt   = 10.0_wp           ! initial time step (sec)
            tf   = 1000.0_wp         ! final time (sec)

            !s%num_rejected_steps = 0
            fevals = 0
            first = .true.
            call s%integrate(t0,x0,dt,tf,xf)    !forward
            call s%status(ierr)
            write(*,*) ''
            write(*,*) 'ierr = ', ierr
            write(*,'(A/,*(F15.6/))') 'Final state:',xf
            write(*,'(A,I5)') 'Function evaluations:', fevals

            !s%num_rejected_steps = 0
            fevals = 0
            !s%report => null()    !disable reporting
            call s%integrate(tf,xf,-dt,t0,x02)  !backwards
            call s%status(ierr)

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
                           stepsize_method=sz)

        fevals = 0
        first = .true.
        x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
              1.0_wp,2.0_wp,3.0_wp]
        t0 = 0.0_wp     !initial time (sec)
        dt = 10.0_wp    !time step (sec)
        tf = 1000.0_wp  !final time (sec)

        call s2%integrate_to_event(t0,x0,dt,tf,tol,tf_actual,xf,gf)
        call s2%status(ierr)

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
