
!*****************************************************************************************
!>
!  Unit test of the [[rk_module]].

    program rk_performance

    use rklib_module, wp => rk_module_rk
    use test_support
    use pyplot_module
    use iso_fortran_env

    implicit none

    real(wp),parameter :: mu = 3.9860043543609593E+05_wp  !! central body gravitational parameter (km^3/s^2) - Earth
    integer,parameter  :: n = 6            !! number of state variables

    class(rk_class),allocatable :: s, s2
    integer :: fevals
    type(pyplot) :: plt
    integer :: istat
    character(len=:),allocatable :: case
    integer :: exp_start
    integer :: exp_stop
    integer :: factor

    integer,parameter :: font_size = 35
    integer,parameter :: legend_fontsize = font_size

    ! defaults:
    case = ''
    exp_start = 8
    exp_stop  = 15
    factor    = 2

    if (wp==real64) then
        case = ' [REAL64]'
    else if (wp==real128) then
        case = ' [REAL128]'
        exp_start = 16
        exp_stop  = 25
        factor    = 2
    end if

    ! initialize plot
    call plt%initialize(grid=.true.,xlabel='Relative Error',&
                        ylabel='Number of Function Evaluations',&
                        figsize=[30,15],font_size=font_size,axes_labelsize=font_size,&
                        xtick_labelsize=font_size, ytick_labelsize=font_size,&
                        legend_fontsize=legend_fontsize,&
                        title='Problem 1: Variable-Step Runge Kutta Methods'//case,legend=.true.)

    ! test all the methods:
    call init()
    allocate(rkbs32_class  :: s);  allocate(s2, source=s); call run('rkbs32',  [0,0,255]);  call finish()
    allocate(rkssp43_class :: s);  allocate(s2, source=s); call run('rkssp43', [255,0,0]);  call finish()
    call done(3)

    call init()
    allocate(rkf45_class   :: s);  allocate(s2, source=s); call run('rkf45',   [0,0,255]);  call finish()
    call done(4)

    call init()
    allocate(rkck54_class  :: s);  allocate(s2, source=s); call run('rkck54',  [255, 102, 0]); call finish()
    allocate(rkdp54_class  :: s);  allocate(s2, source=s); call run('rkdp54',  [189, 90, 25]); call finish()
    allocate(rkt54_class   :: s);  allocate(s2, source=s); call run('rkt54',   [143, 78, 36]); call finish()
    allocate(rks54_class   :: s);  allocate(s2, source=s); call run('rks54',   [243, 78, 36],':'); call finish()
    call done(5)

    call init()
    allocate(rkdp65_class  :: s);  allocate(s2, source=s); call run('rkdp65',  [251, 255, 0]);  call finish()
    allocate(rkc65_class   :: s);  allocate(s2, source=s); call run('rkc65',   [207, 0, 0]);call finish()
    allocate(rktp64_class  :: s);  allocate(s2, source=s); call run('rktp64',  [187, 189, 49]); call finish()
    allocate(rkv65e_class  :: s);  allocate(s2, source=s); call run('rkv65e',  [0,0,255]); call finish()
    allocate(rkv65r_class  :: s);  allocate(s2, source=s); call run('rkv65r',  [65, 71, 41]); call finish()
    allocate(rkv65_class   :: s);  allocate(s2, source=s); call run('rkv65',   [75, 81, 51],':'); call finish()
    allocate(rktf65_class  :: s);  allocate(s2, source=s); call run('rktf65',  [0,0,0],'--'); call finish()
    call done(5)

    call init()
    allocate(rktp75_class  :: s);  allocate(s2, source=s); call run('rktp75',  [0, 255, 38]);    call finish()
    allocate(rktmy7_class  :: s);  allocate(s2, source=s); call run('rktmy7',  [102, 247, 255]); call finish()
    allocate(rkv76e_class  :: s);  allocate(s2, source=s); call run('rkv76e',  [38, 189, 60]);   call finish()
    allocate(rkv76r_class  :: s);  allocate(s2, source=s); call run('rkv76r',  [149, 163, 93],':');   call finish()
    allocate(rkf78_class   :: s);  allocate(s2, source=s); call run('rkf78',   [66, 143, 77],'--');   call finish()
    allocate(rkv78_class   :: s);  allocate(s2, source=s); call run('rkv78',   [77, 105, 81]);   call finish()
    call done(7)

    call init()
    allocate(rktp86_class  :: s);  allocate(s2, source=s); call run('rktp86',  [0, 0, 0]);         call finish()
    allocate(rkdp87_class  :: s);  allocate(s2, source=s); call run('rkdp87',  [0, 0, 255]);        call finish()
    allocate(rkv87e_class  :: s);  allocate(s2, source=s); call run('rkv87e',  [90, 116, 230]);     call finish()
    allocate(rkv87r_class  :: s);  allocate(s2, source=s); call run('rkv87r',  [0,0,0],':');        call finish()
    allocate(rkk87_class  :: s);  allocate(s2, source=s); call run('rkk87',    [255,0,0]);          call finish()
    allocate(rkf89_class   :: s);  allocate(s2, source=s); call run('rkf89',   [116, 133, 207],'--'); call finish()
    allocate(rkv89_class   :: s);  allocate(s2, source=s); call run('rkv89',   [169, 176, 219]);    call finish()
    call done(8)

    call init()
    allocate(rkt98a_class  :: s);  allocate(s2, source=s); call run('rkt98a',  [195, 0, 255]);   call finish()
    allocate(rkv98e_class  :: s);  allocate(s2, source=s); call run('rkv98e',  [192, 52, 235]);  call finish()
    allocate(rkv98r_class  :: s);  allocate(s2, source=s); call run('rkv98r',  [79, 5, 153]);    call finish()
    call done(9)

    call init()
    allocate(rkf108_class  :: s);  allocate(s2, source=s); call run('rkf108',  [0, 0, 0]);  call finish()
    allocate(rkc108_class  :: s);  allocate(s2, source=s); call run('rkc108',  [0, 0, 255],'--'); call finish()
    allocate(rkb109_class  :: s);  allocate(s2, source=s); call run('rkb109',  [232, 0, 0],':'); call finish()
    call done(10)

    call init()
    allocate(rks1110a_class :: s); allocate(s2, source=s); call run('rks1110a',[145, 145, 145]); call finish()
    call done(11)

    call init()
    allocate(rkf1210_class :: s);  allocate(s2, source=s); call run('rkf1210', [94,94,94]);            call finish()
    allocate(rko129_class :: s);   allocate(s2, source=s); call run('rko129',  [145, 145, 145],'--');  call finish()
    call done(12)

    call init()
    allocate(rkf1412_class :: s);  allocate(s2, source=s); call run('rkf1412', [0, 0, 0]);  call finish()
    call done(14)

    contains
!*****************************************************************************************

!*****************************************************************************************
    subroutine init()
        !! initialize plot
        call plt%initialize(grid=.true.,xlabel='Relative Error',&
            ylabel='Number of Function Evaluations',&
            figsize=[30,15],font_size=font_size,axes_labelsize=font_size,&
            xtick_labelsize=font_size, ytick_labelsize=font_size,&
            legend_fontsize=legend_fontsize,&
            title='Problem 1: Variable-Step Runge Kutta Methods'//case,legend=.true.)
    end subroutine init
!*****************************************************************************************

!*****************************************************************************************
    subroutine done(iorder)
        !! save plot
        integer,intent(in) :: iorder
        character(len=3) :: rstr,istr
        write(rstr,'(I3)') wp
        write(istr,'(I3)') iorder
        call plt%savefig(figfile='rk_performance_test_R'//trim(adjustl(rstr))//&
                                 '_ORDER='//trim(adjustl(istr))//'.png',istat=istat)
    end subroutine done
!*****************************************************************************************

!*****************************************************************************************
    subroutine finish(); deallocate(s); deallocate(s2); end subroutine finish
!*****************************************************************************************

!*****************************************************************************************
    subroutine run(method,color,linestyle)
        !! generate a performance plot for all the methods
        character(len=*),intent(in) :: method !! name of the RK method to use
        integer,dimension(3),intent(in) :: color !! color for the plot
        character(len=*),intent(in),optional :: linestyle !! plot line style (e.g,. '.-')

        character(len=:),allocatable :: linestyle_
        integer :: i !! counter
        real(wp),dimension(:),allocatable :: r_error, v_error
        integer,dimension(:),allocatable :: feval
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

        allocate(r_error(factor*exp_start:factor*exp_stop))
        allocate(v_error(factor*exp_start:factor*exp_stop))
        allocate(feval(factor*exp_start:factor*exp_stop))

        write(*,*) trim(method)

        if (present(linestyle)) then
            linestyle_ = trim(linestyle)
        else
            linestyle_ = '-'
        end if

        !initial conditions:
        x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
              1.0_wp,2.0_wp,3.0_wp]
        t0 = 0.0_wp      ! initial time (sec)
        dt = 1.0_wp     ! initial time step (sec)
        tf = 10000.0_wp  ! final time (sec)

        do i = exp_start*factor, exp_stop*factor

            rtol = 10.0_wp ** (-i/real(factor,wp))
            atol = 10.0_wp ** (-i/real(factor,wp))

            !write(*,*) '   rtol', rtol

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
                fevals = 0
                call s%integrate(t0,x0,dt,tf,xf)     !forward
                feval(i) = fevals
                call s%integrate(tf,xf,dt,t0,x02)    !reverse

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
                            linestyle=linestyle_,color=real(color/255.0_wp,wp),&
                            markersize=5,linewidth=4,istat=istat,&
                            xscale='log',yscale='log')

    end subroutine run
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

        r = x(1:3)
        v = x(4:6)
        rmag = norm2(r)
        a_grav = -mu/rmag**3 * r !acceleration due to gravity

        xdot(1:3) = v
        xdot(4:6) = a_grav

        fevals = fevals + 1

        end subroutine twobody
    !*********************************************************

    end program rk_performance
!*****************************************************************************************
