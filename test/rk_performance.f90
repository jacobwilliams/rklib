
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
    allocate(rkbs32_class  :: s); call run([0,0,255])
    allocate(rkssp43_class :: s); call run([255,0,0])
    call done(3)

    call init()
    allocate(rkf45_class   :: s); call run([0,0,255])
    call done(4)

    call init()
    allocate(rkck54_class  :: s); call run([255, 102, 0])
    allocate(rkdp54_class  :: s); call run([189, 90, 25])
    allocate(rkt54_class   :: s); call run([143, 78, 36])
    allocate(rks54_class   :: s); call run([243, 78, 36],':')
    allocate(rkpp54_class  :: s); call run([243, 78, 255],'--')
    call done(5)

    call init()
    allocate(rkdp65_class  :: s); call run([251, 255, 0])
    allocate(rkc65_class   :: s); call run([207, 0, 0])
    allocate(rktp64_class  :: s); call run([187, 189, 49])
    allocate(rkv65e_class  :: s); call run([0,0,255])
    allocate(rkv65r_class  :: s); call run([65, 71, 41])
    allocate(rkv65_class   :: s); call run([75, 81, 51],':')
    allocate(rktf65_class  :: s); call run([0,0,0],'--')
    call done(6)

    call init()
    allocate(rktp75_class  :: s); call run([0, 255, 38])
    allocate(rktmy7_class  :: s); call run([102, 247, 255])
    allocate(rkv76e_class  :: s); call run([38, 189, 60])
    allocate(rkv76r_class  :: s); call run([149, 163, 93],':')
    allocate(rkf78_class   :: s); call run([66, 143, 77],'--')
    allocate(rkv78_class   :: s); call run([77, 105, 81])
    allocate(dverk78_class   :: s); call run([102, 200, 81])
    call done(7)

    call init()
    allocate(rktp86_class  :: s); call run([0, 0, 0])
    allocate(rkdp87_class  :: s); call run([0, 0, 255])
    allocate(rkv87e_class  :: s); call run([90, 116, 230])
    allocate(rkv87r_class  :: s); call run([0,0,0],':')
    allocate(rkk87_class   :: s); call run([255,0,0])
    allocate(rkf89_class   :: s); call run([116, 133, 207],'--')
    allocate(rkv89_class   :: s); call run([169, 176, 219])
    call done(8)

    call init()
    allocate(rkt98a_class  :: s); call run([195, 0, 255])
    allocate(rkv98e_class  :: s); call run([192, 52, 235])
    allocate(rkv98r_class  :: s); call run([79, 5, 153])
    call done(9)

    call init()
    allocate(rkf108_class  :: s); call run([0, 0, 0])
    allocate(rkc108_class  :: s); call run([0, 0, 255],'--')
    allocate(rkb109_class  :: s); call run([232, 0, 0],':')
    call done(10)

    call init()
    allocate(rks1110a_class :: s);  call run([145, 145, 145])
    call done(11)

    call init()
    allocate(rkf1210_class :: s); call run([94,94,94])
    allocate(rko129_class  :: s); call run([145, 145, 145],'--')
    call done(12)

    call init()
    allocate(rkf1412_class :: s); call run([0, 0, 0])
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
    subroutine run(color,linestyle)
        !! generate a performance plot for all the methods
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
        character(len=:),allocatable :: method !! name of the RK method to use
        type(rklib_properties) :: props

        allocate(s2, source=s)
        props = s%properties()
        method = props%short_name

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

        call finish()

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
