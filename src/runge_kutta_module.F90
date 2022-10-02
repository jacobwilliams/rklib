!*****************************************************************************************
!> author: Jacob Williams
!
!  Runge-Kutta integration.
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:
#ifdef REAL32
!      `real(kind=real32)` [4 bytes]
#elif REAL64
!      `real(kind=real64)` [8 bytes]
#elif REAL128
!      `real(kind=real128)` [16 bytes]
#else
!      `real(kind=real64)` [8 bytes]
#endif

    module runge_kutta_module

    use iso_fortran_env
    use root_module

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: rk_module_rk = real32   !! real kind used by this module [4 bytes]
#elif REAL64
    integer,parameter,public :: rk_module_rk = real64   !! real kind used by this module [8 bytes]
#elif REAL128
    integer,parameter,public :: rk_module_rk = real128  !! real kind used by this module [16 bytes]
#else
    integer,parameter,public :: rk_module_rk = real64   !! real kind used by this module [8 bytes]
#endif

    integer,parameter :: wp = rk_module_rk  !! local copy of `rk_module_rk` with a shorter name
    real(wp),parameter :: zero = 0.0_wp

    type,public :: stepsize_class

        !! Algorithms for adjusting the step size for variable-step
        !! Runge-Kutta integrators.

        private

        real(wp) :: hmax           = huge(1.0_wp)           !! maximum allowed step size
        real(wp) :: hmin           = 2.0_wp*epsilon(1.0_wp) !! minimum allowed step size
        real(wp) :: hfactor_reject = 1.0e-3_wp              !! minimum allowed factor for decreasing step size after rejected step
        real(wp) :: hfactor_accept = 100.0_wp               !! maximum allowed factor for increasing step size after accepted step
        integer  :: accept_mode    = 1                      !! method to determine if step is accepted [1,2]
        integer  :: max_attempts   = 100                    !! maximum number of attempts to decrease step size before giving up

        ! the `hfactor` equation is:
        !
        ! if (relative_err) then
        !     hfactor = safety_factor * abs(tol*h/err)**(one/real(p+p_exponent_offset,wp))
        ! else
        !     hfactor = safety_factor * abs(tol/err)**(one/real(p+p_exponent_offset,wp))
        ! end if

        logical  :: relative_err      = .false. !! to use `tol*h` in the `hfactor` equation
        real(wp) :: safety_factor     = 0.9_wp  !! for `hfactor` equation (>0)
        integer  :: p_exponent_offset = 0       !! p + this value in the exponent (0 or 1)

        procedure(norm_func),nopass,pointer :: norm => maxval_func
            !! routine for computing the norm of the state

        contains

        private

        procedure,public :: initialize => stepsize_class_constructor
        procedure,public :: compute_stepsize
        procedure,public :: destroy => destroy_stepsize_class

    end type stepsize_class

    type,abstract,public :: rk_class

        !! main integration class:

        integer :: n = 0  !! user specified number of variables
        procedure(deriv_func),pointer  :: f      => null()  !! user-specified derivative function
        procedure(report_func),pointer :: report => null()  !! user-specified report function
        procedure(event_func),pointer  :: g      => null()  !! event function (stop when this is zero)

        contains

        procedure,public :: destroy !! destructor

    end type rk_class

    type,extends(rk_class),abstract,public :: rk_fixed_step_class

        !! fixed step size class

        contains

        procedure(step_func_fixed),deferred :: step !! the step routine for the rk method
        procedure,public :: initialize => initialize_fixed_step       !! initialize the class (set n,f, and report)

        procedure,public :: integrate => integrate_fixed_step
        procedure,public :: integrate_to_event => integrate_to_event_fixed_step

    end type rk_fixed_step_class

    type,extends(rk_class),abstract,public :: rk_variable_step_class

        !! Main integration class for variable step size Runge-Kutta methods

        private

        class(stepsize_class),allocatable :: stepsize_method  !! the method for varying the step size

        real(wp),dimension(:),allocatable :: rtol  !! relative tolerance (`size(n)`)
        real(wp),dimension(:),allocatable :: atol  !! absolute tolerance (`size(n)`)

        integer :: p = 0 !! order of the method

        integer :: hinit_method = 1 !! if automatically computing the inital step size, which
                                    !! method to use. 1 = `hstart`, 2 = `hinit`.

        integer :: num_rejected_steps = 0 !! number of rejected steps

        contains

        private

        procedure(step_func_variable),deferred :: step !! the step routine for the rk method
        procedure,public :: initialize => initialize_variable_step  !! initialize the class (set n,f, and report)
        procedure,public :: integrate => integrate_variable_step
        procedure,public :: integrate_to_event => integrate_to_event_variable_step

        procedure(order_func),deferred   :: order              !! returns `p`, the order of the method
        procedure :: hstart  !! for automatically computing the initial step size [this is from DDEABM]
        procedure :: hinit   !! for automatically computing the initial step size [this is from DOP853]

    end type rk_variable_step_class

    ! Fixed step methods:
    type,extends(rk_fixed_step_class),public :: rk4_class
        !! 4th order Runge-Kutta method.
        contains
        procedure :: step => rk4
    end type rk4_class
    type,extends(rk_fixed_step_class),public :: rk7_class
        !! 7th order Runge-Kutta method.
        contains
        procedure :: step => rk7
    end type rk7_class
    type,extends(rk_fixed_step_class),public :: rk8_10_class
        !! 8th order Runge-Kutta method.
        contains
        procedure :: step => rk8_10
    end type rk8_10_class

    ! Variable step methods:
    type,extends(rk_variable_step_class),public :: rkf78_class
        !! Runga-Kutta Fehlberg 7(8) method.
        contains
        procedure :: step  => rkf78
        procedure :: order => rkf78_order
    end type rkf78_class
    type,extends(rk_variable_step_class),public :: rkf89_class
        !! Runga-Kutta Fehlberg 8(9) method.
        contains
        procedure :: step  => rkf89
        procedure :: order => rkf89_order
    end type rkf89_class
    type,extends(rk_variable_step_class),public :: rkv89_class
        !! Runga-Kutta Verner 8(9) method.
        contains
        procedure :: step  => rkv89
        procedure :: order => rkv89_order
    end type rkv89_class
    type,extends(rk_variable_step_class),public :: rkf108_class
        !! Runga-Kutta Feagin 8(10) method.
        contains
        procedure :: step  => rkf108
        procedure :: order => rkf108_order
    end type rkf108_class
    type,extends(rk_variable_step_class),public :: rkf1210_class
        !! Runga-Kutta Feagin 12(10) method.
        contains
        procedure :: step  => rkf1210
        procedure :: order => rkf1210_order
    end type rkf1210_class
    type,extends(rk_variable_step_class),public :: rkf1412_class
        !! Runga-Kutta Feagin 14(12) method.
        contains
        procedure :: step  => rkf1412
        procedure :: order => rkf1412_order
    end type rkf1412_class

    interface

        subroutine integrate_func(me,t0,x0,h,tf,xf,ierr)
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)     :: me
            real(wp),intent(in)               :: t0    !! initial time
            real(wp),dimension(:),intent(in)  :: x0    !! initial state
            real(wp),intent(in)               :: h     !! initial abs(time step)
            real(wp),intent(in)               :: tf    !! final time
            real(wp),dimension(:),intent(out) :: xf    !! final state
            integer,intent(out),optional      :: ierr
        end subroutine integrate_func

        subroutine integrate_to_event_func(me,t0,x0,h,tmax,tol,tf,xf,gf,ierr)
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)     :: me
            real(wp),intent(in)               :: t0      !! initial time
            real(wp),dimension(:),intent(in)  :: x0      !! initial state
            real(wp),intent(in)               :: h       !! abs(time step)
            real(wp),intent(in)               :: tmax    !! max final time if event not located
            real(wp),intent(in)               :: tol     !! function tolerance for root finding
            real(wp),intent(out)              :: tf      !! actual final time reached
            real(wp),dimension(:),intent(out) :: xf      !! final state (at tf)
            real(wp),intent(out)              :: gf      !! g value at tf
            integer,intent(out),optional      :: ierr
        end subroutine integrate_to_event_func

        subroutine deriv_func(me,t,x,xdot)  !! derivative function
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)     :: me
            real(wp),intent(in)               :: t    !! time
            real(wp),dimension(:),intent(in)  :: x    !! state vector
            real(wp),dimension(:),intent(out) :: xdot !! derivative of state vector
        end subroutine deriv_func

        subroutine event_func(me,t,x,g)  !! event function
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)    :: me
            real(wp),intent(in)              :: t !! time
            real(wp),dimension(:),intent(in) :: x !! state vector
            real(wp),intent(out)             :: g !! g(t,x). The goal is to stop the integration when g=0.
        end subroutine event_func

        subroutine report_func(me,t,x)  !! report function
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)    :: me
            real(wp),intent(in)              :: t !! time
            real(wp),dimension(:),intent(in) :: x !! state vector
        end subroutine report_func

        subroutine step_func_fixed(me,t,x,h,xf)   !! rk step function
        import :: rk_fixed_step_class,wp
        implicit none
            class(rk_fixed_step_class),intent(inout) :: me
            real(wp),intent(in)                      :: t  !! initial time
            real(wp),dimension(me%n),intent(in)      :: x  !! initial state vector
            real(wp),intent(in)                      :: h  !! time step \( |\Delta t| \)
            real(wp),dimension(me%n),intent(out)     :: xf !! final state vector
        end subroutine step_func_fixed

        subroutine step_func_variable(me,t,x,h,xf,terr)   !! rk step function
        import :: rk_variable_step_class,wp
        implicit none
            class(rk_variable_step_class),intent(inout) :: me
            real(wp),intent(in)                         :: t    !! initial time
            real(wp),dimension(me%n),intent(in)         :: x    !! initial state vector
            real(wp),intent(in)                         :: h    !! time step \( |\Delta t| \)
            real(wp),dimension(me%n),intent(out)        :: xf   !! final state vector
            real(wp),dimension(me%n),intent(out)        :: terr !! truncation error estimate
        end subroutine step_func_variable

        pure function norm_func(x) result(xmag)
        !! Vector norm function. Must return a value \( \ge 0 \).
        import :: wp
        implicit none
            real(wp),dimension(:),intent(in) :: x    !! a vector
            real(wp)                         :: xmag !! the magnitude of the vector
        end function norm_func

        pure function order_func(me) result(p)
        import :: rk_variable_step_class
        implicit none
            class(rk_variable_step_class),intent(in) :: me
            integer :: p !! order of the method
        end function order_func

    end interface

    ! submodule procedures:
    interface
        module subroutine rk4(me,t,x,h,xf)
            implicit none
            class(rk4_class),intent(inout)       :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rk4
        module subroutine rk7(me,t,x,h,xf)
            implicit none
            class(rk7_class),intent(inout)       :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rk7
        module subroutine rk8_10(me,t,x,h,xf)
            implicit none
            class(rk8_10_class),intent(inout)    :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rk8_10
        module subroutine rkf78(me,t,x,h,xf,terr)
            implicit none
            class(rkf78_class),intent(inout)     :: me
            real(wp),intent(in)                  :: t
            real(wp),dimension(me%n),intent(in)  :: x
            real(wp),intent(in)                  :: h
            real(wp),dimension(me%n),intent(out) :: xf
            real(wp),dimension(me%n),intent(out) :: terr
        end subroutine rkf78
        module subroutine rkf89(me,t,x,h,xf,terr)
            implicit none
            class(rkf89_class),intent(inout)     :: me
            real(wp),intent(in)                  :: t
            real(wp),dimension(me%n),intent(in)  :: x
            real(wp),intent(in)                  :: h
            real(wp),dimension(me%n),intent(out) :: xf
            real(wp),dimension(me%n),intent(out) :: terr
        end subroutine rkf89
        module subroutine rkv89(me,t,x,h,xf,terr)
            implicit none
            class(rkv89_class),intent(inout)     :: me
            real(wp),intent(in)                  :: t
            real(wp),dimension(me%n),intent(in)  :: x
            real(wp),intent(in)                  :: h
            real(wp),dimension(me%n),intent(out) :: xf
            real(wp),dimension(me%n),intent(out) :: terr
        end subroutine rkv89
        module subroutine rkf108(me,t,x,h,xf,terr)
            implicit none
            class(rkf108_class),intent(inout)     :: me
            real(wp),intent(in)                  :: t
            real(wp),dimension(me%n),intent(in)  :: x
            real(wp),intent(in)                  :: h
            real(wp),dimension(me%n),intent(out) :: xf
            real(wp),dimension(me%n),intent(out) :: terr
        end subroutine rkf108
        module subroutine rkf1210(me,t,x,h,xf,terr)
            implicit none
            class(rkf1210_class),intent(inout)     :: me
            real(wp),intent(in)                  :: t
            real(wp),dimension(me%n),intent(in)  :: x
            real(wp),intent(in)                  :: h
            real(wp),dimension(me%n),intent(out) :: xf
            real(wp),dimension(me%n),intent(out) :: terr
        end subroutine rkf1210
        module subroutine rkf1412(me,t,x,h,xf,terr)
            implicit none
            class(rkf1412_class),intent(inout)     :: me
            real(wp),intent(in)                  :: t
            real(wp),dimension(me%n),intent(in)  :: x
            real(wp),intent(in)                  :: h
            real(wp),dimension(me%n),intent(out) :: xf
            real(wp),dimension(me%n),intent(out) :: terr
        end subroutine rkf1412
        pure module function rkf78_order(me) result(p)
            implicit none
            class(rkf78_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkf78_order
        pure module function rkf89_order(me) result(p)
            implicit none
            class(rkf89_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkf89_order
        pure module function rkv89_order(me) result(p)
            implicit none
            class(rkv89_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkv89_order
        pure module function rkf108_order(me) result(p)
            implicit none
            class(rkf108_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkf108_order
        pure module function rkf1210_order(me) result(p)
            implicit none
            class(rkf1210_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkf1210_order
        pure module function rkf1412_order(me) result(p)
            implicit none
            class(rkf1412_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkf1412_order
    end interface

    ! public routines:
    public :: norm2_func,maxval_func

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[rk_class]].

    subroutine destroy(me)

    implicit none

    class(rk_class),intent(out)   :: me

    end subroutine destroy
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the [[rk_class]].

    subroutine initialize_fixed_step(me,n,f,report,g)

    implicit none

    class(rk_fixed_step_class),intent(inout)   :: me
    integer,intent(in)              :: n       !! number of variables
    procedure(deriv_func)           :: f       !! derivative function
    procedure(report_func),optional :: report  !! for reporting the steps
    procedure(event_func),optional  :: g       !! for stopping at an event

    call me%destroy()

    me%n = n
    me%f => f
    if (present(report)) me%report => report
    if (present(g))      me%g      => g

    end subroutine initialize_fixed_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Main integration routine for the [[rk_class]].

    subroutine integrate_fixed_step(me,t0,x0,h,tf,xf,ierr)

    implicit none

    class(rk_fixed_step_class),intent(inout)     :: me
    real(wp),intent(in)               :: t0    !! initial time
    real(wp),dimension(:),intent(in)  :: x0    !! initial state
    real(wp),intent(in)               :: h     !! abs(time step)
    real(wp),intent(in)               :: tf    !! final time
    real(wp),dimension(:),intent(out) :: xf    !! final state
    integer,intent(out),optional      :: ierr  !! 0 = no errors,
                                               !! <0 = error.
                                               !! if not present, an error will stop program.

    real(wp) :: t,dt,t2
    real(wp),dimension(me%n) :: x
    logical :: last,export

    if (.not. associated(me%f)) then
        if (present(ierr)) then
            ierr = -1
            return
        else
            error stop 'Error in integrate: f is not associated.'
        end if
    end if

    export = associated(me%report)

    if (export) call me%report(t0,x0)  !first point

    if (h==zero) then
        xf = x0
    else

        t = t0
        x = x0
        dt = sign(h,tf-t0)  !time step (correct sign)
        do
            t2 = t + dt
            last = ((dt>=zero .and. t2>=tf) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tf))         !
            if (last) dt = tf-t                     !
            call me%step(t,x,dt,xf)
            if (last) exit
            if (export) call me%report(t2,xf)   !intermediate point
            x = xf
            t = t2
        end do

    end if

    if (export) call me%report(tf,xf)   !last point

    end subroutine integrate_fixed_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Event-finding integration routine for the [[rk_class]].
!  Integrates until g(t,x)=0, or until t=tf (whichever happens first).
!
!@note There are some efficiency improvements that could be made here.
!      This is a work in progress.

    subroutine integrate_to_event_fixed_step(me,t0,x0,h,tmax,tol,tf,xf,gf)

    implicit none

    class(rk_fixed_step_class),intent(inout)     :: me
    real(wp),intent(in)               :: t0      !! initial time
    real(wp),dimension(:),intent(in)  :: x0      !! initial state
    real(wp),intent(in)               :: h       !! abs(time step)
    real(wp),intent(in)               :: tmax    !! max final time if event not located
    real(wp),intent(in)               :: tol     !! function tolerance for root finding
    real(wp),intent(out)              :: tf      !! actual final time reached
    real(wp),dimension(:),intent(out) :: xf      !! final state (at tf)
    real(wp),intent(out)              :: gf      !! g value at tf

    !local variables:
    real(wp) :: t,dt,t2,ga,gb,dt_root,dum
    real(wp),dimension(me%n) :: x,g_xf
    logical :: first,last,export
    procedure(report_func),pointer :: report
    type(brent_solver) :: solver
    integer :: iflag

    if (.not. associated(me%f)) error stop 'Error in integrate_to_event: f is not associated.'
    if (.not. associated(me%g)) error stop 'Error in integrate_to_event: g is not associated.'
    if (h==zero) error stop 'Error in integrate_to_event: h must not be zero.'

    !If the points are being exported:
    export = associated(me%report)

    !first point:
    if (export) call me%report(t0,x0)

    if (t0==tmax) then
        xf = x0
        tf = t0
        call me%g(t0,x0,gf)
    else

        first = .true.
        t = t0
        x = x0
        call me%g(t0,x0,ga)     !evaluate event function
        dt = sign(h,tmax-t0)    !time step (correct sign)

        do

            t2 = t + dt
            last = ((dt>=zero .and. t2>=tmax) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tmax))         !
            if (last) then
                dt = tmax-t
                t2 = tmax
            end if
            call me%step(t,x,dt,xf)
            call me%g(t2,xf,gb)     !evaluate event function

            if (first .and. abs(ga)<=tol) then

                !we ignore a root at t0 after the first step
                if (abs(gb)<=tol) then !check this one since it could have landed on a root
                    gf = gb
                    tf = t2
                    exit
                else
                    if (last) then  !exiting without having found a root
                        tf = t2
                        gf = gb
                        exit
                    end if
                    if (export) call me%report(t2,xf)   !intermediate point
                    x = xf
                    t = t2
                    ga = gb
                end if

            elseif (ga*gb<=zero) then !there is a root somewhere on [t,t+dt]

                !find the root:
                call solver%initialize(solver_func)
                call solver%solve(zero,dt,dt_root,dum,iflag,fax=ga,fbx=gb)
                t2 = t + dt_root
                gf = solver_func(solver,dt_root)
                tf = t2
                xf = g_xf !computed in the solver function
                exit

            else  !no root yet, continue

                if (last) then  !exiting without having found a root
                    tf = t2
                    gf = gb
                    exit
                end if
                if (export) call me%report(t2,xf)   !intermediate point
                x = xf
                t = t2
                ga = gb

            end if

            if (first) first = .false.

        end do

    end if

    if (export) call me%report(t2,xf)   !last point

    contains

        function solver_func(this,delt) result(g)

        !! root solver function. The input is the dt offset from time t.

        implicit none

        class(root_solver),intent(inout) :: this
        real(wp),intent(in) :: delt  !! from [0 to dt]
        real(wp) :: g

        !take a step from t to t+delt and evaluate g function:
        call me%step(t,x,delt,g_xf)
        call me%g(t+delt,g_xf,g)

        end function solver_func

    end subroutine integrate_to_event_fixed_step
!*****************************************************************************************



!*****************************************************************************************
!>
!  Use intrinsic `norm2(x)` for computing the vector norm.

    pure function norm2_func(x) result(xmag)

    implicit none

    real(wp),dimension(:),intent(in) :: x
    real(wp) :: xmag

    xmag = norm2(x)

    end function norm2_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  Use `maxval(abs(x))` for computing the vector norm.

    pure function maxval_func(x) result(xmag)

    implicit none

    real(wp),dimension(:),intent(in) :: x
    real(wp) :: xmag

    xmag = maxval(abs(x))

    end function maxval_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  Constructor for a [[stepsize_class]].
!
!@warning The `norm` and `compute_h_factor` options aren't public in the module.
!         Need to fix this.

    pure subroutine stepsize_class_constructor(me,hmin,hmax,hfactor_reject,&
                        hfactor_accept,norm,accept_mode,relative_err,&
                        safety_factor,p_exponent_offset,max_attempts)


    implicit none

    class(stepsize_class),intent(inout)       :: me
    real(wp),intent(in),optional              :: hmin             !! minimum allowed step size (>0)
    real(wp),intent(in),optional              :: hmax             !! maximum allowed step size (>0)
    real(wp),intent(in),optional              :: hfactor_reject   !! minimum allowed factor for
                                                                  !! decreasing step size after
                                                                  !! rejected step (>0)
    real(wp),intent(in),optional              :: hfactor_accept   !! maximum allowed factor for
                                                                  !! decreasing step size after
                                                                  !! accepted step (>0)
    procedure(norm_func),optional             :: norm             !! the user-specified \( ||x|| \)
                                                                  !! function
    integer,intent(in),optional               :: accept_mode      !! method to determine if step
                                                                  !! is accepted [1,2]
    integer,intent(in),optional               :: max_attempts     !! max step size change attempts
                                                                  !! after rejected step
    logical,intent(in),optional   :: relative_err       !! to use `tol*h` in the `hfactor` equation
    real(wp),intent(in),optional  :: safety_factor      !! for `hfactor` equation (>0)
    integer,intent(in),optional   :: p_exponent_offset  !! p + this value in the exponent (0 or 1)

    if (present(hmin))             me%hmin             = abs(hmin)
    if (present(hmax))             me%hmax             = abs(hmax)
    if (present(hfactor_reject))   me%hfactor_reject   = abs(hfactor_reject)
    if (present(hfactor_accept))   me%hfactor_accept   = abs(hfactor_accept)
    if (present(norm))             me%norm             => norm
    if (present(accept_mode))      me%accept_mode      = accept_mode
    if (present(max_attempts))     me%max_attempts     = max_attempts

    !if (present(compute_h_factor)) me%compute_h_factor => compute_h_factor
    if (present(relative_err     )) me%relative_err      = relative_err
    if (present(safety_factor    )) me%safety_factor     = abs(safety_factor    )
    if (present(p_exponent_offset)) me%p_exponent_offset = abs(p_exponent_offset)

    end subroutine stepsize_class_constructor
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[stepsize_class]].

    subroutine destroy_stepsize_class(me)

    implicit none

    class(stepsize_class),intent(out) :: me

    end subroutine destroy_stepsize_class
!*****************************************************************************************

!*****************************************************************************************
    pure subroutine compute_stepsize(me,h,tol,err,p,hnew,accept)

    !! Compute the new step size using the specific method.

    implicit none

    class(stepsize_class),intent(in) :: me
    real(wp),intent(in)              :: h      !! current step size (<>0)
    real(wp),intent(in)              :: tol    !! abs error tolerance (>0)
    real(wp),intent(in)              :: err    !! truncation error estimate (>0)
    integer,intent(in)               :: p      !! order of the method
    real(wp),intent(out)             :: hnew   !! new step size (<>0)
    logical,intent(out)              :: accept !! if the step is accepted

    real(wp) :: hfactor  !! step size factor (>0)
    real(wp),parameter :: small = 10.0_wp * epsilon(1.0_wp) !! small error value

    if (err<=small) then ! the error is extremely small

        hfactor = me%hfactor_accept
        accept = .true.

    else

        ! compute base factor based on the selected formula:
        !hfactor = abs(me%compute_h_factor(h,tol,err,p))
        if (me%relative_err) then
            hfactor = abs( me%safety_factor*abs(tol*h/err)**(1.0_wp/real(p+me%p_exponent_offset,wp)) )
        else
            hfactor = abs( me%safety_factor*abs(tol/err)**(1.0_wp/real(p+me%p_exponent_offset,wp)) )
        end if

        ! if the step is to be accepted:
        select case (me%accept_mode)
        case(1) !algorithm 17.12
            accept = (hfactor>=1.0_wp)
        case(2) !algorithm 17.13
            accept = (err<=tol)
        end select

        !...notes:
        ! see: L. Shampine "Some Practical Runge-Kutta Formulas",
        !      Mathematics of Computation, 46(173), Jan 1986.
        ! different conditions for satisfying error conditions:
        !  ||err|| <= tol   -- Error per step (EPS)
        !  ||err|| <= h*tol -- Error per unit step (EPUS)

        !compute the actual hfactor based on the limits:
        if (accept) then
            hfactor = min(me%hfactor_accept, hfactor)
        else
            hfactor = max(me%hfactor_reject, hfactor)
        end if

    end if

    ! compute the new step size (enforce min/max bounds & add sign):
    hnew = sign(max(me%hmin,min(me%hmax,abs(h)*hfactor)),h)

    end subroutine compute_stepsize
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the [[rk_variable_step_class]].

    subroutine initialize_variable_step(me,n,f,rtol,atol,stepsize_method,hinit_method,report,g)

    implicit none

    class(rk_variable_step_class),intent(inout) :: me
    integer,intent(in)                          :: n               !! number of equations
    procedure(deriv_func)                       :: f               !! derivative function
    real(wp),dimension(:),intent(in)            :: rtol            !! relative tolerance (if size=1,
                                                                   !! then same tol used for all
                                                                   !! equations)
    real(wp),dimension(:),intent(in)            :: atol            !! absolute tolerance (if size=1,
                                                                   !! then same tol used for all
                                                                   !! equations)
    class(stepsize_class),intent(in)            :: stepsize_method !! method for varying the step size
    integer,intent(in),optional                 :: hinit_method    !! which method to use for
                                                                   !! automatic initial step size
                                                                   !! computation.
                                                                   !! 1 = use `hstart`, 2 = use `hinit`.
    procedure(report_func),optional             :: report          !! for reporting the steps
    procedure(event_func),optional              :: g               !! for stopping at an event

    call me%destroy()

    me%n = n
    me%f => f

    allocate(me%rtol(n))
    allocate(me%atol(n))
    if (size(rtol)==1) then
        me%rtol = rtol(1) !use this for all equations
    else if (size(rtol)==n) then
        me%rtol = rtol
    else
        error stop 'invalid size for rtol array.'
    end if
    if (size(atol)==1) then
        me%atol = atol(1) !use this for all equations
    else if (size(atol)==n) then
        me%atol = atol
    else
        error stop 'invalid size for atol array.'
    end if

    if (present(hinit_method)) me%hinit_method = hinit_method

    if (present(report)) me%report => report
    if (present(g))      me%g      => g

    allocate(me%stepsize_method, source=stepsize_method)

    me%num_rejected_steps = 0

    end subroutine initialize_variable_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Main integration routine for the [[rk_variable_step_class]].

    subroutine integrate_variable_step(me,t0,x0,h,tf,xf,ierr)

    implicit none

    class(rk_variable_step_class),intent(inout)     :: me
    real(wp),intent(in)               :: t0    !! initial time
    real(wp),dimension(:),intent(in)  :: x0    !! initial state
    real(wp),intent(in)               :: h     !! initial abs(time step)
    real(wp),intent(in)               :: tf    !! final time
    real(wp),dimension(:),intent(out) :: xf    !! final state
    integer,intent(out),optional      :: ierr  !! 0 = no errors,
                                               !! <0 = error.
                                               !! if not present, an error will stop program.

    real(wp) :: t,dt,t2,err,tol,dt_new
    real(wp),dimension(me%n) :: x,terr,etol,xp0
    logical :: last,export,accept
    integer :: i,p

    if (present(ierr)) ierr = 0

    if (.not. associated(me%f)) then
        if (present(ierr)) then
            ierr = -1
            return
        else
            error stop 'Error in integrate: f is not associated.'
        end if
    end if

    me%num_rejected_steps = 0
    export = associated(me%report)

    if (export) call me%report(t0,x0)  !first point

    if (t0==tf) then
        xf = x0
    else

        t = t0
        x = x0

        if (h==zero) then
            ! compute an appropriate initial step size:
            ! WARNING: this may not be working in all cases .....
            etol = me%rtol * me%stepsize_method%norm(x0) + me%atol
            call me%f(t0,x0,xp0)  ! get initial dx/dt
            select case (me%hinit_method)
            case(1)
                call me%hstart(t0,tf,x0,xp0,etol,dt)
            case(2)
                dt = me%hinit(t0,x0,sign(1.0_wp,tf-t0),xp0,me%stepsize_method%hmax,me%atol,me%rtol)
            case default
                if (present(ierr)) then
                    ierr = -2
                    return
                else
                    error stop 'invalid hinit_method selection'
                end if
            end select
            !write(*,*) 'inital step size: ',dt
        else
            ! user-specified initial step size:
            dt = sign(h,tf-t0)  ! (correct sign)
        end if

        p = me%order()     !order of the method
        do
            t2 = t + dt
            last = ((dt>=zero .and. t2>=tf) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tf))         !
            if (last) dt = tf-t                     !

            do i=0,me%stepsize_method%max_attempts

                ! take a step:
                call me%step(t,x,dt,xf,terr)

                ! evaluate error and compute new step size:
                err = me%stepsize_method%norm(terr)
                tol = me%stepsize_method%norm( me%rtol * xf + me%atol )
                call me%stepsize_method%compute_stepsize(dt,tol,err,p,dt_new,accept)
                dt = dt_new

                if (accept) then
                    !accept this step
                    exit
                else
                    !step is rejected, repeat step with new dt
                    me%num_rejected_steps = me%num_rejected_steps + 1

                    !note: if we have reached the min step size, and the error
                    !is still too large, we can't proceed.
                    if (i>=me%stepsize_method%max_attempts) then
                        if (present(ierr)) then
                            ierr = -3
                            return
                        else
                            error stop 'error: too many attempts to reduce step size.'
                        end if
                    end if
                    if (abs(dt) <= abs(me%stepsize_method%hmin)) then
                        if (present(ierr)) then
                            ierr = -4
                            return
                        else
                            error stop 'warning: min step size.'
                        end if
                    end if

                    !......
                    !... if we have two rejected steps and the step size hasn't changed..
                    !    then we need to abort, since no progress is being made...
                    !......

                    last = ((dt>=zero .and. t2>=tf) .or. &  !adjust last time step
                            (dt<zero .and. t2<=tf))         !
                    if (last) dt = tf-t                     !
                    t2 = t + dt

                end if

            end do

            if (last) exit
            if (export) call me%report(t2,xf)   !intermediate point
            x = xf
            t = t2
        end do

    end if

    if (export) call me%report(tf,xf)   !last point

    end subroutine integrate_variable_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Event-finding integration routine for the [[rk_variable_step_class]].
!  Integrates until g(t,x)=0, or until t=tf (whichever happens first).
!
!@note There are some efficiency improvements that could be made here.
!      This is a work in progress.

    subroutine integrate_to_event_variable_step(me,t0,x0,h,tmax,tol,tf,xf,gf,ierr)

    implicit none

    class(rk_variable_step_class),intent(inout) :: me
    real(wp),intent(in)                  :: t0      !! initial time
    real(wp),dimension(me%n),intent(in)  :: x0      !! initial state
    real(wp),intent(in)                  :: h       !! abs(time step)
    real(wp),intent(in)                  :: tmax    !! max final time if event not located
    real(wp),intent(in)                  :: tol     !! function tolerance for root finding
    real(wp),intent(out)                 :: tf      !! actual final time reached
    real(wp),dimension(me%n),intent(out) :: xf      !! final state (at tf)
    real(wp),intent(out)                 :: gf      !! g value at tf
    integer,intent(out),optional      :: ierr  !! 0 = no errors,
                                               !! <0 = error.
                                               !! if not present, an error will stop program.

    real(wp),dimension(me%n) :: etol,xp0
    real(wp),dimension(me%n) :: x,g_xf
    real(wp),dimension(me%n) :: terr !! truncation error estimate
    integer :: i,p,iflag
    real(wp) :: t,dt,t2,ga,gb,dt_root,dum,err,dt_new,stol
    logical :: first,last,export,accept
    procedure(report_func),pointer :: report
    type(brent_solver) :: solver

    if (present(ierr)) ierr = 0

    if (.not. associated(me%f)) then
        if (present(ierr)) then
            ierr = -1
            return
        else
            error stop 'Error in integrate_to_event: f is not associated.'
        end if
    end if
    if (.not. associated(me%g)) then
        if (present(ierr)) then
            ierr = -2
            return
        else
            error stop 'Error in integrate_to_event: g is not associated.'
        end if
    end if

    me%num_rejected_steps = 0
    export = associated(me%report)

    if (export) call me%report(t0,x0)  !first point

    if (t0==tmax) then
        xf = x0
        tf = t0
        call me%g(t0,x0,gf)
    else

        first = .true.
        t = t0
        x = x0
        call me%g(t,x,ga)     !evaluate event function

        if (h==zero) then
            ! compute an appropriate initial step size:
            ! WARNING: this may not be working in all cases .....
            etol = me%rtol * me%stepsize_method%norm(x0) + me%atol
            call me%f(t0,x0,xp0)  ! get initial dx/dt
            select case (me%hinit_method)
            case(1)
                call me%hstart(t0,tmax,x0,xp0,etol,dt)
            case(2)
                dt = me%hinit(t0,x0,sign(1.0_wp,tmax-t0),xp0,me%stepsize_method%hmax,me%atol,me%rtol)
            case default
                if (present(ierr)) then
                    ierr = -3
                    return
                else
                    error stop 'invalid hinit_method selection'
                end if
            end select
        else
            ! user-specified initial step size:
            dt = sign(h,tmax-t0)  ! (correct sign)
        end if

        p = me%order()     !order of the method
        do

            t2 = t + dt
            last = ((dt>=zero .and. t2>=tmax) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tmax))         !
            if (last) then
                dt = tmax-t
                t2 = tmax
            end if

            do i=0,me%stepsize_method%max_attempts

                ! take a step:
                call me%step(t,x,dt,xf,terr)

                ! evaluate error and compute new step size:
                err = me%stepsize_method%norm(terr)
                stol = me%stepsize_method%norm( me%rtol * xf + me%atol )
                call me%stepsize_method%compute_stepsize(dt,stol,err,p,dt_new,accept)
                dt = dt_new

                if (accept) then
                    !accept this step
                    exit
                else
                    !step is rejected, repeat step with new dt
                    me%num_rejected_steps = me%num_rejected_steps + 1

                    !note: if we have reached the min step size, and the error
                    !is still too large, we can't proceed.
                    if (i>=me%stepsize_method%max_attempts) then
                        if (present(ierr)) then
                            ierr = -4
                            return
                        else
                            error stop 'error: too many attempts to reduce step size.'
                        end if
                    end if
                    if (abs(dt) <= abs(me%stepsize_method%hmin)) then
                        if (present(ierr)) then
                            ierr = -5
                            return
                        else
                            error stop 'warning: min step size.'
                        end if
                    end if

                    !......
                    !... if we have two rejected steps and the step size hasn't changed..
                    !    then we need to abort, since no progress is being made...
                    !......

                    last = ((dt>=zero .and. t2>=tmax) .or. &  !adjust last time step
                            (dt<zero .and. t2<=tmax))         !
                    if (last) then
                        dt = tmax-t
                        t2 = tmax
                    else
                        t2 = t + dt
                    end if

                end if

            end do

            call me%g(t2,xf,gb)     !evaluate event function

            if (first .and. abs(ga)<=tol) then

                !we ignore a root at t0 after the first step
                if (abs(gb)<=tol) then !check this one since it could have landed on a root
                    gf = gb
                    tf = t2
                    exit
                else
                    if (last) then  !exiting without having found a root
                        tf = t2
                        gf = gb
                        exit
                    end if
                    if (export) call me%report(t2,xf)   !intermediate point
                    x = xf
                    t = t2
                    ga = gb
                end if

            elseif (ga*gb<=zero) then !there is a root somewhere on [t,t+dt]

                !find the root:
                call solver%initialize(solver_func, rtol=tol, atol=tol)
                call solver%solve(zero,dt,dt_root,dum,iflag,fax=ga,fbx=gb)
                t2 = t + dt_root
                gf = solver_func(solver,dt_root)
                tf = t2
                xf = g_xf !computed in the solver function
                exit

            else  !no root yet, continue

                if (last) then  !exiting without having found a root
                    tf = t2
                    gf = gb
                    exit
                end if
                if (export) call me%report(t2,xf)   !intermediate point
                x = xf
                t = t2
                ga = gb

            end if

            if (first) first = .false.
            if (last) exit
            x = xf
            t = t2
        end do

    end if

    if (export) call me%report(tf,xf)   !last point

    contains

        function solver_func(this,delt) result(g)

        !! root solver function. The input is the dt offset from time t.

        implicit none

        class(root_solver),intent(inout) :: this
        real(wp),intent(in) :: delt  !! from [0 to dt]
        real(wp) :: g

        real(wp),dimension(me%n) :: terr !! truncation error estimate

        !take a step from t to t+delt and evaluate g function:
        ! [we don't check the error because we are within a
        !  step that was already accepted, so it should be ok]
        call me%step(t,x,delt,g_xf,terr)
        call me%g(t+delt,g_xf,g)

        end function solver_func

    end subroutine integrate_to_event_variable_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes a starting step size to be used in solving initial
!  value problems in ordinary differential equations.
!
!  It is based on an estimate of the local lipschitz constant for the
!  differential equation (lower bound on a norm of the jacobian) ,
!  a bound on the differential equation (first derivative), and
!  a bound on the partial derivative of the equation with respect to
!  the independent variable. (all approximated near the initial point a)
!
!@note Subroutine hstart also uses the `me%stepsize_method%norm`
!      function for computing vector norms
!
!@note This routine is from [DDEABM](https://github.com/jacobwilliams/ddeabm).
!
!# History
!   * 820301  date written -- watts, h. a., (snla)
!   * 890531  changed all specific intrinsics to generic.  (wrb)
!   * 890831  modified array declarations.  (wrb)
!   * 890911  removed unnecessary intrinsics.  (wrb)
!   * 891024  changed references from dvnorm to dhvnrm.  (wrb)
!   * 891214  prologue converted to version 4.0 format.  (bab)
!   * 900328  added type section.  (wrb)
!   * 910722  updated author section.  (als)
!   * December, 2015 : Refactored this routine (jw)
!   * April 2016 : Some modifications for the variable-step RK module (jw)

    subroutine hstart(me,a,b,y,yprime,etol,h)

    implicit none

    class(rk_variable_step_class),intent(inout) :: me
    real(wp),intent(in)                :: a         !! the initial point of integration.
    real(wp),intent(in)                :: b         !! a value of the independent variable used to define
                                                    !! the direction of integration. a reasonable choice is to
                                                    !! set `b` to the first point at which a solution is desired.
                                                    !! you can also use `b`, if necessary, to restrict the length
                                                    !! of the first integration step because the algorithm will
                                                    !! not compute a starting step length which is bigger than
                                                    !! `abs(b-a)`, unless `b` has been chosen too close to `a`.
                                                    !! (it is presumed that hstart has been called with `b`
                                                    !! different from `a` on the machine being used. also see the
                                                    !! discussion about the parameter `small`.)
    real(wp),dimension(me%n),intent(in) :: y        !! the vector of initial values of the `neq` solution
                                                    !! components at the initial point `a`.
    real(wp),dimension(me%n),intent(in) :: yprime   !! the vector of derivatives of the `neq`
                                                    !! solution components at the initial point `a`.
                                                    !! (defined by the differential equations in subroutine `me%f`)
    real(wp),dimension(me%n),intent(in) :: etol     !! the vector of error tolerances corresponding to
                                                    !! the `neq` solution components. it is assumed that all
                                                    !! elements are positive. following the first integration
                                                    !! step, the tolerances are expected to be used by the
                                                    !! integrator in an error test which roughly requires that
                                                    !! `abs(local error) <= etol` for each vector component.
    real(wp),intent(out) :: h       !! appropriate starting step size to be attempted by the
                                    !! differential equation method.

    real(wp),dimension(me%n) :: spy  !! work array which provide the routine with needed storage space.
    real(wp),dimension(me%n) :: pv   !! work array which provide the routine with needed storage space.
    real(wp),dimension(me%n) :: yp   !! work array which provide the routine with needed storage space.
    real(wp),dimension(me%n) :: sf   !! work array which provide the routine with needed storage space.

    real(wp),parameter :: small  = epsilon(1.0_wp)
    real(wp),parameter :: big    = huge(1.0_wp)
    real(wp),parameter :: relper = small**0.375_wp

    integer :: j, k, lk
    real(wp) :: absdx, da, delf, dely,&
                dfdub, dfdxb,&
                dx, dy, fbnd,&
                srydpb, tolexp, tolmin, tolp, tolsum, ydpb
    integer :: morder    !! the order of the formula which will be used by
                         !! the initial value method for taking the first integration
                         !! step.

    morder = me%p
    dx = b - a
    absdx = abs(dx)

    ! compute an approximate bound (dfdxb) on the partial
    ! derivative of the equation with respect to the
    ! independent variable. protect against an overflow.
    ! also compute a bound (fbnd) on the first derivative
    ! locally.

    da = sign(max(min(relper*abs(a),absdx),100.0_wp*small*abs(a)),dx)
    if (da == zero) da = relper*dx
    call me%f(a+da,y,sf)
    yp = sf - yprime
    delf = me%stepsize_method%norm(yp)
    dfdxb = big
    if (delf < big*abs(da)) dfdxb = delf/abs(da)
    fbnd = me%stepsize_method%norm(sf)

    ! compute an estimate (dfdub) of the local lipschitz
    ! constant for the system of differential equations. this
    ! also represents an estimate of the norm of the jacobian
    ! locally.  three iterations (two when neq=1) are used to
    ! estimate the lipschitz constant by numerical differences.
    ! the first perturbation vector is based on the initial
    ! derivatives and direction of integration. the second
    ! perturbation vector is formed using another evaluation of
    ! the differential equation.  the third perturbation vector
    ! is formed using perturbations based only on the initial
    ! values. components that are zero are always changed to
    ! non-zero values (except on the first iteration). when
    ! information is available, care is taken to ensure that
    ! components of the perturbation vector have signs which are
    ! consistent with the slopes of local solution curves.
    ! also choose the largest bound (fbnd) for the first
    ! derivative.
    !
    ! perturbation vector size is held
    ! constant for all iterations. compute
    ! this change from the
    ! size of the vector of initial
    ! values.

    dely = relper*me%stepsize_method%norm(y)
    if (dely == zero) dely = relper
    dely = sign(dely,dx)
    delf = me%stepsize_method%norm(yprime)
    fbnd = max(fbnd,delf)
    if (delf == zero) then
        ! cannot have a null perturbation vector
        spy  = zero
        yp   = 1.0_wp
        delf = me%stepsize_method%norm(yp)
    else
        ! use initial derivatives for first perturbation
        spy = yprime
        yp  = yprime
    end if

    dfdub = zero
    lk = min(me%n+1,3)
    do k = 1, lk
        ! define perturbed vector of initial values
        pv = y + yp * (dely/delf)
        if (k == 2) then
            ! use a shifted value of the independent variable
            ! in computing one estimate
            call me%f(a+da,pv,yp)
            pv = yp - sf
        else
            ! evaluate derivatives associated with perturbed
            ! vector and compute corresponding differences
            call me%f(a,pv,yp)
            pv = yp - yprime
        end if
        ! choose largest bounds on the first derivative
        ! and a local lipschitz constant
        fbnd = max(fbnd,me%stepsize_method%norm(yp))
        delf = me%stepsize_method%norm(pv)
        if (delf >= big*abs(dely)) then
            ! protect against an overflow
            dfdub = big
            exit
        end if
        dfdub = max(dfdub,delf/abs(dely))
        if (k == lk) exit

        ! choose next perturbation vector
        if (delf == zero) delf = 1.0_wp
        do j = 1, me%n
            if (k == 2) then
                dy = y(j)
                if (dy == zero) dy = dely/relper
            else
                dy = abs(pv(j))
                if (dy == zero) dy = delf
            end if
            if (spy(j) == zero) spy(j) = yp(j)
            if (spy(j) /= zero) dy = sign(dy,spy(j))
            yp(j) = dy
        end do
        delf = me%stepsize_method%norm(yp)
    end do

    ! compute a bound (ydpb) on the norm of the second derivative

    ydpb = dfdxb + dfdub*fbnd

    ! define the tolerance parameter upon which the starting step
    ! size is to be based.  a value in the middle of the error
    ! tolerance range is selected.

    tolmin = big
    tolsum = zero
    do k = 1, me%n
        tolexp = log10(etol(k))
        tolmin = min(tolmin,tolexp)
        tolsum = tolsum + tolexp
    end do
    tolp = 10.0_wp**(0.5_wp*(tolsum/me%n + tolmin)/(morder+1))

    ! compute a starting step size based on the above first and
    ! second derivative information
    !
    ! restrict the step length to be not bigger
    ! than abs(b-a). (unless b is too close to a)

    h = absdx

    if (ydpb == zero .and. fbnd == zero) then
        ! both first derivative term (fbnd) and second
        ! derivative term (ydpb) are zero
        if (tolp < 1.0_wp) h = absdx*tolp
    elseif (ydpb == zero) then
        ! only second derivative term (ydpb) is zero
        if (tolp < fbnd*absdx) h = tolp/fbnd
    else
        ! second derivative term (ydpb) is non-zero
        srydpb = sqrt(0.5_wp*ydpb)
        if (tolp < srydpb*absdx) h = tolp/srydpb
    end if

    ! further restrict the step length to be not bigger than  1/dfdub
    if (h*dfdub > 1.0_wp) h = 1.0_wp/dfdub

    ! finally, restrict the step length to be not
    ! smaller than 100*small*abs(a). however, if
    ! a=0. and the computed h underflowed to zero,
    ! the algorithm returns small*abs(b) for the step length.
    h = max(h,100.0_wp*small*abs(a))
    if (h == zero) h = small*abs(b)

    ! now set direction of integration
    h = sign(h,dx)

    end subroutine hstart
!*****************************************************************************************

!*****************************************************************************************
!>
!  computation of an initial step size guess
!
!@note This routine is from dop853. It was modified for this module.

    function hinit(me,x,y,posneg,f0,hmax,atol,rtol)

    implicit none

    class(rk_variable_step_class),intent(inout) :: me
    real(wp),intent(in)               :: x
    real(wp),dimension(:),intent(in)  :: y       !! dimension(n)
    real(wp),intent(in)               :: posneg  !! posneg = sign(1.0_wp,xend-x)
    real(wp),dimension(:),intent(in)  :: f0      !! dimension(n)
    real(wp),intent(in)               :: hmax
    real(wp),dimension(:),intent(in)  :: atol
    real(wp),dimension(:),intent(in)  :: rtol

    real(wp) :: der12,der2,dnf,dny,h,h1,hinit,sk
    integer :: i
    integer :: iord  !! order of the method
    real(wp),dimension(me%n) :: f1,y1

    iord = me%p

    ! compute a first guess for explicit euler as
    !   h = 0.01 * norm (y0) / norm (f0)
    ! the increment for explicit euler is small
    ! compared to the solution
    dnf = zero
    dny = zero
    do i = 1 , me%n
        sk = atol(i) + rtol(i)*abs(y(i))
        dnf = dnf + (f0(i)/sk)**2
        dny = dny + (y(i)/sk)**2
    end do
    if ( dnf<=1.0e-10_wp .or. dny<=1.0e-10_wp ) then
        h = 1.0e-6_wp
    else
        h = sqrt(dny/dnf)*0.01_wp
    end if
    h = min(h,hmax)
    h = sign(h,posneg)
    ! perform an explicit euler step
    do i = 1 , me%n
        y1(i) = y(i) + h*f0(i)
    end do
    call me%f(x+h,y1,f1)
    ! estimate the second derivative of the solution
    der2 = zero
    do i = 1 , me%n
        sk = atol(i) + rtol(i)*abs(y(i))
        der2 = der2 + ((f1(i)-f0(i))/sk)**2
    end do
    der2 = sqrt(der2)/h
    ! step size is computed such that
    !  h**iord * max ( norm (f0), norm (der2)) = 0.01
    der12 = max(abs(der2),sqrt(dnf))
    if ( der12<=1.0e-15_wp ) then
        h1 = max(1.0e-6_wp,abs(h)*1.0e-3_wp)
    else
        h1 = (0.01_wp/der12)**(1.0_wp/iord)
    end if

    h = min(100.0_wp*abs(h),h1,hmax)
    hinit = sign(h,posneg)

    end function hinit
!*****************************************************************************************




!*****************************************************************************************
    end module runge_kutta_module
!*****************************************************************************************
