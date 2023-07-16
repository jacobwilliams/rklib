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

    module rklib_module

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

    integer,parameter :: max_error_len = 100 !! max size of error message strings
    integer,parameter,public :: RKLIB_ERROR_TOO_MANY_STEPS       = -10
    integer,parameter,public :: RKLIB_ERROR_INVALID_RTOL_SIZE    = -9
    integer,parameter,public :: RKLIB_ERROR_INVALID_ATOL_SIZE    = -8
    integer,parameter,public :: RKLIB_ERROR_INVALID_H            = -7
    integer,parameter,public :: RKLIB_ERROR_USER_STOPPED         = -6
    integer,parameter,public :: RKLIB_ERROR_MIN_STEP_SIZE        = -5
    integer,parameter,public :: RKLIB_ERROR_TOO_MANY_REDUCTIONS  = -4
    integer,parameter,public :: RKLIB_ERROR_INVALID_HINIT_METHOD = -3
    integer,parameter,public :: RKLIB_ERROR_G_NOT_ASSOCIATED     = -2
    integer,parameter,public :: RKLIB_ERROR_F_NOT_ASSOCIATED     = -1
    integer,parameter,public :: RKLIB_ERROR_NONE                 =  0
    character(len=max_error_len),dimension(RKLIB_ERROR_TOO_MANY_STEPS:RKLIB_ERROR_NONE),parameter :: &
        rklib_error_messages = [&
            'Too many steps                              ', & ! -10
            'Invalid size for rtol array                 ', & ! -9
            'Invalid size for atol array                 ', & ! -8
            'Step size cannot be zero                    ', & ! -7
            'User stopped the integration                ', & ! -6
            'Too many attempts to reduce step size       ', & ! -5
            'Invalid initial step size estimation method ', & ! -4
            'The function procedure f is not associated  ', & ! -3
            'The event procedure g is not associated     ', & ! -2
            'Minimum step size reached                   ', & ! -1
            'Success                                     ' ]  !  0
            !! Status message strings that go with the status codes.
            !! The index in this array is the `istatus` code.

    type,public :: rklib_properties
        !! Properties of an RK method.
        integer :: order = 0 !! order of the method
        integer :: number_of_stages = 0 !! number of stages
        logical :: fsal = .false. !! if it is a FSAL method
        logical :: low_storage = .false. !! if it is a LS method
        logical :: strong_stability_preserving = .false. !! if it is a SSP method
        integer :: number_of_registers = 0 !! number of `f` vectors used
        real(wp) :: cfl = zero !! Courant-Friedrichs-Lewy number
        character(len=:),allocatable :: short_name !! short version of the method name
        character(len=:),allocatable :: long_name !! longer description of the method
    end type rklib_properties

    type,public :: stepsize_class

        !! Algorithms for adjusting the step size for variable-step
        !! Runge-Kutta integrators.

        private

        logical :: fixed_step_mode = .false. !! if true, then the method runs in
                                             !! fixed step mode with not error estimation

        real(wp) :: hmax           = 1.0e+6_wp  !! maximum allowed step size
        real(wp) :: hmin           = 1.0e-6_wp  !! minimum allowed step size
        real(wp) :: hfactor_reject = 0.5_wp     !! minimum allowed factor for decreasing step size after rejected step
        real(wp) :: hfactor_accept = 2.0_wp     !! maximum allowed factor for increasing step size after accepted step
        integer  :: accept_mode    = 2          !! method to determine if step is accepted [1,2]
        integer  :: max_attempts   = 10000      !! maximum number of attempts to decrease step size before giving up

        ! for the `hfactor` equation:
        logical  :: relative_err      = .false. !! to use `tol*h` in the `hfactor` equation
        real(wp) :: safety_factor     = 0.9_wp  !! for `hfactor` equation (>0)
        integer  :: p_exponent_offset = 1       !! `p` + this value in the exponent (0 or 1)

        procedure(norm_func),nopass,pointer :: norm => maxval_func
            !! routine for computing the norm of the state

        contains

        private

        procedure,public :: initialize => stepsize_class_constructor
        procedure,public :: compute_stepsize
        procedure,public :: destroy => destroy_stepsize_class

    end type stepsize_class

    type,abstract,public :: rk_class

        !! main integration class

        private

        integer :: istatus = 0 !! status code
        logical :: stopped = .false. !! if user has stopped the integration in `f` or `report`.
        integer :: num_steps = 0 !! number of accepted steps taken
        integer :: max_number_of_steps  = huge(1) !! maximum number of steps to take
        integer :: report_rate = 1 !! how often to call report function:
                                   !! `0` : no reporting (same as not associating `report`),
                                   !! `1` : report every point,
                                   !! `2` : report every other point, etc.
                                   !! The first and last point are always reported.

        logical :: stop_on_errors = .false. !! if true, then errors will stop the program
        integer :: n = 0  !! user specified number of variables
        procedure(deriv_func),pointer  :: f      => null()  !! user-specified derivative function
        procedure(report_func),pointer :: report => null()  !! user-specified report function
        procedure(event_func),pointer  :: g      => null()  !! event function (stop when this is zero)

        real(wp),dimension(:,:),allocatable :: funcs !! matrix to store the function
                                                     !! evalutaions in the step function.
                                                     !! this will be size (`n` x `number_of_registers`)

        contains

        private

        procedure,public :: destroy !! destructor
        procedure,public :: stop => rk_class_stop !! user-callable method to stop the integration
        procedure,public :: status => rk_class_status !! get status code and message
        procedure,public :: failed

        procedure :: init => initialize_rk_class
        procedure :: begin => begin_integration_rk_class
        procedure :: raise_exception
        procedure :: clear_exception
        procedure :: export_point
        procedure(begin_func),deferred :: begin_integration
        procedure(properties_func),deferred,public :: properties

    end type rk_class

    type,extends(rk_class),abstract,public :: rk_fixed_step_class

        !! fixed step size class

        private

        contains

        private

        procedure(step_func_fixed),deferred :: step !! the step routine for the rk method
        procedure,public :: initialize => initialize_fixed_step !! initialize the class (set n,f, and report)

        procedure,public :: integrate => integrate_fixed_step
        procedure,public :: integrate_to_event => integrate_to_event_fixed_step
        procedure :: begin_integration => begin_integration_rk_fixed_step_class

    end type rk_fixed_step_class

    type,extends(rk_class),abstract,public :: rk_variable_step_class

        !! Main integration class for variable step size Runge-Kutta methods

        private

        type(stepsize_class) :: stepsize_method  !! the method for varying the step size

        real(wp),dimension(:),allocatable :: rtol  !! relative tolerance (`size(n)`)
        real(wp),dimension(:),allocatable :: atol  !! absolute tolerance (`size(n)`)

        integer :: hinit_method = 1 !! if automatically computing the inital step size, which
                                    !! method to use. 1 = `hstart`, 2 = `hinit`.

        integer :: num_rejected_steps = 0 !! number of rejected steps
        real(wp) :: last_accepted_step_size = zero !! the last accepted step size `dt` from the integration
                                                   !! (positive or negative)

        contains

        private

        procedure(step_func_variable),deferred :: step !! the step routine for the rk method
        procedure,public :: initialize => initialize_variable_step  !! initialize the class (set n,f, and report)
        procedure,public :: integrate => integrate_variable_step
        procedure,public :: integrate_to_event => integrate_to_event_variable_step
        procedure,public :: info => info_variable_step

        procedure :: hstart  !! for automatically computing the initial step size [this is from DDEABM]
        procedure :: hinit   !! for automatically computing the initial step size [this is from DOP853]
        procedure :: begin_integration => begin_integration_rk_variable_step_class
        procedure :: compute_initial_step
        procedure :: order !! returns `p`, the order of the method

    end type rk_variable_step_class

    type,extends(rk_variable_step_class),abstract,public :: rk_variable_step_fsal_class
        !! a variable step method with the "first same as last" (FSAL) property.
        !! Cache the last `f` and `x` vectors to use for the next step.
        !!
        !! The assumption is that the nature of the
        !! function has not changed since the last step.
        !! If it has, the user would need to manually call [[destroy_fsal_cache]]
        !! so that the previous point was not reused.
        private
        real(wp),allocatable :: t_saved !! cached `t`
        real(wp),dimension(:),allocatable :: x_saved  !! cached `x`
        real(wp),dimension(:),allocatable :: f_saved !! cached `f`
        contains
        private
        procedure,public :: destroy_fsal_cache
        procedure,public :: check_fsal_cache
        procedure,public :: set_fsal_cache
    end type rk_variable_step_fsal_class

#include "rklib_fixed_classes.i90"
#include "rklib_variable_classes.i90"

    abstract interface

        subroutine begin_func(me)
            !! routine called before integration begins
            !! to set up internal variables.
            import :: rk_class
            class(rk_class),intent(inout) :: me
        end subroutine begin_func

        subroutine deriv_func(me,t,x,xdot)
        !! derivative function
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)     :: me
            real(wp),intent(in)               :: t    !! time
            real(wp),dimension(:),intent(in)  :: x    !! state vector
            real(wp),dimension(:),intent(out) :: xdot !! derivative of state vector
        end subroutine deriv_func

        subroutine event_func(me,t,x,g)
        !! event function
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)    :: me
            real(wp),intent(in)              :: t !! time
            real(wp),dimension(:),intent(in) :: x !! state vector
            real(wp),intent(out)             :: g !! g(t,x). The goal is to stop the integration when g=0.
        end subroutine event_func

        subroutine report_func(me,t,x)
        !! report function
        import :: rk_class,wp
        implicit none
            class(rk_class),intent(inout)    :: me
            real(wp),intent(in)              :: t !! time
            real(wp),dimension(:),intent(in) :: x !! state vector
        end subroutine report_func

        subroutine step_func_fixed(me,t,x,h,xf)
        !! rk step function for the fixed-step methods.
        import :: rk_fixed_step_class,wp
        implicit none
            class(rk_fixed_step_class),intent(inout) :: me
            real(wp),intent(in)                      :: t  !! initial time
            real(wp),dimension(me%n),intent(in)      :: x  !! initial state vector
            real(wp),intent(in)                      :: h  !! time step \( |\Delta t| \)
            real(wp),dimension(me%n),intent(out)     :: xf !! final state vector
        end subroutine step_func_fixed

        subroutine step_func_variable(me,t,x,h,xf,xerr)
        !! rk step function for the variable-step methods.
        import :: rk_variable_step_class,wp
        implicit none
            class(rk_variable_step_class),intent(inout) :: me
            real(wp),intent(in)                         :: t    !! initial time
            real(wp),dimension(me%n),intent(in)         :: x    !! initial state vector
            real(wp),intent(in)                         :: h    !! time step \( |\Delta t| \)
            real(wp),dimension(me%n),intent(out)        :: xf   !! final state vector
            real(wp),dimension(me%n),intent(out)        :: xerr !! truncation error estimate
        end subroutine step_func_variable

        pure function norm_func(x) result(xmag)
        !! Vector norm function. Must return a value \( \ge 0 \).
        import :: wp
        implicit none
            real(wp),dimension(:),intent(in) :: x    !! a vector
            real(wp)                         :: xmag !! the magnitude of the vector
        end function norm_func

        pure function properties_func(me) result(p)
        !! Returns the properties of the method.
        import :: rk_class,rklib_properties
        implicit none
            class(rk_class),intent(in) :: me
            type(rklib_properties) :: p !! properties of the method
        end function properties_func

    end interface

    ! submodule procedures:
    interface
#include "rklib_fixed_step_interfaces.i90"
#include "rklib_variable_step_interfaces.i90"
#include "rklib_fixed_property_interfaces.i90"
#include "rklib_variable_property_interfaces.i90"
    end interface

    ! public routines:
    public :: norm2_func,maxval_func

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the order of the RK method

    pure function order(me) result(p)
        class(rk_variable_step_class),intent(in) :: me
        integer :: p !! order of the method
        type(rklib_properties) :: properties
        properties = me%properties()
        p = properties%order
    end function order
!*****************************************************************************************

!*****************************************************************************************
!>
!  Clear any exception.

    subroutine clear_exception(me)
        class(rk_class),intent(inout) :: me
        me%istatus = RKLIB_ERROR_NONE
    end subroutine clear_exception
!*****************************************************************************************

!*****************************************************************************************
!>
!  Raise an exception.

    subroutine raise_exception(me, error_code)
        class(rk_class),intent(inout) :: me
        integer,intent(in) :: error_code !! the error to raise

        me%istatus = error_code

        if (error_code<0 .and. me%stop_on_errors) then
            error stop trim(rklib_error_messages(error_code))
        end if
    end subroutine raise_exception
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for [[rk_class]].

    subroutine destroy(me)
        class(rk_class),intent(out) :: me
    end subroutine destroy
!*****************************************************************************************

!*****************************************************************************************
!>
!  User-callable method to stop the integration.

    subroutine rk_class_stop(me)
        class(rk_class),intent(inout) :: me
        me%stopped = .true.
    end subroutine rk_class_stop
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns true if there was an error.
!  Can use [[rk_class_status]] to get more info.

    logical function failed(me)
        class(rk_class),intent(in) :: me
        failed = me%istatus < 0
    end function failed
!*****************************************************************************************

!*****************************************************************************************
!>
!  Get the status of an integration.

    subroutine rk_class_status(me,istatus,message)
        class(rk_class),intent(in) :: me
        integer,intent(out),optional :: istatus !! status code (`<0` means an error)
        character(len=:),allocatable,intent(out),optional :: message !! status message
        if (present(istatus)) istatus = me%istatus
        if (present(message)) message = trim(rklib_error_messages(me%istatus))
    end subroutine rk_class_status
!*****************************************************************************************

!*****************************************************************************************
!>
!  Wrapper for exporting points during integration.

    subroutine export_point(me,t,x,first_or_last)
        class(rk_class),intent(inout) :: me
        real(wp),intent(in) :: t
        real(wp),dimension(:),intent(in) :: x
        logical,intent(in),optional :: first_or_last  !! if this is the first or
                                                      !! last point (always reported)

        logical :: export !! if the point is to be exported

        if (associated(me%report) .and. me%report_rate > 0) then

            export = .false.
            if (present(first_or_last)) then
                ! always report first and last step
                if (first_or_last) export = .true.
            end if

            if (.not. export) then
                ! report steps at user-specified rate
                export = modulo(me%num_steps, me%report_rate) == 0
            end if

            if (export) call me%report(t,x)

        end if

    end subroutine export_point
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destructor for the FSAL variables.

    subroutine destroy_fsal_cache(me)
        class(rk_variable_step_fsal_class),intent(inout) :: me
        if (allocated(me%t_saved)) deallocate(me%t_saved)
        if (allocated(me%x_saved)) deallocate(me%x_saved)
        if (allocated(me%f_saved)) deallocate(me%f_saved)
    end subroutine destroy_fsal_cache
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the FSAL cache.

    subroutine check_fsal_cache(me,t,x,f)
        class(rk_variable_step_fsal_class),intent(inout) :: me
        real(wp),intent(in) :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f

        logical :: fsal !! if we can avoid a step due to first-same-as-last

        fsal = .false.
        if (allocated(me%x_saved)) then
            if (size(x) == size(me%x_saved)) then
                fsal = all(x==me%x_saved) .and. t==me%t_saved
            end if
        end if
        if (fsal) then
            f = me%f_saved
        else
            call me%f(t, x, f)
        end if

    end subroutine check_fsal_cache
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the function and add it to the FSAL cache.

    subroutine set_fsal_cache(me,t,x,f)
        class(rk_variable_step_fsal_class),intent(inout) :: me
        real(wp),intent(in) :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: f

        call me%f(t,x,f)

        me%t_saved = t
        me%x_saved = x
        me%f_saved = f

    end subroutine set_fsal_cache
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the [[rk_class]].

    subroutine initialize_rk_class(me,n,f,report,g,stop_on_errors,&
                                   max_number_of_steps,report_rate)

    implicit none

    class(rk_class),intent(inout)   :: me
    integer,intent(in)              :: n       !! number of variables
    procedure(deriv_func)           :: f       !! derivative function
    procedure(event_func),optional  :: g       !! for stopping at an event
    procedure(report_func),optional :: report  !! for reporting the steps
    logical,intent(in),optional     :: stop_on_errors !! stop the program for
                                                      !! any errors (default is False)
    integer,intent(in),optional     :: max_number_of_steps !! max number of steps allowed
    integer,intent(in),optional     :: report_rate !! how often to call report function:
                                                   !! `0` : no reporting (same as not associating `report`),
                                                   !! `1` : report every point,
                                                   !! `2` : report every other point, etc.
                                                   !! The first and last point are always reported.

    type(rklib_properties) :: props !! to get the method properties

    call me%destroy()

    me%n = n
    me%f => f
    if (present(report)) me%report => report
    if (present(g))      me%g      => g
    if (present(stop_on_errors)) me%stop_on_errors = stop_on_errors
    if (present(max_number_of_steps)) me%max_number_of_steps = abs(max_number_of_steps)
    if (present(report_rate)) me%report_rate = abs(report_rate)

    ! allocate the registers:
    props = me%properties()
    if (allocated(me%funcs)) deallocate(me%funcs)
    allocate(me%funcs(n, props%number_of_registers))
    me%funcs = zero

    ! reset internal variables:
    me%num_steps = 0
    me%stopped = .false.

    end subroutine initialize_rk_class
!*****************************************************************************************

!*****************************************************************************************
!>
!  Begin an integration.

    subroutine begin_integration_rk_class(me)
        class(rk_class),intent(inout) :: me
        call me%clear_exception()
        me%num_steps = 0
        me%stopped = .false.
    end subroutine begin_integration_rk_class
!*****************************************************************************************

!*****************************************************************************************
!>
!  Begin a [[rk_fixed_step_class]] integration.

    subroutine begin_integration_rk_fixed_step_class(me)
    class(rk_fixed_step_class),intent(inout) :: me
    call me%begin() ! all we need is base method here.
    end subroutine begin_integration_rk_fixed_step_class
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the [[rk_fixed_step_class]].

    subroutine initialize_fixed_step(me,n,f,report,g,stop_on_errors,&
                                     max_number_of_steps,report_rate)

    implicit none

    class(rk_fixed_step_class),intent(inout)   :: me
    integer,intent(in)              :: n       !! number of variables
    procedure(deriv_func)           :: f       !! derivative function
    procedure(report_func),optional :: report  !! for reporting the steps
    procedure(event_func),optional  :: g       !! for stopping at an event
    logical,intent(in),optional     :: stop_on_errors !! stop the program for
                                                      !! any errors (default is False)
    integer,intent(in),optional     :: max_number_of_steps !! max number of steps allowed
    integer,intent(in),optional     :: report_rate !! how often to call report function:
                                                   !! `0` : no reporting (same as not associating `report`),
                                                   !! `1` : report every point,
                                                   !! `2` : report every other point, etc.
                                                   !! The first and last point are always reported.

    ! base init all we need here:
    call me%init(n,f,report,g,stop_on_errors,max_number_of_steps,report_rate)

    end subroutine initialize_fixed_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Main integration routine for the [[rk_class]].

    subroutine integrate_fixed_step(me,t0,x0,h,tf,xf)

    implicit none

    class(rk_fixed_step_class),intent(inout) :: me
    real(wp),intent(in)               :: t0    !! initial time
    real(wp),dimension(:),intent(in)  :: x0    !! initial state
    real(wp),intent(in)               :: h     !! abs(time step)
    real(wp),intent(in)               :: tf    !! final time
    real(wp),dimension(:),intent(out) :: xf    !! final state

    real(wp) :: t  !! current time value
    real(wp) :: dt !! time step from `t` to `t2`
    real(wp) :: t2 !! time to step to from `t`
    real(wp),dimension(me%n) :: x !! state vector
    logical :: last !! if it is the last step

    if (.not. associated(me%f)) then
        call me%raise_exception(RKLIB_ERROR_F_NOT_ASSOCIATED)
        return
    end if
    if (abs(h)<=zero) then
        call me%raise_exception(RKLIB_ERROR_INVALID_H)
        return
    end if

    call me%begin_integration()

    call me%export_point(t0,x0,.true.)  !first point

    if (abs(h)>zero) then

        t = t0
        x = x0
        dt = sign(h,tf-t0)  !time step (correct sign)
        do
            t2 = t + dt
            last = ((dt>=zero .and. t2>=tf) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tf))         !
            if (last) dt = tf-t                     !
            call me%step(t,x,dt,xf)
            if (me%stopped) return
            me%num_steps = me%num_steps + 1
            if (me%num_steps > me%max_number_of_steps) then
                call me%raise_exception(RKLIB_ERROR_TOO_MANY_STEPS)
                return
            end if
            if (last) exit
            call me%export_point(t2,xf)   !intermediate point
            x = xf
            t = t2
        end do

    else
        xf = x0
    end if

    call me%export_point(tf,xf,.true.)   !last point

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

    class(rk_fixed_step_class),intent(inout) :: me
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
    real(wp),dimension(me%n) :: x !! state vector
    real(wp),dimension(me%n) :: g_xf !! state vector from the root finder
    logical :: first !! it is the first step
    logical :: last  !! it is the last step
    type(brent_solver) :: solver !! for the root finding problem
    integer :: iflag !! return flag from `solver`

    if (.not. associated(me%f)) then
        call me%raise_exception(RKLIB_ERROR_F_NOT_ASSOCIATED)
        return
    end if
    if (.not. associated(me%g)) then
        call me%raise_exception(RKLIB_ERROR_G_NOT_ASSOCIATED)
        return
    end if
    if (abs(h)<=zero) then
        call me%raise_exception(RKLIB_ERROR_INVALID_H)
        return
    end if

    call me%begin_integration()

    call me%export_point(t0,x0,.true.) !first point

    if (abs(t0-tmax)<=zero) then
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
            if (me%stopped) return
            me%num_steps = me%num_steps + 1
            if (me%num_steps > me%max_number_of_steps) then
                call me%raise_exception(RKLIB_ERROR_TOO_MANY_STEPS)
                return
            end if
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
                    call me%export_point(t2,xf)   !intermediate point
                    x = xf
                    t = t2
                    ga = gb
                end if

            elseif (ga*gb<=zero) then !there is a root somewhere on [t,t+dt]

                !find the root:
                call solver%initialize(solver_func)
                call solver%solve(zero,dt,dt_root,dum,iflag,fax=ga,fbx=gb)
                if (me%stopped) return
                t2 = t + dt_root
                gf = solver_func(solver,dt_root)
                if (me%stopped) return
                tf = t2
                xf = g_xf !computed in the solver function
                exit

            else  !no root yet, continue

                if (last) then  !exiting without having found a root
                    tf = t2
                    gf = gb
                    exit
                end if
                call me%export_point(t2,xf)   !intermediate point
                x = xf
                t = t2
                ga = gb

            end if

            if (first) first = .false.

        end do

    end if

    call me%export_point(t2,xf,.true.)   !last point

    contains

        function solver_func(this,delt) result(g)

        !! root solver function. The input is the `dt` offset from time `t`.

        implicit none

        class(root_solver),intent(inout) :: this
        real(wp),intent(in) :: delt  !! from [0 to `dt`]
        real(wp) :: g

        !take a step from t to t+delt and evaluate g function:
        call me%step(t,x,delt,g_xf)
        if (me%stopped) return
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

    pure subroutine stepsize_class_constructor(me,hmin,hmax,hfactor_reject,&
                        hfactor_accept,norm,accept_mode,relative_err,&
                        safety_factor,p_exponent_offset,max_attempts,&
                        fixed_step_mode)

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
    logical,intent(in),optional   :: fixed_step_mode    !! if true, then the method runs in
                                                        !! fixed step mode with not error estimation.
                                                        !! All the other inputs are ignored. Note that
                                                        !! this requires a `dt /= 0` input for the integrator.

    if (present(hmin))                me%hmin                = abs(hmin)
    if (present(hmax))                me%hmax                = abs(hmax)
    if (present(hfactor_reject))      me%hfactor_reject      = abs(hfactor_reject)
    if (present(hfactor_accept))      me%hfactor_accept      = abs(hfactor_accept)
    if (present(norm))                me%norm                => norm
    if (present(accept_mode))         me%accept_mode         = accept_mode
    if (present(max_attempts))        me%max_attempts        = max_attempts

    !if (present(compute_h_factor)) me%compute_h_factor => compute_h_factor
    if (present(relative_err     )) me%relative_err      = relative_err
    if (present(safety_factor    )) me%safety_factor     = abs(safety_factor    )
    if (present(p_exponent_offset)) me%p_exponent_offset = abs(p_exponent_offset)

    if (present(fixed_step_mode)) me%fixed_step_mode = fixed_step_mode

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
!>
!  Compute the new step size using the specific method.

    subroutine compute_stepsize(me,n,h,tol,err,p,hnew,accept)

    implicit none

    class(stepsize_class),intent(in) :: me
    integer,intent(in)               :: n      !! number of variables
    real(wp),intent(in)              :: h      !! current step size (<>0)
    real(wp),dimension(n),intent(in) :: tol    !! abs error tolerance (>0)
    real(wp),dimension(n),intent(in) :: err    !! truncation error estimate (>0)
    integer,intent(in)               :: p      !! order of the method
    real(wp),intent(out)             :: hnew   !! new step size (<>0)
    logical,intent(out)              :: accept !! if the step is accepted

    real(wp) :: e        !! exponent
    real(wp) :: hfactor  !! step size factor (>0)
    real(wp) :: max_err  !! max error for all the elements

    real(wp),parameter :: small = 10.0_wp * epsilon(1.0_wp) !! small error value

    if (me%fixed_step_mode) then
        ! do not adjust the step size
        accept = .true.
        hnew = h
    else

        if (all(err<=small)) then ! the error is extremely small

            hfactor = me%hfactor_accept
            accept = .true.

        else

            ! compute base factor based on the selected formula:
            if (me%relative_err) then
                max_err = me%norm(err/tol*abs(h))
            else
                max_err = me%norm(err/tol)
            end if

            e = 1.0_wp / real(p+me%p_exponent_offset,wp)
            hfactor = abs( me%safety_factor * (1.0_wp/max_err)**e )

            ! if the step is to be accepted:
            select case (me%accept_mode)
            case(1) !algorithm 17.12
                accept = hfactor >= 1.0_wp
            case(2) !algorithm 17.13
                accept = me%norm(err/tol) <= 1.0_wp
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

    end if

    end subroutine compute_stepsize
!*****************************************************************************************

!*****************************************************************************************
!>
!  Begin a [[rk_variable_step_class]] integration.

    subroutine begin_integration_rk_variable_step_class(me)
    class(rk_variable_step_class),intent(inout) :: me

    call me%begin() ! base

    ! variable step params:
    me%num_rejected_steps = 0
    me%last_accepted_step_size = zero
    select type (me)
    class is (rk_variable_step_fsal_class)
        call me%destroy_fsal_cache()
    end select

    end subroutine begin_integration_rk_variable_step_class
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the [[rk_variable_step_class]].

    subroutine initialize_variable_step(me,n,f,rtol,atol,stepsize_method,&
                                        hinit_method,report,g,stop_on_errors,&
                                        max_number_of_steps,report_rate)

    implicit none

    class(rk_variable_step_class),intent(inout) :: me
    integer,intent(in)                          :: n               !! number of equations
    procedure(deriv_func)                       :: f               !! derivative function
    real(wp),dimension(:),intent(in),optional   :: rtol            !! relative tolerance (if size=1,
                                                                   !! then same tol used for all
                                                                   !! equations). If not present, a default
                                                                   !! of `100*eps` is used
    real(wp),dimension(:),intent(in),optional   :: atol            !! absolute tolerance (if size=1,
                                                                   !! then same tol used for all
                                                                   !! equations). If not present, a default
                                                                   !! of `100*eps` is used
    type(stepsize_class),intent(in),optional    :: stepsize_method !! method for varying the step size
    integer,intent(in),optional                 :: hinit_method    !! which method (1 or 2) to use for
                                                                   !! automatic initial step size
                                                                   !! computation.
                                                                   !! 1 = use `hstart`, 2 = use `hinit`.
    procedure(report_func),optional             :: report          !! for reporting the steps
    procedure(event_func),optional              :: g               !! for stopping at an event
    logical,intent(in),optional :: stop_on_errors !! stop the program for
                                                  !! any errors (default is False)
    integer,intent(in),optional :: max_number_of_steps !! max number of steps allowed
    integer,intent(in),optional :: report_rate !! how often to call report function:
                                               !! `0` : no reporting (same as not associating `report`),
                                               !! `1` : report every point,
                                               !! `2` : report every other point, etc.
                                               !! The first and last point are always reported.

    real(wp),parameter :: default_tol = 100*epsilon(1.0_wp) !! if tols not specified

    ! base init:
    call me%init(n,f,report,g,stop_on_errors,max_number_of_steps,report_rate)

    ! variable-step specific inputs:
    if (allocated(me%rtol)) deallocate(me%rtol)
    if (allocated(me%atol)) deallocate(me%atol)
    allocate(me%rtol(n))
    allocate(me%atol(n))

    if (present(rtol)) then
        if (size(rtol)==1) then
            me%rtol = abs(rtol(1)) !use this for all equations
        else if (size(rtol)==n) then
            me%rtol = abs(rtol)
        else
            call me%raise_exception(RKLIB_ERROR_INVALID_RTOL_SIZE)
        end if
    else
        me%rtol = default_tol
    end if

    if (present(atol)) then
        if (size(atol)==1) then
            me%atol = abs(atol(1)) !use this for all equations
        else if (size(atol)==n) then
            me%atol = abs(atol)
        else
            call me%raise_exception(RKLIB_ERROR_INVALID_ATOL_SIZE)
        end if
    else
        me%atol = default_tol
    end if

    if (present(hinit_method)) then
        if (any(hinit_method == [1,2])) then
            me%hinit_method = hinit_method
        else
            call me%raise_exception(RKLIB_ERROR_INVALID_HINIT_METHOD)
            return
        end if
    end if
    if (present(stepsize_method)) me%stepsize_method = stepsize_method

    ! reset internal variables:
    me%num_rejected_steps = 0
    me%last_accepted_step_size = zero

    end subroutine initialize_variable_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Return some info about the integration.

    subroutine info_variable_step(me,num_steps,num_rejected_steps,last_accepted_step_size)

    implicit none

    class(rk_variable_step_class),intent(in) :: me
    integer,intent(out),optional :: num_steps !! number of steps taken
    integer,intent(out),optional :: num_rejected_steps  !! number of rejected steps
    real(wp),intent(out),optional  :: last_accepted_step_size !! the last accepted step size
                                                              !! `dt` from the integration
                                                              !! (positive or negative)

    if (present(num_steps)) num_steps = me%num_steps
    if (present(num_rejected_steps)) num_rejected_steps = me%num_rejected_steps
    if (present(last_accepted_step_size)) last_accepted_step_size = me%last_accepted_step_size

    end subroutine info_variable_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the initial step size.

    function compute_initial_step(me,t0,tf,x0,h0) result(dt)

        class(rk_variable_step_class),intent(inout) :: me
        real(wp),intent(in) :: h0 !! user-input initial step size (if zero, then one is computed)
        real(wp),intent(in) :: t0 !! initial time
        real(wp),intent(in) :: tf !! final time
        real(wp) :: dt !! step size to use

        real(wp),dimension(me%n) :: x0 !! initial state
        real(wp),dimension(me%n) :: etol !! tolerance vector
        real(wp),dimension(me%n) :: f0 !! initial derivative

        if (abs(h0)<=zero) then
            ! compute an appropriate initial step size:
            etol = me%rtol * me%stepsize_method%norm(x0) + me%atol
            call me%f(t0,x0,f0)  ! get initial dx/dt
            select case (me%hinit_method) ! value was checked in initialize_variable_step
            case(1); call me%hstart(t0,tf,x0,f0,etol,dt)
            case(2); dt = me%hinit(t0,x0,sign(1.0_wp,tf-t0),f0,&
                                me%stepsize_method%hmax,&
                                me%atol,me%rtol)
            end select
        else
            ! user-specified initial step size:
            dt = sign(h0,tf-t0)  ! (correct sign)
        end if

    end function compute_initial_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Main integration routine for the [[rk_variable_step_class]].

    subroutine integrate_variable_step(me,t0,x0,h,tf,xf)

    implicit none

    class(rk_variable_step_class),intent(inout)     :: me
    real(wp),intent(in)               :: t0    !! initial time
    real(wp),dimension(:),intent(in)  :: x0    !! initial state
    real(wp),intent(in)               :: h     !! initial abs(time step)
    real(wp),intent(in)               :: tf    !! final time
    real(wp),dimension(:),intent(out) :: xf    !! final state

    real(wp) :: t,dt,t2,dt_new
    real(wp),dimension(me%n) :: x,xerr,tol
    logical :: last !! it is the last step
    logical :: accept !! the step is accepted
    integer :: i !! max step size reduction attempts counter
    integer :: p !! order of the method

    if (.not. associated(me%f)) then
        call me%raise_exception(RKLIB_ERROR_F_NOT_ASSOCIATED)
        return
    end if

    call me%begin_integration()

    call me%export_point(t0,x0,.true.)  !first point

    if (abs(t0-tf)<=zero) then
        xf = x0
    else

        t = t0
        x = x0
        dt = me%compute_initial_step(t0,tf,x0,h)
        p = me%order()     !order of the method

        do
            t2 = t + dt
            last = ((dt>=zero .and. t2>=tf) .or. &  !adjust last time step
                    (dt<zero .and. t2<=tf))         !
            if (last) dt = tf-t                     !

            do i=0,me%stepsize_method%max_attempts

                ! take a step:
                call me%step(t,x,dt,xf,xerr)
                if (me%stopped) return

                if (me%stepsize_method%fixed_step_mode) then
                    ! don't adjust the step size
                    accept = .true.
                    me%last_accepted_step_size = dt ! save it [really only needs to be done once]
                else
                    ! evaluate error and compute new step size:
                    xerr = abs(xerr)
                    tol = me%rtol * abs(xf) + me%atol
                    call me%stepsize_method%compute_stepsize(me%n,dt,tol,xerr,p,dt_new,accept)
                    if (accept) me%last_accepted_step_size = dt ! save it
                    dt = dt_new
                end if

                if (accept) then
                    !accept this step
                    me%num_steps = me%num_steps + 1
                    if (me%num_steps > me%max_number_of_steps) then
                        call me%raise_exception(RKLIB_ERROR_TOO_MANY_STEPS)
                        return
                    end if
                    exit
                else
                    !step is rejected, repeat step with new dt
                    me%num_rejected_steps = me%num_rejected_steps + 1

                    !note: if we have reached the min step size, and the error
                    !is still too large, we can't proceed.
                    if (i>=me%stepsize_method%max_attempts) then
                        call me%raise_exception(RKLIB_ERROR_TOO_MANY_REDUCTIONS)
                        return
                    end if
                    if (abs(dt) < abs(me%stepsize_method%hmin)) then
                        call me%raise_exception(RKLIB_ERROR_MIN_STEP_SIZE)
                        return
                    end if

                    last = ((dt>=zero .and. (t+dt)>=tf) .or. &  !adjust last time step
                            (dt<zero .and. (t+dt)<=tf))         !
                    if (last) dt = tf-t                         !
                    t2 = t + dt

                end if

            end do

            if (last) exit
            call me%export_point(t2,xf)   !intermediate point
            x = xf
            t = t2
        end do

    end if

    call me%export_point(tf,xf,.true.)   !last point

    end subroutine integrate_variable_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Event-finding integration routine for the [[rk_variable_step_class]].
!  Integrates until g(t,x)=0, or until t=tf (whichever happens first).
!
!@note There are some efficiency improvements that could be made here.
!      This is a work in progress.

    subroutine integrate_to_event_variable_step(me,t0,x0,h,tmax,tol,tf,xf,gf)

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

    real(wp),dimension(me%n) :: x,g_xf
    real(wp),dimension(me%n) :: xerr !! truncation error estimate
    real(wp),dimension(me%n) :: stol
    integer :: i,p,iflag
    real(wp) :: t,dt,t2,ga,gb,dt_root,dum,dt_new
    logical :: first,last,accept
    type(brent_solver) :: solver

    if (.not. associated(me%f)) then
        call me%raise_exception(RKLIB_ERROR_F_NOT_ASSOCIATED)
        return
    end if
    if (.not. associated(me%g)) then
        call me%raise_exception(RKLIB_ERROR_G_NOT_ASSOCIATED)
        return
    end if

    call me%begin_integration()

    call me%export_point(t0,x0,.true.)  !first point

    if (abs(t0-tmax)<=zero) then
        xf = x0
        tf = t0
        call me%g(t0,x0,gf)
    else

        first = .true.
        t = t0
        x = x0
        call me%g(t,x,ga)  !evaluate event function
        dt = me%compute_initial_step(t0,tmax,x0,h)
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
                call me%step(t,x,dt,xf,xerr)
                if (me%stopped) return

                if (me%stepsize_method%fixed_step_mode) then
                    ! don't adjust the step size
                    accept = .true.
                    me%last_accepted_step_size = dt ! save it [really only needs to be done once]
                else
                    ! evaluate error and compute new step size:
                    xerr = abs(xerr)
                    stol = me%rtol * abs(xf) + me%atol
                    call me%stepsize_method%compute_stepsize(me%n,dt,stol,xerr,p,dt_new,accept)
                    if (accept) me%last_accepted_step_size = dt ! save it
                    dt = dt_new
                end if

                if (accept) then
                    !accept this step
                    me%num_steps = me%num_steps + 1
                    if (me%num_steps > me%max_number_of_steps) then
                        call me%raise_exception(RKLIB_ERROR_TOO_MANY_STEPS)
                        return
                    end if
                    exit
                else
                    !step is rejected, repeat step with new dt
                    me%num_rejected_steps = me%num_rejected_steps + 1

                    !note: if we have reached the min step size, and the error
                    !is still too large, we can't proceed.
                    if (i>=me%stepsize_method%max_attempts) then
                        call me%raise_exception(RKLIB_ERROR_TOO_MANY_REDUCTIONS)
                        return
                    end if
                    if (abs(dt) < abs(me%stepsize_method%hmin)) then
                        call me%raise_exception(RKLIB_ERROR_MIN_STEP_SIZE)
                        return
                    end if

                    last = ((dt>=zero .and. (t+dt)>=tmax) .or. &  !adjust last time step
                            (dt<zero .and. (t+dt)<=tmax))         !
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
                    call me%export_point(t2,xf)   !intermediate point
                    x = xf
                    t = t2
                    ga = gb
                end if

            elseif (ga*gb<=zero) then !there is a root somewhere on [t,t+dt]

                !find the root:
                call solver%initialize(solver_func, rtol=tol, atol=tol)
                call solver%solve(zero,dt,dt_root,dum,iflag,fax=ga,fbx=gb)
                if (me%stopped) return
                t2 = t + dt_root
                gf = solver_func(solver,dt_root)
                if (me%stopped) return
                tf = t2
                xf = g_xf !computed in the solver function
                exit

            else  !no root yet, continue

                if (last) then  !exiting without having found a root
                    tf = t2
                    gf = gb
                    exit
                end if
                call me%export_point(t2,xf)   !intermediate point
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

    call me%export_point(tf,xf,.true.)   !last point

    contains

        function solver_func(this,delt) result(g)

        !! root solver function. The input is the `dt` offset from time `t`.

        implicit none

        class(root_solver),intent(inout) :: this
        real(wp),intent(in) :: delt  !! from [0 to `dt`]
        real(wp) :: g

        real(wp),dimension(me%n) :: xerr !! truncation error estimate

        !take a step from t to t+delt and evaluate g function:
        ! [we don't check the error because we are within a
        !  step that was already accepted, so it should be ok]
        call me%step(t,x,delt,g_xf,xerr)
        if (me%stopped) return
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

    morder = me%order()
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

    iord = me%order()

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
    end module rklib_module
!*****************************************************************************************