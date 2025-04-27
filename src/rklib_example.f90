program rklib_example

  use rklib_module, wp => rk_module_rk
  use iso_fortran_env, only: output_unit

  implicit none

  integer,parameter :: n = 2 !! dimension of the system
  real(wp),parameter :: tol = 1.0e-12_wp !! integration tolerance
  real(wp),parameter :: t0 = 0.0_wp !! initial t value
  real(wp),parameter :: dt = 1.0_wp !! initial step size
  real(wp),parameter :: tf = 100.0_wp !! endpoint of integration
  real(wp),dimension(n),parameter :: x0 = [0.0_wp,0.1_wp] !! initial x value

  real(wp),dimension(n) :: xf !! final x value
  type(rktp86_class) :: prop
  character(len=:),allocatable :: message

  call prop%initialize(n=n,f=fvpol,rtol=[tol],atol=[tol])
  call prop%integrate(t0,x0,dt,tf,xf)
  call prop%status(message=message)

  write (output_unit,'(A)') message
  write (output_unit,'(A,F7.2/,A,2E18.10)') &
              'tf =',tf ,'xf =',xf(1),xf(2)

contains

  subroutine fvpol(me,t,x,f)
    !! Right-hand side of van der Pol equation

    class(rk_class),intent(inout)     :: me
    real(wp),intent(in)               :: t
    real(wp),dimension(:),intent(in)  :: x
    real(wp),dimension(:),intent(out) :: f

    f(1) = x(2)
    f(2) = 0.2_wp*(1.0_wp-x(1)**2)*x(2) - x(1)

  end subroutine fvpol

end program rklib_example