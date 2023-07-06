  program rklib_example

  use rklib_module, wp => rk_module_rk
  use iso_fortran_env, only: output_unit

  implicit none

  integer,parameter :: n = 2 !! dimension of the system
  real(wp),parameter :: tol = 1.0e-12_wp !! integration tolerance
  real(wp),parameter :: x0 = 0.0_wp !! initial x value
  real(wp),parameter :: dx = 1.0_wp !! initial step size
  real(wp),parameter :: xf = 100.0_wp !! endpoint of integration
  real(wp),dimension(n),parameter :: y0 = [0.0_wp,0.1_wp] !! initial y value

  type(rktp86_class) :: prop
  real(wp),dimension(n) :: yf
  character(len=:),allocatable :: message

  call prop%initialize(n=n,f=fvpol,rtol=[tol],atol=[tol])
  call prop%integrate(x0,y0,dx,xf,yf)
  call prop%status(message=message)

  write (output_unit,'(A)') message
  write (output_unit,'(A,F7.2/,A,2E18.10)') &
              'xf =',xf ,'yf =',yf(1),yf(2)

contains

  subroutine fvpol(me,x,y,f)
  !! Right-hand side of van der Pol's equation

  implicit none

  class(rk_class),intent(inout)     :: me
  real(wp),intent(in)               :: x
  real(wp),dimension(:),intent(in)  :: y
  real(wp),dimension(:),intent(out) :: f

  real(wp),parameter :: mu  = 0.2_wp

  f(1) = y(2)
  f(2) = mu*(1.0_wp-y(1)**2)*y(2) - y(1)

  end subroutine fvpol

end program rklib_example
