
    module rk78_module

        use runge_kutta_module, wp => rk_module_rk

        implicit none

        private

        public :: rk78

        abstract interface
            subroutine func(t,x,xdot)
                import :: wp
                implicit none
                real(wp),intent(in)               :: t
                real(wp),dimension(:),intent(in)  :: x
                real(wp),dimension(*),intent(out) :: xdot
            end subroutine func
        end interface

    contains

    !!  runge-kutta 7(8) method as given by the fehlberg coefficients
    !!
    !!  Note: user must call rk78cn before calling rk78.
    !!  Based on: https://oceancolor.gsfc.nasa.gov/docs/ocssw/rk78_8f_source.html
    !!  References: jpl em 312/87-153, 20 april 1987
    !!              nasa tr r-287, october 1968

    subroutine rk78(der,t,tout,n,x,h,relerr,abserr)

    implicit none

    procedure(func) :: der
    real(wp),intent(inout) :: t
    real(wp),intent(in)    :: tout
    integer,intent(in)     :: n
    real(wp),intent(inout) :: x(n)
    real(wp),intent(in)    :: h
    real(wp),intent(in)    :: relerr
    real(wp),intent(in)    :: abserr

    real(wp) :: alph(13) , beta(13,12) , ch(13) , f(n,13) , xdum(n)
    real(wp) :: const , dt , err , sdt , tdum , temp, ter , tol
    integer :: i , j , k , kk , l , nrej , nstp

    real(wp),parameter :: zero  = 0.0_wp
    real(wp),parameter :: one   = 1.0_wp
    real(wp),parameter :: two   = 2.0_wp
    real(wp),parameter :: three = 3.0_wp
    integer,parameter :: norder = 8
    integer,parameter :: ntimes = 13
    real(wp),parameter :: ordrcp = one/norder

    real(wp),parameter :: fact = 0.7_wp !! fact is a scaling factor to reduce the estimated stepsize
                                        !! to avoid excessive number of rejection.
    real(wp),parameter :: small =  epsilon(1.0_wp) !1.0e-8_wp !! small is a small number compared to the size of t

    logical,save :: initialized = .false.

    if (.not. initialized) then
        call rk78cn()
        initialized = .true.
    end if

    ! initialization

    nstp = 0
    nrej = 0
    sdt = sign(one,tout-t)
    dt = abs(h)*sdt

    do

        ! start integration and save initial info

        nstp = nstp + 1
        tdum = t
        do i = 1 , n
            xdum(i) = x(i)
        enddo

        ! check for t+dt passing tout

        if ( abs(dt)>abs(tout-t) ) dt = tout - t

        ! check for reaching tout

        if ( abs(t-tout)<small ) return

        ! start function evaluations

        call der(t,x,f(1,1))
        do k = 2 , ntimes
            kk = k - 1
            do i = 1 , n
                temp = zero
                do j = 1 , kk
                    temp = temp + beta(k,j)*f(i,j)
                enddo
                x(i) = xdum(i) + dt*temp
            enddo
            t = tdum + alph(k)*dt
            call der(t,x,f(1,k))
        enddo

        ! compute new state

        do i = 1 , n
            temp = zero
            do l = 1 , ntimes
                temp = temp + ch(l)*f(i,l)
            enddo
            x(i) = xdum(i) + dt*temp
        enddo

        ! start local truncation error computation

        err = relerr
        do i = 1 , n
            ter = abs((f(i,1)+f(i,11)-f(i,12)-f(i,13))*ch(12)*dt)
            tol = abs(x(i))*relerr + abserr
            const = ter/tol
            if ( const>err ) err = const
        enddo

        ! estimate new step-size
        dt = fact*dt*(one/err)**ordrcp

        ! reject last step if er<one

        if ( err>one ) then  ! ter > tol
            t = tdum
            do i = 1 , n
                x(i) = xdum(i)
            enddo
            nrej = nrej + 1
            nstp = nstp - 1
        endif

    end do

    contains

        subroutine rk78cn()

        !! computes the fehlberg coefficients for a runge-kutta 78 integrator
        !! this routine must be called before calling routine [[rk78]]

        ! initialize
        beta = zero
        alph = zero
        ch = zero

        ch(6)  = 34.0_wp/105.0_wp
        ch(7)  = 9.0_wp/35.0_wp
        ch(8)  = ch(7)
        ch(9)  = 9.0_wp/280.0_wp
        ch(10) = ch(9)
        ch(12) = 41.0_wp/840.0_wp
        ch(13) = ch(12)

        alph(2)  = two/27.0_wp
        alph(3)  = one/9.0_wp
        alph(4)  = one/6.0_wp
        alph(5)  = 5.0_wp/12.0_wp
        alph(6)  = 0.5_wp
        alph(7)  = 5.0_wp/6.0_wp
        alph(8)  = one/6.0_wp
        alph(9)  = two/three
        alph(10) = one/three
        alph(11) = one
        alph(13) = one

        beta(2,1)   = two/27.0_wp
        beta(3,1)   = one/36.0_wp
        beta(4,1)   = one/24.0_wp
        beta(5,1)   = 5.0_wp/12.0_wp
        beta(6,1)   = 0.5e-1_wp
        beta(7,1)   = -25.0_wp/108.0_wp
        beta(8,1)   = 31.0_wp/300.0_wp
        beta(9,1)   = two
        beta(10,1)  = -91.0_wp/108.0_wp
        beta(11,1)  = 2383.0_wp/4100.0_wp
        beta(12,1)  = three/205.0_wp
        beta(13,1)  = -1777.0_wp/4100.0_wp
        beta(3,2)   = one/12.0_wp
        beta(4,3)   = one/8.0_wp
        beta(5,3)   = -25.0_wp/16.0_wp
        beta(5,4)   = -beta(5,3)
        beta(6,4)   = 0.25_wp
        beta(7,4)   = 125.0_wp/108.0_wp
        beta(9,4)   = -53.0_wp/6.0_wp
        beta(10,4)  = 23.0_wp/108.0_wp
        beta(11,4)  = -341.0_wp/164.0_wp
        beta(13,4)  = beta(11,4)
        beta(6,5)   = 0.2_wp
        beta(7,5)   = -65.0_wp/27.0_wp
        beta(8,5)   = 61.0_wp/225.0_wp
        beta(9,5)   = 704.0_wp/45.0_wp
        beta(10,5)  = -976.0_wp/135.0_wp
        beta(11,5)  = 4496.0_wp/1025.0_wp
        beta(13,5)  = beta(11,5)
        beta(7,6)   = 125.0_wp/54.0_wp
        beta(8,6)   = -two/9.0_wp
        beta(9,6)   = -107.0_wp/9.0_wp
        beta(10,6)  = 311.0_wp/54.0_wp
        beta(11,6)  = -301.0_wp/82.0_wp
        beta(12,6)  = -6.0_wp/41.0_wp
        beta(13,6)  = -289.0_wp/82.0_wp
        beta(8,7)   = 13.0_wp/900.0_wp
        beta(9,7)   = 67.0_wp/90.0_wp
        beta(10,7)  = -19.0_wp/60.0_wp
        beta(11,7)  = 2133.0_wp/4100.0_wp
        beta(12,7)  = -three/205.0_wp
        beta(13,7)  = 2193.0_wp/4100.0_wp
        beta(9,8)   = three
        beta(10,8)  = 17.0_wp/6.0_wp
        beta(11,8)  = 45.0_wp/82.0_wp
        beta(12,8)  = -three/41.0_wp
        beta(13,8)  = 51.0_wp/82.0_wp
        beta(10,9)  = -one/12.0_wp
        beta(11,9)  = 45.0_wp/164.0_wp
        beta(12,9)  = three/41.0_wp
        beta(13,9)  = 33.0_wp/164.0_wp
        beta(11,10) = 18.0_wp/41.0_wp
        beta(12,10) = 6.0_wp/41.0_wp
        beta(13,10) = 12.0_wp/41.0_wp
        beta(13,12) = one

        end subroutine rk78cn

    end subroutine rk78

    end module rk78_module
