
    module test_support

    use rklib_module, wp => rk_module_rk

    implicit none

    private

    real(wp),parameter,public :: deg2rad = acos(-1.0_wp) / 180.0_wp

    public :: orbital_elements_to_rv

    contains
!*****************************************************************************************

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
            aop_tmp = 0.0_wp
        else
            aop_tmp = aop
        end if

        if (equatorial) then   ! node undefined
            raan_tmp = 0.0_wp
        else
            raan_tmp = raan
        end if

        ! perifocal position and velocity:
        ctru   = cos(tru)
        stru   = sin(tru)
        r_pqw  = [ctru, stru] * p/(1.0_wp+ecc*ctru)
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
        equatorial = (1.0_wp - abs(cos(inc))) < tol  ! 0 or 180 deg

        end subroutine orbit_check
    !*****************************************************************************************

    end module test_support