
    module test_support

    use rklib_module, wp => rk_module_rk

    implicit none

    private

    real(wp),parameter,public :: deg2rad = acos(-1.0_wp) / 180.0_wp

    public :: orbital_elements_to_rv
    public :: hslToRgb

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

    !*****************************************************************************************
    !>
    !  Converts an HSL color value to RGB.
    !
    !  See: https://stackoverflow.com/questions/2353211/hsl-to-rgb-color-conversion

       function hslToRgb(hsl) result(rgb)

        real(wp),dimension(3),intent(in) :: hsl !! [h,s,l] in range [0.0, 1.0]
        integer,dimension(3) :: rgb !! [r,g,b] in range [0, 255]

        real(wp) :: h,s,l,r,g,b,p,q

        h = hsl(1)
        s = hsl(2)
        l = hsl(3)

        if (l < 0.5_wp) then
            q = l * (1 + s)
        else
            q = l + s - l * s
        end if

        p = 2.0_wp * l - q

        if (s == 0.0_wp) then
            r = l; g = l; b = l  ! achromatic
        else
            r = hueToRgb(p, q, h + 1.0_wp/3.0_wp)
            g = hueToRgb(p, q, h)
            b = hueToRgb(p, q, h - 1.0_wp/3.0_wp)
        end if
        rgb = [to255(r), to255(g), to255(b)]

        contains
            integer function to255(v)
                !! Helper method that converts hue to rgb
                real(wp),intent(in) :: v
                to255 = int(min(255.0_wp,256.0_wp*v))
            end function to255

            function hueToRgb(p, q, t) result(r)
                real(wp),value :: p,q,t
                real(wp) :: r
                if (t < 0.0_wp) t = t + 1.0_wp
                if (t > 1.0_wp) t = t - 1.0_wp
                if (t < 1.0_wp/6.0_wp) then
                    r = p + (q - p) * 6.0_wp * t
                else if (t < 1.0_wp/2.0_wp) then
                    r = q
                else if (t < 2.0_wp/3.0_wp) then
                    r = p + (q - p) * (2.0_wp/3.0_wp - t) * 6.0_wp
                else
                    r = p
                end if
            end function hueToRgb

        end function hslToRgb

    end module test_support