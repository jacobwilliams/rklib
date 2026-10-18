!*****************************************************************************************
!>
!  Dense output coefficients for the variable-step RK formulas.
!
!  Each `*_dense_coefficients` function returns the matrix `b` (`number_of_stages` x `degree`)
!  defining the interpolant over the step \( [t, t+h] \) in the form used in `DOPRI5`
!  and `DOP853` (see [[dense_coefficients]]). Coefficients published in the polynomial form
!  \( x(t + \theta h) = x(t) + h \sum_{i} b_i(\theta) f_i \) can be
!  converted using [[poly_to_dopri_form]].
!
!  If the dense output needs function evaluations in addition to the ones from the step,
!  these are computed in the `*_dense_extra_stages` routine, and stored in `me%funcs`
!  after the ones from the step.

    submodule(rklib_module) rklib_variable_dense

    implicit none

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dense output for the Dormand-Prince 5(4) method [[rkdp54]].
!  This is the 4th order continuous extension due to Shampine (as used in `DOPRI5`).
!  The coefficients are given here in the polynomial form \( b_i(\theta) \)
!  (row `i` contains the coefficients of \( \theta^1 ... \theta^4 \)).
!
!### References
!  * L. F. Shampine, "Some Practical Runge-Kutta Formulas",
!    Mathematics of Computation, 46(173), 1986, p. 135-150.
!  * E. Hairer, S. P. Norsett, G. Wanner, "Solving Ordinary Differential Equations I:
!    Nonstiff Problems", 2nd ed., Springer, 1993, p. 191-192.

    module procedure rkdp54_dense_coefficients

    real(wp),dimension(7,4) :: bp !! coefficients in the polynomial form

    bp(1,:) = [ 1.0_wp, &
               -8048581381.0_wp / 2820520608.0_wp, &
                8663915743.0_wp / 2820520608.0_wp, &
               -12715105075.0_wp / 11282082432.0_wp ]
    bp(2,:) = 0.0_wp
    bp(3,:) = [ 0.0_wp, &
                131558114200.0_wp / 32700410799.0_wp, &
               -68118460800.0_wp / 10900136933.0_wp, &
                87487479700.0_wp / 32700410799.0_wp ]
    bp(4,:) = [ 0.0_wp, &
               -1754552775.0_wp / 470086768.0_wp, &
                14199869525.0_wp / 1410260304.0_wp, &
               -10690763975.0_wp / 1880347072.0_wp ]
    bp(5,:) = [ 0.0_wp, &
                127303824393.0_wp / 49829197408.0_wp, &
               -318862633887.0_wp / 49829197408.0_wp, &
                701980252875.0_wp / 199316789632.0_wp ]
    bp(6,:) = [ 0.0_wp, &
               -282668133.0_wp / 205662961.0_wp, &
                2019193451.0_wp / 616988883.0_wp, &
               -1453857185.0_wp / 822651844.0_wp ]
    bp(7,:) = [ 0.0_wp, &
                40617522.0_wp / 29380423.0_wp, &
               -110615467.0_wp / 29380423.0_wp, &
                69997945.0_wp / 29380423.0_wp ]

    b = poly_to_dopri_form(bp)

    end procedure rkdp54_dense_coefficients
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dense output for the Dormand-Prince 8(5) method [[rkdp85]].
!  This is the 7th order continuous extension from `DOP853`. It requires 4 additional
!  function evaluations (see [[rkdp85_dense_extra_stages]]).
!
!  The coefficients are given here in the form used in `DOP853`:
!
!  \( x(t + \theta h) = r_0 + \theta(r_1 + \theta_1(r_2 + \theta(r_3 + \theta_1(r_4 +
!     \theta(r_5 + \theta_1(r_6 + \theta r_7)))))) \),  where \( \theta_1 = 1 - \theta \)
!
!  where \( r_0 = x(t) \) and \( r_j = h \sum_i b_{ij} f_i \).
!
!### References
!  * E. Hairer, S. P. Norsett, G. Wanner, "Solving Ordinary Differential Equations I:
!    Nonstiff Problems", 2nd ed., Springer, 1993.
!  * [DOP853](https://github.com/jacobwilliams/dop853)

    module procedure rkdp85_dense_coefficients

    ! weights of the 8th order solution (see [[rkdp85]]):
    real(wp),dimension(12),parameter :: w = [ 5.42937341165687622380535766363e-2_wp, &
                                              0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                                              4.45031289275240888144113950566_wp, &
                                              1.89151789931450038304281599044_wp, &
                                             -5.8012039600105847814672114227_wp, &
                                              3.1116436695781989440891606237e-1_wp, &
                                             -1.52160949662516078556178806805e-1_wp, &
                                              2.01365400804030348374776537501e-1_wp, &
                                              4.47106157277725905176885569043e-2_wp ]

    allocate(b(16,7))
    b = 0.0_wp

    ! r1 = x(t+h) - x(t)
    b(1:12,1) = w
    ! r2 = h*f1 - r1
    b(1:12,2) = -w
    b(1,2)    = b(1,2) + 1.0_wp
    ! r3 = r1 - h*f13 - r2
    b(1:12,3) = 2.0_wp*w
    b(1,3)    = b(1,3) - 1.0_wp
    b(13,3)   = -1.0_wp

    ! r4 - r7:
    b([1,6,7,8,9,10,11,12,13,14,15,16],4) = [ &
        -0.84289382761090128651353491142e+01_wp,  0.56671495351937776962531783590e+00_wp, &
        -0.30689499459498916912797304727e+01_wp,  0.23846676565120698287728149680e+01_wp, &
         0.21170345824450282767155149946e+01_wp, -0.87139158377797299206789907490e+00_wp, &
         0.22404374302607882758541771650e+01_wp,  0.63157877876946881815570249290e+00_wp, &
        -0.88990336451333310820698117400e-01_wp,  0.18148505520854727256656404962e+02_wp, &
        -0.91946323924783554000451984436e+01_wp, -0.44360363875948939664310572000e+01_wp ]
    b([1,6,7,8,9,10,11,12,13,14,15,16],5) = [ &
         0.10427508642579134603413151009e+02_wp,  0.24228349177525818288430175319e+03_wp, &
         0.16520045171727028198505394887e+03_wp, -0.37454675472269020279518312152e+03_wp, &
        -0.22113666853125306036270938578e+02_wp,  0.77334326684722638389603898808e+01_wp, &
        -0.30674084731089398182061213626e+02_wp, -0.93321305264302278729567221706e+01_wp, &
         0.15697238121770843886131091075e+02_wp, -0.31139403219565177677282850411e+02_wp, &
        -0.93529243588444783865713862664e+01_wp,  0.35816841486394083752465898540e+02_wp ]
    b([1,6,7,8,9,10,11,12,13,14,15,16],6) = [ &
         0.19985053242002433820987653617e+02_wp, -0.38703730874935176555105901742e+03_wp, &
        -0.18917813819516756882830838328e+03_wp,  0.52780815920542364900561016686e+03_wp, &
        -0.11573902539959630126141871134e+02_wp,  0.68812326946963000169666922661e+01_wp, &
        -0.10006050966910838403183860980e+01_wp,  0.77771377980534432092869265740e+00_wp, &
        -0.27782057523535084065932004339e+01_wp, -0.60196695231264120758267380846e+02_wp, &
         0.84320405506677161018159903784e+02_wp,  0.11992291136182789328035130030e+02_wp ]
    b([1,6,7,8,9,10,11,12,13,14,15,16],7) = [ &
        -0.25693933462703749003312586129e+02_wp, -0.15418974869023643374053993627e+03_wp, &
        -0.23152937917604549567536039109e+03_wp,  0.35763911791061412378285349910e+03_wp, &
         0.93405324183624310003907691704e+02_wp, -0.37458323136451633156875139351e+02_wp, &
         0.10409964950896230045147246184e+03_wp,  0.29840293426660503123344363579e+02_wp, &
        -0.43533456590011143754432175058e+02_wp,  0.96324553959188282948394950600e+02_wp, &
        -0.39177261675615439165231486172e+02_wp, -0.14972683625798562581422125276e+03_wp ]

    end procedure rkdp85_dense_coefficients
!*****************************************************************************************

!*****************************************************************************************
!>
!  The additional function evaluations for the dense output of [[rkdp85]].
!  These are stored in `me%funcs(:,13:16)`.

    module procedure rkdp85_dense_extra_stages

    real(wp),parameter :: c14 = 0.1_wp
    real(wp),parameter :: c15 = 0.2_wp
    real(wp),parameter :: c16 = 7.0_wp / 9.0_wp

    real(wp),parameter :: a141  =  5.61675022830479523392909219681e-2_wp
    real(wp),parameter :: a147  =  2.53500210216624811088794765333e-1_wp
    real(wp),parameter :: a148  = -2.46239037470802489917441475441e-1_wp
    real(wp),parameter :: a149  = -1.24191423263816360469010140626e-1_wp
    real(wp),parameter :: a1410 =  1.5329179827876569731206322685e-1_wp
    real(wp),parameter :: a1411 =  8.20105229563468988491666602057e-3_wp
    real(wp),parameter :: a1412 =  7.56789766054569976138603589584e-3_wp
    real(wp),parameter :: a1413 = -8.298e-3_wp

    real(wp),parameter :: a151  =  3.18346481635021405060768473261e-2_wp
    real(wp),parameter :: a156  =  2.83009096723667755288322961402e-2_wp
    real(wp),parameter :: a157  =  5.35419883074385676223797384372e-2_wp
    real(wp),parameter :: a158  = -5.49237485713909884646569340306e-2_wp
    real(wp),parameter :: a1511 = -1.08347328697249322858509316994e-4_wp
    real(wp),parameter :: a1512 =  3.82571090835658412954920192323e-4_wp
    real(wp),parameter :: a1513 = -3.40465008687404560802977114492e-4_wp
    real(wp),parameter :: a1514 =  1.41312443674632500278074618366e-1_wp

    real(wp),parameter :: a161  = -4.28896301583791923408573538692e-1_wp
    real(wp),parameter :: a166  = -4.69762141536116384314449447206e0_wp
    real(wp),parameter :: a167  =  7.68342119606259904184240953878e0_wp
    real(wp),parameter :: a168  =  4.06898981839711007970213554331e0_wp
    real(wp),parameter :: a169  =  3.56727187455281109270669543021e-1_wp
    real(wp),parameter :: a1613 = -1.39902416515901462129418009734e-3_wp
    real(wp),parameter :: a1614 =  2.9475147891527723389556272149e0_wp
    real(wp),parameter :: a1615 = -9.15095847217987001081870187138e0_wp

    associate (t   => me%dense_t,      &
               h   => me%dense_h,      &
               x   => me%dense_x,      &
               f1  => me%funcs(:,1),   &
               f6  => me%funcs(:,6),   &
               f7  => me%funcs(:,7),   &
               f8  => me%funcs(:,8),   &
               f9  => me%funcs(:,9),   &
               f10 => me%funcs(:,10),  &
               f11 => me%funcs(:,11),  &
               f12 => me%funcs(:,12),  &
               f13 => me%funcs(:,13),  &
               f14 => me%funcs(:,14),  &
               f15 => me%funcs(:,15),  &
               f16 => me%funcs(:,16))

        call me%f(t+h,     me%dense_xf, f13)
        if (me%stopped) return
        call me%f(t+c14*h, x+h*(a141*f1+a147*f7+a148*f8+a149*f9+a1410*f10+&
                                a1411*f11+a1412*f12+a1413*f13), f14)
        if (me%stopped) return
        call me%f(t+c15*h, x+h*(a151*f1+a156*f6+a157*f7+a158*f8+a1511*f11+&
                                a1512*f12+a1513*f13+a1514*f14), f15)
        if (me%stopped) return
        call me%f(t+c16*h, x+h*(a161*f1+a166*f6+a167*f7+a168*f8+a169*f9+&
                                a1613*f13+a1614*f14+a1615*f15), f16)

    end associate

    end procedure rkdp85_dense_extra_stages
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert dense output coefficients from the polynomial form:
!
!  \( x(t + \theta h) = x(t) + h \sum_i \sum_{k=1}^{degree} b_{ik} \theta^k f_i \)
!
!  to the form used in `DOPRI5` and `DOP853` (see [[dense_coefficients]]):
!
!  \( x(t + \theta h) = x(t) + \theta(r_1 + \theta_1(r_2 + \theta(r_3 + \theta_1(r_4 + ...)))) \),
!  where \( \theta_1 = 1 - \theta \) and \( r_j = h \sum_i r_{ij} f_i \).
!
!  The \( r_j \) term is multiplied by \( P_j(\theta) =
!  \theta^{\lceil j/2 \rceil} \theta_1^{\lfloor j/2 \rfloor} \), which has degree `j`
!  and leading coefficient \( (-1)^{\lfloor j/2 \rfloor} \), so the conversion is
!  done by eliminating the highest power of \( \theta \) first.

    pure function poly_to_dopri_form(b) result(r)

    real(wp),dimension(:,:),intent(in) :: b !! coefficients in the polynomial form (`stages` x `degree`)
    real(wp),dimension(size(b,1),size(b,2)) :: r !! coefficients in the `DOPRI` form (`stages` x `degree`)

    real(wp),dimension(size(b,1),size(b,2)) :: bt !! remaining polynomial coefficients
    real(wp),dimension(size(b,2)) :: pj !! coefficients of \( \theta^1 ... \theta^{degree} \) in \( P_j \)
    integer :: j !! counter for the `r` terms
    integer :: k !! counter for the binomial expansion of \( (1-\theta)^{\lfloor j/2 \rfloor} \)
    integer :: p !! power of \( \theta \) in \( P_j \)
    integer :: q !! power of \( 1-\theta \) in \( P_j \)
    real(wp) :: binom !! signed binomial coefficient

    bt = b
    do j = size(b,2), 1, -1
        ! expand P_j:
        p = (j+1)/2
        q = j/2
        pj = 0.0_wp
        binom = 1.0_wp
        do k = 0, q
            pj(p+k) = binom
            binom = -binom * real(q-k,wp) / real(k+1,wp)
        end do
        ! eliminate the theta^j term:
        r(:,j) = bt(:,j) / pj(j)
        do k = 1, j
            bt(:,k) = bt(:,k) - r(:,j) * pj(k)
        end do
    end do

    end function poly_to_dopri_form
!*****************************************************************************************

!*****************************************************************************************
    end submodule rklib_variable_dense
!*****************************************************************************************
