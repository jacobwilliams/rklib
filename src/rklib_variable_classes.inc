    ! Variable step methods:
    type,extends(rk_variable_step_fsal_class),public :: rkbs32_class
        !! Cash & Karp 5(4) order method.
        contains
        procedure :: step  => rkbs32
        procedure :: order => rkbs32_order
    end type rkbs32_class
    type,extends(rk_variable_step_class),public :: rkf45_class
        !! Fehlberg's 4(5) method
        contains
        procedure :: step  => rkf45
        procedure :: order => rkf45_order
    end type rkf45_class
    type,extends(rk_variable_step_class),public :: rkck54_class
        !! Cash & Karp 5(4) order method.
        contains
        procedure :: step  => rkck54
        procedure :: order => rkck54_order
    end type rkck54_class
    type,extends(rk_variable_step_fsal_class),public :: rkdp54_class
        !! Dormand Prince 5(4) order method.
        contains
        procedure :: step  => rkdp54
        procedure :: order => rkdp54_order
    end type rkdp54_class
    type,extends(rk_variable_step_fsal_class),public :: rkt54_class
        !! Tsitouras 5(4)
        contains
        procedure :: step  => rkt54
        procedure :: order => rkt54_order
    end type rkt54_class
    type,extends(rk_variable_step_class),public :: rkdp65_class
        !! Dormand-Prince 6(5) method.
        contains
        procedure :: step  => rkdp65
        procedure :: order => rkdp65_order
    end type rkdp65_class
    type,extends(rk_variable_step_class),public :: rkc65_class
        !! Calvo 6(5) method.
        contains
        procedure :: step  => rkc65
        procedure :: order => rkc65_order
    end type rkc65_class
    type,extends(rk_variable_step_class),public :: rktp64_class
        !! Tsitouras & Papakostas NEW6(4).
        contains
        procedure :: step  => rktp64
        procedure :: order => rktp64_order
    end type rktp64_class
    type,extends(rk_variable_step_fsal_class),public :: rkv65e_class
        !! Verner's 'most efficient' Runge-Kutta (9,6(5))
        contains
        procedure :: step  => rkv65e
        procedure :: order => rkv65e_order
    end type rkv65e_class
    type,extends(rk_variable_step_fsal_class),public :: rktf65_class
        !! Tsitouras & Famelis 6(5)
        contains
        procedure :: step  => rktf65
        procedure :: order => rktf65_order
    end type rktf65_class
    type,extends(rk_variable_step_fsal_class),public :: rkv65r_class
        !! Verner's 'most robust' Runge-Kutta (9,6(5))
        contains
        procedure :: step  => rkv65r
        procedure :: order => rkv65r_order
    end type rkv65r_class
    type,extends(rk_variable_step_class),public :: rkv76e_class
        !! Verner's 'most efficient' Runge-Kutta (10:7(6))
        contains
        procedure :: step  => rkv76e
        procedure :: order => rkv76e_order
    end type rkv76e_class
    type,extends(rk_variable_step_class),public :: rkv76r_class
        !! Verner's 'most robust' Runge-Kutta (10:7(6))
        contains
        procedure :: step  => rkv76r
        procedure :: order => rkv76r_order
    end type rkv76r_class
    type,extends(rk_variable_step_class),public :: rkf78_class
        !! Runga-Kutta Fehlberg 7(8) method.
        contains
        procedure :: step  => rkf78
        procedure :: order => rkf78_order
    end type rkf78_class
    type,extends(rk_variable_step_class),public :: rkdp87_class
        !! Dormand & Prince RK8(7)13M method.
        contains
        procedure :: step  => rkdp87
        procedure :: order => rkdp87_order
    end type rkdp87_class
    type,extends(rk_variable_step_class),public :: rkv87e_class
        !! Verner's "most efficient" Runge-Kutta (8)7 method.
        contains
        procedure :: step  => rkv87e
        procedure :: order => rkv87e_order
    end type rkv87e_class
    type,extends(rk_variable_step_class),public :: rkv87r_class
        !! Verner's "most robust" Runge-Kutta (8)7 method.
        contains
        procedure :: step  => rkv87r
        procedure :: order => rkv87r_order
    end type rkv87r_class
    type,extends(rk_variable_step_class),public :: rkk87_class
        !! Kovalnogov, Fedorov, Karpukhina, Simos, Tsitouras 8(7) method.
        contains
        procedure :: step  => rkk87
        procedure :: order => rkk87_order
    end type rkk87_class
    type,extends(rk_variable_step_class),public :: rkv78_class
        !! Runga-Kutta Verner 7(8) method.
        contains
        procedure :: step  => rkv78
        procedure :: order => rkv78_order
    end type rkv78_class
    type,extends(rk_variable_step_class),public :: rktp75_class
        !! Tsitouras & Papakostas NEW7(5).
        contains
        procedure :: step  => rktp75
        procedure :: order => rktp75_order
    end type rktp75_class
    type,extends(rk_variable_step_class),public :: rktmy7_class
        !! Seventh-order Runge-Kutta Method, by M. Tanaka, S. Muramatsu and S. Yamashita
        contains
        procedure :: step  => rktmy7
        procedure :: order => rktmy7_order
    end type rktmy7_class
    type,extends(rk_variable_step_class),public :: rktp86_class
        !! Tsitouras & Papakostas NEW8(6).
        contains
        procedure :: step  => rktp86
        procedure :: order => rktp86_order
    end type rktp86_class
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
    type,extends(rk_variable_step_class),public :: rkt98a_class
        !! Tsitouras combined order 9(8) Runge-Kutta scheme A.
        contains
        procedure :: step  => rkt98a
        procedure :: order => rkt98a_order
    end type rkt98a_class
    type,extends(rk_variable_step_class),public :: rkv98e_class
        !! Verner's "most efficient" Runge-Kutta (16:9(8)) method.
        contains
        procedure :: step  => rkv98e
        procedure :: order => rkv98e_order
    end type rkv98e_class
    type,extends(rk_variable_step_class),public :: rkv98r_class
        !! Verner's "most robust" Runge-Kutta (16:9(8)) method.
        contains
        procedure :: step  => rkv98r
        procedure :: order => rkv98r_order
    end type rkv98r_class
    type,extends(rk_variable_step_class),public :: rkf108_class
        !! Runga-Kutta Feagin 8(10) method.
        contains
        procedure :: step  => rkf108
        procedure :: order => rkf108_order
    end type rkf108_class
    type,extends(rk_variable_step_class),public :: rkc108_class
        !! Curtis 10(8) method.
        contains
        procedure :: step  => rkc108
        procedure :: order => rkc108_order
    end type rkc108_class
    type,extends(rk_variable_step_class),public :: rks1110a_class
        !! Stone 11(10) method.
        contains
        procedure :: step  => rks1110a
        procedure :: order => rks1110a_order
    end type rks1110a_class
    type,extends(rk_variable_step_class),public :: rkf1210_class
        !! Runga-Kutta Feagin 12(10) method.
        contains
        procedure :: step  => rkf1210
        procedure :: order => rkf1210_order
    end type rkf1210_class
    type,extends(rk_variable_step_class),public :: rko129_class
        !! Ono 12(9) method
        contains
        procedure :: step  => rko129
        procedure :: order => rko129_order
    end type rko129_class
    type,extends(rk_variable_step_class),public :: rkf1412_class
        !! Runga-Kutta Feagin 14(12) method.
        contains
        procedure :: step  => rkf1412
        procedure :: order => rkf1412_order
    end type rkf1412_class
