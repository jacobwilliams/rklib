    ! Variable step methods:

    type,extends(rk_variable_step_fsal_class),public :: rkbs32_class
        !! Bogacki & Shampine 3(2)
        contains
        procedure :: step => rkbs32
        procedure :: properties => rkbs32_properties
    end type rkbs32_class

    type,extends(rk_variable_step_class),public :: rkssp43_class
        !! 4-stage, 3rd order SSP
        contains
        procedure :: step => rkssp43
        procedure :: properties => rkssp43_properties
    end type rkssp43_class

    type,extends(rk_variable_step_class),public :: rkf45_class
        !! Fehlberg 4(5)
        contains
        procedure :: step => rkf45
        procedure :: properties => rkf45_properties
    end type rkf45_class

    type,extends(rk_variable_step_class),public :: rkck54_class
        !! Cash & Karp 5(4)
        contains
        procedure :: step => rkck54
        procedure :: properties => rkck54_properties
    end type rkck54_class

    type,extends(rk_variable_step_fsal_class),public :: rkdp54_class
        !! Dormand-Prince 5(4)
        contains
        procedure :: step => rkdp54
        procedure :: properties => rkdp54_properties
    end type rkdp54_class

    type,extends(rk_variable_step_fsal_class),public :: rkt54_class
        !! Tsitouras 5(4)
        contains
        procedure :: step => rkt54
        procedure :: properties => rkt54_properties
    end type rkt54_class

    type,extends(rk_variable_step_fsal_class),public :: rks54_class
        !! Stepanov 5(4)
        contains
        procedure :: step => rks54
        procedure :: properties => rks54_properties
    end type rks54_class

    type,extends(rk_variable_step_fsal_class),public :: rkpp54_class
        !! Papakostas-PapaGeorgiou 5(4)
        contains
        procedure :: step => rkpp54
        procedure :: properties => rkpp54_properties
    end type rkpp54_class

    type,extends(rk_variable_step_fsal_class),public :: rkpp54b_class
        !! Papakostas-PapaGeorgiou 5(4) b
        contains
        procedure :: step => rkpp54b
        procedure :: properties => rkpp54b_properties
    end type rkpp54b_class

    type,extends(rk_variable_step_class),public :: rkbs54_class
        !! Bogacki & Shampine 5(4)
        contains
        procedure :: step => rkbs54
        procedure :: properties => rkbs54_properties
    end type rkbs54_class

    type,extends(rk_variable_step_class),public :: rkss54_class
        !! Sharp & Smart 5(4)
        contains
        procedure :: step => rkss54
        procedure :: properties => rkss54_properties
    end type rkss54_class

    type,extends(rk_variable_step_class),public :: rkdp65_class
        !! Dormand-Prince 6(5)
        contains
        procedure :: step => rkdp65
        procedure :: properties => rkdp65_properties
    end type rkdp65_class

    type,extends(rk_variable_step_class),public :: rkc65_class
        !! Calvo 6(5)
        contains
        procedure :: step => rkc65
        procedure :: properties => rkc65_properties
    end type rkc65_class

    type,extends(rk_variable_step_class),public :: rktp64_class
        !! Tsitouras & Papakostas NEW6(4)
        contains
        procedure :: step => rktp64
        procedure :: properties => rktp64_properties
    end type rktp64_class

    type,extends(rk_variable_step_fsal_class),public :: rkv65e_class
        !! Verner efficient (9,6(5))
        contains
        procedure :: step => rkv65e
        procedure :: properties => rkv65e_properties
    end type rkv65e_class

    type,extends(rk_variable_step_fsal_class),public :: rkv65r_class
        !! Verner robust (9,6(5))
        contains
        procedure :: step => rkv65r
        procedure :: properties => rkv65r_properties
    end type rkv65r_class

    type,extends(rk_variable_step_class),public :: rkv65_class
        !! Verner 6(5)
        contains
        procedure :: step => rkv65
        procedure :: properties => rkv65_properties
    end type rkv65_class

    type,extends(rk_variable_step_class),public :: dverk65_class
        !! Verner 6(5) "DVERK"
        contains
        procedure :: step => dverk65
        procedure :: properties => dverk65_properties
    end type dverk65_class

    type,extends(rk_variable_step_fsal_class),public :: rktf65_class
        !! Tsitouras & Famelis 6(5)
        contains
        procedure :: step => rktf65
        procedure :: properties => rktf65_properties
    end type rktf65_class

    type,extends(rk_variable_step_class),public :: rktp75_class
        !! Tsitouras & Papakostas NEW7(5)
        contains
        procedure :: step => rktp75
        procedure :: properties => rktp75_properties
    end type rktp75_class

    type,extends(rk_variable_step_class),public :: rktmy7_class
        !! 7th order Tanaka-Muramatsu-Yamashita
        contains
        procedure :: step => rktmy7
        procedure :: properties => rktmy7_properties
    end type rktmy7_class

    type,extends(rk_variable_step_class),public :: rktmy7s_class
        !! 7th order Stable Tanaka-Muramatsu-Yamashita
        contains
        procedure :: step => rktmy7s
        procedure :: properties => rktmy7s_properties
    end type rktmy7s_class

    type,extends(rk_variable_step_class),public :: rkv76e_class
        !! Verner efficient (10:7(6))
        contains
        procedure :: step => rkv76e
        procedure :: properties => rkv76e_properties
    end type rkv76e_class

    type,extends(rk_variable_step_class),public :: rkv76r_class
        !! Verner robust (10:7(6))
        contains
        procedure :: step => rkv76r
        procedure :: properties => rkv76r_properties
    end type rkv76r_class

    type,extends(rk_variable_step_class),public :: rkss76_class
        !! Sharp & Smart 7(6)
        contains
        procedure :: step => rkss76
        procedure :: properties => rkss76_properties
    end type rkss76_class

    type,extends(rk_variable_step_class),public :: rkf78_class
        !! Fehlberg 7(8)
        contains
        procedure :: step => rkf78
        procedure :: properties => rkf78_properties
    end type rkf78_class

    type,extends(rk_variable_step_class),public :: rkv78_class
        !! Verner 7(8)
        contains
        procedure :: step => rkv78
        procedure :: properties => rkv78_properties
    end type rkv78_class

    type,extends(rk_variable_step_class),public :: dverk78_class
        !! Verner "Maple" 7(8)
        contains
        procedure :: step => dverk78
        procedure :: properties => dverk78_properties
    end type dverk78_class

    type,extends(rk_variable_step_class),public :: rkdp85_class
        !! Dormand-Prince 8(5)
        contains
        procedure :: step => rkdp85
        procedure :: properties => rkdp85_properties
    end type rkdp85_class

    type,extends(rk_variable_step_class),public :: rktp86_class
        !! Tsitouras & Papakostas NEW8(6)
        contains
        procedure :: step => rktp86
        procedure :: properties => rktp86_properties
    end type rktp86_class

    type,extends(rk_variable_step_class),public :: rkdp87_class
        !! Dormand & Prince RK8(7)13M
        contains
        procedure :: step => rkdp87
        procedure :: properties => rkdp87_properties
    end type rkdp87_class

    type,extends(rk_variable_step_class),public :: rkv87e_class
        !! Verner efficient (8)7
        contains
        procedure :: step => rkv87e
        procedure :: properties => rkv87e_properties
    end type rkv87e_class

    type,extends(rk_variable_step_class),public :: rkv87r_class
        !! Verner robust (8)7
        contains
        procedure :: step => rkv87r
        procedure :: properties => rkv87r_properties
    end type rkv87r_class

    type,extends(rk_variable_step_class),public :: rkev87_class
        !! Enright-Verner (8)7
        contains
        procedure :: step => rkev87
        procedure :: properties => rkev87_properties
    end type rkev87_class

    type,extends(rk_variable_step_class),public :: rkk87_class
        !! Kovalnogov-Fedorov-Karpukhina-Simos-Tsitouras 8(7)
        contains
        procedure :: step => rkk87
        procedure :: properties => rkk87_properties
    end type rkk87_class

    type,extends(rk_variable_step_class),public :: rkf89_class
        !! Fehlberg 8(9)
        contains
        procedure :: step => rkf89
        procedure :: properties => rkf89_properties
    end type rkf89_class

    type,extends(rk_variable_step_class),public :: rkv89_class
        !! Verner 8(9)
        contains
        procedure :: step => rkv89
        procedure :: properties => rkv89_properties
    end type rkv89_class

    type,extends(rk_variable_step_class),public :: rkt98a_class
        !! Tsitouras 9(8) A
        contains
        procedure :: step => rkt98a
        procedure :: properties => rkt98a_properties
    end type rkt98a_class

    type,extends(rk_variable_step_class),public :: rkv98e_class
        !! Verner efficient (16:9(8))
        contains
        procedure :: step => rkv98e
        procedure :: properties => rkv98e_properties
    end type rkv98e_class

    type,extends(rk_variable_step_class),public :: rkv98r_class
        !! Verner robust (16:9(8))
        contains
        procedure :: step => rkv98r
        procedure :: properties => rkv98r_properties
    end type rkv98r_class

    type,extends(rk_variable_step_class),public :: rks98_class
        !! Sharp 9(8)
        contains
        procedure :: step => rks98
        procedure :: properties => rks98_properties
    end type rks98_class

    type,extends(rk_variable_step_class),public :: rkf108_class
        !! Feagin 8(10)
        contains
        procedure :: step => rkf108
        procedure :: properties => rkf108_properties
    end type rkf108_class

    type,extends(rk_variable_step_class),public :: rkc108_class
        !! Curtis 10(8)
        contains
        procedure :: step => rkc108
        procedure :: properties => rkc108_properties
    end type rkc108_class

    type,extends(rk_variable_step_class),public :: rkb109_class
        !! Baker 10(9)
        contains
        procedure :: step => rkb109
        procedure :: properties => rkb109_properties
    end type rkb109_class

    type,extends(rk_variable_step_class),public :: rks1110a_class
        !! Stone 11(10)
        contains
        procedure :: step => rks1110a
        procedure :: properties => rks1110a_properties
    end type rks1110a_class

    type,extends(rk_variable_step_class),public :: rkf1210_class
        !! Feagin 12(10)
        contains
        procedure :: step => rkf1210
        procedure :: properties => rkf1210_properties
    end type rkf1210_class

    type,extends(rk_variable_step_class),public :: rko129_class
        !! Ono 12(9)
        contains
        procedure :: step => rko129
        procedure :: properties => rko129_properties
    end type rko129_class

    type,extends(rk_variable_step_class),public :: rkf1412_class
        !! Feagin 14(12)
        contains
        procedure :: step => rkf1412
        procedure :: properties => rkf1412_properties
    end type rkf1412_class

