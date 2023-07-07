!*****************************************************************************************
!>
!  Procedures specifying the order of the various variable-step RK methods.

    submodule(rklib_module) rklib_variable_orders

    implicit none

    contains
!*****************************************************************************************

    module procedure rkbs32_order
        !! Returns the order of the [[rkbs32]] method
        p = 3
    end procedure rkbs32_order

    module procedure rkf45_order
        !! Returns the order of the [[rkf45]] method
        p = 4
    end procedure rkf45_order

    module procedure rkdp54_order
        !! Returns the order of the [[rkdp54]] method
        p = 5
    end procedure rkdp54_order

    module procedure rkt54_order
        !! Returns the order of the [[rkt54]] method
        p = 5
    end procedure rkt54_order

    module procedure rks54_order
        !! Returns the order of the [[rks54]] method
        p = 5
    end procedure rks54_order

    module procedure rkck54_order
        !! Returns the order of the [[rkck54]] method
        p = 5
    end procedure rkck54_order

    module procedure rkdp65_order
        !! Returns the order of the [[rkdp65]] method
        p = 6
    end procedure rkdp65_order

    module procedure rktp64_order
        !! Returns the order of the [[rktp64]] method
        p = 6
    end procedure rktp64_order

    module procedure rkc65_order
        !! Returns the order of the [[rkc65]] method
        p = 6
    end procedure rkc65_order

    module procedure rkv65e_order
        !! Returns the order of the [[rkv65e]] method
        p = 6
    end procedure rkv65e_order

    module procedure rkv65r_order
        !! Returns the order of the [[rkv65r]] method
        p = 6
    end procedure rkv65r_order

    module procedure rktf65_order
        !! Returns the order of the [[rktf65]] method
        p = 6
    end procedure rktf65_order

    module procedure rktp75_order
        !! Returns the order of the [[rktp75]] method
        p = 7
    end procedure rktp75_order

    module procedure rktmy7_order
        !! Returns the order of the [[rktmy7]] method
        p = 7
    end procedure rktmy7_order

    module procedure rkv76e_order
        !! Returns the order of the [[rkv76e]] method
        p = 7
    end procedure rkv76e_order

    module procedure rkv76r_order
        !! Returns the order of the [[rkv76r]] method
        p = 7
    end procedure rkv76r_order

    module procedure rkf78_order
        !! Returns the order of the [[rkf78]] method
        p = 7
    end procedure rkf78_order

    module procedure rkv78_order
        !! Returns the order of the [[rkv78]] method
        p = 7
    end procedure rkv78_order

    module procedure rktp86_order
        !! Returns the order of the [[rktp86]] method
        p = 8
    end procedure rktp86_order

    module procedure rkdp87_order
        !! Returns the order of the [[rkdp87]] method
        p = 8
    end procedure rkdp87_order

    module procedure rkv87e_order
        !! Returns the order of the [[rkv87e]] method
        p = 8
    end procedure rkv87e_order

    module procedure rkv87r_order
        !! Returns the order of the [[rkv87r]] method
        p = 8
    end procedure rkv87r_order

    module procedure rkk87_order
        !! Returns the order of the [[rkk87]] method
        p = 8
    end procedure rkk87_order

    module procedure rkf89_order
        !! Returns the order of the [[rkf89]] method
        p = 8
    end procedure rkf89_order

    module procedure rkv89_order
        !! Returns the order of the [[rkv89]] method
        p = 8
    end procedure rkv89_order

    module procedure rkt98a_order
        !! Returns the order of the [[rkt98a]] method
        p = 9
    end procedure rkt98a_order

    module procedure rkv98e_order
        !! Returns the order of the [[rkv98e]] method
        p = 9
    end procedure rkv98e_order

    module procedure rkv98r_order
        !! Returns the order of the [[rkv98r]] method
        p = 9
    end procedure rkv98r_order

    module procedure rkf108_order
        !! Returns the order of the [[rkf108]] method.
        p = 10
    end procedure rkf108_order

    module procedure rkc108_order
        !! Returns the order of the [[rkc108]] method.
        p = 10
    end procedure rkc108_order

    module procedure rks1110a_order
        !! Returns the order of the [[rks1110a]] method.
        p = 11
    end procedure rks1110a_order

    module procedure rkf1210_order
        !! Returns the order of the [[rkf1210]] method.
        p = 12
    end procedure rkf1210_order

    module procedure rko129_order
        !! Returns the order of the [[rko129]] method.
        p = 12
    end procedure rko129_order

    module procedure rkf1412_order
        !! Returns the order of the [[rkf1412]] method.
        p = 14
    end procedure rkf1412_order

!*****************************************************************************************
    end submodule rklib_variable_orders
!*****************************************************************************************