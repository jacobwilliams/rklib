    ! Fixed step methods:

    type,extends(rk_fixed_step_class),public :: euler_class
        !! Euler
        contains
        procedure :: step => euler
        procedure :: properties => euler_properties
    end type euler_class

    type,extends(rk_fixed_step_class),public :: midpoint_class
        !! Midpoint
        contains
        procedure :: step => midpoint
        procedure :: properties => midpoint_properties
    end type midpoint_class

    type,extends(rk_fixed_step_class),public :: heun_class
        !! Heun
        contains
        procedure :: step => heun
        procedure :: properties => heun_properties
    end type heun_class

    type,extends(rk_fixed_step_class),public :: rkssp22_class
        !! 2-stage, 2nd order TVD Runge-Kutta Shu-Osher
        contains
        procedure :: step => rkssp22
        procedure :: properties => rkssp22_properties
    end type rkssp22_class

    type,extends(rk_fixed_step_class),public :: rk3_class
        !! 3th order Runge-Kutta
        contains
        procedure :: step => rk3
        procedure :: properties => rk3_properties
    end type rk3_class

    type,extends(rk_fixed_step_class),public :: rkssp33_class
        !! 3-stage, 3rd order TVD Runge-Kutta Shu-Osher
        contains
        procedure :: step => rkssp33
        procedure :: properties => rkssp33_properties
    end type rkssp33_class

    type,extends(rk_fixed_step_class),public :: rkssp53_class
        !! 5-stage, 3rd order SSP Runge-Kutta Spiteri-Ruuth
        contains
        procedure :: step => rkssp53
        procedure :: properties => rkssp53_properties
    end type rkssp53_class

    type,extends(rk_fixed_step_class),public :: rk4_class
        !! Classic 4th order Runge-Kutta
        contains
        procedure :: step => rk4
        procedure :: properties => rk4_properties
    end type rk4_class

    type,extends(rk_fixed_step_class),public :: rks4_class
        !! 4th order Runge-Kutta Shanks
        contains
        procedure :: step => rks4
        procedure :: properties => rks4_properties
    end type rks4_class

    type,extends(rk_fixed_step_class),public :: rkr4_class
        !! 4th order Runge-Kutta Ralston
        contains
        procedure :: step => rkr4
        procedure :: properties => rkr4_properties
    end type rkr4_class

    type,extends(rk_fixed_step_class),public :: rkls44_class
        !! 4-stage, 4th order low storage non-TVD Runge-Kutta Jiang-Shu
        contains
        procedure :: step => rkls44
        procedure :: properties => rkls44_properties
    end type rkls44_class

    type,extends(rk_fixed_step_class),public :: rkls54_class
        !! 5-stage, 4th order low storage Runge-Kutta Carpenter-Kennedy
        contains
        procedure :: step => rkls54
        procedure :: properties => rkls54_properties
    end type rkls54_class

    type,extends(rk_fixed_step_class),public :: rkssp54_class
        !! 5-stage, 4th order SSP Runge-Kutta Spiteri-Ruuth
        contains
        procedure :: step => rkssp54
        procedure :: properties => rkssp54_properties
    end type rkssp54_class

    type,extends(rk_fixed_step_class),public :: rks5_class
        !! 5th order Runge-Kutta Shanks
        contains
        procedure :: step => rks5
        procedure :: properties => rks5_properties
    end type rks5_class

    type,extends(rk_fixed_step_class),public :: rk5_class
        !! 5th order Runge-Kutta
        contains
        procedure :: step => rk5
        procedure :: properties => rk5_properties
    end type rk5_class

    type,extends(rk_fixed_step_class),public :: rkc5_class
        !! 5th order Runge-Kutta Cassity
        contains
        procedure :: step => rkc5
        procedure :: properties => rkc5_properties
    end type rkc5_class

    type,extends(rk_fixed_step_class),public :: rkl5_class
        !! 5th order Runge-Kutta Lawson
        contains
        procedure :: step => rkl5
        procedure :: properties => rkl5_properties
    end type rkl5_class

    type,extends(rk_fixed_step_class),public :: rklk5a_class
        !! 5th order Runge-Kutta Luther-Konen 1
        contains
        procedure :: step => rklk5a
        procedure :: properties => rklk5a_properties
    end type rklk5a_class

    type,extends(rk_fixed_step_class),public :: rklk5b_class
        !! 5th order Runge-Kutta Luther-Konen 2
        contains
        procedure :: step => rklk5b
        procedure :: properties => rklk5b_properties
    end type rklk5b_class

    type,extends(rk_fixed_step_class),public :: rkb6_class
        !! 6th order Runge-Kutta Butcher
        contains
        procedure :: step => rkb6
        procedure :: properties => rkb6_properties
    end type rkb6_class

    type,extends(rk_fixed_step_class),public :: rk7_class
        !! 7th order Runge-Kutta Shanks
        contains
        procedure :: step => rk7
        procedure :: properties => rk7_properties
    end type rk7_class

    type,extends(rk_fixed_step_class),public :: rk8_10_class
        !! 10-stage, 8th order Runge-Kutta Shanks
        contains
        procedure :: step => rk8_10
        procedure :: properties => rk8_10_properties
    end type rk8_10_class

    type,extends(rk_fixed_step_class),public :: rkcv8_class
        !! 11-stage, 8th order Runge-Kutta Cooper-Verner
        contains
        procedure :: step => rkcv8
        procedure :: properties => rkcv8_properties
    end type rkcv8_class

    type,extends(rk_fixed_step_class),public :: rk8_12_class
        !! 12-stage, 8th order Runge-Kutta Shanks
        contains
        procedure :: step => rk8_12
        procedure :: properties => rk8_12_properties
    end type rk8_12_class

    type,extends(rk_fixed_step_class),public :: rks10_class
        !! 10th order Runge-Kutta Stepanov
        contains
        procedure :: step => rks10
        procedure :: properties => rks10_properties
    end type rks10_class

    type,extends(rk_fixed_step_class),public :: rkz10_class
        !! 10th order Runge-Kutta Zhang
        contains
        procedure :: step => rkz10
        procedure :: properties => rkz10_properties
    end type rkz10_class

    type,extends(rk_fixed_step_class),public :: rko10_class
        !! 10th order Runge-Kutta Ono
        contains
        procedure :: step => rko10
        procedure :: properties => rko10_properties
    end type rko10_class

    type,extends(rk_fixed_step_class),public :: rkh10_class
        !! 10th order Runge-Kutta Hairer
        contains
        procedure :: step => rkh10
        procedure :: properties => rkh10_properties
    end type rkh10_class

