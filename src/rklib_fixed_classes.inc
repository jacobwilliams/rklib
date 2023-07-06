    ! Fixed step methods:
    type,extends(rk_fixed_step_class),public :: euler_class
        !! Basic Euler method
        contains
        procedure :: step => euler
    end type euler_class
    type,extends(rk_fixed_step_class),public :: midpoint_class
        !! Basic Midpoint method
    contains
        procedure :: step => midpoint
    end type midpoint_class
    type,extends(rk_fixed_step_class),public :: heun_class
        !! Basic Heun method
        contains
        procedure :: step => heun
    end type heun_class
    type,extends(rk_fixed_step_class),public :: rk3_class
        !! 3th order Runge-Kutta method.
        contains
        procedure :: step => rk3
    end type rk3_class
    type,extends(rk_fixed_step_class),public :: rk4_class
        !! 4th order Runge-Kutta method.
        contains
        procedure :: step => rk4
    end type rk4_class
    type,extends(rk_fixed_step_class),public :: rks4_class
        !! 4th order Runge-Kutta method.
        contains
        procedure :: step => rks4
    end type rks4_class
    type,extends(rk_fixed_step_class),public :: rks5_class
        !! 5th order Runge-Kutta Shanks method.
        contains
        procedure :: step => rks5
    end type rks5_class
    type,extends(rk_fixed_step_class),public :: rkb6_class
        !! 6th order Runge-Kutta Butcher method.
        contains
        procedure :: step => rkb6
    end type rkb6_class
    type,extends(rk_fixed_step_class),public :: rk7_class
        !! 7th order Runge-Kutta method.
        contains
        procedure :: step => rk7
    end type rk7_class
    type,extends(rk_fixed_step_class),public :: rk8_10_class
        !! 8th order Runge-Kutta method.
        contains
        procedure :: step => rk8_10
    end type rk8_10_class
    type,extends(rk_fixed_step_class),public :: rk8_12_class
        !! 8th order Runge-Kutta method.
        contains
        procedure :: step => rk8_12
    end type rk8_12_class
    type,extends(rk_fixed_step_class),public :: rkcv8_class
        !! Cooper-Verner 11 stage, 8th order Runge-Kutta method.
        contains
        procedure :: step => rkcv8
    end type rkcv8_class
    type,extends(rk_fixed_step_class),public :: rkz10_class
        !! Zhang's 10th order method
        contains
        procedure :: step => rkz10
    end type rkz10_class
