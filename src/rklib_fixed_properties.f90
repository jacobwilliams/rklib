
!*****************************************************************************************
!>
!  Procedures specifying the properties of the various fixed-step RK methods.

    submodule(rklib_module) rklib_fixed_properties

    implicit none

    contains
!*****************************************************************************************

    module procedure euler_properties
        !! Returns the properties of the [[euler]] method
        p%short_name = 'euler'
        p%long_name = 'Euler'
        p%order = 1
        p%number_of_stages = 1
        p%number_of_registers = 1
        p%cfl = 1.0_wp
    end procedure euler_properties

    module procedure midpoint_properties
        !! Returns the properties of the [[midpoint]] method
        p%short_name = 'midpoint'
        p%long_name = 'Midpoint'
        p%order = 2
        p%number_of_stages = 2
        p%number_of_registers = 2
    end procedure midpoint_properties

    module procedure heun_properties
        !! Returns the properties of the [[heun]] method
        p%short_name = 'heun'
        p%long_name = 'Heun'
        p%order = 2
        p%number_of_stages = 2
        p%number_of_registers = 2
    end procedure heun_properties

    module procedure rkssp22_properties
        !! Returns the properties of the [[rkssp22]] method
        p%short_name = 'rkssp22'
        p%long_name = '2-stage, 2nd order TVD Runge-Kutta Shu-Osher'
        p%order = 2
        p%number_of_stages = 2
        p%number_of_registers = 1
        p%strong_stability_preserving = .true.
        p%cfl = 1.0_wp
    end procedure rkssp22_properties

    module procedure rk3_properties
        !! Returns the properties of the [[rk3]] method
        p%short_name = 'rk3'
        p%long_name = '3th order Runge-Kutta'
        p%order = 3
        p%number_of_stages = 3
        p%number_of_registers = 3
    end procedure rk3_properties

    module procedure rkssp33_properties
        !! Returns the properties of the [[rkssp33]] method
        p%short_name = 'rkssp33'
        p%long_name = '3-stage, 3rd order TVD Runge-Kutta Shu-Osher'
        p%order = 3
        p%number_of_stages = 3
        p%number_of_registers = 1
        p%strong_stability_preserving = .true.
        p%cfl = 1.0_wp
    end procedure rkssp33_properties

    module procedure rkssp53_properties
        !! Returns the properties of the [[rkssp53]] method
        p%short_name = 'rkssp53'
        p%long_name = '5-stage, 3rd order SSP Runge-Kutta Spiteri-Ruuth'
        p%order = 3
        p%number_of_stages = 5
        p%number_of_registers = 2
        p%strong_stability_preserving = .true.
        p%cfl = 2.65_wp
    end procedure rkssp53_properties

    module procedure rk4_properties
        !! Returns the properties of the [[rk4]] method
        p%short_name = 'rk4'
        p%long_name = 'Classic 4th order Runge-Kutta'
        p%order = 4
        p%number_of_stages = 4
        p%number_of_registers = 4
    end procedure rk4_properties

    module procedure rks4_properties
        !! Returns the properties of the [[rks4]] method
        p%short_name = 'rks4'
        p%long_name = '4th order Runge-Kutta Shanks'
        p%order = 4
        p%number_of_stages = 4
        p%number_of_registers = 4
    end procedure rks4_properties

    module procedure rkr4_properties
        !! Returns the properties of the [[rkr4]] method
        p%short_name = 'rkr4'
        p%long_name = '4th order Runge-Kutta Ralston'
        p%order = 4
        p%number_of_stages = 4
        p%number_of_registers = 4
    end procedure rkr4_properties

    module procedure rkls44_properties
        !! Returns the properties of the [[rkls44]] method
        p%short_name = 'rkls44'
        p%long_name = '4-stage, 4th order low storage non-TVD Runge-Kutta Jiang-Shu'
        p%order = 4
        p%number_of_stages = 4
        p%number_of_registers = 2
        p%low_storage = .true.
    end procedure rkls44_properties

    module procedure rkls54_properties
        !! Returns the properties of the [[rkls54]] method
        p%short_name = 'rkls54'
        p%long_name = '5-stage, 4th order low storage Runge-Kutta Carpenter-Kennedy'
        p%order = 4
        p%number_of_stages = 5
        p%number_of_registers = 2
        p%low_storage = .true.
        p%cfl = 0.32_wp
    end procedure rkls54_properties

    module procedure rkssp54_properties
        !! Returns the properties of the [[rkssp54]] method
        p%short_name = 'rkssp54'
        p%long_name = '5-stage, 4th order SSP Runge-Kutta Spiteri-Ruuth'
        p%order = 4
        p%number_of_stages = 5
        p%number_of_registers = 4
        p%strong_stability_preserving = .true.
        p%cfl = 1.51_wp
    end procedure rkssp54_properties

    module procedure rks5_properties
        !! Returns the properties of the [[rks5]] method
        p%short_name = 'rks5'
        p%long_name = '5th order Runge-Kutta Shanks'
        p%order = 5
        p%number_of_stages = 5
        p%number_of_registers = 5
    end procedure rks5_properties

    module procedure rk5_properties
        !! Returns the properties of the [[rk5]] method
        p%short_name = 'rk5'
        p%long_name = '5th order Runge-Kutta'
        p%order = 5
        p%number_of_stages = 6
        p%number_of_registers = 6
    end procedure rk5_properties

    module procedure rkc5_properties
        !! Returns the properties of the [[rkc5]] method
        p%short_name = 'rkc5'
        p%long_name = '5th order Runge-Kutta Cassity'
        p%order = 5
        p%number_of_stages = 6
        p%number_of_registers = 6
    end procedure rkc5_properties

    module procedure rkl5_properties
        !! Returns the properties of the [[rkl5]] method
        p%short_name = 'rkl5'
        p%long_name = '5th order Runge-Kutta Lawson'
        p%order = 5
        p%number_of_stages = 6
        p%number_of_registers = 6
    end procedure rkl5_properties

    module procedure rkb6_properties
        !! Returns the properties of the [[rkb6]] method
        p%short_name = 'rkb6'
        p%long_name = '6th order Runge-Kutta Butcher'
        p%order = 6
        p%number_of_stages = 7
        p%number_of_registers = 7
    end procedure rkb6_properties

    module procedure rk7_properties
        !! Returns the properties of the [[rk7]] method
        p%short_name = 'rk7'
        p%long_name = '7th order Runge-Kutta Shanks'
        p%order = 7
        p%number_of_stages = 9
        p%number_of_registers = 9
    end procedure rk7_properties

    module procedure rk8_10_properties
        !! Returns the properties of the [[rk8_10]] method
        p%short_name = 'rk8_10'
        p%long_name = '10-stage, 8th order Runge-Kutta Shanks'
        p%order = 8
        p%number_of_stages = 10
        p%number_of_registers = 10
    end procedure rk8_10_properties

    module procedure rkcv8_properties
        !! Returns the properties of the [[rkcv8]] method
        p%short_name = 'rkcv8'
        p%long_name = '11-stage, 8th order Runge-Kutta Cooper-Verner'
        p%order = 8
        p%number_of_stages = 11
        p%number_of_registers = 11
    end procedure rkcv8_properties

    module procedure rk8_12_properties
        !! Returns the properties of the [[rk8_12]] method
        p%short_name = 'rk8_12'
        p%long_name = '12-stage, 8th order Runge-Kutta Shanks'
        p%order = 8
        p%number_of_stages = 12
        p%number_of_registers = 12
    end procedure rk8_12_properties

    module procedure rkz10_properties
        !! Returns the properties of the [[rkz10]] method
        p%short_name = 'rkz10'
        p%long_name = '10th order Runge-Kutta Zhang'
        p%order = 10
        p%number_of_stages = 16
        p%number_of_registers = 16
    end procedure rkz10_properties

    module procedure rko10_properties
        !! Returns the properties of the [[rko10]] method
        p%short_name = 'rko10'
        p%long_name = '10th order Runge-Kutta Ono'
        p%order = 10
        p%number_of_stages = 17
        p%number_of_registers = 17
    end procedure rko10_properties

    module procedure rkh10_properties
        !! Returns the properties of the [[rkh10]] method
        p%short_name = 'rkh10'
        p%long_name = '10th order Runge-Kutta Hairer'
        p%order = 10
        p%number_of_stages = 17
        p%number_of_registers = 17
    end procedure rkh10_properties

!*****************************************************************************************
    end submodule rklib_fixed_properties
!*****************************************************************************************