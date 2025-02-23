
!*****************************************************************************************
!>
!  Procedures specifying the properties of the various variable-step RK methods.

    submodule(rklib_module) rklib_variable_properties

    implicit none

    contains
!*****************************************************************************************

    module procedure rkbs32_properties
        !! Returns the properties of the [[rkbs32]] method
        p%short_name = 'rkbs32'
        p%long_name = 'Bogacki & Shampine 3(2)'
        p%order = 3
        p%number_of_stages = 4
        p%number_of_registers = 4
        p%fsal = .true.
    end procedure rkbs32_properties

    module procedure rkssp43_properties
        !! Returns the properties of the [[rkssp43]] method
        p%short_name = 'rkssp43'
        p%long_name = '4-stage, 3rd order SSP'
        p%order = 3
        p%number_of_stages = 4
        p%number_of_registers = 2
        p%low_storage = .true.
        p%strong_stability_preserving = .true.
        p%cfl = 2.0_wp
    end procedure rkssp43_properties

    module procedure rkf45_properties
        !! Returns the properties of the [[rkf45]] method
        p%short_name = 'rkf45'
        p%long_name = 'Fehlberg 4(5)'
        p%order = 4
        p%number_of_stages = 6
        p%number_of_registers = 6
    end procedure rkf45_properties

    module procedure rkck54_properties
        !! Returns the properties of the [[rkck54]] method
        p%short_name = 'rkck54'
        p%long_name = 'Cash & Karp 5(4)'
        p%order = 5
        p%number_of_stages = 6
        p%number_of_registers = 6
    end procedure rkck54_properties

    module procedure rkdp54_properties
        !! Returns the properties of the [[rkdp54]] method
        p%short_name = 'rkdp54'
        p%long_name = 'Dormand-Prince 5(4)'
        p%order = 5
        p%number_of_stages = 7
        p%number_of_registers = 7
        p%fsal = .true.
    end procedure rkdp54_properties

    module procedure rkt54_properties
        !! Returns the properties of the [[rkt54]] method
        p%short_name = 'rkt54'
        p%long_name = 'Tsitouras 5(4)'
        p%order = 5
        p%number_of_stages = 7
        p%number_of_registers = 7
        p%fsal = .true.
    end procedure rkt54_properties

    module procedure rks54_properties
        !! Returns the properties of the [[rks54]] method
        p%short_name = 'rks54'
        p%long_name = 'Stepanov 5(4)'
        p%order = 5
        p%number_of_stages = 7
        p%number_of_registers = 7
        p%fsal = .true.
    end procedure rks54_properties

    module procedure rkpp54_properties
        !! Returns the properties of the [[rkpp54]] method
        p%short_name = 'rkpp54'
        p%long_name = 'Papakostas-PapaGeorgiou 5(4)'
        p%order = 5
        p%number_of_stages = 7
        p%number_of_registers = 7
        p%fsal = .true.
    end procedure rkpp54_properties

    module procedure rkpp54b_properties
        !! Returns the properties of the [[rkpp54b]] method
        p%short_name = 'rkpp54b'
        p%long_name = 'Papakostas-PapaGeorgiou 5(4) b'
        p%order = 5
        p%number_of_stages = 7
        p%number_of_registers = 7
        p%fsal = .true.
    end procedure rkpp54b_properties

    module procedure rkbs54_properties
        !! Returns the properties of the [[rkbs54]] method
        p%short_name = 'rkbs54'
        p%long_name = 'Bogacki & Shampine 5(4)'
        p%order = 5
        p%number_of_stages = 8
        p%number_of_registers = 8
    end procedure rkbs54_properties

    module procedure rkss54_properties
        !! Returns the properties of the [[rkss54]] method
        p%short_name = 'rkss54'
        p%long_name = 'Sharp & Smart 5(4)'
        p%order = 5
        p%number_of_stages = 7
        p%number_of_registers = 7
    end procedure rkss54_properties

    module procedure rkdp65_properties
        !! Returns the properties of the [[rkdp65]] method
        p%short_name = 'rkdp65'
        p%long_name = 'Dormand-Prince 6(5)'
        p%order = 6
        p%number_of_stages = 8
        p%number_of_registers = 8
    end procedure rkdp65_properties

    module procedure rkc65_properties
        !! Returns the properties of the [[rkc65]] method
        p%short_name = 'rkc65'
        p%long_name = 'Calvo 6(5)'
        p%order = 6
        p%number_of_stages = 9
        p%number_of_registers = 9
    end procedure rkc65_properties

    module procedure rktp64_properties
        !! Returns the properties of the [[rktp64]] method
        p%short_name = 'rktp64'
        p%long_name = 'Tsitouras & Papakostas NEW6(4)'
        p%order = 6
        p%number_of_stages = 7
        p%number_of_registers = 7
    end procedure rktp64_properties

    module procedure rkv65e_properties
        !! Returns the properties of the [[rkv65e]] method
        p%short_name = 'rkv65e'
        p%long_name = 'Verner efficient (9,6(5))'
        p%order = 6
        p%number_of_stages = 9
        p%number_of_registers = 9
        p%fsal = .true.
    end procedure rkv65e_properties

    module procedure rkv65r_properties
        !! Returns the properties of the [[rkv65r]] method
        p%short_name = 'rkv65r'
        p%long_name = 'Verner robust (9,6(5))'
        p%order = 6
        p%number_of_stages = 9
        p%number_of_registers = 9
        p%fsal = .true.
    end procedure rkv65r_properties

    module procedure rkv65_properties
        !! Returns the properties of the [[rkv65]] method
        p%short_name = 'rkv65'
        p%long_name = 'Verner 6(5)'
        p%order = 6
        p%number_of_stages = 8
        p%number_of_registers = 8
    end procedure rkv65_properties

    module procedure dverk65_properties
        !! Returns the properties of the [[dverk65]] method
        p%short_name = 'dverk65'
        p%long_name = 'Verner 6(5) "DVERK"'
        p%order = 6
        p%number_of_stages = 8
        p%number_of_registers = 8
    end procedure dverk65_properties

    module procedure rktf65_properties
        !! Returns the properties of the [[rktf65]] method
        p%short_name = 'rktf65'
        p%long_name = 'Tsitouras & Famelis 6(5)'
        p%order = 6
        p%number_of_stages = 9
        p%number_of_registers = 9
        p%fsal = .true.
    end procedure rktf65_properties

    module procedure rktp75_properties
        !! Returns the properties of the [[rktp75]] method
        p%short_name = 'rktp75'
        p%long_name = 'Tsitouras & Papakostas NEW7(5)'
        p%order = 7
        p%number_of_stages = 9
        p%number_of_registers = 9
    end procedure rktp75_properties

    module procedure rktmy7_properties
        !! Returns the properties of the [[rktmy7]] method
        p%short_name = 'rktmy7'
        p%long_name = '7th order Tanaka-Muramatsu-Yamashita'
        p%order = 7
        p%number_of_stages = 10
        p%number_of_registers = 10
    end procedure rktmy7_properties

    module procedure rktmy7s_properties
        !! Returns the properties of the [[rktmy7s]] method
        p%short_name = 'rktmy7s'
        p%long_name = '7th order Stable Tanaka-Muramatsu-Yamashita'
        p%order = 7
        p%number_of_stages = 10
        p%number_of_registers = 10
    end procedure rktmy7s_properties

    module procedure rkv76e_properties
        !! Returns the properties of the [[rkv76e]] method
        p%short_name = 'rkv76e'
        p%long_name = 'Verner efficient (10:7(6))'
        p%order = 7
        p%number_of_stages = 10
        p%number_of_registers = 10
    end procedure rkv76e_properties

    module procedure rkv76r_properties
        !! Returns the properties of the [[rkv76r]] method
        p%short_name = 'rkv76r'
        p%long_name = 'Verner robust (10:7(6))'
        p%order = 7
        p%number_of_stages = 10
        p%number_of_registers = 10
    end procedure rkv76r_properties

    module procedure rkss76_properties
        !! Returns the properties of the [[rkss76]] method
        p%short_name = 'rkss76'
        p%long_name = 'Sharp & Smart 7(6)'
        p%order = 7
        p%number_of_stages = 11
        p%number_of_registers = 11
    end procedure rkss76_properties

    module procedure rkf78_properties
        !! Returns the properties of the [[rkf78]] method
        p%short_name = 'rkf78'
        p%long_name = 'Fehlberg 7(8)'
        p%order = 7
        p%number_of_stages = 13
        p%number_of_registers = 13
    end procedure rkf78_properties

    module procedure rkv78_properties
        !! Returns the properties of the [[rkv78]] method
        p%short_name = 'rkv78'
        p%long_name = 'Verner 7(8)'
        p%order = 7
        p%number_of_stages = 13
        p%number_of_registers = 13
    end procedure rkv78_properties

    module procedure dverk78_properties
        !! Returns the properties of the [[dverk78]] method
        p%short_name = 'dverk78'
        p%long_name = 'Verner "Maple" 7(8)'
        p%order = 7
        p%number_of_stages = 13
        p%number_of_registers = 13
    end procedure dverk78_properties

    module procedure rkdp85_properties
        !! Returns the properties of the [[rkdp85]] method
        p%short_name = 'rkdp85'
        p%long_name = 'Dormand-Prince 8(5)'
        p%order = 8
        p%number_of_stages = 12
        p%number_of_registers = 12
    end procedure rkdp85_properties

    module procedure rktp86_properties
        !! Returns the properties of the [[rktp86]] method
        p%short_name = 'rktp86'
        p%long_name = 'Tsitouras & Papakostas NEW8(6)'
        p%order = 8
        p%number_of_stages = 12
        p%number_of_registers = 12
    end procedure rktp86_properties

    module procedure rkdp87_properties
        !! Returns the properties of the [[rkdp87]] method
        p%short_name = 'rkdp87'
        p%long_name = 'Dormand & Prince RK8(7)13M'
        p%order = 8
        p%number_of_stages = 13
        p%number_of_registers = 13
    end procedure rkdp87_properties

    module procedure rkv87e_properties
        !! Returns the properties of the [[rkv87e]] method
        p%short_name = 'rkv87e'
        p%long_name = 'Verner efficient (8)7'
        p%order = 8
        p%number_of_stages = 13
        p%number_of_registers = 13
    end procedure rkv87e_properties

    module procedure rkv87r_properties
        !! Returns the properties of the [[rkv87r]] method
        p%short_name = 'rkv87r'
        p%long_name = 'Verner robust (8)7'
        p%order = 8
        p%number_of_stages = 13
        p%number_of_registers = 13
    end procedure rkv87r_properties

    module procedure rkev87_properties
        !! Returns the properties of the [[rkev87]] method
        p%short_name = 'rkev87'
        p%long_name = 'Enright-Verner (8)7'
        p%order = 8
        p%number_of_stages = 13
        p%number_of_registers = 13
    end procedure rkev87_properties

    module procedure rkk87_properties
        !! Returns the properties of the [[rkk87]] method
        p%short_name = 'rkk87'
        p%long_name = 'Kovalnogov-Fedorov-Karpukhina-Simos-Tsitouras 8(7)'
        p%order = 8
        p%number_of_stages = 13
        p%number_of_registers = 13
    end procedure rkk87_properties

    module procedure rkf89_properties
        !! Returns the properties of the [[rkf89]] method
        p%short_name = 'rkf89'
        p%long_name = 'Fehlberg 8(9)'
        p%order = 8
        p%number_of_stages = 17
        p%number_of_registers = 17
    end procedure rkf89_properties

    module procedure rkv89_properties
        !! Returns the properties of the [[rkv89]] method
        p%short_name = 'rkv89'
        p%long_name = 'Verner 8(9)'
        p%order = 8
        p%number_of_stages = 16
        p%number_of_registers = 16
    end procedure rkv89_properties

    module procedure rkt98a_properties
        !! Returns the properties of the [[rkt98a]] method
        p%short_name = 'rkt98a'
        p%long_name = 'Tsitouras 9(8) A'
        p%order = 9
        p%number_of_stages = 16
        p%number_of_registers = 16
    end procedure rkt98a_properties

    module procedure rkv98e_properties
        !! Returns the properties of the [[rkv98e]] method
        p%short_name = 'rkv98e'
        p%long_name = 'Verner efficient (16:9(8))'
        p%order = 9
        p%number_of_stages = 16
        p%number_of_registers = 16
    end procedure rkv98e_properties

    module procedure rkv98r_properties
        !! Returns the properties of the [[rkv98r]] method
        p%short_name = 'rkv98r'
        p%long_name = 'Verner robust (16:9(8))'
        p%order = 9
        p%number_of_stages = 16
        p%number_of_registers = 16
    end procedure rkv98r_properties

    module procedure rks98_properties
        !! Returns the properties of the [[rks98]] method
        p%short_name = 'rks98'
        p%long_name = 'Sharp 9(8)'
        p%order = 9
        p%number_of_stages = 16
        p%number_of_registers = 16
    end procedure rks98_properties

    module procedure rkf108_properties
        !! Returns the properties of the [[rkf108]] method
        p%short_name = 'rkf108'
        p%long_name = 'Feagin 8(10)'
        p%order = 10
        p%number_of_stages = 17
        p%number_of_registers = 17
    end procedure rkf108_properties

    module procedure rkc108_properties
        !! Returns the properties of the [[rkc108]] method
        p%short_name = 'rkc108'
        p%long_name = 'Curtis 10(8)'
        p%order = 10
        p%number_of_stages = 21
        p%number_of_registers = 21
    end procedure rkc108_properties

    module procedure rkb109_properties
        !! Returns the properties of the [[rkb109]] method
        p%short_name = 'rkb109'
        p%long_name = 'Baker 10(9)'
        p%order = 10
        p%number_of_stages = 21
        p%number_of_registers = 21
    end procedure rkb109_properties

    module procedure rks1110a_properties
        !! Returns the properties of the [[rks1110a]] method
        p%short_name = 'rks1110a'
        p%long_name = 'Stone 11(10)'
        p%order = 11
        p%number_of_stages = 26
        p%number_of_registers = 26
    end procedure rks1110a_properties

    module procedure rkf1210_properties
        !! Returns the properties of the [[rkf1210]] method
        p%short_name = 'rkf1210'
        p%long_name = 'Feagin 12(10)'
        p%order = 12
        p%number_of_stages = 25
        p%number_of_registers = 25
    end procedure rkf1210_properties

    module procedure rko129_properties
        !! Returns the properties of the [[rko129]] method
        p%short_name = 'rko129'
        p%long_name = 'Ono 12(9)'
        p%order = 12
        p%number_of_stages = 29
        p%number_of_registers = 29
    end procedure rko129_properties

    module procedure rkf1412_properties
        !! Returns the properties of the [[rkf1412]] method
        p%short_name = 'rkf1412'
        p%long_name = 'Feagin 14(12)'
        p%order = 14
        p%number_of_stages = 35
        p%number_of_registers = 35
    end procedure rkf1412_properties

!*****************************************************************************************
    end submodule rklib_variable_properties
!*****************************************************************************************