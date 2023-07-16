#
# generate files for the methods. This generates everything but the step function.
#
# also generates the tables in the Readme (has to be manually pasted in the README.md file)
#

# fixed:
#                 Name       |  Description                                               | Properties | Order | Stages   | Registers | CFL  | Reference
fixed_methods = [('euler'    , 'Euler'                                                        , '     ', 1     , 1        , 1         , 1.0  , '[Euler (1768)](https://archive.org/details/institutionescal020326mbp)'),
                 ('midpoint' , 'Midpoint'                                                     , '     ', 2     , 2        , 2         , None , '?'),
                 ('heun'     , 'Heun'                                                         , '     ', 2     , 2        , 2         , None , '?'),
                 ('rkssp22'  , '2-stage, 2nd order TVD Runge-Kutta Shu-Osher'                 , ' SSP ', 2     , 2        , 1         , 1.0  , '[Shu & Oscher (1988)](https://ntrs.nasa.gov/api/citations/19880014833/downloads/19880014833.pdf)'),
                 ('rk3'      , '3th order Runge-Kutta'                                        , '     ', 3     , 3        , 3         , None , '?'),
                 ('rkssp33'  , '3-stage, 3rd order TVD Runge-Kutta Shu-Osher'                 , ' SSP ', 3     , 3        , 1         , 1.0  , '[Shu & Oscher (1988)](https://ntrs.nasa.gov/api/citations/19880014833/downloads/19880014833.pdf)'),
                 ('rkssp53'  , '5-stage, 3rd order SSP Runge-Kutta Spiteri-Ruuth'             , ' SSP ', 3     , 5        , 2         , 2.65 , '[Ruuth (2006)](https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf)'),
                 ('rk4'      , 'Classic 4th order 4th order Runge-Kutta'                      , '     ', 4     , 4        , 4         , None , '[Kutta (1901)](https://archive.org/stream/zeitschriftfrma12runggoog#page/n449/mode/2up)'),
                 ('rks4'     , '4th order Runge-Kutta Shanks'                                 , '     ', 4     , 4        , 4         , None , '[Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)'),
                 ('rkls44'   , '4-stage, 4th order low storage non-TVD Runge-Kutta Jiang-Shu' , ' LS  ', 4     , 4        , 2         , None , '[Jiang and Shu (1988)](https://ntrs.nasa.gov/api/citations/19960007052/downloads/19960007052.pdf)'),
                 ('rkls54'   , '5-stage, 4th order low storage Runge-Kutta Carpenter-Kennedy' , ' LS  ', 4     , 5        , 2         , 0.32 , '[Carpenter & Kennedy (1994)](https://ntrs.nasa.gov/api/citations/19940028444/downloads/19940028444.pdf)'),
                 ('rkssp54'  , '5-stage, 4th order SSP Runge-Kutta Spiteri-Ruuth'             , ' SSP ', 4     , 5        , 4         , 1.51 , '[Ruuth (2006)](https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf)'),
                 ('rks5'     , '5th order Runge-Kutta Shanks'                                 , '     ', 5     , 5        , 5         , None , '[Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)'),
                 ('rk5'      , '5th order Runge-Kutta'                                        , '     ', 5     , 6        , 6         , None , '?'),
                 ('rkc5'     , '5th order Runge-Kutta Cassity'                                , '     ', 5     , 6        , 6         , None , '[Cassity (1966)](https://epubs.siam.org/doi/10.1137/0703052)'),
                 ('rkl5'     , '5th order Runge-Kutta Lawson'                                 , '     ', 5     , 6        , 6         , None , '[Lawson (1966)](https://epubs.siam.org/doi/abs/10.1137/0703051)'),
                 ('rkb6'     , '6th order Runge-Kutta Butcher'                                , '     ', 6     , 7        , 7         , None , '[Butcher (1963)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/40DFE501CAB781C9AAE1439B6B8F481A/S1446788700023387a.pdf/div-class-title-on-runge-kutta-processes-of-high-order-div.pdf)'),
                 ('rk7'      , '7th order Runge-Kutta Shanks'                                 , '     ', 7     , 9        , 9         , None , '[Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)'),
                 ('rk8_10'   , '10-stage, 8th order Runge-Kutta Shanks'                       , '     ', 8     , 10       , 10        , None , '[Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)'),
                 ('rkcv8'    , '11-stage, 8th order Runge-Kutta Cooper-Verner'                , '     ', 8     , 11       , 11        , None , '[Cooper & Verner (1972)](https://epubs.siam.org/doi/abs/10.1137/0709037)'),
                 ('rk8_12'   , '12-stage, 8th order Runge-Kutta Shanks'                       , '     ', 8     , 12       , 12        , None , '[Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)'),
                 ('rkz10'    , '10th order Runge-Kutta Zhang'                                 , '     ', 10    , 16       , 16        , None , '[Zhang (2019)](https://arxiv.org/abs/1911.00318)'),
                 ('rko10'    , '10th order Runge-Kutta Ono'                                   , '     ', 10    , 17       , 17        , None , '[Ono (2003)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10f_1.pdf)'),
                 ('rkh10'    , '10th order Runge-Kutta Hairer'                                , '     ', 10    , 17       , 17        , None , '[Hairer (2003)](https://www.researchgate.net/publication/31221486_A_Runge-Kutta_Method_of_Order_10)')]

#variable:
#                         Name  |  Description                                         | Properties | Order | Stages | Registers | CFL  | Reference
variable_methods = [('rkbs32'   , 'Bogacki & Shampine 3(2)'                            , 'FSAL    ' , 3     , 4      , 4         , None , '[Bogacki & Shampine (1989)](https://www.sciencedirect.com/science/article/pii/0893965989900797)'),
                    ('rkssp43'  , '4-stage, 3rd order SSP'                             , 'SSP, LS ' , 3     , 4      , 2         , 2.0  , '[Kraaijevanger (1991)](https://doi.org/10.1007/BF01933264), [Conde et al. (2018)](https://doi.org/10.48550/arXiv.1806.08693)'),
                    ('rkf45'    , 'Fehlberg 4(5)'                                      , '        ' , 4     , 6      , 6         , None , '[Fehlberg (1969)](https://ntrs.nasa.gov/api/citations/19690021375/downloads/19690021375.pdf)'),
                    ('rkck54'   , 'Cash & Karp 5(4)'                                   , '        ' , 5     , 6      , 6         , None , '[Cash & Karp (1990)](http://www.elegio.it/mc2/rk/doc/p201-cash-karp.pdf)'),
                    ('rkdp54'   , 'Dormand-Prince 5(4)'                                , 'FSAL    ' , 5     , 7      , 7         , None , '[Dormand & Prince (1980)](https://www.sciencedirect.com/science/article/pii/0771050X80900133?via%3Dihub)'),
                    ('rkt54'    , 'Tsitouras 5(4)'                                     , 'FSAL    ' , 5     , 7      , 7         , None , '[Tsitouras (2011)](https://www.sciencedirect.com/science/article/pii/S0898122111004706/pdf)'),
                    ('rks54'    , 'Stepanov 5(4)'                                      , 'FSAL    ' , 5     , 7      , 7         , None , '[Stepanov (2022)](https://arxiv.org/pdf/2108.12590.pdf)'),
                    ('rkpp54'   , 'Papakostas-PapaGeorgiou 5(4)'                       , 'FSAL    ' , 5     , 7      , 7         , None , '[Papakostas & Papageorgiou (1996)](https://www.jstor.org/stable/2153797)'),
                    ('rkdp65'   , 'Dormand-Prince 6(5)'                                , '        ' , 6     , 8      , 8         , None , '[Dormand & Prince (1981)](https://www.sciencedirect.com/science/article/pii/0771050X81900103)'),
                    ('rkc65'    , 'Calvo 6(5)'                                         , '        ' , 6     , 9      , 9         , None , '[Calvo (1990)](https://www.sciencedirect.com/science/article/pii/089812219090064Q)'),
                    ('rktp64'   , 'Tsitouras & Papakostas NEW6(4)'                     , '        ' , 6     , 7      , 7         , None , '[Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)'),
                    ('rkv65e'   , 'Verner efficient (9,6(5))'                          , 'FSAL    ' , 6     , 9      , 9         , None , '[Verner (1994)](https://www.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT)'),
                    ('rkv65r'   , 'Verner robust (9,6(5))'                             , 'FSAL    ' , 6     , 9      , 9         , None , '[Verner (1994)](https://www.sfu.ca/~jverner/RKV65.IIIXb.Robust.00010102836.081204.RATOnWeb)'),
                    ('rkv65'    , 'Verner 6(5)'                                        , '        ' , 6     , 8      , 8         , None , '[Verner (2006)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK6/RKcoeff6e_3.pdf)'),
                    ('rktf65'   , 'Tsitouras & Famelis 6(5)'                           , 'FSAL    ' , 6     , 9      , 9         , None , '[Tsitouras & Famelis (2006)](http://users.uoa.gr/~tsitourasc/ModifiedRK-ICNAAM2006.pdf)'),
                    ('rktp75'   , 'Tsitouras & Papakostas NEW7(5)'                     , '        ' , 7     , 9      , 9         , None , '[Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)'),
                    ('rktmy7'   , '7th order Tanaka-Muramatsu-Yamashita'               , '        ' , 7     , 10     , 10        , None , '[Tanaka, Muramatsu & Yamashita (1992)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK7/RKcoeff7d_4.pdf)'),
                    ('rkv76e'   , 'Verner efficient (10:7(6))'                         , '        ' , 7     , 10     , 10        , None , '[Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)'),
                    ('rkv76r'   , 'Verner robust (10:7(6))'                            , '        ' , 7     , 10     , 10        , None , '[Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)'),
                    ('rkf78'    , 'Fehlberg 7(8)'                                      , '        ' , 7     , 13     , 13        , None , '[Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)'),
                    ('rkv78'    , 'Verner 7(8)'                                        , '        ' , 7     , 13     , 13        , None , '[Verner (1978)](https://www.jstor.org/stable/2156853)'),
                    ('dverk78'  , 'Verner "Maple" 7(8)'                                , '        ' , 7     , 13     , 13        , None , '[Verner (?)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK8/RKcoeff8c_2.pdf)'),
                    ('rktp86'   , 'Tsitouras & Papakostas NEW8(6)'                     , '        ' , 8     , 12     , 12        , None , '[Tsitouras & Papakostas (1999)](https://epubs.siam.org/doi/abs/10.1137/S1064827596302230?journalCode=sjoce3)'),
                    ('rkdp87'   , 'Dormand & Prince RK8(7)13M'                         , '        ' , 8     , 13     , 13        , None , '[Prince & Dormand (1981)](https://www.sciencedirect.com/science/article/pii/0771050X81900103)'),
                    ('rkv87e'   , 'Verner efficient (8)7'                              , '        ' , 8     , 13     , 13        , None , '[Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)'),
                    ('rkv87r'   , 'Verner robust (8)7'                                 , '        ' , 8     , 13     , 13        , None , '[Verner (1978)](https://epubs.siam.org/doi/10.1137/0715051)'),
                    ('rkk87'    , 'Kovalnogov-Fedorov-Karpukhina-Simos-Tsitouras 8(7)' , '        ' , 8     , 13     , 13        , None , '[Kovalnogov, Fedorov, Karpukhina, Simos, Tsitouras (2022)](https://www.researchgate.net/publication/363396601_Runge-Kutta_Embedded_Methods_of_Orders_87_for_Use_in_Quadruple_Precision_Computations)'),
                    ('rkf89'    , 'Fehlberg 8(9)'                                      , '        ' , 8     , 17     , 17        , None , '[Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281)'),
                    ('rkv89'    , 'Verner 8(9)'                                        , '        ' , 8     , 16     , 16        , None , '[Verner (1978)](https://www.jstor.org/stable/2156853)'),
                    ('rkt98a'   , 'Tsitouras 9(8) A'                                   , '        ' , 9     , 16     , 16        , None , '[Tsitouras (2001)](https://www.sciencedirect.com/science/article/abs/pii/S0168927401000253)'),
                    ('rkv98e'   , 'Verner efficient (16:9(8))'                         , '        ' , 9     , 16     , 16        , None , '[Verner (1978)](https://www.jstor.org/stable/2156853)'),
                    ('rkv98r'   , 'Verner robust (16:9(8))'                            , '        ' , 9     , 16     , 16        , None , '[Verner (1978)](https://www.jstor.org/stable/2156853)'),
                    ('rkf108'   , 'Feagin 8(10)'                                       , '        ' , 10    , 17     , 17        , None , '[Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk108.txt)'),
                    ('rkc108'   , 'Curtis 10(8)'                                       , '        ' , 10    , 21     , 21        , None , '[Curtis (1975)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10a(8)_2.pdf)'),
                    ('rkb109'   , 'Baker 10(9)'                                        , '        ' , 10    , 21     , 21        , None , '[Baker (?)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10c_1.pdf)'),
                    ('rks1110a' , 'Stone 11(10)'                                       , '        ' , 11    , 26     , 26        , None , '[Stone (2015)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK11/RKcoeff11_a.pdf)'),
                    ('rkf1210'  , 'Feagin 12(10)'                                      , '        ' , 12    , 25     , 25        , None , '[Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1210.txt)'),
                    ('rko129'   , 'Ono 12(9)'                                          , '        ' , 12    , 29     , 29        , None , '[Ono (2006)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK12/RKcoeff12h(9)_1.pdf)'),
                    ('rkf1412'  , 'Feagin 14(12)'                                      , '        ' , 14    , 35     , 35        , None , '[Feagin (2006)](https://sce.uhcl.edu/rungekutta/rk1412.txt)') ]

################################################################################################
def header(fixed_or_variable : str):
    return f"""
!*****************************************************************************************
!>
!  Procedures specifying the properties of the various {fixed_or_variable}-step RK methods.

    submodule(rklib_module) rklib_{fixed_or_variable}_properties

    implicit none

    contains
!*****************************************************************************************
"""

################################################################################################
def footer(fixed_or_variable : str):
    return f"""!*****************************************************************************************
    end submodule rklib_{fixed_or_variable}_properties
!*****************************************************************************************"""

################################################################################################
def write_property_file(fixed_or_variable : str, methods : list):
    """Creates the property methods in a submodule"""
    with open(f'./src/rklib_{fixed_or_variable}_properties.f90', 'w') as f:
        f.write(header(fixed_or_variable)+'\n')
        for m in methods:
            short_name, long_name, props, order, stages, registers, cfl, reference = m
            f.write(f'    module procedure {short_name}_properties\n')
            f.write(f'        !! Returns the properties of the [[{short_name}]] method\n')
            f.write(f'        p%short_name = \'{short_name}\'\n')
            f.write(f'        p%long_name = \'{long_name}\'\n')
            f.write(f'        p%order = {order}\n')
            f.write(f'        p%number_of_stages = {stages}\n')
            if registers:
                f.write(f'        p%number_of_registers = {registers}\n')
            if 'FSAL' in props:
                f.write(f'        p%fsal = .true.\n')
            if 'LS' in props:
                f.write(f'        p%low_storage = .true.\n')
            if 'SSP' in props:
                f.write(f'        p%strong_stability_preserving = .true.\n')
            if cfl:
                f.write(f'        p%cfl = {float(cfl)}_wp\n')
            f.write(f'    end procedure {short_name}_properties\n\n')
        f.write(footer(fixed_or_variable))

################################################################################################
def write_property_interface_file(fixed_or_variable : str, methods : list):
    """Interfaces for the property methods (creates an include file)"""
    with open(f'./src/rklib_{fixed_or_variable}_property_interfaces.i90', 'w') as f:
        for m in methods:
            short_name, long_name, props, order, stages, registers, cfl, reference = m
            f.write(f'pure module function {short_name}_properties(me) result(p)\n')
            f.write(f'    implicit none\n')
            f.write(f'    class({short_name}_class),intent(in) :: me\n')
            f.write(f'    type(rklib_properties) :: p !! properties of the method\n')
            f.write(f'end function {short_name}_properties\n\n')

################################################################################################
def write_class_file(fixed_or_variable : str, methods : list):
    """Defines the integrator classes (creates an include file)"""
    with open(f'./src/rklib_{fixed_or_variable}_classes.i90', 'w') as f:
        f.write(f'    ! {fixed_or_variable.capitalize()} step methods:\n\n')
        for m in methods:
            short_name, long_name, props, order, stages, registers, cfl, reference = m
            fsal = 'FSAL' in props
            if fsal:
                f.write(f'    type,extends(rk_{fixed_or_variable}_step_fsal_class),public :: {short_name}_class\n')
            else:
                f.write(f'    type,extends(rk_{fixed_or_variable}_step_class),public :: {short_name}_class\n')
            f.write(f'        !! {long_name}\n')
            f.write(f'        contains\n')
            f.write(f'        procedure :: step => {short_name}\n')
            f.write(f'        procedure :: properties => {short_name}_properties\n')
            f.write(f'    end type {short_name}_class\n\n')

def write_step_interface_file(fixed_or_variable : str, methods : list):
    """Interfaces for the step methods (creates an include file)"""
    with open(f'./src/rklib_{fixed_or_variable}_step_interfaces.i90', 'w') as f:
        f.write(f'    ! {fixed_or_variable} step interfaces\n\n')
        for m in methods:
            short_name, long_name, props, order, stages, registers, cfl, reference = m
            if (fixed_or_variable=='variable'):
                f.write(f'    module subroutine {short_name}(me,t,x,h,xf,xerr)\n')
            else:
                f.write(f'    module subroutine {short_name}(me,t,x,h,xf)\n')
            f.write(f'        implicit none\n')
            f.write(f'        class({short_name}_class),intent(inout) :: me\n')
            f.write(f'        real(wp),intent(in) :: t !! initial time\n')
            f.write(f'        real(wp),dimension(me%n),intent(in) :: x !! initial state\n')
            f.write(f'        real(wp),intent(in) :: h !! time step\n')
            f.write(f'        real(wp),dimension(me%n),intent(out) :: xf !! state at time `t+h`\n')
            if (fixed_or_variable=='variable'):
                f.write(f'        real(wp),dimension(me%n),intent(out) :: xerr !! truncation error estimate for `x`\n')
            f.write(f'    end subroutine {short_name}\n\n')

################################################################################################
def write_readme_tables(fixed_or_variable : str, methods : list):
    """generate the tables in the readme"""

    print(f'### {fixed_or_variable.capitalize()}-step methods:\n')
    print(f'Name       | Description| Properties | Order | Stages   | Registers | CFL  | Reference')
    print(f'---        | ---        | ---        | ---   | ---      | ---       | ---  | ---')
    for m in methods:
        short_name, long_name, props, order, stages, registers, cfl, reference = m
        print(f'`{short_name}` | {long_name} | {props.strip()} | {order} | {stages} | {registers} | {cfl} | {reference}'.replace('None',''))
    print('')

################################################################################################
def run_all(fixed_or_variable : str, methods : list):
    """Generate all the files"""
    write_property_file(fixed_or_variable, methods)
    write_property_interface_file(fixed_or_variable, methods)
    write_class_file(fixed_or_variable, methods)
    write_step_interface_file(fixed_or_variable, methods)
    write_readme_tables(fixed_or_variable, methods)  # this one just prints to the console

################################################################################################

run_all('fixed',    fixed_methods)
run_all('variable', variable_methods)
