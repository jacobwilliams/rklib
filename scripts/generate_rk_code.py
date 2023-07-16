#
#  Generate a template Runge-Kutta Fortran procedure for a p(p-1) type method.
#
#  Author: Jacob Williams, 2023
#

import sys

def generate(n_steps : int, name : str, associate : bool = True):
    """Generate the code"""

    print(f'module procedure {name}\n')

    for i in range(2,n_steps+1):
        print(f'real(wp),parameter :: a{i} = 0')
    print('')

    for i in range(2,n_steps+1):
        for j in range(1,i):
            print(f'real(wp),parameter :: b{i}{j} = 0')
    print('')

    for i in range(1,n_steps+1):
        print(f'real(wp),parameter :: c{i} = 0')
    print('')

    for i in range(1,n_steps+1):
        print(f'real(wp),parameter :: d{i} = 0')
    print('')

    for i in range(1,n_steps+1):
        print(f'real(wp),parameter :: e{i}  = c{i}  - d{i}')
    print('')

    if not associate:
        s = 'real(wp),dimension(me%n) :: '
        for i in range(1,n_steps+1):
            s = f'{s}f{i},'
        s = s.strip(',')
        print(s)
        # print('')

    else:
        s = 'associate ('
        for i in range(1,n_steps+1):
            if i>1:
                indent = '           '
            else:
                indent = ''
            s = f'{s}{indent}f{i} => me%funcs(:,{i})'
            if i<n_steps:
                s = f'{s}, &\n'
            else:
                s = f'{s})\n'
        print(s)

    # print('if (h==zero) then')
    # print('    xf = x')
    # print('    xerr = zero')
    # print('    return')
    # print('end if')
    # print('')

    for i in range(n_steps):
        j = i + 1

        t = f't+a{j}*h'
        f = f'f{j}'

        if j==1:
            x = 'x'
        else:
            x = f'x+h*('
            for k in range(j-1):
                x = f'{x}b{j}{k+1}*f{k+1}+'
            x = x.strip('+')
            x = f'{x})'

        s = f'call me%f({t},{x},{f})'

        print(s)
    print('')

    s = 'xf = x+h*('
    for i in range(1,n_steps+1):
        s = f'{s}c{i}*f{i}+'
    s = s.strip('+')
    s = f'{s})'
    print(s)
    print('')

    s = 'xerr = h*('
    for i in range(1,n_steps+1):
        s = f'{s}e{i}*f{i}+'
    s = s.strip('+')
    s = f'{s})'
    print(s)

    if associate:
        print('\nend associate')

    print(f'\nend procedure {name}')

############################################
if __name__ == "__main__":
    """
    Usage: generate_rk_code n_steps name
    """

    n_steps = int(sys.argv[1])  # number of steps
    name    = sys.argv[2]

    generate(n_steps, name)