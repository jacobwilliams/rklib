#
# generate RK Fortran code.
#

def generate(n_steps : int):
    """Generate the code"""

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
    for i in range(n_steps):
        s = f'{s}c{i+1}*f{i+1}+'
    s = s.strip('+')
    s = f'{s})'
    print(s)
    print('')

    s = 'terr = h*('
    for i in range(n_steps):
        s = f'{s}e{i+1}*f{i+1}+'
    s = s.strip('+')
    s = f'{s})'
    print(s)

############################################
# generate(29)
generate(21)
