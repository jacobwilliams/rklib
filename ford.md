project: rklib
src_dir: ./src
         ./example
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/rklib
summary: Fixed and variable-step Runge-Kutta solvers in Modern Fortran
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
source: true
graph: true
search: true
preprocessor: gfortran -E
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README.md!}