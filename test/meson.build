executable('load', 'load.cpp', include_directories : incdir)
executable('chebyshev', 'chebyshev.cpp', include_directories : incdir)
executable('interpolations', 'interpolations.cpp', include_directories : incdir)

copy = find_program('copy.py')
run_command(copy, 'matrix_from_efit.dat', 'matrix_from_efit.dat', check : true)