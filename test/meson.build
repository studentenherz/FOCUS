files = executable('files', 'files.cpp', include_directories : incdir)
# executable('interpolations', 'interpolations.cpp', include_directories : incdir)
# executable('chebyshev_expansion', 'chebyshev_expansion.cpp', include_directories : incdir)
# derivatives = executable('chebyshev_derivatives', 'chebyshev_derivatives.cpp', include_directories : incdir)

copy = find_program('copy.py')
run_command(copy, 'matrix_from_efit.dat', 'matrix_from_efit.dat', check : true)

test('Load and dump files', files, args : 'matrix_from_efit.dat')
# test('Derivatives of Chebyshev', derivatives)