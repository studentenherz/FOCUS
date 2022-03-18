files = executable('files', 'files.cpp', include_directories : incdir)
executable('interpolations', 'interpolations.cpp', include_directories : incdir)
# executable('chebyshev_expansion', 'chebyshev_expansion.cpp', include_directories : incdir)
# derivatives = executable('chebyshev_derivatives', 'chebyshev_derivatives.cpp', include_directories : incdir)

copy = find_program('copy.py')
cmp = find_program('cmp')

# Copy data file for testing
run_command(copy, 'matrix_from_efit.dat', 'matrix_from_efit.dat', check : true)

# TESTS
test('Load and dump files', files, 
	args : [
		join_paths(meson.build_root(), meson.current_build_dir(), 'matrix_from_efit.dat'), 
		join_paths(meson.build_root(), meson.current_build_dir(), 'dump_matrix_from_efit.dat')], 
	is_parallel: false)

test('Compare read and dumped', cmp, 
	args : [
		join_paths(meson.build_root(), meson.current_build_dir(), 'matrix_from_efit.dat'), 
		join_paths(meson.build_root(), meson.current_build_dir(), 'dump_matrix_from_efit.dat')], 
	is_parallel: false)
# test('Derivatives of Chebyshev', derivatives)