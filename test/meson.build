files = executable('files', 'files.cpp', include_directories : incdir)
interpol = executable('interpolations', 'interpolations.cpp', include_directories : incdir)
cheby_exp = executable('chebyshev_expansion', 'chebyshev_expansion.cpp', include_directories : incdir)
cheby_der = executable('chebyshev_derivatives', 'chebyshev_derivatives.cpp', include_directories : incdir)
cheby_der_exp = executable('chebyshev_derivative_expansion', 'chebyshev_derivative_expansion.cpp', include_directories: incdir)
geqdsk = executable('geqdsk', 'geqdsk.cpp', include_directories: incdir)
# vector = executable('vector', 'vector.cpp', include_directories: incdir)
odeint = executable('odeint', 'odeint.cpp', include_directories: incdir)
magnetic_field = executable('magnetic_field', 'magnetic_field.cpp', include_directories: incdir, cpp_args: ['-g', '-O0'])


copy = find_program('scripts/copy.py')
cmp = find_program('cmp')

# Copy data file for testing
run_command(copy, 'data', 'data', check : true)

# TESTS
test('Load and dump files', files, 
	args : [
		join_paths(meson.build_root(), meson.current_build_dir(), 'data/matrix_from_efit.dat'), 
		join_paths(meson.build_root(), meson.current_build_dir(), 'dump_matrix_from_efit.dat')], 
	is_parallel: false)

test('Compare read and dumped', cmp, 
	args : [
		join_paths(meson.build_root(), meson.current_build_dir(), 'data/matrix_from_efit.dat'), 
		join_paths(meson.build_root(), meson.current_build_dir(), 'dump_matrix_from_efit.dat')], 
	is_parallel: false)

# test('Intarpolations', interpol, 
# 	args : [
# 		join_paths(meson.build_root(), meson.current_build_dir(), 'matrix_from_efit.dat')]
# 	)

test('Chebyshev Expansion', cheby_exp, 
	args : [
		join_paths(meson.build_root(), meson.current_build_dir(), 'data/matrix_from_efit.dat'),
		'23',
		'0.02']
	)

test('Derivatives of Chebyshev', cheby_der)