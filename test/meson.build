git_version = run_command('git', 'describe', '--dirty=*', '--always',  '--match', '"NOT A TAG"', capture: true, check:true)

files = executable('files', 'files.cpp', include_directories : incdir, implicit_include_directories: false)
interpol = executable('interpolations', 'interpolations.cpp', include_directories : incdir, implicit_include_directories: false)
cheby_exp = executable('chebyshev_expansion', 'chebyshev_expansion.cpp', include_directories : incdir, implicit_include_directories: false)
cheby_der = executable('chebyshev_derivatives', 'chebyshev_derivatives.cpp', include_directories : incdir, implicit_include_directories: false)
cheby_der_exp = executable('chebyshev_derivative_expansion', 'chebyshev_derivative_expansion.cpp', include_directories: incdir)
geqdsk = executable('geqdsk', 'geqdsk.cpp', include_directories: incdir, implicit_include_directories: false)
vector = executable('vector', 'vector.cpp', include_directories: incdir, implicit_include_directories: false)
odeint = executable('odeint', 'odeint.cpp', include_directories: incdir, implicit_include_directories: false)
magnetic_field = executable('magnetic_field', 'magnetic_field.cpp', include_directories: incdir, implicit_include_directories: false)
random = executable('random', 'random.cpp', include_directories: incdir, implicit_include_directories: false)
collisions = executable('collisions', 'collisions.cpp', include_directories: incdir, implicit_include_directories: false)
input_gacode = executable('input_gacode', 'input_gacode.cpp', include_directories: incdir, implicit_include_directories: false)
particle_states = executable('particle_states', 'particle_states.cpp', include_directories: incdir, implicit_include_directories: false, cpp_args: '-DGIT_VERSION="@0@"'.format(git_version.stdout().strip()))


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