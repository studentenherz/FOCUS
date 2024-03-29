git_version = run_command('git', 'describe', '--dirty=*', '--always',  '--match', '"NOT A TAG"', capture: true, check:true)

array = executable('failin_array', 'array.cpp', include_directories: incdir)
matrix_2d = executable('failin_matrix_2d', 'matrix_2d.cpp', include_directories: incdir)
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
slow_down = executable('slow_down', 'slow_down.cpp', include_directories: incdir, implicit_include_directories: false, cpp_args: ['-O3', '-DFOCUS_TEST_SLOW_DOWN'])
input_gacode = executable('input_gacode', 'input_gacode.cpp', include_directories: incdir, implicit_include_directories: false)
particle_states = executable('particle_states', 'particle_states.cpp', include_directories: incdir, implicit_include_directories: false, cpp_args: '-DGIT_VERSION="@0@"'.format(git_version.stdout().strip()))
pitch_E_2_state = executable('pitch_E_2_state', 'pitch_E_2_state.cpp', include_directories: incdir, implicit_include_directories: false, cpp_args: '-DGIT_VERSION="@0@"'.format(git_version.stdout().strip()))

subdir('cuda')

copy = find_program('scripts/copy.py')
cmp = find_program('cmp')

# Copy data file for testing
run_command(copy, 'data', 'data', check : true)

# TESTS

test('Out of bound array fails', array,
	should_fail: true
)

test('Out of bound matrix fails', matrix_2d,
	should_fail: true
)

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