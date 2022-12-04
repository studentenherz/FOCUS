options = [
	'-Xcompiler', 
	'-Wno-pedantic',
  '-maxrregcount=24',
  '-O3',
]

# Collect the source files.
c = run_command('grabber.sh', '*.cu', check: true)
sources = c.stdout().strip().split('\n')

# Create a build target for each source file.
foreach source:sources
	exe_name = source.substring(0, -3)
  executable(exe_name, source, 
		cuda_args: options, 
		include_directories: incdir, 
		implicit_include_directories: false,
	)
endforeach
