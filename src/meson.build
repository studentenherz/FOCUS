git_version = run_command('git', 'describe', '--dirty=*', '--always',  '--match', '"NOT A TAG"', capture: true, check:true)
focus_serial = executable('focus', 'focus.cpp', include_directories: incdir, implicit_include_directories: false, cpp_args: '-DGIT_VERSION="@0@"'.format(git_version.stdout().strip()))
