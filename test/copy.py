#!/usr/bin/env python3 
'''
	Copy file from source to build directory.
	Usage: copy <input_file> <output_file>
'''
import os, sys, shutil

# Get absolute input and output paths
input_path = os.path.join(os.getenv('MESON_SOURCE_ROOT'), os.getenv('MESON_SUBDIR'), sys.argv[1])
output_path = os.path.join(os.getenv('MESON_BUILD_ROOT'), os.getenv('MESON_SUBDIR'), sys.argv[2])

# Make sure destination directory exists
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# Finally copy the file
shutil.copyfile(input_path, output_path)