#!/usr/bin/env python3 
'''
	Converts a matrix file to a `x y value`
	file, where x and y are the indexes 
	Usage: copy <input_file> <output_file>
'''
import sys

fname_in = sys.argv[1]
fname_out = sys.argv[2]

with open(fname_in, 'r') as file_in:
	with open(fname_out, 'w') as file_out:
		lines = file_in.readlines()
		for i, line in enumerate(lines):
			values = line.strip().split(' ')
			for j, value in enumerate(values):
				file_out.write(f'{i} {j} {value}\n')
			file_out.write('\n')