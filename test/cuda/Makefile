OPTIONS = -I../include/ -I../extern/include -maxrregcount=24 -O3
SOURCES_PATH = ../test/cuda/

%.exe: $(SOURCES_PATH)%.cu
	nvcc $(OPTIONS) $< -o $@

