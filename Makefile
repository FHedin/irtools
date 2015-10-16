all:
	make cpu
	make gpu

cpu:
	gcc -O3 -march=native -ffast-math -std=c99 src/ir2dspectra.c -o ir2d -lm -lfftw3

gpu:
	nvcc -O3 --use_fast_math cuda/irspec.cu -o ir2d_cuda -lcufft

clean:
	rm -f ir2d ir2d_cuda
