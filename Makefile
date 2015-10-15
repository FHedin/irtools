all:
	make cpu
	make gpu

cpu:
# 	gcc -O2 -march=native -std=c99 src/ir2dspectral2.c -o ir2d -lm -lfftw3
	gcc -O3 -march=native -std=c99 src/ir2dspectral_many.c -o ir2d_many -lm -lfftw3

gpu:
	nvcc -O3 cuda/irspec.cu -o ir2d_cuda -lcufft

clean:
	rm -f ir2d ir2d_many ir2d_cuda
