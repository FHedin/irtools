all:
	make cpu
# 	make gpu

cpu:
	gcc -O2 src/ir2dspectral2.c -o ir2d -lm -lfftw3
	gcc -O2 src/ir2dspectral_many.c -o ird2d_many -lm -lfftw3

# gpu:
# 	nvcc -O2 cuda/irspec.cu -o ir2d_cuda -lcufft

clean:
	rm -f ir2d ird2d_many ir2d_cuda
