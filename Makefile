
CC=gcc
CUCC=nvcc

TARGET=ir2d
MKDIR=mkdir -p ./obj

all:
	@$(MKDIR)
	gcc  -O2 -std=c99 -I. -c src/ir2d.c        -o obj/ir2d.o
	nvcc -O2 -D USE_CUDA -I. -c src/ir2dcuda.cu    -o obj/ir2dcuda.o
	gcc  -O2 -std=c99 -D USE_CUDA -I. -c src/main.c -o obj/main.o
	nvcc -O2 obj/*.o -o $(TARGET) -lfftw3 -lcufft -lm

clean:
	rm -f $(TARGET) obj/*.o
