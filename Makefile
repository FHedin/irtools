# Copyright (c) 2013-2015, Florent Hédin, Pierre-André Cazade, Markus Meuwly, 
# and the University of Basel, Switzerland. All rights reserved.
# 
# The 3-clause BSD license is applied to this software.
# 
# See LICENSE.txt
#

# Please be sure to have you env. variables properly configured if you want to build a CUDA version
# your $PATH should contain the path to the nvcc compiler
# and the $LD_LIBRARY_PATH should contain the path to the cuda libraries

CUDA=OFF

COPT= -std=c99 -Wall -Wextra -O2 
LOPT= -lm -lfftw3

CUOPT= -O2
LCUOPT= -lcufft

ifeq ($(CUDA),ON)
	CC=gcc
	CUCC=nvcc
	DEFINES= USE_CUDA
else
	CC=gcc
	CUCC=gcc
	DEFINES= NO_CUDA
endif

TARGET=ir2d
MKDIR=mkdir -p ./obj

all:
ifeq ($(CUDA),ON)
	make cuda
else
	make nocuda
endif
	
nocuda:
	@$(MKDIR)
	gcc $(COPT) -I. -c src/ir2d.c -o obj/ir2d.o
	gcc $(COPT) -D $(DEFINES) -I. -c src/main.c -o obj/main.o
	gcc $(COPT) obj/*.o -o $(TARGET) $(LOPT)
cuda:
	@$(MKDIR)
	gcc $(COPT) -I. -c src/ir2d.c -o obj/ir2d.o
	nvcc -O2 -D $(DEFINES) -I. -c src/ir2dcuda.cu -o obj/ir2dcuda.o
	gcc $(COPT) -D $(DEFINES) -I. -c src/main.c -o obj/main.o
	nvcc -O2 obj/*.o -o $(TARGET) $(LOPT) $(LCUOPT)

clean:
	rm -f $(TARGET) obj/*.o
