# IR tools : tools for 1 Dimension and 2 Dimensions Infra Red (IR & 2DIR) spectroscopy modelling
------------------------------------------------------------------------------------------------

This repository provides source code and analysis code for producing 2DIR spectra, 
following the approach detailed by Peter Hamm et al in their book : "Concepts and Methods of 2D Infrared Spectroscopy" ISBN: 9781107000056

Written in C99 and CUDA : you need the CUDA toolkit (https://developer.nvidia.com/cuda-toolkit)
for compiling a version using CUDA capable GPUs.

CUDA build is **optional** (see below).

## How to compile
-----------------

* Dependency : 
  * FFTW Version >= 3 : http://www.fftw.org/ : Tested 3.3.4 on 16 Oct 2015 on Fedora linux 21 x86_64
  * If a CUDA version is wanted :
    * CUDA TOOLKIT : Tested version 7 on 16 Oct 2015 on Fedora linux 21 x86_64

## How to build
---------------

* Use the provided makefile : 
  * run the command **make** if you want a version **WITHOUT CUDA capabilities**
  * run the command **make CUDA=ON** if you want a version **WITH CUDA capabilities**


## How to run 
-------------
Once compiled you should get the executable **ir2d**

Use **ir2d -h** for more instructions

You will find 2 input files for both CPU or GPU versions : everything after a ** \# ** is considered as comments

## Licensing and Copyright
--------------------------

Copyright (c) 2013-2015, Florent Hédin, Pierre-André Cazade, Markus Meuwly, 
and the University of Basel, Switzerland. All rights reserved.

The 3-clause BSD license is applied to this software.

See LICENSE.txt

