#!/usr/bin/env bash

home=/Users/snively


export FC=gfortran
export F77=gfortran
export CC=cc
export CXX=c++
export FCFLAGS="-O3"
export FFLAGS="-O3"
CMAKE_C_COMPILER=cc
CMAKE_C_FLAGS="-O3 -Wall"
CMAKE_CXX_COMPILER=c++
CMAKE_CXX_FLAGS="-O3 -Wall"
CMAKE_Fortran_COMPILER=gfortran
CMAKE_Fortran_FLAGS="-O3 -Wall -Wno-unused-dummy-argument"
CMAKE_EXE_LINKER_FLAGS="-Wl,-no_compact_unwind"


# For now, only gcc/g++/gfortran seem to work.

cmake configure    \
    -Dmpi=on \
    -Dclawpack=on \
    ../forestclaw
