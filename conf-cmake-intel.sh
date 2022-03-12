#!/usr/bin/env bash

home=/Users/snively

export FC=ifort
export F77=ifort
export CC=icc
export CXX=icpc
export CFLAGS="-std=c11"
export CXXFLAGS="-std=c++11"
export FCFLAGS="-O3"
export FFLAGS="-O3 -assume buffered_io"
export CMAKE_C_COMPILER=icc
export CMAKE_C_FLAGS="-O3 -Wall -std=c11"
export CMAKE_CXX_COMPILER=icpc
export CMAKE_CXX_FLAGS="-O3 -Wall -std=c++11"
export CMAKE_Fortran_COMPILER=ifort 
export CMAKE_Fortran_FLAGS="-O3 -Wall -assume buffered_io"
export CMAKE_EXE_LINKER_FLAGS="-Wl,-no_compact_unwind"


# For now, only gcc/g++/gfortran seem to work.

cmake configure    \
    -Dmpi=on \
    -Dclawpack=on \
    ../forestclaw
