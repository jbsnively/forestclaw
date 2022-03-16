#!/usr/bin/env bash

home=/Users/snively


export FC=gfortran
export F77=gfortran
export CC=cc
export CXX=c++
export FCFLAGS="-O3"
export FFLAGS="-O3"
export LDFLAGS=' -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -F/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/'
export LIBRARY_PATH=$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
export CMAKE_C_COMPILER=cc
export CMAKE_C_FLAGS="-O3 -Wall"
export CMAKE_CXX_COMPILER=c++
export CMAKE_CXX_FLAGS="-O3 -Wall"
export CMAKE_Fortran_COMPILER=gfortran
export CMAKE_Fortran_FLAGS="-O3 -Wall -Wno-unused-dummy-argument"
export CMAKE_EXE_LINKER_FLAGS="-Wl,-no_compact_unwind"


cmake configure    \
    -Dmpi=on \
    -Dclawpack=on \
    ../forestclaw
