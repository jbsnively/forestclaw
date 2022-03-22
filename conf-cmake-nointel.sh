#!/usr/bin/env bash

export FC=gfortran
export F77=gfortran
export CC=cc
export CXX=c++
export CFLAGS="-O3 -std=gnu11 -Wno-error=implicit-function-declaration"
export CXXFLAGS="-O3 -Wno-error=implicit-function-declaration"
export FCFLAGS="-O3 -Wall -Wno-unused-dummy-argument  -fdefault-double-8 -funsafe-math-optimizations"
export FFLAGS="-O3 -Wall -Wno-unused-dummy-argument  -fdefault-double-8 -funsafe-math-optimizations"
export LDFLAGS=' -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -F/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/'
export LIBRARY_PATH=$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
export CMAKE_C_COMPILER=cc
export CMAKE_C_FLAGS="-O3 -Wall -std=gnu11 -Wno-error=implicit-function-declaration"
export CMAKE_CXX_COMPILER=c++
export CMAKE_CXX_FLAGS="-O3 -Wall -Wno-error=implicit-function-declaration"
export CMAKE_Fortran_COMPILER=gfortran
export CMAKE_Fortran_FLAGS="-O3 -Wall -Wno-unused-dummy-argument -fdefault-double-8 -funsafe-math-optimizations"
export CMAKE_EXE_LINKER_FLAGS="-Wl,-no_compact_unwind"


cmake configure    \
    -Dmpi=on \
    -Dclawpack=on \
    ../forestclaw
