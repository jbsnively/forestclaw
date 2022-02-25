#!/usr/bin/env bash

home=/Users/snively

# For now, only gcc/g++/gfortran seem to work.

cmake configure    \
    -DCMAKE_C_COMPILER=cc \
    -DCMAKE_C_FLAGS="-O3 -Wall -DFCLAW2D_PATCHDIM=2 -DFCLAW2D_REFINEDIM=2" \
    -DCMAKE_CXX_COMPILER=c++ \
    -DCMAKE_CXX_FLAGS="-O3 -Wall -DFCLAW2D_PATCHDIM=2 -DFCLAW2D_REFINEDIM=2"\
    -DCMAKE_Fortran_COMPILER=gfortran  \
    -DCMAKE_Fortran_FLAGS="-O3 -Wall -Wno-unused-dummy-argument" \
    -DCMAKE_EXE_LINKER_FLAGS="-Wl,-no_compact_unwind" \
    -Dmpi=on \
    -Dclawpack=on \
    ../forestclaw
