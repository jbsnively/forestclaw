#!/usr/bin/env bash

home=/Users/snively

# For now, only gcc/g++/gfortran seem to work.

cmake configure    \
    -DCMAKE_C_COMPILER=icc \
    -DCMAKE_C_FLAGS="-O3 -Wall -std=c11" \
    -DCMAKE_CXX_COMPILER=icpc \
    -DCMAKE_CXX_FLAGS="-O3 -Wall -std=c11"\
    -DCMAKE_Fortran_COMPILER=ifort  \
    -DCMAKE_Fortran_FLAGS="-O3 -assume buffered_io" \
    -DCMAKE_EXE_LINKER_FLAGS="-Wl" \
    -Dmpi=on \
    -Dclawpack=on \
    ../forestclaw
