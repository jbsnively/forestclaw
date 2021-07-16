export LIBS=""
export FCLAW=/Users/snively/projects/forestclaw
export FCLAW_BIN=$FCLAW/local/bin

./configure \
    --enable-mpi \
    --enable-clawpack \
    --enable-magic \
    --disable-shared \
    --without-blas \
    F77="mpif77" \
    CC="mpicc" \
    CXX="mpicxx" \
    CFLAGS="-Wall -O3" \
    CXXFLAGS="-O3 -Wall" \
    FFLAGS="-fast -assume buffered_io" 
#    LIBS="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential"


