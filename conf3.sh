export LIBS=""
export FCLAW=/Users/snively/projects/forestclaw
export FCLAW_BIN=$FCLAW/local/bin

./configure \
    --enable-mpi \
    --enable-clawpack3 \
    --enable-magic3d \
    --disable-shared \
    --without-blas \
    F77="mpif77" \
    CC="mpicc" \
    CXX="mpicxx" \
    CFLAGS="-Wall -O3 -DFCLAW2D_PATCHDIM=3 -DFCLAW2D_REFINEDIM=2" \
    CXXFLAGS="-O3 -Wall -DFCLAW2D_PATCHDIM=3 -DFCLAW2D_REFINEDIM=2" \
    FFLAGS="-fast -assume buffered_io" 
#    LIBS="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential"
#    --enable-clawpack \


