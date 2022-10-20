export FCLAW_LOCAL_PATH=/Users/snively/projects/forestclaw/local/
export CMAKE_BUILD_TYPE=Release
cd zlib
cmake -B build -DCMAKE_PREFIX_PATH=$FCLAW_LOCAL_PATH -DCMAKE_INSTALL_PREFIX=$FCLAW_LOCAL_PATH -DZLIB_COMPAT=TRUE -DZLIB_ENABLE_TESTS=FALSE 
-DCMAKE_POSITION_INDEPENDENT_CODE=TRUE
cmake --build build
cmake --install build
cd ../sc
cmake -B build -Dmpi=TRUE -DCMAKE_PREFIX_PATH=$FCLAW_LOCAL_PATH -DCMAKE_INSTALL_PREFIX=$FCLAW_LOCAL_PATH
cmake --build build
cmake --install build
cd ../p4est
cmake -B build -Dmpi=TRUE -DCMAKE_PREFIX_PATH=$FCLAW_LOCAL_PATH -DCMAKE_INSTALL_PREFIX=$FCLAW_LOCAL_PATH
cmake --build build
cmake --install build
cd ../
cmake -B build -Dmpi=TRUE -DCMAKE_PREFIX_PATH=$FCLAW_LOCAL_PATH -DCMAKE_INSTALL_PREFIX=$FCLAW_LOCAL_PATH -Dclawpatch=TRUE -Dclawpack=TRUE -Dapplications=TRUE
cmake --build build
cmake --install build

