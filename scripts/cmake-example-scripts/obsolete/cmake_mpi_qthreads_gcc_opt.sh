rm -f CMakeCache.txt

cmake \
-D CMAKE_BUILD_TYPE:STRING=Release \
-D CMAKE_C_COMPILER:STRING=/path/to/mpicc \
-D CMAKE_CXX_COMPILER:STRING=/path/to/mpicxx \
-D CMAKE_CXX_FLAGS:STRING="-O2 -std=c++1y" \
-D EXODUS_INCLUDE_DIR:PATH=/path/to/seacas/include \
-D EXODUS_LIB_DIR:PATH=/path/to/seacas/lib \
-D USE_PURE_MPI:BOOL=ON \
-D HAVE_QTHREADS=ON \
-D QTHREADS_INCLUDE_DIR:PATH=/path/to/qthreads/include \
-D QTHREADS_LIB_DIR:PATH=/path/to/qthreads/lib \
../../Nimble/
