rm -f CMakeCache.txt

cmake \
-D CMAKE_BUILD_TYPE:STRING=Release \
-D CMAKE_C_COMPILER:STRING=/path/to/mpicc \
-D CMAKE_CXX_COMPILER:STRING=/path/to/mpicxx \
-D CMAKE_CXX_FLAGS:STRING="-O3 -std=c++11" \
-D CMAKE_CXX_EXTENSIONS:BOOL=OFF \
-D EXODUS_INCLUDE_DIR:PATH=/path/to/seacas/include \
-D EXODUS_LIB_DIR:PATH=/path/to/seacas/lib \
-D HAVE_KOKKOS:BOOL=ON \
-D KOKKOS_ROOT_DIR:PATH=/path/to/kokkos \
-D KOKKOS_ENABLE_LIBRT:BOOL=OFF \
-D USE_PURE_MPI:BOOL=ON \
../../Nimble/



