rm -f CMakeCache.txt

cmake \
-D CMAKE_BUILD_TYPE:STRING=Release \
-D CMAKE_C_COMPILER:STRING=/path/to/clang \
-D CMAKE_CXX_COMPILER:STRING=/path/to/clang++ \
-D CMAKE_Fortran_COMPILER:STRING=/path/to/gfortran \
-D CMAKE_CXX_FLAGS:STRING="-O2 -std=c++1y" \
-D EXODUS_INCLUDE_DIR:PATH=/path/to/seacas/include \
-D EXODUS_LIB_DIR:PATH=/path/to/seacas/lib \
-D HAVE_BVH=ON \
-D TinyMath_DIR:PATH=/path/to/tinymath \
-D BVH_DIR:PATH=/path/to/bvh \
-D HAVE_DARMA=OFF \
../../Nimble/
