
rm -f CMakeCache.txt

# DARMA_FRONTEND_INSTALL_PATH=`spack location --install-dir darma-frontend`
# DARMA_THREADS_INSTALL_PATH=`spack location --install-dir darma-threads`

cmake \
-D CMAKE_BUILD_TYPE:STRING=Release \
-D EXODUS_INCLUDE_DIR:PATH=/path/to/seacas/include \
-D EXODUS_LIB_DIR:PATH=/path/to/seacas/lib \
-D HAVE_DARMA=ON \
-D DARMA_BACKEND_PKG=Threads \
../../Nimble/
