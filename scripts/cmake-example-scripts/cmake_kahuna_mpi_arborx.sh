#
# Pseudo-script to build NimbleSM with MPI, Kokkos, and ArborX
# Note that ArborX requires Kokkos
#
# This example was last updated on July 22, 2021 for kahuna.ca.sandia.gov
# Author: U. Hetmaniuk (ulhetma@sandia.gov)
#

#
# Load modules
#
module load gcc/9.3.0
module load gcc10-support
module load cmake/3.19.1
module load openmpi/3.1.6
module load exodusii/2016-08-09
module load seacas/2020-08-13
#
# Build Kokkos 
#
cd /path/to/kokkos/build
rm -f CMakeCache.txt
cmake \
-DCMAKE_CXX_FLAGS="-fopenmp" \
-DKokkos_ENABLE_OPENMP=TRUE \
-DKokkos_ENABLE_SERIAL=TRUE \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_INSTALL_PREFIX="/projects/wg-nimble/users/ulhetma/local/" \
/path/to/kokkos/CMakeLists.txt
#
make -j4 install
#
# Build ArborX
#
cd /path/to/arborx/build
rm -f CMakeCache.txt
cmake \
-DARBORX_ENABLE_MPI=ON \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_INSTALL_PREFIX="/projects/wg-nimble/users/ulhetma/local/" \
-DCMAKE_CXX_EXTENSIONS=OFF \
-DCMAKE_PREFIX_PATH="/projects/wg-nimble/users/ulhetma/local/" \
/path/to/arborx/CMakeLists.txt
#
make -j4 install
#
# Build NimbleSM
#
cd /path/to/nimblesm/build
rm -f CMakeCache.txt
#
cmake \
-D CMAKE_BUILD_TYPE:STRING=Release \
-DNimbleSM_ENABLE_MPI=ON \
-DNimbleSM_ENABLE_ARBORX=ON \
-DArborX_ROOT="/projects/wg-nimble/users/ulhetma/local" \
-DNimbleSM_ENABLE_KOKKOS=ON \
-DKokkos_ROOT="/projects/wg-nimble/users/ulhetma/local/" \
/path/to/nimblesm/CMakeLists.txt
#
make -j4
#
# Succesful build will create one executable "build/src/NimbleSM"
#
