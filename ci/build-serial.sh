set -e
set -x
. /opt/spack/share/spack/setup-env.sh
spack env activate nimble

cmake -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_C_COMPILER=gcc-11 \
     -DCMAKE_CXX_COMPILER=g++-11 \
     -DCMAKE_CXX_FLAGS="-Werror" \
     -DNimbleSM_ENABLE_UNIT_TESTS=ON \
     -DNimbleSM_ENABLE_MPI=OFF \
     -DNimbleSM_ENABLE_KOKKOS=OFF \
     -DNimbleSM_ENABLE_TRILINOS=OFF \
     -DNimbleSM_ENABLE_UQ=OFF \
     -DNimbleSM_ENABLE_ARBORX=OFF \
     -S /opt/src/NimbleSM -B /opt/build/NimbleSM
cmake --build /opt/build/NimbleSM --parallel $(nproc)
