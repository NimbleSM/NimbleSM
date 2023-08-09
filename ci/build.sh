set -e
set -x
. /opt/spack/share/spack/setup-env.sh
spack env activate nimble

cmake -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_C_COMPILER=gcc-11 \
     -DCMAKE_CXX_COMPILER=g++-11 \
     -DCMAKE_CXX_FLAGS="-Werror" \
     -DCMAKE_PREFIX_PATH=/opt/local \
     -DNimbleSM_ENABLE_UNIT_TESTS=ON \
     -DNimbleSM_ENABLE_MPI=$NimbleSM_ENABLE_MPI \
     -DNimbleSM_ENABLE_KOKKOS=$NimbleSM_ENABLE_KOKKOS \
     -DNimbleSM_ENABLE_TRILINOS=$NimbleSM_ENABLE_TRILINOS \
     -DNimbleSM_ENABLE_UQ=$NimbleSM_ENABLE_UQ \
     -DNimbleSM_ENABLE_ARBORX=$NimbleSM_ENABLE_ARBORX \
     -S /opt/src/NimbleSM -B /opt/build/NimbleSM
cmake --build /opt/build/NimbleSM --parallel $(nproc)
