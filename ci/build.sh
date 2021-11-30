. /opt/spack/share/spack/setup-env.sh
spack env activate /opt/spack-environment

cmake -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_C_COMPILER=gcc-11 \
     -DCMAKE_CXX_COMPILER=g++-11 \
     -DCMAKE_CXX_FLAGS="-Werror" \
     -DNimbleSM_ENABLE_UNIT_TESTS=ON \
     -DNimbleSM_ENABLE_MPI=$NimbleSM_ENABLE_MPI \
     -DNimbleSM_ENABLE_KOKKOS=$NimbleSM_ENABLE_KOKKOS \
     -DNimbleSM_ENABLE_TRILINOS=$NimbleSM_ENABLE_TRILINOS \
     -S /opt/src/NimbleSM -B /opt/build/NimbleSM
cmake --build /opt/build/NimbleSM --parallel $(nproc)