set -e
set -x

if [ "$NimbleSM_CONFIGURATION_NAME" = "NimbleSMMPI+Kokkos+ArborX" ]; then \
     echo "Using MPI+KOKKOS+ARBORX NimbleSM configuration."; \
     export NimbleSM_SPACK_ENV_NAME="nimble-mpi-kokkos-arborx"; \

elif [ "$NimbleSM_CONFIGURATION_NAME" = "NimbleSMMPI+Kokkos" ]; then \
     echo "Using MPI+KOKKOS NimbleSM configuration."; \
     export NimbleSM_SPACK_ENV_NAME="nimble-mpi-kokkos"; \

elif [ "$NimbleSM_CONFIGURATION_NAME" = "NimbleSMMPI+Trilinos" ]; then \
     echo "Using MPI+TRILINOS NimbleSM configuration."; \
     export NimbleSM_SPACK_ENV_NAME="nimble-mpi-trilinos"; \

elif [ "$NimbleSM_CONFIGURATION_NAME" = "NimbleSMMPI" ]; then \
     echo "Using MPI NimbleSM configuration."; \
     export NimbleSM_SPACK_ENV_NAME="nimble-mpi"; \

else \
     echo "No existing NimbleSM configuration with provided arguments. Exiting."; \
     exit 1; \
fi
. /opt/spack/share/spack/setup-env.sh
spack env activate $NimbleSM_SPACK_ENV_NAME

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
     -DNimbleSM_ENABLE_P3A=ON \
     -S /opt/src/NimbleSM -B /opt/build/NimbleSM
cmake --build /opt/build/NimbleSM --parallel $(nproc)
