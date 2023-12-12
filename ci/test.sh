#!/usr/bin/env bash

set -x
set -e

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

# Ensure that OpenMP environment is set up correctly
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

pushd /opt/build/NimbleSM
ret_code=0
ctest --output-on-failure || ret_code=$?
# We collect the test logs for exporting
echo "ctest returned: $ret_code"
mkdir -p /tmp/artifacts/
cp /opt/build/NimbleSM/Testing/Temporary/LastTest.log /tmp/artifacts/
cp /opt/build/NimbleSM/test/*/*/*.log /tmp/artifacts/
echo ${ret_code} > /tmp/artifacts/success_flag.txt
ls /tmp/artifacts
popd
