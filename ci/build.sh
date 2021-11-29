. /opt/spack/share/spack/setup-env.sh
spack env activate /opt/spack-environment

case $NIMBLE_MPI_CONFIGURATION in
  serial)
    NIMBLESM_ENABLE_MPI=OFF ;;
  mpi)
    NIMBLESM_ENABLE_MPI=ON ;;
esac

cmake -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_C_COMPILER=gcc-11 \
     -DCMAKE_CXX_COMPILER=g++-11 \
     -DCMAKE_CXX_FLAGS="-Werror" \
     -DNimbleSM_ENABLE_UNIT_TESTS=ON \
     -DNimbleSM_ENABLE_MPI=$NIMBLESM_ENABLE_MPI
     -S /opt/src/NimbleSM -B /opt/builds/NimbleSM
cmake --build /opt/builds/NimbleSM --parallel $(nproc)