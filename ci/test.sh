. /opt/spack/share/spack/setup-env.sh
spack env activate /opt/spack-environment

pushd /opt/build/NimbleSM

ctest --output-on-failure

popd