#!/usr/bin/env bash

set -e
. /opt/spack/share/spack/setup-env.sh
spack env activate nimble

pushd /opt/build/NimbleSM
ctest --output-on-failure

cp /opt/build/NimbleSM/Testing/Temporary/LastTest.log /tmp/artifacts
cp /opt/build/NimbleSM/test/*/*/*.log /tmp/artifacts

popd
