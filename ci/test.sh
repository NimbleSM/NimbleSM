#!/usr/bin/env bash

set -e
. /opt/spack/share/spack/setup-env.sh
spack env activate nimble

pushd /opt/build/NimbleSM
ctest --output-on-failure
popd
