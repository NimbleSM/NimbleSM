set -x

. /opt/spack/share/spack/setup-env.sh
spack env activate nimble

git clone https://github.com/sandialabs/p3a.git /opt/src/p3a
cd /opt/src/p3a && git checkout ee33b596b16a3baba629859096bb42de00bb8bd3 && cd -

cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=/opt/local -D CMAKE_PREFIX_PATH=/opt/local -S /opt/src/p3a -B /opt/build/p3a
cmake --build /opt/build/p3a --parallel $(nproc)
cmake --build /opt/build/p3a --target install
