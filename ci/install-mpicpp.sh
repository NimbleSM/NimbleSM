set -x

. /opt/spack/share/spack/setup-env.sh
spack env activate nimble

git clone https://github.com/sandialabs/mpicpp.git /opt/src/mpicpp
cd /opt/src/mpicpp && git checkout db872e297a311a5ad361868d1f74e054b7877b68 && cd -

cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=/opt/local -S /opt/src/mpicpp -B /opt/build/mpicpp
cmake --build /opt/build/mpicpp --parallel $(nproc)
cmake --build /opt/build/mpicpp --target install
