FROM ubuntu:20.04

ARG NIMBLE_MPI_CONFIGURATION

RUN apt-get update \
  && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
      git \
      python3 \
      python3-distutils \
      xz-utils \
      bzip2 \
      zip \
      gpg \
      wget \
      gpgconf \
      software-properties-common \
      libsigsegv2 \
      libsigsegv-dev \
      pkg-config \
      zlib1g \
      zlib1g-dev \
      m4 \
  && rm -rf /var/lib/apt/lists/*

# Cmake ppa
RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null
RUN echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ focal main' | tee /etc/apt/sources.list.d/kitware.list >/dev/null

# gcc ppa
RUN add-apt-repository ppa:ubuntu-toolchain-r/test

RUN apt-get update \
  && apt-get install -y \
     gcc-11 \
     g++-11 \
     gfortran-11 \
     cmake-data=3.21.3-0kitware1ubuntu20.04.1 \
     cmake=3.21.3-0kitware1ubuntu20.04.1 \
  && rm -rf /var/lib/apt/lists/*

# Now we install spack and find compilers/externals
RUN mkdir -p /opt/ && cd /opt/ && git clone https://github.com/spack/spack.git
RUN . /opt/spack/share/spack/setup-env.sh && spack compiler find
RUN . /opt/spack/share/spack/setup-env.sh && spack external find --not-buildable && spack external list
RUN mkdir -p /opt/spack-environment
ADD ./ci/spack-depends-$NIMBLE_MPI_CONFIGURATION.yml /opt/spack-environment/spack.yaml
RUN cd /opt/spack-environment \
  && . /opt/spack/share/spack/setup-env.sh \
  && spack env activate . \
  && spack install --fail-fast \
  && spack gc -y

# Add current source dir into the image
ADD . /opt/src/NimbleSM
RUN mkdir -p /opt/build/NimbleSM

# Build using the spack environment we created
RUN . /opt/spack/share/spack/setup-env.sh \
  && spack env activate /opt/spack-environment \
  && cmake -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_C_COMPILER=gcc-11 \
     -DCMAKE_CXX_COMPILER=g++-11 \
     -DCMAKE_CXX_FLAGS="-Werror" \
     -DNimbleSM_ENABLE_UNIT_TESTS=ON \
     -S /opt/src/NimbleSM -B /opt/builds/NimbleSM \
  && cmake --build /opt/builds/NimbleSM --parallel $(nproc)

RUN . /opt/spack/share/spack/setup-env.sh \
  && spack env activate /opt/spack-environment \
  && cd /opt/builds/NimbleSM && ctest --output-on-failure