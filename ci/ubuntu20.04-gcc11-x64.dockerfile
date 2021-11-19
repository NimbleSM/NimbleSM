FROM ubuntu:20.04

RUN apt-get update \
  && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
      git \
      python3 \
      xz-utils \
      bzip2 \
      zip \
      gpg \
      wget \
      gpgconf \
      software-properties-common \
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

# mpich
RUN apt-get update \
  && apt-get install -y \
     mpich \
     libmpich-dev \
  && rm -rf /var/lib/apt/lists/*

# Now we install spack and find compilers/externals
RUN mkdir -p /opt/ && cd /opt/ && git clone https://github.com/spack/spack.git
RUN . /opt/spack/share/spack/setup-env.sh && spack compiler find
RUN . /opt/spack/share/spack/setup-env.sh && spack external find --not-buildable && spack external list
RUN . /opt/spack/share/spack/setup-env.sh && spack install --fail-fast exodusii