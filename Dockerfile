#
# Dockerfile for building C++ project
#
# https://github.com/tatsy/dockerfile/ubuntu/cxx
#

# OS image
FROM ubuntu:16.04

MAINTAINER gmedders "https://github.com/gmedders"

# Install
RUN \
  apt-get update -y &&  \
  apt-get upgrade -y && \
  apt-get dist-upgrade -y && \
  apt-get install build-essential software-properties-common -y && \
  add-apt-repository ppa:ubuntu-toolchain-r/test -y && \
  apt-get update -y && \
  apt-get install gcc-8 g++-8 -y && \
  update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 60 --slave /usr/bin/g++ g++ /usr/bin/g++-8 && \
  update-alternatives --config gcc

RUN apt-get install git -y
RUN apt-get install cmake -y

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

# Define default command.
CMD ["bash"]

## Update CMake to v3.5.1
RUN git clone --depth 1 -b v3.5.1 https://github.com/Kitware/CMake.git
RUN \
  cd CMake && \
  mkdir build && \
  cd build && \
  cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr && \
  make -j4 && \
  make install && \
  ldconfig && \
  cd ../.. && \
  cmake --version

# Install armadillo dependencies
RUN apt-get install gfortran-8 -y
RUN apt-get install liblapack-dev -y

RUN git clone -q --branch=v0.2.20 --depth=1 https://github.com/xianyi/OpenBLAS.git
RUN cd OpenBLAS && \
    make -j4 && \
    make install
RUN ldconfig

# Install armadillo
RUN git clone --depth 1 -b 8.400.x https://gitlab.com/conradsnicta/armadillo-code.git && cd armadillo-code \
    && cmake . \
    && make -j4 \
    && make install

# Install fftw3
RUN apt-get install wget
RUN wget http://www.fftw.org/fftw-3.3.8.tar.gz \
    && tar -xf fftw-3.3.8.tar.gz \
    && cd fftw-3.3.8 \
    && ./configure --prefix=/usr --enable-shared=yes \
    && make -j4 \
    && make install


## Clean working directory
#RUN rm -rf CMake

# Show environments
RUN echo "--- Build Enviroment ---"
RUN gcc --version | grep gcc
RUN g++ --version | grep g++
RUN cmake --version | grep version
RUN echo "------------------------"



#RUN apt-get install libarmadillo-dev libfftw3-dev -y

WORKDIR /app
ADD . /app

RUN mkdir build && cd build && cmake .. && cmake --build . && ctest -VV
