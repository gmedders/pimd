# OS image
FROM ubuntu:16.04

MAINTAINER gmedders "https://github.com/gmedders"

# Install g++-8, from https://gist.github.com/application2000/73fd6f4bf1be6600a2cf9f56315a2d91
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

# Install remaining git, cmake, armadillo, and FFTW-3
RUN apt-get install git -y
RUN apt-get install cmake -y
RUN apt-get install libarmadillo-dev libfftw3-dev -y

# Define default command.
CMD ["bash"]

WORKDIR /app
ADD . /app

# Build the packaged
RUN mkdir build && cd build && cmake .. && cmake --build .

# When the image is run, perform the tests
WORKDIR /app/build
ENTRYPOINT ["ctest", "-VV"]
