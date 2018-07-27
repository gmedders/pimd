# OS image
FROM gmedders/cpp_scientific_libraries:v0.1.2

MAINTAINER gmedders "https://github.com/gmedders"

# Define default command.
CMD ["bash"]

WORKDIR /app
ADD . /app

# Build the packaged
RUN mkdir build && cd build && cmake .. && cmake --build .

# When the image is run, go directly to the executables
WORKDIR /app/build

# Define default command.
CMD ["bash"]
