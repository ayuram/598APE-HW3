# Use official lightweight image
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install only C++ build tools and OpenMP
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ make libomp-dev \
    && rm -rf /var/lib/apt/lists/*

# Set work directory
WORKDIR /app

# Copy local source code into container
COPY . .

# Build the project
RUN make

# Set OMP max threads to 128 via command line
ENV OMP_NUM_THREADS=128

# Default command
CMD ["./main.exe", "1000", "5000", "0.5"]
