# First stage: Build the SURVIVOR binary
FROM python:3.10-slim AS build

# Install build dependencies (gcc, make, etc.) and git to clone repository
RUN apt-get update && apt-get install -y \
    libbz2-dev \
    build-essential \
    git \
    python3-dev \
    procps \
    gcc \
 && rm -rf /var/lib/apt/lists/*

# Clone the SURVIVOR repository and build the binary
RUN git clone https://github.com/fritzsedlazeck/SURVIVOR.git
WORKDIR /SURVIVOR/Debug
RUN make

# Second stage: Final minimal image
FROM python:3.10-slim

# Install runtime dependencies (without build tools)
RUN apt-get update && apt-get install -y \
    libbz2-dev \
    procps \
 && rm -rf /var/lib/apt/lists/*

# Copy the compiled binary from the build stage
COPY --from=build /SURVIVOR/Debug /SURVIVOR/Debug

# Set environment variable for SURVIVOR
ENV PATH="/SURVIVOR/Debug:${PATH}"