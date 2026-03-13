# First stage: Build the SURVIVOR binary
FROM python:3.10-slim AS build

RUN apt-get update && apt-get install -y \
    libbz2-dev \
    build-essential \
    git \
    python3-dev \
    procps \
    gcc \
 && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/fritzsedlazeck/SURVIVOR.git
WORKDIR /SURVIVOR/Debug
RUN make

# Second stage: Final minimal image
FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    libbz2-dev \
    procps \
 && rm -rf /var/lib/apt/lists/*

COPY --from=build /SURVIVOR/Debug /SURVIVOR/Debug

ENV PATH="/SURVIVOR/Debug:${PATH}"