# First stage: Build the SURVIVOR binary
FROM python:3.10-slim AS build

ARG SURVIVOR_REF=07404b74d3fe42b20f1362f4d8625f284b9d536c

RUN apt-get update && apt-get install -y \
    libbz2-dev \
    build-essential \
    git \
    python3-dev \
    procps \
    gcc \
 && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/fritzsedlazeck/SURVIVOR.git && \
    cd SURVIVOR && \
    git checkout ${SURVIVOR_REF}
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