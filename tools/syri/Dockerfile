# First stage
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4 AS builder

ARG SYRI_REF=v1.7.1

RUN apt-get update && apt-get install -y \
    git \
    wget \
    unzip \
    build-essential \
    gcc \
    python3-dev \
    libc-dev \
    libbz2-dev \
    libxml2-dev \
    libz-dev \
    libffi-dev \
    libglpk-dev \
    libpng-dev \
    libblas-dev \
    liblapack-dev \
    pkg-config \
    procps \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

WORKDIR /build
RUN git clone --branch ${SYRI_REF} --depth 1 https://github.com/schneebergerlab/syri.git


WORKDIR /build/syri

RUN pip install --upgrade pip && \
    pip install cython numpy

RUN pip install \
    pandas \
    "pysam>=0.22.0" \
    "scipy>=1.11.0" \
    "cython>=3.0.0" \
    "python-igraph>=0.11.3" \
    matplotlib \
    psutil



# Second stage: Final minimal image
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4

RUN apt-get update && apt-get install -y \
    libbz2-dev \
    libxml2 \
    libz1 \
    libffi8 \
    libglpk40 \
    libpng16-16 \
    libblas3 \
    liblapack3 \
    procps \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/syri/syri /usr/local/bin/syri

RUN pip install --no-cache-dir \
    pandas \
    "pysam>=0.22.0" \
    "scipy>=1.11.0" \
    "cython>=3.0.0" \
    "python-igraph>=0.11.3" \
    matplotlib \
    psutil

RUN chmod +x /usr/local/bin/syri