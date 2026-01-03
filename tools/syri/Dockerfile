# First stage
FROM python:3.10-slim AS builder

RUN apt-get update && apt-get install -y \
    git \
    wget \
    unzip \
    build-essential \
    gcc \
    python3-dev \
    libbz2-dev \
    procps \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

WORKDIR /build
RUN git clone https://github.com/schneebergerlab/syri.git

WORKDIR /build/syri

RUN pip install \
    pandas \
    numpy \
    pysam \
    scipy \
    cython \
    igraph \
    matplotlib \
    psutil

RUN pip wheel . --no-deps -w /wheels

# Second stage: Final minimal image
FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    libbz2-dev \
    procps \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /wheels /wheels

RUN pip install --no-cache-dir /wheels/*

RUN chmod +x /usr/local/bin/syri