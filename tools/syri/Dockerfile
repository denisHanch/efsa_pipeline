
FROM python:3.10-alpine3.23@sha256:c99275d6bc0c37d8e98b388d4c404861fda7dad5ff87e2995fe7b7bb33898aed

ARG SYRI_REF=v1.7.1

RUN echo "http://dl-cdn.alpinelinux.org/alpine/v3.23/main" > /etc/apk/repositories && \
    echo "http://dl-cdn.alpinelinux.org/alpine/v3.23/community" >> /etc/apk/repositories && \apk add --no-cache \
    bash \
    git \
    wget \
    unzip \
    build-base \
    bzip2-dev \
    libxml2-dev \
    zlib-dev \
    libffi-dev \
    glpk \
    libpng-dev \
    openblas \
    procps-ng

RUN pip install --upgrade pip


WORKDIR /syri
RUN git clone --branch ${SYRI_REF} --depth 1 https://github.com/schneebergerlab/syri.git .

# Install dependencies
RUN pip install pandas numpy pysam scipy cython igraph matplotlib psutil

# Install SyRI
RUN pip install .


# Ensure syri is executable
RUN chmod +x /usr/local/bin/syri