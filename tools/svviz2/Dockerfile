FROM python:3.8-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends git build-essential zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Required by seqlib
RUN pip install Cython

# Install seqlib from GitHub
# seqlib0-0.0.6 (in PIP) is not compatible with the newer GCC
# Similar issue https://github.com/lh3/bwa/issues/314
RUN pip install git+https://github.com/kepbod/seqlib.git

RUN pip install git+https://github.com/nspies/svviz2.git