# First stage: Build the binary
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4 AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    gcc \
    python3-dev \
    libbz2-dev \
    git \
    procps \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip setuptools wheel

# Install cigar module directly (workaround for broken PyPI tarball)
RUN curl -fsSL https://raw.githubusercontent.com/brentp/cigar/master/cigar.py \
    -o "$(python3 -c 'import sysconfig; print(sysconfig.get_paths()[\"purelib\"])')/cigar.py"

# Clone cuteSV source and patch out Cigar from install_requires
RUN git clone --branch cuteSV-v2.1.3 --depth 1 https://github.com/tjiangHIT/cuteSV.git /tmp/cutesv && \
    sed -i "/'Cigar',/d" /tmp/cutesv/setup.py

# Build wheels for all dependencies and patched cuteSV
RUN pip wheel pysam scipy biopython -w /tmp/wheels && \
    pip wheel /tmp/cutesv -w /tmp/wheels

# Second stage: Final minimal image
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4

RUN apt-get update && apt-get install -y \
    libbz2-1.0 \
    procps \
    curl \
 && rm -rf /var/lib/apt/lists/*

# Install cigar module directly
RUN curl -fsSL https://raw.githubusercontent.com/brentp/cigar/master/cigar.py \
    -o "$(python3 -c 'import sysconfig; print(sysconfig.get_paths()[\"purelib\"])')/cigar.py"

COPY --from=builder /tmp/wheels /wheels
RUN ls -lh /wheels && pip install --no-cache-dir /wheels/*