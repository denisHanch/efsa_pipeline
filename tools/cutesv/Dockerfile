# First stage: Build the binary
FROM python:3.10-slim AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    gcc \
    python3-dev \
    libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip setuptools wheel

# cigar's PyPI tarball has a broken ez_setup.py; grab the single module file
RUN curl -fsSL https://raw.githubusercontent.com/brentp/cigar/master/cigar.py \
    -o "$(python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')/cigar.py"

# Build cuteSV wheel without pulling cigar (already installed above)
RUN pip wheel cuteSV==2.1.3 --no-deps --no-cache-dir -w /wheels

# Also grab wheels for cuteSV's other deps (pysam etc.), excluding cigar
RUN pip wheel pysam --no-cache-dir -w /wheels

# Second stage: Final minimal image
FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    libbz2-1.0 \
    procps \
 && rm -rf /var/lib/apt/lists/*

# Install cigar module directly
COPY --from=builder /usr/local/lib/python3.10/site-packages/cigar.py /usr/local/lib/python3.10/site-packages/cigar.py

COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir /wheels/*