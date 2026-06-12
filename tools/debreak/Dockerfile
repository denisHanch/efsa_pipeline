# First stage: Build the debreak binary
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4 AS builder

ARG PYSAM_VERSION=0.23.3
ARG DEBREAK_REF=v1.2
ARG MINIMAP2_VERSION=2.30

RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    python3-dev \
    git \
    libbz2-dev \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

RUN pip wheel pysam==${PYSAM_VERSION} --no-cache-dir -w /wheels

WORKDIR /build
RUN git clone --branch ${DEBREAK_REF} --depth 1 https://github.com/Maggi-Chen/DeBreak.git
RUN git clone --branch v${MINIMAP2_VERSION} --depth 1 https://github.com/lh3/minimap2.git \
 && make -C minimap2

# Second stage: Final minimal image
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4

RUN apt-get update && apt-get install -y \
    libbz2-1.0 \
    procps \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir /wheels/*

COPY --from=builder /build/DeBreak /DeBreak
COPY --from=builder /build/minimap2/minimap2 /usr/local/bin/minimap2

RUN chmod +x /DeBreak/debreak

ENV PATH="/DeBreak:${PATH}"

WORKDIR /DeBreak
