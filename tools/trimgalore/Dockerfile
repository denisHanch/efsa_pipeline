# First stage: Install cutadapt
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4 AS builder

ARG CUTADAPT_VERSION=4.9

RUN pip install --upgrade pip \
 && pip wheel cutadapt==${CUTADAPT_VERSION} --no-cache-dir -w /wheels

# Second stage: Final minimal image
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4

ARG TRIMGALORE_VERSION=0.6.10

RUN apt-get update && apt-get install -y --no-install-recommends \
    perl \
    pigz \
    procps \
    curl \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir /wheels/* && rm -rf /wheels

RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/${TRIMGALORE_VERSION}.tar.gz \
    | tar -xzf - -C /opt \
 && ln -s /opt/TrimGalore-${TRIMGALORE_VERSION}/trim_galore /usr/local/bin/trim_galore \
 && chmod +x /usr/local/bin/trim_galore
