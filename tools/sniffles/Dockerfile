# Builder stage: build Sniffles wheel
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4 AS builder

ARG SNIFFLES_VERSION=2.7.3

RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    python3-dev \
    libbz2-dev \
    procps \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip wheel

RUN pip wheel --no-cache-dir sniffles==${SNIFFLES_VERSION} -w /wheels

# Final stage: minimal runtime image
FROM python:3.12-slim@sha256:3d5ed973e45820f5ba5e46bd065bd88b3a504ff0724d85980dcd05eab361fcf4

RUN apt-get update && apt-get install -y \
    libbz2-dev \
    procps \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip wheel

COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir /wheels/*