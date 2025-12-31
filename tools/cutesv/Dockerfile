# First stage: Build the SURVIVOR binary
FROM python:3.10-slim AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    python3-dev \
    libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

RUN pip wheel cuteSV --no-cache-dir -w /wheels

# Second stage: Final minimal image
FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    libbz2-1.0 \
    procps \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /wheels /wheels

RUN pip install --no-cache-dir /wheels/*