# First stage
FROM python:3.10-slim AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    python3-dev \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip

RUN pip wheel Truvari --no-cache-dir -w /wheels

# Second stage: Final minimal image
FROM python:3.10-slim

RUN apt-get update && apt-get install -y \
    procps \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /wheels /wheels

RUN pip install --no-cache-dir /wheels/*