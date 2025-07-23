# Base image
FROM ubuntu:24.04

# Metadata
LABEL description="Docker image with Nextflow and OpenJDK 21 on Ubuntu 22.04"

# Environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64 \
    NEXTFLOW_VERSION=23.10.1

# Install dependencies and OpenJDK 21
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    unzip \
    ca-certificates \
    git \
    pigz \
    && apt-get update && apt-get install -y openjdk-21-jdk \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME explicitly and ensure it's on PATH
ENV PATH="$JAVA_HOME/bin:$PATH"
COPY .inputrc /root/

# Install Nextflow
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow -version

RUN apt-get update && apt-get install -y --no-install-recommends \
    docker-cli sudo \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# RUN usermod -aG sudo $LOCAL_USER

# Set working directory
WORKDIR /efsa_pipeline