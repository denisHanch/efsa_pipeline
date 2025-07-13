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

# # WIP DIND
# ENV DOCKER_TLS_CERTDIR=""

# # Install Docker (Docker-ce)
# RUN curl -fsSL https://get.docker.com | sh \
#     && apt-get update && apt-get install -y --no-install-recommends \
#     docker-ce docker-ce-cli containerd.io docker-compose-plugin \
#     && apt-get clean && rm -rf /var/lib/apt/lists/*

# # Set up Docker Daemon (DinD)
# # Use dockerd directly for the Docker daemon
# CMD ["dockerd", "--host=unix:///var/run/docker.sock", "--host=tcp://0.0.0.0:2375"]

# Define build arguments for UID and GID
# ARG LOCAL_UID=501
# ARG LOCAL_USER=esimonova
# ARG LOCAL_GID=20
# ARG LOCAL_GROUP=staff
# ARG CUSTOM_HOMEDIR=/Users/esimonova

# # Create group and user with the same UID/GID as local
# RUN useradd -m -u $LOCAL_UID -g $LOCAL_GID -d $CUSTOM_HOMEDIR -s /bin/bash $LOCAL_USER
# RUN usermod -aG root $LOCAL_USER

RUN apt-get update && apt-get install -y --no-install-recommends \
    docker-cli sudo \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# RUN usermod -aG sudo $LOCAL_USER

# Set working directory
WORKDIR /efsa_pipeline