FROM python:3.10-slim

# Install required packages

 RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      openssh-client \
      git \
      default-jdk \
      curl \
      docker.io \
    && rm -rf /var/lib/apt/lists/*

# Install the Nextflow launcher
RUN curl -fsSL https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/ \
 && chmod +x /usr/local/bin/nextflow

# COPY .inputrc /root/

WORKDIR /EFSA_workspace
