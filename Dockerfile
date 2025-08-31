FROM python:3.10-slim-bookworm

# Install necessary packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    openssh-client \
    git \
    curl \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Install the Nextflow launcher
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow

# Copy shell configuration for better user experience
COPY .devcontainer/.inputrc /root/

# Set working directory
WORKDIR /EFSA_workspace

COPY . /EFSA_workspace/

# Default command to bash for interactive use
CMD ["/bin/bash"]