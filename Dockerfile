FROM python:3.10-slim-bookworm
# Install necessary packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    docker.io \
    vim less \
    openssh-client \
    git \
    curl \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Copy the Nextflow binary from VM and make it executable
COPY nextflow /usr/local/bin/nextflow
RUN chmod +x /usr/local/bin/nextflow

# Copy shell configuration for better user experience
COPY .devcontainer/.inputrc /root/

# Set working directory
WORKDIR /EFSA_workspace
COPY . /EFSA_workspace/

# Default command to bash for interactive use
CMD ["/bin/bash"]