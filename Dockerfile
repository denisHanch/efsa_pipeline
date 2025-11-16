FROM docker:28.5.0-dind-alpine3.22

RUN echo "http://dl-cdn.alpinelinux.org/alpine/v3.22/main" > /etc/apk/repositories && \
    echo "http://dl-cdn.alpinelinux.org/alpine/v3.22/community" >> /etc/apk/repositories && \
    apk update --no-cache --allow-untrusted


# Install necessary packages
RUN apk update && \
    apk add --no-cache \
    procps \
    bash \
    emacs-nox \
    mc \
    python3 \
    py3-pip \
    openssh-client \
    git \
    curl \
    tree \
    openjdk17-jre-headless \
    && rm -rf /var/cache/apk/*


    
# Copy the Nextflow binary from VM and make it executable
COPY nextflow /usr/local/bin/nextflow
RUN chmod +x /usr/local/bin/nextflow

WORKDIR /EFSA_workspace

# Copy shell configuration for better user experience
COPY .devcontainer/.inputrc /root/