FROM docker:28.5.0-dind-alpine3.22
RUN echo "http://dl-cdn.alpinelinux.org/alpine/v3.22/main" > /etc/apk/repositories && \
    echo "http://dl-cdn.alpinelinux.org/alpine/v3.22/community" >> /etc/apk/repositories && \
    apk update --no-cache --allow-untrusted

# Install necessary packages including compression tools
RUN apk update && \
    apk add --no-cache \
        bash \
        emacs-nox \
        mc \
        python3 \
        py3-pip \
        py3-pandas \
        py3-openpyxl \
        openssh-client \
        git \
        curl \
        ca-certificates \
        tree \
        openjdk17-jre-headless \
        gzip \
        bzip2 \
        pigz \
    && update-ca-certificates \
    && rm -rf /var/cache/apk/*

# Install build dependencies and build pbzip2 from Debian source as pbzip2 is not available in Alpine
RUN apk add --no-cache build-base bzip2-dev && \
    cd /tmp && \
    wget http://ftp.debian.org/debian/pool/main/p/pbzip2/pbzip2_1.1.13.orig.tar.gz -O pbzip2.tar.gz && \
    tar xzf pbzip2.tar.gz && \
    cd pbzip2-1.1.13 && \
    make && \
    make install && \
    cd / && \
    rm -rf /tmp/pbzip2* && \
    apk del build-base bzip2-dev

# Build and install gffread from source (not available in Alpine repos)
RUN apk add --no-cache --virtual .gffread-build-deps \
        build-base \
    && cd /tmp \
    && git clone https://github.com/gpertea/gffread.git \
    && cd gffread \
    && make release \
    && cp gffread /usr/local/bin/ \
    && chmod +x /usr/local/bin/gffread \
    && cd / \
    && rm -rf /tmp/gffread \
    && apk del .gffread-build-deps

# Copy and install validation package
COPY modules/validation/ /tmp/validation/
# Install build dependencies for Python packages
RUN apk add --no-cache --virtual .build-deps \
        gcc \
        g++ \
        python3-dev \
        musl-dev \
        linux-headers \
    && pip3 install --no-cache-dir --break-system-packages /tmp/validation/validation-pkg/ \
    && rm -rf /tmp/validation/ \
    && apk del .build-deps

# Copy the Nextflow binary from VM and make it executable
COPY nextflow /usr/local/bin/nextflow
RUN chmod +x /usr/local/bin/nextflow

# Copy the validation wrapper script
COPY modules/validation/validation.sh /usr/local/bin/validation.sh
RUN chmod +x /usr/local/bin/validation.sh

WORKDIR /EFSA_workspace

# Copy shell configuration for better user experience
COPY .devcontainer/.inputrc /root/

# Setup bash aliases
RUN echo 'alias validate="validation.sh"' >> /root/.bashrc

# Set bash as default shell
ENV SHELL=/bin/bash