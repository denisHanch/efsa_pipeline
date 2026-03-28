FROM docker:29.3.1-dind-alpine3.23@sha256:4d90f1f6c400315c2dba96d3ec93c01e64198395cbba04f79d12adce4f737029
ARG GFFREAD_REF=v0.12.7
RUN echo "http://dl-cdn.alpinelinux.org/alpine/v3.23/main" > /etc/apk/repositories && \
    echo "http://dl-cdn.alpinelinux.org/alpine/v3.23/community" >> /etc/apk/repositories && \
    apk update --no-cache --allow-untrusted

# Install necessary packages including compression tools
RUN apk update && \
    apk add --no-cache \
        bash \
        emacs-nox \
        mc \
        python3 \
        py3-pip \
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

 # Copy pre-built pbzip2 binary from dedicated image (built from tools/minimap2/Dockerfile)
 COPY --from=ecomolegmo/pbzip2:v1.1.13@sha256:4a308661071e1ba5e111199067353a8885a964bf93cacb00b5bb149097bacdcb /usr/bin/pbzip2 /usr/bin/pbzip2

# Copy pre-built minimap2 binary from dedicated image (built from tools/minimap2/Dockerfile)
COPY --from=ecomolegmo/minimap2:v2.30@sha256:50d38b713d7d68e105aa3870950492407d82128aa9f3c7c20307632edcab50a5 /usr/local/bin/minimap2 /usr/local/bin/minimap2

# Build and install gffread from source (not available in Alpine repos)
RUN apk add --no-cache --virtual .gffread-build-deps \
        build-base \
    && cd /tmp \
    && git clone --branch ${GFFREAD_REF} --depth 1 https://github.com/gpertea/gffread.git \
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

RUN apk add --no-cache py3-pandas

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