FROM docker:29.3.1-dind-alpine3.23@sha256:4d90f1f6c400315c2dba96d3ec93c01e64198395cbba04f79d12adce4f737029

RUN printf '%s\n' \
      'http://dl-cdn.alpinelinux.org/alpine/v3.23/main' \
      'http://dl-cdn.alpinelinux.org/alpine/v3.23/community' \
      > /etc/apk/repositories

# Install runtime packages in the main image, with explicit version pinning.
RUN apk add --no-cache \
      bash=5.3.3-r1 \
      emacs-nox=30.2-r0 \
      mc=4.8.33-r2 \
      python3=3.12.12-r0 \
      py3-pip=25.1.1-r1 \
      openssh-client-default=10.2_p1-r0 \
      git=2.52.0-r0 \
      curl=8.17.0-r1 \
      ca-certificates=20251003-r0 \
      tree=2.2.1-r0 \
      openjdk17-jre-headless=17.0.18_p8-r0 \
      gzip=1.14-r2 \
      bzip2=1.0.8-r6 \
      pigz=2.8-r1 \
      py3-pandas=2.3.3-r0 \
    && update-ca-certificates


# Copy the Nextflow binary from VM and make it executable
COPY nextflow /usr/local/bin/nextflow
RUN chmod +x /usr/local/bin/nextflow

WORKDIR /EFSA_workspace

# Copy shell configuration for better user experience
COPY .devcontainer/.inputrc /root/

ENV SHELL=/bin/bash