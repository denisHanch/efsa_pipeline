FROM docker:29.3.1-dind-alpine3.23@sha256:4d90f1f6c400315c2dba96d3ec93c01e64198395cbba04f79d12adce4f737029
ARG GFFREAD_REF=v0.12.7

RUN printf '%s\n' \
      'http://dl-cdn.alpinelinux.org/alpine/v3.23/main' \
      'http://dl-cdn.alpinelinux.org/alpine/v3.23/community' \
      > /etc/apk/repositories

# Install runtime packages in the main image, with explicit version pinning.
# Note: these versions should be kept aligned with Alpine v3.23 for your target arch.
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

# Copy pre-built pbzip2 binary from dedicated image
COPY --from=ecomolegmo/pbzip2:v1.1.13@sha256:4a308661071e1ba5e111199067353a8885a964bf93cacb00b5bb149097bacdcb /usr/bin/pbzip2 /usr/bin/pbzip2

# Copy pre-built minimap2 binary from dedicated image
COPY --from=ecomolegmo/minimap2:v2.30@sha256:50d38b713d7d68e105aa3870950492407d82128aa9f3c7c20307632edcab50a5 /usr/local/bin/minimap2 /usr/local/bin/minimap2

# Copy pre-built gffread binary from dedicated image (built from tools/gffread/Dockerfile)
COPY --from=ecomolegmo/gffread:v0.12.7@sha256:dad98757a1b8dcfae49f452678d59599bfb81d10c1a5272e418a352d1521b9aa /usr/local/bin/gffread /usr/local/bin/gffread

# Copy validation package source
COPY modules/validation/ /tmp/validation/

# Build deps only while creating an isolated venv for the validation package.
RUN apk add --no-cache --virtual .build-deps \
      gcc \
      g++ \
      python3-dev \
      musl-dev \
      linux-headers \
    && python3 -m venv /opt/validation-venv \
    && /opt/validation-venv/bin/pip install --upgrade pip setuptools wheel \
    && /opt/validation-venv/bin/pip install --no-cache-dir /tmp/validation/validation-pkg/ \
    && rm -rf /tmp/validation \
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

ENV SHELL=/bin/bash