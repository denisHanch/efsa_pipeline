FROM docker:29.3.1-alpine3.23@sha256:4d90f1f6c400315c2dba96d3ec93c01e64198395cbba04f79d12adce4f737029

ARG MOSDEPTH_VERSION=0.3.11

RUN printf '%s\n' \
      'http://dl-cdn.alpinelinux.org/alpine/v3.23/main' \
      'http://dl-cdn.alpinelinux.org/alpine/v3.23/community' \
      > /etc/apk/repositories \
    && apk add --no-cache \
      bash \
      ca-certificates \
      curl \
      gcompat \
      libgcc \
      libstdc++ \
      zlib \
    && update-ca-certificates \
    && curl -kfsSL "https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth" \
         -o /usr/local/bin/mosdepth \
    && chmod +x /usr/local/bin/mosdepth
