FROM docker:29.3.1-alpine3.23@sha256:4d90f1f6c400315c2dba96d3ec93c01e64198395cbba04f79d12adce4f737029
ARG GFFREAD_REF=v0.12.7

RUN echo "http://dl-cdn.alpinelinux.org/alpine/v3.23/main" > /etc/apk/repositories && \
    echo "http://dl-cdn.alpinelinux.org/alpine/v3.23/community" >> /etc/apk/repositories && \
    apk add --no-cache --virtual .gffread-build-deps \
        build-base \
        git \
    && cd /tmp \
    && git clone --branch ${GFFREAD_REF} --depth 1 https://github.com/gpertea/gffread.git \
    && cd gffread \
    && make release \
    && cp gffread /usr/local/bin/ \
    && chmod +x /usr/local/bin/gffread \
    && cd / \
    && rm -rf /tmp/gffread \
    && apk del .gffread-build-deps
