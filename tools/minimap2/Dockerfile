FROM alpine:3.22
ARG MINIMAP2_VERSION=2.28

RUN echo "http://dl-cdn.alpinelinux.org/alpine/v3.22/main" > /etc/apk/repositories && \
    echo "http://dl-cdn.alpinelinux.org/alpine/v3.22/community" >> /etc/apk/repositories && \
    apk add --no-cache --virtual .minimap2-build-deps \
        build-base \
        zlib-dev \
        curl \
    && cd /tmp \
    && curl -L https://github.com/lh3/minimap2/archive/refs/tags/v${MINIMAP2_VERSION}.tar.gz | tar -xzf - \
    && cd minimap2-${MINIMAP2_VERSION} \
    && make \
    && cp minimap2 /usr/local/bin/ \
    && chmod +x /usr/local/bin/minimap2 \
    && cd / \
    && rm -rf /tmp/minimap2-${MINIMAP2_VERSION} \
    && apk del .minimap2-build-deps
