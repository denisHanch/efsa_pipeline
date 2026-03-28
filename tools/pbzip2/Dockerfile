FROM docker:29.3.1-alpine3.23@sha256:4d90f1f6c400315c2dba96d3ec93c01e64198395cbba04f79d12adce4f737029

RUN echo "http://dl-cdn.alpinelinux.org/alpine/v3.23/main" > /etc/apk/repositories && \
    echo "http://dl-cdn.alpinelinux.org/alpine/v3.23/community" >> /etc/apk/repositories && \
    apk add --no-cache build-base bzip2-dev && \
    cd /tmp && \
    wget http://ftp.debian.org/debian/pool/main/p/pbzip2/pbzip2_1.1.13.orig.tar.gz -O pbzip2.tar.gz && \
    tar xzf pbzip2.tar.gz && \
    cd pbzip2-1.1.13 && \
    make && \
    make install && \
    cd / && \
    rm -rf /tmp/pbzip2* && \
    apk del build-base bzip2-dev