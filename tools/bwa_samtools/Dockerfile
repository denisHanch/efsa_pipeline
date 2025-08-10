FROM ubuntu:24.04

# Install basic tools
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    autoconf \
    automake \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    xz-utils \
    pkg-config \
    git \
    ca-certificates && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install BWA
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    cp bwa /usr/local/bin/ && \
    cd .. && rm -rf bwa

# Install Samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar -xjf samtools-1.17.tar.bz2 && cd samtools-1.17 && \
    ./configure --prefix=/usr/local && make && make install && \
    cd .. && rm -rf samtools-1.17*

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
