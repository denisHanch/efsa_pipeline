# Use a base image with Python
FROM python:3.10-slim

# Install dependencies (including wget and unzip for example genome files)
RUN apt-get update && apt-get install -y \
    git \
    wget \
    unzip \
    libbz2-dev \
    build-essential \
    python3-dev \
    procps \
    gcc \
    libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

# Clone the repo
RUN git clone https://github.com/schneebergerlab/syri.git

# Set working directory inside the container
WORKDIR /syri

# # Install Python dependencies
RUN pip install pandas numpy pysam scipy cython igraph matplotlib psutil

RUN pip install .

# # Make the 'syri' script executable
RUN chmod +x syri


