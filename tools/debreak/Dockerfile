# Use a base image with Python
FROM python:3.9-slim

# Install dependencies (including wget and unzip for example genome files)
RUN apt-get update && apt-get install -y \
    libbz2-dev \
    git \
    build-essential \
    python3-dev \
    gcc \
    libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

# Install requirments
RUN pip install pysam

RUN git clone https://github.com/Maggi-Chen/DeBreak.git

RUN chmod +x /DeBreak/debreak

ENV PATH="/DeBreak:${PATH}"

WORKDIR /DeBreak
