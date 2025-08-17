# Use a base image with Python
FROM python:3.10-slim

# Install dependencies (including wget and unzip for example genome files)
RUN apt-get update && apt-get install -y \
    libbz2-dev \
    build-essential \
    python3-dev \
    gcc \
    libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

# Install cuteSV
RUN pip install --upgrade pip
RUN pip install cuteSV