# Use a base image with Python
FROM python:3.10-slim

# Install dependencies (including wget and unzip for example genome files)
RUN apt-get update && apt-get install -y \
    build-essential \
    procps \
    python3-dev \
 && rm -rf /var/lib/apt/lists/*

# Install requirments
RUN pip install Truvari