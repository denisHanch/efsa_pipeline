# Use an official Python runtime as a parent image
FROM python:3.11-slim

# Set working directory inside the container
WORKDIR /app

# Install required system dependencies (git for cloning)
RUN apt-get update && \
    apt-get install -y --no-install-recommends git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Clone the repository
RUN git clone https://github.com/schneebergerlab/plotsr.git

# Change working directory to plotsr
WORKDIR /app/plotsr

# Install the package
RUN python setup.py install