#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Customize this prefix to match your pipeline’s image namespace
# e.g. your DockerHub username or org: pipeline-validation:dev becomes myuser/pipeline-validation:dev
IMAGE_PREFIX="pipeline"

# The directory where your module Dockerfiles live
DOCKER_DIR="configs/docker"

echo "Building all custom module images from ${DOCKER_DIR}/…"

# Loop over every Dockerfile in configs/docker
for df in "${dockerfiles[@]}"; do
  module="$(basename "${df}" .Dockerfile)"
  image_tag="${IMAGE_PREFIX}-${module}:dev"

  echo "Building image '${image_tag}' from '${df}'…"
  docker build \
    --pull \
    --file "${df}" \
    --tag "${image_tag}" \
    "${DOCKER_DIR}"
  echo "Built ${image_tag}"
done

echo "ll images built!"