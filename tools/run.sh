#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
version="v1.0.4"

for tool_dir in "$script_dir"/*; do
    [[ -d "$tool_dir" ]] || continue
    tool="$(basename "$tool_dir")"

    # Build only tool directories that have a Dockerfile
    if [[ ! -f "$tool_dir/Dockerfile" ]]; then
        echo "Skipping $tool (no Dockerfile)"
        continue
    fi

    tag="ecomolegmo/$tool:$version"
    latest="ecomolegmo/$tool:latest"

    echo "Processing $tool, tag: $tag, latest: $latest"
    if ! docker build -f "$tool_dir/Dockerfile" "$repo_root" -t "$tag" -t "$latest"; then
        echo "ERROR: Build failed for $tool, removing image and continuing..."
        docker rmi "$tag" "$latest" 2>/dev/null || true
        continue
    fi

    docker push "$tag"
    docker push "$latest"
done
