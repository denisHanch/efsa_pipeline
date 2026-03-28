#!/usr/bin/env bash

base=$(pwd)
version="v1.0.3"
for tool in *; do 
    if [[ -f $tool ]]; then
        echo "Skipping $tool"
        continue
    fi

    tag="ecomolegmo/$tool:$version"
    latest="ecomolegmo/$tool:latest"
    
    echo "Processing $tool, tag: $tag, latest: $latest"
    cd "$base/$tool"
    if ! docker build . -t "$tag" -t "$latest"; then
        echo "ERROR: Build failed for $tool, removing image and continuing..."
        docker rmi "$tag" "$latest" 2>/dev/null
        cd "$base"
        continue
    fi
    docker push $tag
    docker push $latest
    cd "$base"
done
