#!/usr/bin/env bash

base=$(pwd)
version="v1.0.0"
for tool in *; do 
    if [[ -f $tool ]]; then
        echo "Skipping $tool"
        continue
    fi

    tag="ecomolegmo/$tool:$version"
    latest="ecomolegmo/$tool:latest"
    
    echo "Processing $tool, tag: $tag, latest: $latest"
    cd "$base/$tool"
    docker build . -t $tag -t $latest
    docker push $tag
    docker push $latest
    cd "$base"
done
