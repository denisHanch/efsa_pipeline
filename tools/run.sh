#!/usr/bin/env bash

version="v1.0.0"
for tool in *; do 
    tag="ecomolegmo/$tool:$version"
    latest="ecomolegmo/$tool:latest"
    
    echo "Processing $tool, tag: $tag, latest: $latest"
    cd $tool
    docker build . -t $tag -t $latest
    docker push $tag
    docker push $latest
    cd .. 
done
