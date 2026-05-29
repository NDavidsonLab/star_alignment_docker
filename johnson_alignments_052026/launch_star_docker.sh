#!/usr/bin/env bash

# IMAGE_NAME="mgibio/star:2.7.0f"
IMAGE_NAME="mgibio-star-plus-qcs:2.7.0f"

docker build -t ${IMAGE_NAME} - < ../Dockerfile && \
docker run --rm -it \
    --name star-aligner \
    -v ${PWD}:/app/ \
    -v /data/projects/johnson_rnaseq_05062026/:/app/data/ \
    -w /app \
    ${IMAGE_NAME} /bin/bash
