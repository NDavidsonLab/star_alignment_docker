#!/usr/bin/env bash

# IMAGE_NAME="mgibio/star:2.7.0f"
IMAGE_NAME="mgibio-star-plus-qcs:2.7.0f"

docker build -t ${IMAGE_NAME} - < ../Dockerfile && \
docker run --rm -it \
    --name star-aligner \
    -v ${PWD}:/app/ \
    -v /data/projects/woodruff_bitler_rnaseq_05062026/Bitler/20260416_LH00407_0240_A23JNVYLT3/Bitler_03272026/:/app/data/ \
    -w /app \
    ${IMAGE_NAME} /bin/bash
