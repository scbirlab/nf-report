#!/usr/bin/env bash

set -exuo pipefail

MODE=${1:-no}
script_dir="$(readlink -f $(dirname $0))"
start_dir="$(pwd)"

if [ "$MODE" == "gh" ]
then
    export NXF_CONTAINER_ENGINE=docker
    # export NXF_DOCKER_OPTS="-u $(id -u):$(id -g)"
    docker_flag='-profile gh'
elif [ "$MODE" == "crick" ]
then 
    export SINGULARITY_FAKEROOT=1
    export NXF_CONTAINER_ENGINE=singularity
    docker_flag='-profile standard'
    module load Singularity Nextflow
else
    docker_flag=''
fi

# Examples without sample sheet
nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    -work-dir "$script_dir"/work \
    --test \
    --outputs "$script_dir"/output-83332 \
    --organism_id 83332

# Examples with sample sheet
cd "$script_dir"/sheet
nextflow run ../.. \
    -resume $docker_flag \
    -work-dir "$script_dir"/work \
    --test
cd "$start_dir"
