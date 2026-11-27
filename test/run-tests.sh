#!/usr/bin/env bash

set -exuo pipefail

MODE=${1:-no}
script_dir="$(readlink -f $(dirname $0))"
start_dir="$(pwd)"

if [ "$MODE" == "gh" ]
then
    docker_flag='-profile gh'
elif [ "$MODE" == "crick" ]
then
    docker_flag='-profile standard'
    module load Singularity Nextflow Java
else
    docker_flag='-profile local'
fi

# Examples without sample sheet
nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    -work-dir "$script_dir"/work \
    -with-dag "$script_dir"/dag.html\
    --outputs "$script_dir"/output-83332 \
    --organism_id 83332 --test

# Examples with sample sheet
cd "$script_dir"/sheet
nextflow run ../.. \
    -resume $docker_flag \
    -work-dir "$script_dir"/work \
    -with-dag "$script_dir"/dag2.html\
    --test
cd "$start_dir"
