#!/bin/bash

if [[ $# -eq 0 ]]
then
    echo "Usage: run_nplinker_docker.sh <path to nplinker.toml>"
    echo ""
    echo "Note: the image is configured to use the directory containing the"
    echo ".toml file as the location to store all downloaded data if you are"
    echo "using the paired omics platform"
    exit
fi

set -o errexit
set -o nounset

DATADIR=`dirname $1`
echo "NPLinker configuration file: $1"
echo "Local working directory: ${DATADIR}"

docker run --rm -p 5006:5006/tcp -v ${DATADIR}:/data:rw andrewramsay/nplinker:latest
