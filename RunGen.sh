#!/bin/bash
set -e

RUN_DIRECTORY=$PWD
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $DIR


START=$(date +%s)
if [ ! -d "OutputFluxes" ]
then
    mkdir OutputFluxes
fi

echo $PWD



if (($# != 0))
then
    for arg in "$@"
    do
        ./HNLFluxGen $RUN_DIRECTORY/$arg
    done
else
    ./HNLFluxGen
fi

#END=$(date +%s)
#DIFF=$(( $END - $START ))
#echo "Execution time = $DIFF s"