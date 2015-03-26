#!/bin/sh

WORKING_DIR="/gpfs01/bethge/home/ffranzen/criticality"
MATLAB_BIN="/usr/local/MATLAB/R2014a/bin/matlab"
STIMULUS="cb2"

NET_SIZE=1
IDX_DATA=1

qsub -q agbethge -S /bin/bash -o "$WORKING_DIR/logs/" -j oe -l vmem=5G,walltime=3:0 -N "hJV${STIMULUS}_${NET_SIZE}_${IDX_DATA}" -d "${WORKING_DIR}/code" -v MATLAB_BIN=${MATLAB_BIN},NET_SIZE=${NET_SIZE},IDX_DATA=${IDX_DATA},STIMULUS=${STIMULUS} ${WORKING_DIR}/code/shell/runMaxEnt.pbs
