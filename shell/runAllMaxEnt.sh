#!/bin/sh
WORKING_DIR="/gpfs01/bethge/home/ffranzen/criticality"
MATLAB_BIN="/usr/local/MATLAB/R2014a/bin/matlab"
STIMULUS="cb2"

for NET_SIZE in $(seq 2 2 12); do
	for IDX_DATA in $(seq 1 5); do
        qsub -S /bin/bash \
        	 -o "$WORKING_DIR/logs/" \
        	 -j oe -l nodes=1,vmem=5G,walltime=3:0:0 \
        	 -N "hJV${STIMULUS}_${NET_SIZE}_${IDX_DATA}" \
        	 -d "${WORKING_DIR}/code" \
        	 -v MATLAB_BIN=${MATLAB_BIN},NET_SIZE=${NET_SIZE},IDX_DATA=${IDX_DATA},STIMULUS=${STIMULUS} \
        	 ${WORKING_DIR}/code/shell/runMaxEnt.pbs  
	done
done
