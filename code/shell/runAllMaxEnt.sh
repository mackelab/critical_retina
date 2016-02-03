#!/bin/sh
if [ "$#" -ne 1 ]; then
	echo "Usage: $0 stimulus"
	exit 1
fi
STIMULUS="$1"

WORKING_DIR="/gpfs01/bethge/home/ffranzen/criticality"
MATLAB_BIN="/usr/local/MATLAB/R2014a/bin/matlab"

for NET_SIZE in $(seq 2 2 12); do
	for IDX_DATA in $(seq 1 5); do
        qsub -S /bin/bash \
        	 -o "$WORKING_DIR/logs/" \
        	 -j oe \
		 -l nodes=1:ppn=1,vmem=5G,walltime=73:0:0 \
        	 -N "hJV${STIMULUS}_${NET_SIZE}_${IDX_DATA}" \
        	 -d "${WORKING_DIR}/code" \
        	 -v MATLAB_BIN=${MATLAB_BIN},NET_SIZE=${NET_SIZE},IDX_DATA=${IDX_DATA},STIMULUS=${STIMULUS} \
        	 ${WORKING_DIR}/code/shell/runMaxEnt.pbs  
	done
done
