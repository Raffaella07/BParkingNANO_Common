#!/bin/bash

#--------------------
# This script launches the nanoAOD production of a given file on slurm
#--------------------

workdir="/scratch/"${3}"/"${4}"/job_nj"${5}
echo "creating workdir "$workdir
mkdir -p $workdir

echo "copying driver to workdir"
cp run_nano_hnl_cfg.py $workdir
cd $workdir

echo "going to run nano step"
DATE_START=`date +%s`
cmsRun run_nano_hnl_cfg.py inputFile=${1} outputFile="bparknano.root"
DATE_END=`date +%s`
echo "finished running nano step"

echo "copying the file"
xrdcp bparknano.root root://t3dcachedb.psi.ch:1094/${2} 

echo "clearing the workdir"
rm -r $workdir

cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

runtime=$((DATE_END-DATE_START))
echo "Wallclock running time: $runtime s"

