#!/bin/bash

#--------------------
# This script launches the nanoAOD production of a given file on slurm
#--------------------

workdir="/scratch/"${3}"/"${4}"/job_nj"${5}
echo "creating workdir "$workdir
mkdir -p $workdir

echo "copying driver to workdir"
cp run_nano_hnl_cfg.py $workdir

# copy the ntupliser 
if [ ${7} == 1 ] ; then
  echo "copying ntupliser to workdir"
  cp simple_tester.py $workdir 
  cp flat_tree_branches.py $workdir
fi

cd $workdir

echo "going to run nano step"
DATE_START=`date +%s`
cmsRun run_nano_hnl_cfg.py inputFile=${1} outputFile="bparknano.root"
DATE_END=`date +%s`
echo "finished running nano step"

echo "copying the file"
if [ ${6} == 0 ] ; then
  xrdcp bparknano.root root://t3dcachedb.psi.ch:1094/${2}/bparknano_nj${5}.root 
else
  xrdcp bparknano.root root://t3dcachedb.psi.ch:1094/${2}/bparknano_${6}_nj${5}.root 
fi


# run the ntupliser on top of the nanofile
if [ ${7} == 1 ] ; then
  echo "creating directory for flat ntuples"
  mkdir ${2}/flat

  echo "running the ntupliser on top of the nanofile"
  python simple_tester.py  

  echo "copying the file"
  if [ ${6} == 0 ] ; then
    xrdcp flat_bparknano.root root://t3dcachedb.psi.ch:1094/${2}/flat/flat_bparknano_nj${5}.root
  else
    xrdcp flat_bparknano.root root://t3dcachedb.psi.ch:1094/${2}/flat/flat_bparknano_${6}_nj${5}.root
  fi
fi

echo "clearing the workdir"
rm -r $workdir

cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

runtime=$((DATE_END-DATE_START))
echo "Wallclock running time: $runtime s"

