#!/bin/bash

#--------------------
# This script launches the nanoAOD production of a given file on slurm
# ${1}: outdir
# ${2}: usr
# ${3}: pl
# ${4}: tag
# ${5}: isMC
# ${6}: isRemote
# ${7}: flt
# ${8}: filelist
#--------------------


workdir="/scratch/"${2}"/"${3}"/job_nj"$SLURM_ARRAY_TASK_ID
echo "creating workdir "$workdir
mkdir -p $workdir

echo "copying driver to workdir"
cp run_nano_hnl_cfg.py $workdir

echo "copying the file list to workdir"
cp ${8} $workdir/filelist.txt

# copy the ntupliser 
if [ ${7} == 1 ] ; then
  echo "copying ntupliser to workdir"
  cp simple_tester.py $workdir 
  cp flat_tree_branches.py $workdir
fi

cd $workdir

if [ ${5} == 1 ] ; then #isMC

  if [ ${6} == 0 ] ; then  #private MC
    echo "going to run nano step on "$(sed $SLURM_ARRAY_TASK_ID'!d' filelist.txt)
    DATE_START=`date +%s`
    cmsRun run_nano_hnl_cfg.py inputFile=$(sed $SLURM_ARRAY_TASK_ID'!d' filelist.txt) outputFile="bparknano.root" isMC=True
    DATE_END=`date +%s`
    echo "finished running nano step"
  else
    echo "going to run nano step on "$(sed $SLURM_ARRAY_TASK_ID'!d' filelist.txt)
    DATE_START=`date +%s`
    cmsRun run_nano_hnl_cfg.py inputFiles=$(sed $SLURM_ARRAY_TASK_ID'!d' filelist.txt) outputFile="bparknano.root" isMC=True
    DATE_END=`date +%s`
    echo "finished running nano step"
  fi

  echo "copying the file"
  if [ ${4} == 0 ] ; then
    xrdcp bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_nj$SLURM_ARRAY_TASK_ID.root 
  else
    xrdcp bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_${4}_nj$SLURM_ARRAY_TASK_ID.root 
  fi


  # run the ntupliser on top of the nanofile
  if [ ${7} == 1 ] ; then
    echo "creating directory for flat ntuples"
    mkdir ${1}/flat

    echo "running the ntupliser on top of the nanofile"
    python simple_tester.py  

    echo "copying the file"
    if [ ${4} == 0 ] ; then
      xrdcp flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano_nj$SLURM_ARRAY_TASK_ID.root
    else
      xrdcp flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano_${4}_nj$SLURM_ARRAY_TASK_ID.root
    fi
  fi

else #isData

  echo "going to run nano step on "$(sed $SLURM_ARRAY_TASK_ID'!d' filelist.txt)
  DATE_START=`date +%s`
  cmsRun run_nano_hnl_cfg.py inputFiles=$(sed $SLURM_ARRAY_TASK_ID'!d' filelist.txt) outputFile="bparknano.root" isMC=False
  DATE_END=`date +%s`
  echo "finished running nano step"

  echo "copying the file"
  if [ ${4} == 0 ] ; then
    xrdcp bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_nj$SLURM_ARRAY_TASK_ID.root 
  else
    xrdcp bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_${4}_nj$SLURM_ARRAY_TASK_ID.root 
  fi
fi

echo "content of the workdir"
ls -l

echo "clearing the workdir"
rm -r $workdir

cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

runtime=$((DATE_END-DATE_START))
echo "Wallclock running time: $runtime s"

