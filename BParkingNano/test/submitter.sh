#!/bin/bash

#--------------------
# This script launches the nanoAOD production of a given file on slurm
# ${1}:  outdir
# ${2}:  usr
# ${3}:  pl
# ${4}:  tag
# ${5}:  isMC
# ${6}:  isRemote
# ${7}:  filelist
# ${8}:  isResubmission (false if first launch)
# ${9}:  doTagAndProbe
#--------------------


workdir="/scratch/"${2}"/"${3}"/job_nj"${SLURM_JOB_ID}"_"${SLURM_ARRAY_TASK_ID}
echo "creating workdir "$workdir
mkdir -p $workdir

echo "copying driver to workdir"
if [ ${9} == 0 ] ; then
  cp run_nano_hnl_cfg.py $workdir
else
  cp tag_and_probe_cfg.py $workdir
fi

echo "copying the file list(-s) to workdir"
if [ ${8} == 0 ] ; then
  cp ${7} $workdir/filelist.txt
else # different treatment in case of resubmission
  cp -r ${7}* $workdir
  # copy filelist to pnfs
  #xrdcp -r ${7}*$SLURM_ARRAY_TASK_ID* root://t3dcachedb.psi.ch:1094/${1}
  #rm ${7}*$SLURM_ARRAY_TASK_ID*
  xrdcp -r ${7}"_"$SLURM_ARRAY_TASK_ID".txt" root://t3dcachedb.psi.ch:1094/${1}
  rm ${7}"_"$SLURM_ARRAY_TASK_ID".txt"
fi

cd $workdir

inputFilename=''
if [ ${8} == 0 ] ; then
  inputFilename=$(sed $SLURM_ARRAY_TASK_ID'!d' filelist.txt)
else # different treatment in case of resubmission
  inputFilename=$(sed '1!d' *nj$SLURM_ARRAY_TASK_ID.txt)
fi
echo "inputfilename: "$inputFilename

# index of the output file
if [ ${5} == 1 ] ; then #isMC
  if [ ${6} == 0 ] ; then  #private MC
    search='_nj'
    start=$inputFilename
    end=${start##*$search}
    outIdx=$((${#start} - ${#end}))
    outSuffix=${inputFilename:outIdx}
  else
    outSuffix=$SLURM_ARRAY_TASK_ID".root" 
  fi
else # isData
  outSuffix=$SLURM_ARRAY_TASK_ID".root" 
fi
echo "suffix of the outputfile: "$outSuffix


if [ ${5} == 1 ] ; then #isMC

  if [ ${6} == 0 ] ; then  #private MC
    echo "going to run nano step on "$inputFilename 
    DATE_START=`date +%s`
    if [ ${9} == 0 ] ; then
      cmsRun run_nano_hnl_cfg.py inputFile=$inputFilename outputFile="bparknano.root" isMC=True
    else
      cmsRun tag_and_probe_cfg.py inputFile=$inputFilename outputFile="bparknano.root" isMC=True
    fi
    DATE_END=`date +%s`
    echo "finished running nano step"
  else # central MC
    echo "going to run nano step on "$inputFilename
    DATE_START=`date +%s`
    if [ ${9} == 0 ] ; then
      cmsRun run_nano_hnl_cfg.py inputFiles=$inputFilename outputFile="bparknano.root" isMC=True
    else
      cmsRun tag_and_probe_cfg.py inputFiles=$inputFilename outputFile="bparknano.root" isMC=True
    fi
    DATE_END=`date +%s`
    echo "finished running nano step"
  fi

else #isData

  echo "going to run nano step on "$inputFilename
  DATE_START=`date +%s`
  if [ ${9} == 0 ] ; then
    cmsRun run_nano_hnl_cfg.py inputFiles=$inputFilename outputFile="bparknano.root" isMC=False
  else
    cmsRun tag_and_probe_cfg.py inputFiles=$inputFilename outputFile="bparknano.root" isMC=False
  fi
  DATE_END=`date +%s`
  echo "finished running nano step"

fi

echo "copying the file"
if [ ${4} == 0 ] ; then
  xrdcp -f bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_nj$outSuffix 
else
  echo "xrdcp -f bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_${4}_nj$outSuffix"
  xrdcp -f bparknano.root root://t3dcachedb.psi.ch:1094/${1}/bparknano_${4}_nj$outSuffix 
fi

# copy filelist to pnfs
xrdcp -f filelist.txt root://t3dcachedb.psi.ch:1094/${1}/filelist_${3}.txt

echo "content of the workdir"
ls -l

echo "clearing the workdir"
rm -r $workdir

cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

runtime=$((DATE_END-DATE_START))
echo "Wallclock running time: $runtime s"

