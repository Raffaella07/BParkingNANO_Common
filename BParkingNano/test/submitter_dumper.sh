#!/bin/bash

#--------------------
# This script launches the ntuplising tool 
# ${1}:  outdir  -> keep
# ${2}:  usr -> keep
# ${3}:  pl -> keep
# ${4}:  tag -> keep
# ${5}:  isMC --> keep
# ${6}:  isRemote -> remove
# ${7}:  doflat -> remove
# ${8}:  filelist -> remove 
# ${9}:  isResubmission (false if first launch) -> remove
#--------------------


#workdir="/scratch/"${2}"/"${3}"/job_nj"$SLURM_ARRAY_TASK_ID
workdir="/scratch/"${2}"/"${3}"/dumperjob"
echo "creating workdir "$workdir
mkdir -p $workdir

echo "copying ntupliser to workdir"
#cp starter_${3}.C $workdir 
cp nanoTools.py $workdir
cp ../plugins/dumper/NanoDumper.C $workdir 
cp ../plugins/dumper/NanoDumper.h $workdir 
cp ../plugins/dumper/NanoRunDumper.C $workdir 
cp ../plugins/dumper/NanoRunDumper.h $workdir 

cd $workdir

#echo "creating directory for flat ntuples"
#mkdir ${1}/flat


echo "creating the starter"
echo "python nanoTools.py --writestarter --outdir ${1} --tag ${4} --ismc ${5}"
python nanoTools.py --writestarter --outdir ${1} --tag ${4} --ismc ${5} # --myoptions

# for test 
cp starter.C $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

echo "running the ntupliser on top of the nanofile"
#if [ ${5} == 1 ] ; then #isMC

DATE_START_DUMP=`date +%s`
#root -l -q -b "starter_${3}.C+(\"true\")" 
#root -l -q -b "starter.C+(\"true\")" 
root -l -q -b "starter.C+" 
DATE_END_DUMP=`date +%s`

#else #isData

#  DATE_START_DUMP=`date +%s`
#  #root -l -q -b "starter_${3}.C+(\"false\")" 
#  root -l -q -b "starter.C+(\"false\")" 
#  DATE_END_DUMP=`date +%s`
#fi

echo "copying the file"
if [ ${4} == 0 ] ; then
  xrdcp -f flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano.root
else
  xrdcp -f flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano_${4}.root
fi

echo "content of the workdir"
ls -l

echo "clearing the workdir"
rm -r $workdir

cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

runtime_dump=$((DATE_END_DUMP-DATE_START_DUMP))
echo "Wallclock running time: $runtime_dump s"
