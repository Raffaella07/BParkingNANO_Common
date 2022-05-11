#!/bin/bash

#--------------------
# This script launches the ntuplising tool 
# ${1}:  outdir  
# ${2}:  usr 
# ${3}:  pl 
# ${4}:  tag 
# ${5}:  isMC
# ${6}:  dosignal
# ${7}:  docontrol
# ${8}:  dohnl
# ${9}:  doTagAndProbe
#--------------------

if [ ${5} == 1 ] ; then #isMC
  tag=${4}
else
  tag="0"
fi

workdir="/scratch/"${2}"/"${3}"/dumperjob_"${SLURM_JOB_ID}
echo "creating workdir "$workdir
mkdir -p $workdir

echo "copying ntupliser to workdir"
cp ./files/starter_${3}.C $workdir/starter.C
cp ../data/json/golden_2018.json $workdir
cp ../plugins/dumper/utils.C $workdir 
if [ ${5} == 1 ] ; then
  cp ../plugins/dumper/NanoRunDumper.C $workdir 
  cp ../plugins/dumper/NanoRunDumper.h $workdir 
fi
if [ ${6} == 1 ] ; then
  cp ../plugins/dumper/BToMuMuPiDumper.C $workdir 
  cp ../plugins/dumper/BToMuMuPiDumper.h $workdir 
fi
if [ ${7} == 1 ] ; then
  cp ../plugins/dumper/BToKMuMuDumper.C $workdir 
  cp ../plugins/dumper/BToKMuMuDumper.h $workdir 
fi
if [ ${8} == 1 ] ; then
  cp ../plugins/dumper/HNLToMuPiDumper.C $workdir 
  cp ../plugins/dumper/HNLToMuPiDumper.h $workdir 
fi
if [ ${9} == 1 ] ; then
  cp ../plugins/dumper/TagAndProbeDumper.C $workdir 
  cp ../plugins/dumper/TagAndProbeDumper.h $workdir 
fi

echo "copying starter"
if [ ${4} == 0 ] ; then
  xrdcp -f ./files/starter_${3}.C root://t3dcachedb.psi.ch:1094/${1}/flat/starter.C
else
  xrdcp -f ./files/starter_${3}.C root://t3dcachedb.psi.ch:1094/${1}/flat/starter_${4}.C
fi
rm ./files/starter_${3}.C

cd $workdir

echo "running the ntupliser on top of the nanofile"
DATE_START_DUMP=`date +%s`
root -l -q -b "starter.C+" 
DATE_END_DUMP=`date +%s`

echo "copying the file"
if [ ${4} == 0 ] ; then
  echo "xrdcp -f flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano.root"
  xrdcp -f flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano.root
else
  echo "xrdcp -f flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano_${4}.root"
  xrdcp -f flat_bparknano.root root://t3dcachedb.psi.ch:1094/${1}/flat/flat_bparknano_${4}.root
fi

echo "content of the workdir"
ls -l

echo "clearing the workdir"
rm -r $workdir

cd $CMSSW_BASE/src/PhysicsTools/BParkingNano/test

runtime_dump=$((DATE_END_DUMP-DATE_START_DUMP))
echo "Wallclock running time: $runtime_dump s"
