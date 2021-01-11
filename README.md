# nanoAOD producer customized for BParking analysis 

## Getting started

```shell
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

### Add energy regression and July20-depth13-700trees model for LPT electron ID

```shell
scp <username>@lxplus.cern.ch:/afs/cern.ch/user/c/crovelli/public/4BParking/sparse-checkout .git/info/sparse-checkout
git remote add crovelli git@github.com:crovelli/cmssw.git
git fetch crovelli
git checkout -b from-CMSSW_10_2_15__ID-2020Jul26-depth13-700__WithFinalReg crovelli/from-CMSSW_10_2_15__ID-2020Jul26-depth13-700__WithFinalReg
```

### Add the modification needed to use post-fit quantities for electrons  

```shell
git cms-addpkg TrackingTools/TransientTrack
git cms-merge-topic -u CMSBParking:GsfTransientTracks
```

### Add the modification needed to use the KinematicParticleVertexFitter  

```shell
git cms-merge-topic -u CMSBParking:fixKinParticleVtxFitter
```

### Add the BParkingNano package and build everything

```shell
git clone git@github.com:BParkHNLs/BParkingNANO.git ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b -j 8
```

### After first installation

```shell
cd CMSSW_10_2_15/src/
cmsenv 
```


### To run on a test file

```shell
cd PhysicsTools/BParkingNano/test/
cmsRun run_nano_cfg.py
```

## Nano samples production

```shell
cd PhysicsTools/BParkingNano/test/
```

### Locally
Modify the inputFiles in run_nano_hnl_cfg.py, make sure to compile and do

```
cmsRun run_nano_hnl_cfg.py 
```

Runs by default on MC and over -1 events. The outputfile will be saved locally in the current working directory.

### On the batch
Runs on a slurm-based engine. 

Since the outputfiles will be saved on the Storage Element, do not forget to activate your proxy

```
voms-proxy-init --voms cms --valid 186:00
```

Then do

```
python nanoLauncher.py <options>
```
The options are:
```
*   --mcprivate          : runs the tool on signal samples
OR  --mccentral          : runs the tool on central QCD MC samples
OR  --data               : runs the tool on BParking data
* --pl <prodLabel>       : 
  * with --mcprivate     : must correspond to the production label under which 
                           the miniAOD sample was produced (e.g V11_inclB_n4200000_njt200) 
  * with --mccentral/data: any production label of your choice
* --ds <dataset>         : to be used with --data or --mccentral only 
(e.g /ParkingBPH4/Run2018B-05May2019-v2/MINIAOD, /QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM)
* --user <user>          : with --mcprivate only, optional; username that produced the miniAOD sample
* --tag <tag>            : optional, tag to be appended to the rootfile name 
* --maxfiles <maxfiles>  : optional, maximum number of files to process
```

Once ready, merge the different nano steps by doing

```
python nanoMerger.py --pl <prodLabel> --tag <tag> --<mcprivate/mccentral/data>
```
Note that the production label and tag have to be consistent with those of the nanoAOD production. For data, the production label is typically the name of the corresponding pnfs dir.



Note:

To make contributions to the central code, see intructions in https://github.com/CMSBParking/BParkingNANO

