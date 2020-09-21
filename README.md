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
scram b
```

### After first installation

```shell
cd CMSSW_10_2_15/src/
cmsenv 
```


### To run on a test file

```shell
cd PhysicsTools/BParkingNano/test/
cmsenv 
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
python nanoLauncher.py --pl <prodLabel> --tag <tag>
```

The tag is optional, and would be appended to the outputfile name.

Once ready, merge the different nano steps by doing

```
python nanoMerger.py --pl <prodLabel> --tag <tag>
```



Note:

To make contributions to the central code, see intructions in https://github.com/CMSBParking/BParkingNANO
