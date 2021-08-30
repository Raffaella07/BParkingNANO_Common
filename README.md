# NanoAOD producer customized for BParking analysis 

## Installation

Setup the environment
```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```

Add the BParkingNano tool

```
git clone git@github.com:BParkHNLs/BParkingNANO.git ./PhysicsTools
```

Import the BParking modifications on the TransientTracks, KinematicVertexFitter, ElectronRegression and GBRForest
```
git cms-merge-topic -u amlyon:BHNLNano
```

Clone CMSBParking branch for ElectronIdentification
```
git clone --single-branch --branch from-CMSSW_10_2_15_2020Sept15 git@github.com:CMSBParking/RecoEgamma-ElectronIdentification.git $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data
```

To run on CRAB, do
```
git cms-addpkg RecoEgamma/ElectronIdentification
mkdir -p $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data/LowPtElectrons
cp $CMSSW_BASE/external/$SCRAM_ARCH/data/RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/data/LowPtElectrons
```

Build everything

```
git clone https://github.com/Raffaella07/BParkingNANO_Common.git ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b -j 8
```

### After first installation

```shell
cd CMSSW_10_2_15/src/
cmsenv 
```


## Nano samples production

```
cd PhysicsTools/BParkingNano/test/
```

### Locally
Modify the inputFiles in run_nano_hnl_cfg.py, make sure to compile and do

```
cmsRun run_nano_hnl_cfg.py 
```

Runs by default with isMC=True and maxEvents=-1. The outputfile will be saved locally in the current working directory.

### On the batch
The outputfiles will be saved on the Tier3 Storage Element, do not forget to activate your proxy

```
voms-proxy-init --voms cms --valid 186:00
```

The nanoLauncher is a tool that submits a production via slurm arrays of jobs. It allows to launch the nano and/or the ntuplising steps. To run it, do

```
python nanoLauncher.py <options>
```
Options:

* Indicate whether to run on data or mc (central or private)
  * --mcprivate or   
  * --mccentral or
  *  --data       
* --pl `<prodLabel>` 
  * with --mcprivate:  must correspond to the production label of the miniAOD sample (e.g V15_full) 
  * with --mccentral/data: any production label of your choice
* --ds `<dataset>`:  to be used with --data or --mccentral only. Datasets listed in data/samples 
* Indicate which process(es) to run
  * --dosignal: B->mumupi
  * --docontrol: B->Kmumu
  * --dohnl: HNL->mupi
  * --dotageprobe: Jpsi->mumu
* Indicate on which steps to run
  * --donano: launch the nano step
  * --doflat: launch the ntuplising step
  * --domergenano: launch the merging tool automatically after the nano step. Not recommended as doubles the storage space needed
* --user `<user>`: with --mcprivate only; username where the miniAOD samples are stored
* --tagnano `<tagnano>`: optional, tag to be appended to the nano rootfile name 
* --tagflat `<tagflat>`: optional, tag to be appended to the flat rootfile name 
* --maxfiles `<maxfiles>`: optional, maximum number of files to process
* --doquick: optional run on slurm quick partition (time/job < 1h)
* --docompile: optional, compiles the BParkingNano tool before launching

Examples of usage:
```
python nanoLauncher.py --pl V15_full --user mratti --dosignal --donano --doflat --mcprivate
```
```
python nanoLauncher.py --pl V01 --ds /ParkingBPH1/Run2018A-05May2019-v1/MINIAOD --dosignal --donano --doflat --data
```
```
python nanoLauncher.py --pl V01 --ds /QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM --dosignal --donano --doflat --mccentral
```

If not done at launching, you can merge a posteriori the different nano steps by doing

```
python nanoMerger.py --pl <prodLabel> --ds <dataset> --tagnano <tagnano> --donano --<mcprivate/mccentral/data>
```

Note that the production label and tag have to be consistent with those of the nanoAOD production.

Examples of usage:
```
python nanoMerger.py --pl V15_full --donano --mcprivate
```
```
python nanoMerger.py --pl V01 --ds /ParkingBPH1/Run2018A-05May2019-v1/MINIAOD --donano --data
```
```
python nanoMerger.py --pl V01 --ds /QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM --donano --mccentral
```

The nanoProdManager reports on the status of an on-going BParking data production, and serves as a resubmission tool
```
python nanoProdManager.py <options>
```
Options:
* --data 
* --pl `<prodLabel>`
* --ds `<dataset>`: optional, if not set, will loop on all datasets produced under `<prodLabel>`
* --dofullreport: expands the status report with per-chunk information and details on reasons of jobs failure
* --dofetchtime: returns how much time the jobs took to run
* --docheckfile: checking if process branches are in the sample
* --dosignal/--docontrol/--dohnl/--dotageprobe
* --doresubmit: resubmits the failed jobs

Example of usage:
```
python nanoProdManager.py --pl F1 --dofullreport --dosignal --doresubmit --data
```

#### Note

To make contributions to the central code, see intructions in https://github.com/CMSBParking/BParkingNANO


# Tool for Tag&Probe studies

```
cd TagAndProbe/test
```

## Compute efficiencies
Configure tagProbeFitTreeAnalyzer_JPsiMuMu_cfg.py and run it
```
cmsRun tagProbeFitTreeAnalyzer_JPsiMuMu_cfg.py
```

## Plot fits and efficiencies
Adapt the global variables in savePlots.C and run the script
```
root -l -b savePlots.C+
```

## Compute the scale factors
Adapt the global variables in getScaleFactor.C and run the script
```
root -l -b getScaleFactor.C+
```

