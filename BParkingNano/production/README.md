## Crab submission instuctions

There are two main scripts to submit mini -> nano step through BParkingNANO in this directory.

### For unpublished dataset
Either stored at CERN or on local tiers:
```
crab submit nano_on_crab.py
```

or (dryrun usage, better):
```

crab submit --dryrun nano_on_crab.py #runs a job as it would be on crab, but in the tmp dir of your current working machine
  #if command line suggest so, do
crab proceed # which sends jobs on crab batch

```

The script is standard crab submission exploiting xrootd: input files need to be written in a .txt files in xrootd accessible paths (https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookXrootdService). Examples of .txt files are included int the directory.
Site of storage and path to storage directory can be set using the flags `config.Site.storageSite` and `config.Data.outLFNDirBase`


### For CMSDAS published datasets
use:
```

python submit_on_crab.py
```
The datasets to run on can be specified in samples.ylm, together with the lumimask and if the dataset has to be run as MC samples or not. Number of files per job can be also specified.



