import os
import glob
import ROOT


class NanoTools(object):
  def getPointDirs(self, location):
    return [f for f in glob.glob(location+'/*')]


  def getNanoDirectories(self, location, prodlabel, dataset):
    if dataset == None:
      dirs = [f for f in glob.glob('{loc}/{pl}/*'.format(loc=location, pl=prodlabel))]
    else:
      dirs = [f for f in glob.glob('{loc}/{pl}/{ds}'.format(loc=location, pl=prodlabel, ds=dataset))]
    if len(dirs) == 0:
      raise RuntimeError('No samples with the production label "{pl}" were found in {loc}'.format(pl=prodlabel if dataset==None else dataset+'_'+prodlabel, loc=location))
    return dirs


  def getLogDir(self, file_, prodlabel, isData):
   if isData: # probably to be modified for mc
     label = file_[file_.find('/',file_.find(prodlabel))+1:file_.find('Chunk')-1] 
     chunk = file_[file_.find('Chunk'):file_.find('bparknano')-1]
   return '/work/anlyon/logs/{}/{}/{}'.format(label, prodlabel, chunk) # this will have to be modified


  def getFilesLocation(self, isData):  
    location = '/pnfs/psi.ch/cms/trivcat/store/user/{usr}/BHNLsGen'.format(usr=os.environ["USER"])
    if isData:
      location += '/data'
    return location
  

  def getLogFile(self, logdir, file_):
    return '{}/nanostep_{}.log'.format(logdir, file_[file_.rfind('_nj')+1:file_.rfind('.root')])


  def getLocalMiniAODFiles(self, user, prodlabel, point):
    pointdir = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}/{}/'.format(user, prodlabel, point)
    return [f for f in glob.glob(pointdir+'/step4_nj*.root')]


  def checkLocalFile(self, nanofile, cond=True):
    rootFile = ROOT.TNetXNGFile.Open(nanofile, 'r')
    if not rootFile: return False
    if cond and not rootFile.GetListOfKeys().Contains('Events'): return False
    else: return True


  def checkFileExists(self, file_):
    import os.path
    from os import path
    return path.exists(file_)


  def getNFiles(self, file_):
    return sum(1 for line in open(file_))


  def getNExpectedNanoFiles(self, dir_):
    n_exp = dir_[dir_.rfind('_n')+2:len(dir_)]
    return int(n_exp)


  def getJobId(self, job):
    return int(job[job.find('job')+4:])


  def getJobIdsList(self, jobIds):
    listIds = ''
    for jobId in jobIds:
      listIds += '{}:'.format(jobId)
    return listIds[:len(listIds)-1]


  def getDataLabel(self, dataset):
    idx = dataset.find('-')
    return dataset.replace('/', '_')[1:idx]
  

  def getMCLabel(self, dataset):
    idx = dataset[1:].find('/')
    label = dataset[1:idx+1]
    if 'ext' in self.dataset: label += '_ext'
    return label


  def getStep(self, file_): 
    return file_[file_.rfind('_nj')+3:file_.rfind('.root')]

    

