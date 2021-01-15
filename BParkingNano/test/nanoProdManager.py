import sys 
import os
import glob
import subprocess


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to report the status of a nanoAOD production and resubmit failed files', add_help=True)
  parser.add_argument('--pl'        , type=str, dest='pl'          , help='label of the sample file'                                        , default=None)
  parser.add_argument('--ds'        , type=str, dest='ds'          , help='[optional] name of dataset'                                      , default=None)
  parser.add_argument('--mcprivate'           , dest='mcprivate'   , help='run the resubmitter on a private MC sample' , action='store_true', default=False)
  parser.add_argument('--mccentral'           , dest='mccentral'   , help='run the resubmitter on a central MC sample' , action='store_true', default=False)
  parser.add_argument('--data'                , dest='data'        , help='run the resubmitter on a data sample'       , action='store_true', default=False)
  parser.add_argument('--dofullreport'        , dest='dofullreport', help='add to report chunks and failure reason'    , action='store_true', default=False)
  parser.add_argument('--dofetchtime'         , dest='dofetchtime' , help='add to report time fetch'                   , action='store_true', default=False)
  parser.add_argument('--doresubmit'          , dest='doresubmit'  , help='resubmit failed jobs'                       , action='store_true', default=False)
  return parser.parse_args()


def checkParser(opt):

  if opt.pl==None:
    raise RuntimeError('Please indicate the production label of the sample')

  if opt.mcprivate==False and opt.mccentral==False and opt.data==False:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding either --data or--mcprivate or --mccentral to the command line')

  if opt.mcprivate + opt.mccentral + opt.data > 1:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding only --data or --mcprivate or --mccentral to the command line')



class NanoProdManager(object):
  def __init__(self, opt):
    self.prodlabel    = vars(opt)['pl']
    self.dataset      = vars(opt)['ds']
    self.mcprivate    = vars(opt)['mcprivate']
    self.mccentral    = vars(opt)['mccentral']
    self.data         = vars(opt)['data']
    self.dofullreport = vars(opt)['dofullreport'] 
    self.dofetchtime  = vars(opt)['dofetchtime'] 
    self.doresubmit   = vars(opt)['doresubmit'] 


  def getFilesLocation(self):  
    location = '/pnfs/psi.ch/cms/trivcat/store/user/{usr}/BHNLsGen'.format(usr=os.environ["USER"])
    if self.data:
      location += '/data'
    return location


  def getDirectories(self, location):
    if self.dataset == None:
      dirs = [f for f in glob.glob('{loc}/*{pl}'.format(loc=location, pl=self.prodlabel))]
    else:
      dirs = [f for f in glob.glob('{loc}/{ds}_{pl}'.format(loc=location, ds=self.dataset, pl=self.prodlabel))]
    if len(dirs) == 0:
      raise RuntimeError('No samples with the production label "{pl}" were found in {loc}'.format(pl=self.prodlabel if self.dataset==None else self.dataset+'_'+self.prodlabel, loc=location))
    return dirs


  def getNExpectedFiles(self, dir_):
    n_exp = dir_[dir_.rfind('_n')+2:len(dir_)]
    return int(n_exp)


  def getNOutputFiles(self, dir_): # not needed anymore?
    command = 'ls -al {}/bparknano_nj*.root | wc -l'.format(dir_)
    n_out = subprocess.check_output(command, shell=True)
    return int(n_out)


  def checkFileExists(self, file_):
    import os.path
    from os import path
    return path.exists(file_)


  def getStep(self, file_): 
    return file_[file_.rfind('_nj')+3:file_.rfind('.root')]


  def getLogDir(self, file_):
   if self.data: # probably to be modified for mc
     label = file_[file_.find('/',file_.find('data'))+1:file_.find(self.prodlabel)-1] 
     chunk = file_[file_.find('Chunk'):file_.find('bparknano')-1]
   return '/work/anlyon/logs/{}/{}/{}'.format(label, self.prodlabel, chunk) # this will have to be modified


  def getLogFile(self, logdir, file_):
    return '{}/nanostep_{}.log'.format(logdir, file_[file_.rfind('_nj')+1:file_.rfind('.root')])


  def scanLogFile(self, logfile, key):
    with open(logfile) as f:
      if key in f.read():
        return True
      else:
        return False
   

  def isJobFinished(self, logfile):
    return self.scanLogFile(logfile, 'finished running nano step') or self.scanLogFile(logfile, 'slurmstepd: error')
      

  def fetchTime(self, logfile, key):
    if key == 'wallclock':
      command = 'grep  "Wallclock running time:" {}'.format(logfile)
      try:
        output = subprocess.check_output(command, shell=True)
        time = float(output[output.find('time:')+6:len(output)-2])
      except:
        time = 0
    if key == 'cpu':
      command = 'grep  "event loop Real/event" {}'.format(logfile)
      try:
        output = subprocess.check_output(command, shell=True)
        time = float(output[output.find('event =')+8:len(output)])
      except:
        time = 0
    return time


  def checkFailureReason(self, logfile):
    if self.scanLogFile(logfile, "An exception of category 'FallbackFileOpenError'"):
      error_label = 'xrootd'

    elif self.scanLogFile(logfile, "An exception of category 'FileReadError'"):
      error_label = 'readerror'
    
    elif self.scanLogFile(logfile, 'DUE TO TIME LIMIT ***'):
      error_label = 'slurm_timeout'

    elif self.scanLogFile(logfile, 'oom-kill event'):
      error_label = 'slurm_memout'
      
    elif self.scanLogFile(logfile, 'DUE TO NODE FAILURE'):
      error_label = 'slurm_nodefailure'

    else:
      error_label = 'other'

    return error_label


  def writeFileList(self, failed_files):

    logdir = self.getLogDir(failed_files[0])
    label = logdir[logdir.find('logs')+5:].replace('/', '_')

    if self.data:
      filename = './files/resubmit_data_{}'.format(label)
    elif self.mcprivate:
      filename = './files/resubmit_mcprivate_{}'.format(label)
    elif self.mccentral:
      filename = './files/resubmit_mccentral_{}'.format(label)

    for file_ in failed_files:
      # open file list
      filelist = open(filename + '_nj{}.txt'.format(self.getStep(file_)), 'w+')

      # get the file to reprocess
      logfile = self.getLogFile(logdir, file_)
      command = 'grep "going to run nano step on" {}'.format(logfile)
      output = subprocess.check_output(command, shell=True)
      file_toresubmit = output[output.find('step on')+8:len(output)]

      # write to file list
      filelist.write(file_toresubmit + '\n')
      
      # close file list
      filelist.close()
  
      #print 'created {}_nj{}.txt'.format(filename, self.getStep(file_))

    return filename 

    # this is to avoid slurm to deal with too large arrays --> to remove if per-chunk resubmission
    #if len(failed_files) > 1000:
    #  command_split = 'split -l 1000 {fn}.txt {fn} --additional-suffix=.txt'.format(fn=filename)
    #  os.system(command_split)
    #  os.system('rm {fn}.txt'.format(fn=filename))


  def getArray(self, failed_files): 
    idx = []
    for file_ in failed_files:
      idx.append(self.getStep(file_))
    return ','.join(idx)
    

  def resubmit(self, failed_files):
    # strategy: per chunk resubmission
    #           submit job arrays with indices corresponding to the stepId of the failed jobs

    logdir    = self.getLogDir(failed_files[0]) 
    label     = logdir[logdir.find('logs')+5:].replace('/', '_')
    array     = self.getArray(failed_files)
    outputdir = failed_files[0][0:failed_files[0].find('bparknano')]
    filelist  = self.writeFileList(failed_files) 
    

    command = 'sbatch -p wn --account=t3 -o {ld}/nanostep_nj%a.log -e {ld}/nanostep_nj%a.log --job-name=nanostep_nj%a_{pl} --array {ar} --time=03:00:00 submitter.sh {outdir} {usr} {pl} {tag} {isMC} {rmt} {flt} {lst} 1'.format(
      ld      = logdir,
      pl      = label,
      ar      = self.getArray(failed_files), 
      outdir  = outputdir,
      usr     = os.environ["USER"], 
      tag     = 0, #if self.tag == None else self.tag, # adapt
      isMC    = 1 if self.mcprivate or self.mccentral else 0,
      rmt     = 0 if self.mcprivate else 1,
      flt     = 0, #1 if self.doflat == True else 0, # adapt
      lst     =  filelist,
    )

    #print command
    os.system(command)


  def process(self):
    print '---------------------------------'
    print '         Nano Job Manager        '
    print '---------------------------------'

    if self.data and self.dataset==None:
      print '\nINFO: the JobManager Tool is going to be run on all dataset having "{}" as production label\n'.format(opt.pl)

    # pnfs directory where nano samples are located
    location = self.getFilesLocation()

    # get the directories associated to the production label
    dirs = self.getDirectories(location)
    #dirs = [f for f in glob.glob('{loc}/ParkingBPH1_Run2018B*{pl}'.format(loc=location, pl=self.prodlabel))]
    if len(dirs) == 0:
      raise RuntimeError('No samples with the production label "{pl}" were found in {loc}'.format(pl=self.prodlabel, loc=location))

    n_good        = 0
    n_failed      = 0
    n_unfinished  = 0
    n_unprocessed = 0

    n_failure_xrootd  = 0
    n_failure_readerr = 0
    n_failure_timeout = 0
    n_failure_memout  = 0
    n_failure_node    = 0
    n_failure_other   = 0

    time_wallclock = 0
    time_cpu = 0

    for dir_ in dirs:
      print '\n\n {}'.format(dir_)

      # define counters
      n_good_perdir        = 0
      n_failed_perdir      = 0
      n_unfinished_perdir  = 0
      n_unprocessed_perdir = 0
      #n_tot_perdir         = 0
      #n_unstarted_chunk    = 0

      n_failure_xrootd_perdir  = 0
      n_failure_readerr_perdir = 0
      n_failure_timeout_perdir = 0
      n_failure_memout_perdir  = 0
      n_failure_node_perdir    = 0
      n_failure_other_perdir   = 0

      time_wallclock_perdir = 0
      time_cpu_perdir       = 0


      # get different chunks
      chunks = [f for f in glob.glob('{}/Chunk*'.format(dir_))]

      for chunk_ in chunks:
        if self.dofullreport: print '\n -- {} --'.format(chunk_[chunk_.rfind('/')+1:len(chunk_)])

        failed_files = []

        n_good_perchunk        = 0
        n_failed_perchunk      = 0
        n_unfinished_perchunk  = 0
        n_unprocessed_perchunk = 0

        time_wallclock_perchunk = 0
        time_cpu_perchunk       = 0

        n_exp = self.getNExpectedFiles(chunk_)
        #n_out = self.getNOutputFiles(chunk_)

        #n_tot_perdir += n_exp

        #if n_out == 0: # production of this chunk of samples has not started yet
        #  print '-> samples in this chunk have not started to be processed yet'
        #  n_unstarted_chunk += 1
          
        #else:  
        #if n_exp == n_out: # no failed jobs
        #  n_good_perchunk = n_exp

        #else: # there are some failed or unfinished jobs
        files = [chunk_+'/bparknano_nj'+str(nj)+'.root' for nj in range(1, n_exp+1)]

        for file_ in files:
          # get the log file
          logdir = self.getLogDir(file_)
          logfile = self.getLogFile(logdir, file_)
          
          # idle jobs
          if not self.checkFileExists(logfile): 
            n_unprocessed_perchunk += 1
            continue

          if self.checkFileExists(file_): # successfull job
            n_good_perchunk += 1

            if self.dofetchtime:  
              time_wallclock_perchunk += self.fetchTime(logfile, 'wallclock')
              time_cpu_perchunk += self.fetchTime(logfile, 'cpu')

          else: # failed or still running job
           
            if not self.isJobFinished(logfile):
               n_unfinished_perchunk += 1
            else:
               n_failed_perchunk += 1

               if self.dofullreport:
                 try:
                   failure_reason = self.checkFailureReason(logfile)
                 except:
                  pass
                 if failure_reason == 'xrootd': n_failure_xrootd_perdir  += 1
                 if failure_reason == 'readerror': n_failure_readerr_perdir  += 1
                 if failure_reason == 'slurm_timeout': n_failure_timeout_perdir  += 1
                 if failure_reason == 'slurm_memout': n_failure_memout_perdir  += 1
                 if failure_reason == 'slurm_nodefailure': n_failure_node_perdir  += 1
                 if failure_reason == 'other': n_failure_other_perdir  += 1
     
                 #if failure_reason == 'other':
                 # print '{} does not exist'.format(file_)

                 #print '{} does not exist'.format(file_)
                 failed_files.append(file_)
                 #failed_logdir.append(self.getLogDir(file_))
  
       
        if self.dofullreport:
          print 'number of successfull jobs in this chunk: {}'.format(n_good_perchunk)
          print 'number of failed jobs in this chunk     : {}'.format(n_failed_perchunk)
          if n_unfinished_perchunk != 0:
            print 'number of running jobs in this chunk    : {}'.format(n_unfinished_perchunk)
          if n_unprocessed_perchunk != 0:
            print 'number of idle jobs in this chunk       : {}'.format(n_unprocessed_perchunk)
          if self.dofetchtime:
            print 'average wallclock time in this chunk    : {}min '.format(round(time_wallclock_perchunk/(n_good_perchunk*60), 1) if n_good_perchunk != 0 else 0)
            print 'average CPU time/event in this chunk    : {}s/event '.format(round(time_cpu_perchunk/n_good_perchunk, 3) if n_good_perchunk != 0 else 0)

        
        # resubmission
        if self.doresubmit and len(failed_files) != 0:
          print ' --> resubmission of failed files ({})'.format(self.getArray(failed_files))
          self.resubmit(failed_files)


        n_good_perdir        += n_good_perchunk
        n_failed_perdir      += n_failed_perchunk
        n_unfinished_perdir  += n_unfinished_perchunk
        n_unprocessed_perdir += n_unprocessed_perchunk

        time_wallclock_perdir += time_wallclock_perchunk
        time_cpu_perdir       += time_cpu_perchunk


      n_tot_perdir = n_good_perdir + n_failed_perdir + n_unfinished_perdir + n_unprocessed_perdir

      print '\n'
      print ' --> number of successfull jobs      : {}    {}%'.format(n_good_perdir, round(n_good_perdir/float(n_tot_perdir)*100, 2))
      print ' --> number of failed jobs           : {}    {}%'.format(n_failed_perdir, round(n_failed_perdir/float(n_tot_perdir)*100, 2))
      if self.dofullreport:
        print '      - xrootd error:        {}'.format(n_failure_xrootd_perdir)
        print '      - file read error:     {}'.format(n_failure_readerr_perdir)
        print '      - slurm timeout error: {}'.format(n_failure_timeout_perdir)
        print '      - slurm memout error:  {}'.format(n_failure_memout_perdir)
        print '      - slurm node failure:  {}'.format(n_failure_node_perdir)
        print '      - other:               {}'.format(n_failure_other_perdir)
      print ' --> number of running jobs          : {}    {}%'.format(n_unfinished_perdir, round(n_unfinished_perdir/float(n_tot_perdir)*100, 2))
      print ' --> number of idle jobs             : {}    {}%'.format(n_unprocessed_perdir, round(n_unprocessed_perdir/float(n_tot_perdir)*100, 2))
      if self.dofetchtime:
        print ' --> average wallclock time          : {}min '.format(round(time_wallclock_perdir/(n_good_perdir*60), 1) if n_good_perdir != 0 else 0)
        print ' --> average CPU time/event          : {}s/event '.format(round(time_cpu_perdir/(n_good_perdir), 3) if n_good_perdir != 0 else 0)
      #if n_unstarted_chunk != 0:
      #  print ' --> number of yet unprocessed chunks: {}'.format(n_unstarted_chunk)

      n_good        += n_good_perdir
      n_failed      += n_failed_perdir
      n_unfinished  += n_unfinished_perdir
      n_unprocessed += n_unprocessed_perdir

      n_failure_xrootd  += n_failure_xrootd_perdir
      n_failure_readerr += n_failure_readerr_perdir
      n_failure_timeout += n_failure_timeout_perdir
      n_failure_memout  += n_failure_memout_perdir
      n_failure_node    += n_failure_node_perdir
      n_failure_other   += n_failure_other_perdir

      time_wallclock += time_wallclock_perdir
      time_cpu       += time_cpu_perdir

    n_tot = n_good + n_failed + n_unfinished + n_unprocessed

    print '\n'
    print '----------------------------------------------------------'
    print '                      Status Report                       '

    print '\n'
    print ' ---> number of successfull jobs           : {}    {}%'.format(n_good, round(n_good/float(n_tot)*100, 2))
    print ' ---> number of failed jobs                : {}    {}%'.format(n_failed, round(n_failed/float(n_tot)*100, 2))
    print ' ---> number of running jobs               : {}    {}%'.format(n_unfinished, round(n_unfinished/float(n_tot)*100, 2))
    print ' ---> number of idle jobs                  : {}    {}%'.format(n_unprocessed, round(n_unprocessed/float(n_tot)*100, 2))

    if self.dofullreport:
      print '\n'
      print ' Error summary: '
      print '     - xrootd error:        {} ({}%)'.format(n_failure_xrootd, round(n_failure_xrootd/float(n_failed)*100,1) if float(n_failed)!=0 else 0)
      print '     - file read error:     {} ({}%)'.format(n_failure_readerr, round(n_failure_readerr/float(n_failed)*100,1) if float(n_failed)!=0 else 0)
      print '     - slurm timeout error: {} ({}%)'.format(n_failure_timeout, round(n_failure_timeout/float(n_failed)*100,1) if float(n_failed)!=0 else 0)
      print '     - slurm memout error:  {} ({}%)'.format(n_failure_memout, round(n_failure_memout/float(n_failed)*100,1) if float(n_failed)!=0 else 0)
      print '     - slurm node failure:  {} ({}%)'.format(n_failure_node, round(n_failure_node/float(n_failed)*100,1) if float(n_failed)!=0 else 0)
      print '     - other:               {} ({}%)'.format(n_failure_other, round(n_failure_other/float(n_failed)*100,1) if float(n_failed)!=0 else 0)

    if self.dofetchtime:
      print '\n'
      print ' Time summary: '
      print '     - average wallclock time: {}min'.format(round(time_wallclock/(n_good*60), 1) if n_good != 0 else 0)
      print '     - average CPU time/event: {}s/event'.format(round(time_cpu/n_good, 3) if n_good != 0 else 0)

    print '\n'
    print '----------------------------------------------------------'




if __name__ == "__main__":

  opt = getOptions()

  checkParser(opt)

  NanoProdManager(opt).process()


