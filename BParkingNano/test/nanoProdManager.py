import sys 
import os
import glob
import subprocess

from nanoTools import NanoTools

def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='Script to report the status of a nanoAOD production and resubmit failed files', add_help=True)
  parser.add_argument('--pl'        , type=str, dest='pl'          , help='label of the sample file'                                             , default=None)
  parser.add_argument('--ds'        , type=str, dest='ds'          , help='[optional] name of dataset'                                           , default=None)
  parser.add_argument('--tag'       , type=str, dest='tag'         , help='[optional] tag name'                                                  , default=None)
  parser.add_argument('--dosignal'            , dest='dosignal'    , help='run the BToMuMuPi process'                       , action='store_true', default=False)
  parser.add_argument('--docontrol'           , dest='docontrol'   , help='run the BToKMuMu process'                        , action='store_true', default=False)
  parser.add_argument('--dohnl'               , dest='dohnl'       , help='run the HNLToMuMuPi process'                     , action='store_true', default=False)
  parser.add_argument('--dotageprobe'         , dest='dotageprobe' , help='run the JpsiToMuMu process (tag and probe study)', action='store_true', default=False)
  parser.add_argument('--mcprivate'           , dest='mcprivate'   , help='run the resubmitter on a private MC sample'      , action='store_true', default=False)
  parser.add_argument('--mccentral'           , dest='mccentral'   , help='run the resubmitter on a central MC sample'      , action='store_true', default=False)
  parser.add_argument('--data'                , dest='data'        , help='run the resubmitter on a data sample'            , action='store_true', default=False)
  parser.add_argument('--dofullreport'        , dest='dofullreport', help='add to report chunks and failure reason'         , action='store_true', default=False)
  parser.add_argument('--dofetchtime'         , dest='dofetchtime' , help='add to report time fetch'                        , action='store_true', default=False)
  parser.add_argument('--docheckfile'         , dest='docheckfile' , help='check the content of the nano files'             , action='store_true', default=False)
  parser.add_argument('--doresubmit'          , dest='doresubmit'  , help='resubmit failed jobs'                            , action='store_true', default=False)
  return parser.parse_args()


def checkParser(opt):

  if opt.pl==None:
    raise RuntimeError('Please indicate the production label of the sample')

  if opt.mcprivate==False and opt.mccentral==False and opt.data==False:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding either --data or--mcprivate or --mccentral to the command line')

  if opt.mcprivate + opt.mccentral + opt.data > 1:
    raise RuntimeError('Please indicate if you want to run on data or MC by adding only --data or --mcprivate or --mccentral to the command line')

  if opt.docheckfile and not (opt.dosignal or opt.docontrol or opt.dohnl or opt.dotageprobe):
    raise RuntimeError('Please indicate with branch you would like to check on (--dosignal or --docontrol or --dohnl or --dotageprobe)')



class NanoProdManager(NanoTools):
  def __init__(self, opt):
    self.prodlabel    = vars(opt)['pl']
    self.dataset      = vars(opt)['ds']
    self.tag          = vars(opt)['tag']
    self.dosignal     = vars(opt)["dosignal"]
    self.docontrol    = vars(opt)["docontrol"]
    self.dohnl        = vars(opt)["dohnl"]
    self.dotageprobe  = vars(opt)["dotageprobe"]
    self.mcprivate    = vars(opt)['mcprivate']
    self.mccentral    = vars(opt)['mccentral']
    self.data         = vars(opt)['data']
    self.dofullreport = vars(opt)['dofullreport'] 
    self.dofetchtime  = vars(opt)['dofetchtime'] 
    self.docheckfile  = vars(opt)['docheckfile'] 
    self.doresubmit   = vars(opt)['doresubmit'] 


  def isJobFinished(self, logfile):
    return self.scanLogFile(logfile, 'finished running nano step') or self.scanLogFile(logfile, 'slurmstepd: error')


  def scanLogFile(self, logfile, key):
    with open(logfile) as f:
      if key in f.read():
        return True
      else:
        return False
   

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


  def writeFileList(self, chunk, failed_files, label):

    logdir = NanoTools.getLogDir(self, failed_files[0], self.prodlabel, self.tag, self.data, 'mcprivate' if self.mcprivate else '')
    label = logdir[logdir.find('logs')+5:].replace('/', '_')

    if self.data:
      filename = './files/resubmit_data_{}'.format(label)
    elif self.mcprivate:
      filename = './files/resubmit_mcprivate_{}'.format(label)
    elif self.mccentral:
      filename = './files/resubmit_mccentral_{}'.format(label)

    for file_ in failed_files:
      # get the file to reprocess
      filelist = '{}/filelist_{}.txt'.format(chunk, label)
      command = 'head -n {line} {fl} | tail -1'.format(line=NanoTools.getStep(self, file_), fl=filelist)  
      file_toresubmit = subprocess.check_output(command, shell=True)

      # open resubmit file
      resubmit_file = open(filename + '_nj{}.txt'.format(NanoTools.getStep(self, file_)), 'w+')

      # write to file list
      resubmit_file.write(file_toresubmit + '\n')
      
      # close file list
      resubmit_file.close()

      #print 'created {}_nj{}.txt'.format(filename, NanoTools.getStep(self, file_))

    return filename 


  def getArray(self, failed_files): 
    idx = []
    for file_ in failed_files:
      idx.append(NanoTools.getStep(self, file_))
    return ','.join(idx)
    

  def resubmit(self, chunk, failed_files):
    # strategy: per chunk resubmission
    #           submit job arrays with indices corresponding to the stepId of the failed jobs

    logdir    = NanoTools.getLogDir(self, failed_files[0], self.prodlabel, self.tag, self.data, 'mcprivate' if self.mcprivate else '') 
    label     = logdir[logdir.find('logs')+5:].replace('/', '_')
    array     = self.getArray(failed_files)
    outputdir = failed_files[0][0:failed_files[0].find('bparknano')]
    filelist  = self.writeFileList(chunk, failed_files, label) 

    command = 'sbatch -p standard --account=t3 -o {ld}/nanostep_nj%a.log -e {ld}/nanostep_nj%a.log --job-name=nanostep_nj%a_{pl} --array {ar} --time=03:00:00 submitter.sh {outdir} {usr} {pl} {tag} {isMC} {rmt} {lst} 1 {dosig} {doctrl} {dohnl} {dotep}'.format(
      ld      = logdir,
      pl      = label,
      ar      = self.getArray(failed_files), 
      outdir  = outputdir,
      usr     = os.environ["USER"], 
      tag     = 0 if self.tag == None else self.tag,
      isMC    = 1 if self.mcprivate or self.mccentral else 0,
      rmt     = 0 if self.mcprivate else 1,
      lst     =  filelist,
      dosig     = 1 if self.dosignal else 0, 
      doctrl    = 1 if self.docontrol else 0, 
      dohnl     = 1 if self.dohnl else 0, 
      dotep     = 1 if self.dotageprobe else 0, 
    )

    os.system(command)


  def process(self):
    print '---------------------------------'
    print '         Nano Job Manager        '
    print '---------------------------------'

    if self.data and self.dataset==None:
      print '\nINFO: the JobManager Tool is going to be run on all dataset having "{}" as production label\n'.format(opt.pl)

    # pnfs directory where nano samples are located
    if self.data: dirtag = 'data'
    elif self.mccentral: dirtag = 'mccentral'
    else:  dirtag = ''
    location = NanoTools.getFilesLocation(self, dirtag)

    # get the directories associated to the production label
    dirs = NanoTools.getNanoDirectories(self, location, self.prodlabel, self.dataset, 'mcprivate' if self.mcprivate else '')
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

        n_exp = self.getNExpectedNanoFiles(chunk_)
        
        if not self.mcprivate:
          if self.tag == None:
            files = [chunk_+'/bparknano_nj'+str(nj)+'.root' for nj in range(1, n_exp+1)]
          else:
            files = [chunk_+'/bparknano_{}_nj'.format(self.tag) +str(nj)+'.root' for nj in range(1, n_exp+1)]

          logdir = NanoTools.getLogDir(self, files[0], self.prodlabel, self.tag, self.data, 'mcprivate' if self.mcprivate else '')
          logfiles = [NanoTools.getLogFile(self, logdir, file_) for file_ in files]

        else:
          point = chunk_[chunk_.find(self.prodlabel)+len(self.prodlabel)+1:chunk_.find('/nanoFiles')]
          if self.tag == None:
            filelist = chunk_ + '/filelist_' + self.prodlabel + '_' + point + '_' + chunk_[chunk_.find('Chunk'):] + '.txt'
          else:
            filelist = chunk_ + '/filelist_' + self.prodlabel + '_' + point + '_' + self.tag + '_' + chunk_[chunk_.find('Chunk'):] + '.txt'
          try: f = open(filelist)
          except:
            print ' -> no files found in this chunk'
            print ' --> skipping'
            continue
          lines = f.readlines()

          files = [chunk_+'/bparknano_{}_nj'.format(self.tag) +str(NanoTools.getStep(self, lines[nj-1]))+'.root' for nj in range(1, n_exp+1)]
          files_nlog = [chunk_+'/bparknano_{}_nj'.format(self.tag) +str(nj)+'.root' for nj in range(1, n_exp+1)]

          logdir = NanoTools.getLogDir(self, files[0], self.prodlabel, self.tag, self.data, 'mcprivate' if self.mcprivate else '')
          logfiles = [NanoTools.getLogFile(self, logdir, file_) for file_ in files_nlog]

        for ifile, file_ in enumerate(files):
          # get the log file
          logfile = logfiles[ifile]
          
          # idle jobs
          if not NanoTools.checkFileExists(self, logfile): 
            n_unprocessed_perchunk += 1
            continue

          branchname = ''
          if self.dosignal: branchname = 'nBToMuMuPi'
          elif self.docontrol: branchname = 'nBToKMuMu'
          elif self.dohnl: branchname = 'nHNLToMuPi'
          elif self.dotageprobe: branchname = 'nJPsiToMuMu'
          extra_cond = NanoTools.checkLocalFile(self, file_, cond=True, branch_check=True, branchname=branchname) if self.docheckfile else 'True'
          #if NanoTools.checkFileExists(self, file_) and NanoTools.checkLocalFile(self, file_, cond=True): # successfull job
          if NanoTools.checkFileExists(self, file_) and extra_cond: # successfull job
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
                # print '{} does not exist'.format(logfile)

                #print '{} does not exist'.format(file_)
                failed_files.append(file_)
  
       
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
          self.resubmit(chunk_, failed_files)


        n_good_perdir        += n_good_perchunk
        n_failed_perdir      += n_failed_perchunk
        n_unfinished_perdir  += n_unfinished_perchunk
        n_unprocessed_perdir += n_unprocessed_perchunk

        time_wallclock_perdir += time_wallclock_perchunk
        time_cpu_perdir       += time_cpu_perchunk


      n_tot_perdir = n_good_perdir + n_failed_perdir + n_unfinished_perdir + n_unprocessed_perdir
      if n_tot_perdir == 0: continue

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


