#!/usr/bin/python

'''
**********************
sort_and_index_bams.py
**********************

Given a directory, sorts and then indexes every bam file in the directory
using a `DRMAA <http://www.drmaa.org/>`_ enabled cluster.
 
.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:version: 1.3

======================
Command-line Arguments
======================

**usage\:** 
    sort_and_index_bams.pl 
    :option:`-d|--datapath` *<path>*
    :option:`-l|--logfile` *<file>* 
    [:option:`--samtoolspath` *<path>*]
    [:option:`--clustq` *<str>*] 
    [:option:`--tmpdir` *<path>*]
    [:option:`-v|--verbose`] 
    [:option:`-h|--help`]

-------------------
Required Parameters
-------------------

.. option:: -d <path>, --datapath <path>

    Path to the data directory with the bam files to be processed

.. option:: -l <filename>, --logfile <filename>

    The name (inc. path) of the log file to write to.

Other options
-------------

.. option:: --samtoolspath <path> (Default: /local/bin/samtools)

  Specify the path to the samtools executable. 

.. option:: --clustq <str> (Default: 64bit-pri.q)

  Set the q to use to run cluster jobs.

.. option:: --tmpdir <path> (Default: ./.temp)

  The full path of the temporary directory you with to use for processing
  intermediate files.

.. option:: -v, --verbose

    Turn on verbose logging.

.. option:: -h, --help

  Print a basic description of the tool and its options to STDOUT.

======
Output
======

Sorted bam file (affixed with 'sorted') and a bam file index (appended with 
the extension .bai) for every bam file in the directory.

'''

import os, sys, re, warnings, drmaa
from optparse import OptionParser, OptionGroup

from Modules.housekeeping import custom_formatwarning
from Modules.housekeeping import parse_required_options, parse_allowed_values
from Modules.housekeeping import parse_option_type
from Modules.housekeeping import timeStr
from Modules.housekeeping import write_log_header

imported_modules = ["os", "sys", "re", "warnings", "optparse", "drmaa"]

version_string="1.3"

# force sequential output to the screen, rather than buffered output.
def printf(fmt, *args): sys.stdout.write(fmt % args)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

# override warning messages with a custom format
warnings.formatwarning = custom_formatwarning


def sortOrIndexBams(files, options, mode="sort"):
    
    ''' Sort or index a set of bam files.
    
    This function wraps :func:`sortOrIndexBam`, calling this function 
    for each bam file in a directory.
    
    :param list files: a list of bam files to sort, sorted bam files to index.
    
    :param options: an :py:class:`optparse.Values` instance.
    
    :param str mode: either 'sort' or 'index' (default: *sort*).
    
    :param str outfile: if **mode** is *sort*, specify an output filename.
    
    The **options** must include:
          
        * *log* - specify the logfile to write to.
        
        * *samtoolspath* - the location of samtools to use.
        
        * *tmpdir* - the path of the temporary directory to use for processing.
        
        * *clustq* - the cluster queue to use (can be *None*).
        
        * *verbose* - toggle on/off verbose logging
         
    The function returns a list of cluster output filenames.
    ''' 

    # sanity checks
    modes = ["sort", "index"]
    if mode not in modes:
        raise ValueError("Mode not recognized. Please choose either 'sort' " \
                         "or 'index'.")
    
    # get queue information for logging
    queue = "default"
    if options.clustq is not None:
        queue=options.clustq
    
    log = open(options.log,"a")
    if mode=="sort":
        log.write("\n %s: Cluster enabled: Running samtools sort " \
                  "on a cluster node in the %s queue for each bam file.\n" \
                  "" % (timeStr(),queue))
    elif mode=="index":
        log.write("\n %s: Cluster enabled: Running samtools index " \
              "on a cluster node in the %s queue for each bam file.\n" \
              "" % (timeStr(),queue))
        
    # use drmaa to spread multiple jobs across the cluster.
    # open a drmaa
    c_session = drmaa.Session()
    c_session.initialize()
        
    # check to see if the options.tmpdir path exists: if it doesn't, create it.
    if not os.path.exists(options.tmpdir):
        if options.verbose:#
            log = open(options.log,"a")
            log.write("\t\tCreating temp directory: %s\n" % options.tmpdir)
            log.close()
        os.makedirs(options.tmpdir)
    
    # loop throught he files calling a new cluster job for each
    joblist = []
    finished_files={}
    for bam_file in files:
        outfile = False
        if mode=="sort":
            outfile = "%s.sorted" % re.sub("\.bam","",bam_file)
        
        jobid, finished_file = sortOrIndexBam(c_session, bam_file, options,
                                              mode = mode, outfile = outfile)
                
        joblist.append(jobid)
        
        if mode=="sort":
            finished_files[outfile] = finished_file
        elif mode=="index":
            finished_files[bam_file] = finished_file
    
    # wait for all these jobs to finish before moving on.
    # the final 'false' here deals with cleaning up the account
    # info for each job. Since we want to check that each job
    # exited successfully we set it to false and then call
    # 'wait' on each job to check exit status. 
    c_session.synchronize(joblist, 
                          drmaa.Session.TIMEOUT_WAIT_FOREVER, 
                          False)
    
    for curjob in joblist:
        # collect job accoutning information
        retval = c_session.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
                    
        # if job exited pooly, log it!
        log = open(options.log,"a")
        if retval.exitStatus != 0:
            log.write("ERR: Job %s finished with errors! Check out '%s/" \
                      "samtools.o%s' " % (retval.jobId, options.tmpdir, 
                                     retval.jobId))
            sys.exit("ERR: Job %s finished with errors! Check out '%s/" \
                      "samtools.o%s' " % (retval.jobId, options.tmpdir, 
                                     retval.jobId))
        else:
            log.write("\tJob %s finished sucessfully\n" % retval.jobId)
        
        log.close()        
    
    c_session.exit()
        
    return(finished_files)

def sortOrIndexBam(c_session, bam_file, options, mode="sort", outfile=False):
    
    ''' Calls a cluster run of samtools sort/index on a bam file 
    
    :param c_session: a :py:class:`drmaa.Session` object for talking to the 
                      cluster.
    
    :param str bam_file: filename of bam to sort, or a sorted bam file to index.
    
    :param options: an :py:class:`optparse.Values` instance.
    
    :param str mode: either 'sort' or 'index' (default: *sort*).
    
    :param str outfile: if **mode** is *sort*, specify an output filename.
    
    The **options** must include the keywords:
        
        * *log* - specify the logfile to write to.
        
        * *samtoolspath* - the location of samtools to use.
        
        * *tmpdir* - the path of the temporary directory to use for processing.
        
        * *clustq* - the cluster queue to use (can be *None*).
        
        * *verbose* - toggle on/off verbose logging 

    The function returns the filename of the cluster output file.
    '''
    
    # sanity checks
    modes = ["sort", "index"]
    if mode not in modes:
        raise ValueError("Mode not recognized. Please choose either 'sort' " \
                         "or 'index'.")
    
    if mode=="sort" and not outfile:
        raise ValueError("Please specify an output filename when running in " \
                         "'sort' mode.")
    
    log = open(options.log,"a")
    if mode=="sort":
        log.write("\tSorting file %s...\n" % bam_file)
        args = ["sort", bam_file, outfile]
    elif mode=="index":
        log.write("\tIndexing file %s...\n" % bam_file)
        args = ["index", bam_file]
    
    # check to see if the options.tmpdir path exists: if it doesn't, create it.
    if not os.path.exists(options.tmpdir):
        if options.verbose:
            log = open(options.log,"a")
            log.write("\t\tCreating temp directory: %s\n" % options.tmpdir)
            log.close()
        os.makedirs(options.tmpdir)
    
    # specify the cluster job with a command, its arguments, the output and
    # error paths for the cluster log files and the environment for the job
    c_job = c_session.createJobTemplate()
    c_job.remoteCommand = options.samtoolspath
    c_job.args = args
    c_job.outputPath = ":%s" % options.tmpdir
    c_job.errorPath = ":%s" % options.tmpdir
    c_job.jobEnvironment = os.environ
    
    # add support for different cluster queue specifications
    if options.clustq is not None:
        c_job.nativeSpecification = "-clear -q '%s'" % options.clustq
    
    # verbose logging commands
    if options.verbose:
        log.write("\t\tdrmaa command line:\n\t\t%s %s\n" \
              "" % (options.samtoolspath," ".join(c_job.args)))        
        log.write("\t\tdrmaa output intermediates written too: %s\n" \
                  "" % options.tmpdir)
    
    # run the job
    jobid = c_session.runJob(c_job)
    
    # log the result
    log.write("\t\tJob submitted with id: %s\n" % jobid)
    log.close()
    
    # return the cluster output filename
    return(jobid, "%s/samtools.o%s" % (options.tmpdir, 
                                      jobid))

if __name__ == '__main__':

    # parse command line options - note that because we're using python 2.6.4,
    # and hence optparse rather than argparse, I'm forcing everything to be an
    # 'option' even if it is required.
    
    # define options
    optslist = []
    optslist.append({
                     "short": "-d",
                     "long": "--datapath",
                     "dest": "datapath", 
                     "action": "store",
                     "help": "Path to the data directory of the experiment.",
                     "group": "required",
                     "default": None,
                     "opt_type": "input_path"
                    })
    optslist.append({
                     "short": "-l",
                     "long": "--logfile",
                     "dest": "log", 
                     "action": "store",
                     "help": "The name (inc. path if different from the " \
                             "current working directory) of the log file " \
                             "from running this script.",
                     "group": "required",
                     "default": None,
                     "opt_type": "output_file"
                    })
    optslist.append({
                     "short": "-v",
                     "long": "--verbose",
                     "dest": "verbose", 
                     "action": "store_true",
                     "help": "Turn on verbose logging.",
                     "group": None,
                     "default": False
                    })
    optslist.append({
                     "short": None,
                     "long": "--tmpdir",
                     "dest": "tmpdir", 
                     "action": "store",
                     "help": "The full path of the temp dir you with to use " \
                             "for storing intermediate files. The default " \
                             "location for this is a dir called '.temp' in " \
                             "the current directory.",
                     "group": None,
                     "default": "%s/.temp" % os.getcwd(),
                    })
    optslist.append({
                     "short": None,
                     "long": "--clustq",
                     "dest": "clustq", 
                     "action": "store",
                     "help": "Set the q to use to run cluster jobs. " \
                             "Default is 64bit-pri.q.",
                     "group": None,
                     "default": None,
                    })
    optslist.append({
                     "short": None,
                     "long": "--samtoolspath",
                     "dest": "samtoolspath", 
                     "action": "store",
                     "help": "The full path to the samtools instalation you " \
                             "want to use to manipulate sam/bam files. The " \
                             "default location for samtools is /local/bin/" \
                             "samtools",
                     "group": None,
                     "default": "/local/bin/samtools",
                     "opt_type": "input_path"
                    })
    
    # build options from commandline
    usage = "\n\t%prog -d|--datapath <path> -l|--logfile <file> " \
            "\n\t[--samtoolspath <path>] [--clustq <str>] " \
            "\n\t[--tmpdir <path>] [-v|--verbose] [--help]"
    
    parser = OptionParser(usage, version="%prog v"+version_string)
        
    #define option groups
    reqgroup = OptionGroup(parser, "Required Arguments")
    
    # add options to the groups
    for val in optslist:
        if val["group"]=="required":
            reqgroup.add_option(val["short"], val["long"], 
                                action=val["action"], dest=val["dest"],
                                help=val["help"], default=val["default"])
        else:
            if val["short"] is not None:
                parser.add_option(val["short"], val["long"], 
                                  action=val["action"], dest=val["dest"], 
                                  help=val["help"], default=val["default"])
            else:
                parser.add_option(val["long"], action=val["action"], 
                                  dest=val["dest"], help=val["help"], 
                                  default=val["default"])
        
    # add groups to the parser
    parser.add_option_group(reqgroup)
    
    # parse arguments
    options, arguments = parser.parse_args()
    
    # check that all required options have been provided
    parse_required_options(options, optslist)
    
    # check all values are allowed, or default to, well, the defaults!
    parse_allowed_values(options, optslist)
    
    # check that all options are of the right type (if specified)
    parse_option_type(options, optslist)
    
    # open the logfile and write the basic header information
    log = write_log_header(options.log, options, optslist,
                           version_string=version_string, 
                           imported_modules=imported_modules)
    
    dirlist = os.listdir(options.datapath)
    bamfiles=[]
    for entry in dirlist:
        if re.match("^.+\.bam", entry):
            bamfiles.append("%s/%s" % (options.datapath, entry))
        
    sorted_files = sortOrIndexBams(bamfiles, options, mode="sort")
    indexed_files = sortOrIndexBams(sorted_files.keys(), options,
                                    mode="index")
    
    log = open(options.log,"a")
    log.write("\n %s: finished" % timeStr())
    log.close()

