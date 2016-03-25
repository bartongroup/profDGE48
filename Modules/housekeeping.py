'''

====================
Modules.housekeeping
====================

This is the housekeeping module for the two python scripts that are in The
Great RNA-seq Experiment codebase. These functions include logging, parsing,
custom warning messages and custom formatting for time and date strings.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.1

:created_on: 2015-09-04

'''

__version__ = "1.1"

import os, sys, time, re, subprocess

try:
    import pysam
except:
    pass


def createNewTempdir(options, remove=True, create=True):

    ''' Removes an existing tmpdir and create a new one.

    :param  options: a valid :py:class:`optparse.Values` instance.

    :param bool remove: Remove any existing directory (Default: True)

    :param bool create: Create a new directory (Default True)

    The path for the temporary directory is taken from the *tmpdir* keyword
    of the **options** instance, which must also include a *log* keyword
    containing the filename of the log file to write to..
    '''

    # open log
    log = open(options.log,"a",0)

    # remove any existing tmpdir
    if remove and os.path.exists(options.tmpdir):

        log.write("\n %s: Removing existing tmpdir %s..." \
                  ""% (timeStr(), options.tmpdir))
        try:
            subprocess.check_call(['rm','-rf',options.tmpdir])
        except subprocess.CalledProcessError as e:
            log.write("\n %s: Failed to remove temporary directory '%s'" \
                      ""% (timeStr(),options.tmpdir)
                      )
            log.write("\n %s %s" % (e.cmd,
                                    e.retcode)
                      )
        log.write("done.\n")

    # create new tmpdir
    if create:
        log.write("\n %s: Creating tmpdir %s..." \
                  ""% (timeStr(), options.tmpdir))
        os.makedirs(options.tmpdir)
        log.write("done.\n")

    log.close()

def custom_formatwarning(message, category, filename, lineno, line):

    """ Print a custom warning message format that stands out a bit more!

    Takes in the standard elements of a warning message from the warnings
    module (i.e.,  message, a category, a filename the line number and the
    line) and returns a formatted string encapsulating this information.

    :param str message: the message to print

    :param str category: the category of the warning

    :param str filename: the filename of the script

    :param int lineno: the line number generating the warning

    :param str line: the line that generated the warning

    """

    warn_type = re.match(".+\.(.+)\'\>",str(category)).group(1)
    return "\n---\n%s, line %s %s: %s\n---\n" \
           "" % (filename, lineno, warn_type, message)


def parse_required_options(options, optslist):

    """ Checks that all required options have been specified.

    :param options: an :py:class:`optparse.Values` instance.

    :param list optslist: a list of options.

    The **options** list should be a list of dictionaries with the
    'group' key. If 'group'=='required' then its a required option. If the
    value of this option in the :py:modeule:optparse instance is not specified,
    then exit the script and provide help.

    """

    for val in optslist:
        if val["group"]=="required" and options.__dict__[val["dest"]] is None:
            sys.exit(
                     "A required argument (%s) is missing.\nPlease specify " \
                     " a value using '%s arg'|'%s=arg' or see -h|--help for " \
                     " more details on how to run this script."
                     % (val["long"],
                        val["short"],
                        val["long"])
                    )

def parse_allowed_values(options, optslist):

    """ Checks that all options conform to their allowed values

    :param options: an :py:class:`optparse.Values` instance

    :param list optslist: a list of options.

    The **options** list should be a list of dictionaries with an
    optional 'allowed_vals' key. 'allowed_vals' should contain a list of the
    allowed values. If the value of this option in the :py:module:optparse
    instance is not in this list then warn the user and use the default value.

    """

    for val in optslist:
        if "allowed_vals" in val.keys():
            if options.__dict__[val["dest"]] not in val["allowed_vals"]:
                warnstring = "%s is not an allowed option for %s. " \
                             "Defaulting to %s. See -h|--help for more " \
                             "details on how to run this script." \
                             "" % (options.__dict__[val["dest"]],
                                   val["long"],
                                   val["default"])
                warning = warnings.warn(warnstring, SyntaxWarning)
                options.__dict__[val["dest"]]=val["default"]

def parse_option_type(options, optslist):

    """ Checks that all options conform to their specified types

    :param options: an :py:class:`optparse.Values` instance

    :param list optslist: a list of options.

    The **options** list should be a list of dictionaries with an
    optional 'opt_type' key. 'opt_type' should be a string specifying the type
    of the option; examples are input_path, output_path, output_file, string,
    int, float. all the 'opt_type' values should be valid python types except
    the paths. If the value of this option in the :py:module:optparse instance
    is not of this type then either warn the user and use the default value
    (if there is one) or exit the script and provide help.

    """

    for val in optslist:
        if "opt_type" in val.keys():

            # check all paths to make sure they are valid paths
            # paths are all essential so exit if they fail
            if val["opt_type"] in ["input_path", "output_path", "output_file"]:
                
                options.__dict__[val["dest"]] = os.path.abspath(options.__dict__[val["dest"]])
                
                # trim output file name from path
                if val["opt_type"]=="output_file":
                    check_path = os.path.dirname(options.__dict__[val["dest"]])
                else:
                    check_path = options.__dict__[val["dest"]]

                if re.match("^$",check_path):
                    check_path = os.getcwd()

                if not os.path.exists(check_path):
                    sys.exit(
                             "The path specified for %s, %s, is not valid.\n" \
                             "Please specify a valid path. See -h|--help " \
                             "for more details on how to run this script."
                             % (val["long"],
                                check_path)
                             )

            elif val["opt_type"]=="input_file":
                if not options.__dict__[val["dest"]] is None:
                    if not os.path.exists(options.__dict__[val["dest"]]):
                        sys.exit(
                                "The path specified for %s, %s, is not valid.\n" \
                                "Please specify a valid path. See -h|--help " \
                                "for more details on how to run this script."
                                % (val["long"],
                                   val["dest"])
                                )
                    else:
                        options.__dict__[val["dest"]] = os.path.abspath(options.__dict__[val["dest"]])
            else:
                if val["opt_type"]=="string":
                    try:
                        string(options.__dict__[val["dest"]])
                    except ValueError:
                        if val["default"] is not None:
                            warnstring = "The value specified for %s cannot " \
                                         "be successfully cast as a string. " \
                                         "Defaulting to %s. See -h|--help " \
                                         "for more details on how to run " \
                                         "this script." % (val["long"],
                                                      val["default"])
                            warning = warnings.warn(warnstring, SyntaxWarning)
                            options.__dict__[val["dest"]]=val["default"]
                        else:
                            sys.exit("The value specified for %s cannot be " \
                                     "successfully cast as a string and has " \
                                     "no default value. please specify a " \
                                     "valid value. See -h|--help for more " \
                                     "details on how to run this script." \
                                     "" % val["long"]
                                     )
                elif val["opt_type"]=="int":
                    try:
                        int(options.__dict__[val["dest"]])
                    except ValueError:
                        if val["default"] is not None:
                            warnstring = "The value specified for %s cannot " \
                                         "be successfully cast as a int. " \
                                         "Defaulting to %s. See -h|--help " \
                                         "for more details on how to run " \
                                         "this script." % (val["long"],
                                                      val["default"])
                            warning = warnings.warn(warnstring, SyntaxWarning)
                            options.__dict__[val["dest"]]=val["default"]
                        else:
                            sys.exit("The value specified for %s cannot be " \
                                     "successfully cast as a int and has " \
                                     "no default value. please specify a " \
                                     "valid value. See -h|--help for more " \
                                     "details on how to run this script." \
                                     "" % val["long"]
                                     )
                elif val["opt_type"]=="float":
                    try:
                        float(options.__dict__[val["dest"]])
                    except ValueError:
                        if val["default"] is not None:
                            warnstring = "The value specified for %s cannot " \
                                         "be successfully cast as a float. " \
                                         "Defaulting to %s. See -h|--help " \
                                         "for more details on how to run " \
                                         "this script." % (val["long"],
                                                      val["default"])
                            warning = warnings.warn(warnstring, SyntaxWarning)
                            options.__dict__[val["dest"]]=val["default"]
                        else:
                            sys.exit("The value specified for %s cannot be " \
                                     "successfully cast as a float and has " \
                                     "no default value. please specify a " \
                                     "valid value. See -h|--help for more " \
                                     "details on how to run this script." \
                                     "" % val["long"]
                                     )
    
    return(options)
    
def timeAndDateStr():

    """ A pretty current time and date timestring.

    Looks like 12:34 2015/01/01"""

    ttime = time.gmtime()
    timestring =  "%02.0d:%02.0d %02.0d/%02.0d/%02.0d" \
                  "" % (ttime[3],ttime[4],ttime[0],ttime[1],ttime[2])
    return(timestring)

def timeStr():

    """ A pretty current time timestring

    Looks like 12:34"""

    ttime = time.gmtime()
    timestring =  "%02.0d:%02.0d:%02.0d" % (ttime[3],ttime[4],ttime[5])
    return(timestring)

def write_log_header(options, optslist, version_string="",
                     imported_modules=[]):

    """ Write the standard header of a new log file

    Writes a standard set of details information to the header of a logfile.
    The information includes the script name, the command line, the values of
    all command parameters the script was run with (and whether they are the
    parameter default values). If set to 'verbose' mode, it also outputs
    information about:

        * python version.
        * imported modules (where possible).
        * perl version (where found).
        * R version (where found).
        * Rscript version (where found).
        * environmental variables:
            * PYTHONPATH
            * PERLLIB
            * PERL5LIB
            * DRMAA_LIBRARY_PATH
            * R_HOME
            * LD_LIBRARY_PATH
            * PATH
        * :ref:`group_by_gene.pl <group_by_gene_perldoc>` version (where found).
        * :ref:`add_gene_name_column.pl <group_by_gene_perldoc>` version
          (where found).
        * `samtools <http://samtools.sourceforge.net/>`_ version (where found).

    :param options: an :py:class:`optparse.Values` instance

    :param list optslist: a list of options

    :param str version_string: the version of the script for logging

    :param list imported_modules: a list of string of the names of the modules
                                  imported by the script being logged.

    The options list should be a list of dictionaries each with the format::

        {"short": None,
         "long": "--name",
         "dest": "name",
         "action": "store",
         "help": "Some help",
         "group": optional/required,
         "default": "a default",
         allowed_vals: ["some", "values"]
        })

    """

    # get script name
    script_title = "%s v%s" % (re.sub("^.+\/","",sys.argv[0]), version_string)
    title_lines = "=" * len(script_title)

    log = open(options.log,"w")

    # write the title and log time
    log.write("\n" \
              " %s\n" \
              " %s\n" \
              " %s\n" \
              "\n" \
              " Logfile opened at %s\n" \
              "\n" \
              "" % (title_lines, script_title, title_lines, timeAndDateStr())
              )

    # write the full command line
    log.write(" Command line:\n\t%s\n\n" % " ".join(sys.argv))

    log.write(" Full script options:\n")
    for option in optslist:
        opt_printstr = "%s:" % option["long"]
        optstring = "\t%s%s" % (opt_printstr.ljust(16, " "),
                                options.__dict__[option["dest"]])
        if option["short"] not in sys.argv and option["long"] not in sys.argv:
            optstring = "%s (default)" % optstring

        log.write("%s\n" % optstring)

    if options.verbose:

        # list details of python
        log.write("\n Programs:\n"\
                  "\tPython %s\n" % re.sub("\n","",sys.version)
                  )

        # check for group_by_gene.pl
        try:
            # get version info from 'group_by_gene.pl'
            gbg_call = subprocess.Popen([options.gbgpath, "--version"],
                                        stdout=subprocess.PIPE)
            gbg_version_info = re.split("\n",gbg_call.communicate()[0])
            gbg_version = None
            other_versions = []
            for info in gbg_version_info:
                if re.match("group_by_gene.pl",info):
                    gbg_version = info
                elif info!="":
                    other_versions.append(info)

            gbg_version_string = "%s (uses %s)" \
                                 "" % (gbg_version,", ".join(other_versions))

            log.write("\t%s\n" % gbg_version_string)
        except AttributeError:
            pass

        # check for add_gene_name_column.pl
        try:
            agnc_call = subprocess.Popen([options.agncpath, "--version"],
                                        stdout=subprocess.PIPE)
            agnc_version_info = re.split("\n",agnc_call.communicate()[0])
            agnc_version = None
            ensembl_version = None
            for info in agnc_version_info:
                if re.match("add_gene_name_column.pl",info):
                    agnc_version = info
                elif re.match("ensembl API",info):
                    ensembl_version = info

            agnc_version_string = "%s (uses %s)" \
                                 "" % (agnc_version, ensembl_version)

            log.write("\t%s\n" % agnc_version_string)
        except AttributeError:
            pass

        # check for samtools
        try:
            samtools_call = subprocess.Popen([options.samtoolspath],
                                              stderr=subprocess.PIPE)
            samtool_string = re.sub("\n","",samtools_call.communicate()[1])
            samtools_version_info = re.match(".+Version: (.+?)Usage:.+",
                                             samtool_string).group(1)

            samtools_version_string = "samtools %s" % (samtools_version_info)

            log.write("\t%s\n" % samtools_version_string)
        except AttributeError:
            pass
        except OSError:
            log.write("unable to capture samtools version information.\n" \
                      "Are you sure it exists?!")

        # check for R
        try:
            # get version info from 'R'
            R_call = subprocess.Popen([options.rpath, "--version"],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
            R_response = R_call.communicate()
            if R_response[0]!="":
                # if its the R binary it will use this
                R_string = re.split("\n",R_response[0])[0]
                R_version_string = re.split(" -- ", R_string)[0]
            else:
                # if its the RScript binary it will use this
                R_version_string = R_response[1].strip()

            R_version_string = re.sub("^R","%s" % options.rpath,
                                      R_version_string)
            log.write("\t%s\n" % R_version_string)

        except AttributeError:
            log.write("***Cannot get R version! Something may be wrong.\n")

        # check for perl
        try:
            # get version info from 'perl'
            perl_call = subprocess.Popen([options.perlpath, "--version"],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
            perl_response = perl_call.communicate()
            perl_version_info = re.split("\n", perl_response[0])
            perl_version_string = None
            for line in perl_version_info:
                match = re.match("This is perl.+\((.+?)\) built.+", line)
                if match:
                    perl_version_string = match.group(1)

            perl_version_string = "%s %s" % (options.perlpath,
                                             perl_version_string)

            log.write("\t%s\n" % perl_version_string)
        except AttributeError:
            log.write("***Cannot get perl version! Something may be wrong.\n")

        # check for Rscript
        try:
            RScript_call = subprocess.Popen([options.rpath, options.script_file,
                                       "--version"],
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE
                                      )

            RScript_response = RScript_call.communicate()
            # The info should be here
            RScript_version_info = re.split("\n",RScript_response[0])
            skip=True
            script_version = None
            otherinfo = []
            for info in RScript_version_info:
                match = re.match("%s" % os.path.basename(options.script_file),
                                 info)
                if match and skip:
                    skip=False
                    script_version = re.sub("\t"," v", info)
                elif not skip and info!="":
                    otherinfo.append(re.sub("\t"," v", info))

            script_version = re.sub("^R","%s" % options.rpath, script_version)
            script_version = re.sub("version","", script_version)

            R_version_string = "%s\n\t\t(uses %s)" % (script_version,
                                                      ", ".join(otherinfo))
            log.write("\t%s\n" % R_version_string)
        except:
            try:
                log.write("\n %s: Could not get %s script version info" \
                          "" % (timeStr(),
                                os.path.basename(options.script_file)))
            except AttributeError:
                pass

        # munge module information
        standard_modules = []
        third_party_modules = []
        for module in imported_modules:
            i = __import__(module)
            try:
                module_version = i.__version__
                third_party_modules.append("%s %s" % (module, module_version))
            except AttributeError:
                standard_modules.append(module)

        # describe standard imports
        log.write("\n Standard python modules:\n")
        for module in standard_modules:
            log.write("\t%s\n" % module)

        # describe third party python modules
        log.write("\n Contributed python modules:\n")
        for module in third_party_modules:
            if re.match("pysam", module):
                module = "%s (uses samtools %s)" \
                         "" % (module, pysam.version.__samtools_version__)
            log.write("\t%s\n" % module)

        # print pythonpath environment variable
        log.write("\n PYTHONPATH environment variable:\n")
        for line in sys.path:
            if line!="":
                log.write("\t%s\n" % line)

    log.close()

    return(options.log)
