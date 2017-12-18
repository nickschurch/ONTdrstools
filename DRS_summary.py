#!/usr/bin/python
'''
--------------------------
ONTdrstools.DRS_summary.py
--------------------------

This script parses an aligned set of ONT DRS data, generating several sets of
summary statistics, making several plots, and generating a wig file for the
pileup of 3' ends of the reads.

.. moduleauthor:: Nick Schurch <nschurch@dundee.ac.uk>

:module_version: 1.0
:created_on: 2017-12-18

Command-line Arguments
======================

**usage\:** 
    DRS_summary.py
    :param: <input bam file>
    :option:`-l|--log` *<file>*
    [:option:`-p|--prefix` <str>]
    [:option:`-v|--verbose`] 
    [:option:`--version`] 
    [:option:`--help`]

Required Parameters
-------------------

:para: <input bam file>

  The input bam file

:option:`--logfile|-l`        

  The name (inc. path) of the log file from the wrapper.

Optional Parameter
------------------

:option: `--prefix|-p` <str>

  Prefix string for output filenames. Can optionally include a full path. 
  Defaults to the input filename.

:option:`--help|-h`

  Print a basic description of the tool and its options to STDOUT.

:option:`--version`    

  Show program's version number and exit.
    
:option:`--verbose|-v`     

  Turn on verbose logging (recommended).

Output
======

Undefined, as yet :D
'''

ver=1.30

__scriptname__= "DRS_summary"
__version__ = str(ver)
__usage__ = "\n\t%s <input bam file> -l|--logfile [-p|--prefix <str>]\n\t" \
            "[--version][-v|--verbose][--help]"
__progdesc__ = '''
This script parses an aligned set of ONT DRS data, generating several sets of
summary statistics, making several plots, and generating a wig file for the
pileup of 3' ends of the reads..
'''

__progepi__ = '''
--------------------------
ONTdrstools.DRS_summary.py
--------------------------
'''

import sys, pysam
import script_options.standard_parsers as sp
import script_logging.standard_logging as sl
#from parsing_routines.gff_gtf_tools import annotation

def addScriptOptions(parser, pos_args, kw_args):
    
    """ add script-specific script options """
    
    script_options_group = parser.add_argument_group('Options')
    
    hlpstr = "Prefix string for output filenames. Can optionally include a " \
             "full path. Defaults to the input filename."
    option_short_name = "p"
    option_name = "prefix"
    options_default = None
    
    script_options_group.add_argument('-%s' % option_short_name,
                                      '--%s' % option_name,
                                      action = 'store',
                                      dest = option_name,
                                      type = str,
                                      help = hlpstr,
                                      )
    kw_args.append((option_name, option_name, options_default))
    
    return(parser, pos_args, kw_args)

if __name__ == '__main__':

    # parse command line options
    # Set standard parser
    parser, pos_args, kw_args = sp.standard_parser(__version__,
                                                   prog = __scriptname__, 
                                                   usage = __usage__,
                                                   description = __progdesc__,
                                                   epilog = __progepi__,
                                                   infile = True,
                                                   outfile = False,
                                                   tmpdir = False)
        
    parser, pos_args, kw_args = addScriptOptions(parser, pos_args, kw_args)
    
    args = parser.parse_args()
           
    # setup standard logging information
    script_logger = sl.standard_logger(__version__, sys.argv, args, pos_args, 
                                       kw_args, script_name=__scriptname__)
        
    script_logger.info("Parsing alignment data from %s..." % args.infile)
    
    thisbam = pysam.AlignmentFile(args.infile, "rb")
    refs = thisbam.references
    read_lengths={}
    read_lens=[]
    read_aligned_lengths={}
    aligned_len=[]
    read_end_pos={}
    for read in thisbam.fetch("1",0,10000):
        print read.query_name
        read_lens.append(read.query_length)
        try: 
            read_lengths[read.query_length].append(read.query_name)
        except KeyError:
            read_lengths[read.query_length]=[read.query_name]
        
        aligned_len.append(read.reference_length)
        try: 
            read_aligned_lengths[read.reference_length].append(read.query_name)
        except KeyError:
            read_aligned_lengths[read.reference_length]=[read.query_name]
        
        read_ref = refs[read.reference_id]
        if read_ref not in read_end_pos.keys():
            read_end_pos[read_ref]={}
        
        if read.is_reverse:
            read_end = read.reference_start
        else:
            read_end = read.reference_end
        
        try:
            read_end_pos[read_ref][read_end].append(read.query_name)
        except KeyError:
            read_end_pos[read_ref][read_end]=[read.query_name]
    
    print 
    
    script_logger.info("Finished. Have a nice day! ;)")