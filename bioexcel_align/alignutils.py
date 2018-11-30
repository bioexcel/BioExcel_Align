#!/usr/bin/env python

"""
This script contains several utility functions for use with the SeqQC portion
of the IGMM Cancer Genome Sequencing workflow developed as part of BioExcel
"""

import os
import sys
import argparse
from argparse import HelpFormatter
from datetime import datetime
import shutil

class MyFormatter(HelpFormatter):
    """
        From: https://stackoverflow.com/questions/9642692/argparse-help-without-duplicate-allcaps
    """

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_optional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(option_string)

                return '%s %s' % (', '.join(parts), args_string)

            return ', '.join(parts)

    def _get_default_metavar_for_optional(self, action):
        return action.dest.upper()
    '''
    https://stackoverflow.com/questions/32888815/max-help-position-is-not-works-in-python-argparse-library
    '''
    def add_argument(self, action):
        if action.help is not argparse.SUPPRESS:

            # find all invocations
            get_invocation = self._format_action_invocation
            invocations = [get_invocation(action)]
            current_indent = self._current_indent
            for subaction in self._iter_indented_subactions(action):
                # compensate for the indent that will be added
                indent_chg = self._current_indent - current_indent
                added_indent = 'x'*indent_chg
                invocations.append(added_indent+get_invocation(subaction))
            # print('inv', invocations)

            # update the maximum item length
            invocation_length = max([len(s) for s in invocations])
            action_length = invocation_length + self._current_indent
            self._action_max_length = max(self._action_max_length,
                                          action_length)

            # add the item to the list
            self._add_item(self._format_action, [action])

def parse_command_line(description):
    """
    Parser of command line arguments for BioExcel_Align
    """
    formatter_class = lambda prog: MyFormatter(prog, max_help_position=40,
                                                    width=80)
    parser = argparse.ArgumentParser(description=description,
                                    formatter_class=formatter_class)

    maingroup = parser.add_argument_group('Generic arguments',
                    'Common arguments used when running pipeline.')
    maingroup.add_argument("-f", "--files", nargs='+', metavar=('F1'),
                            help="Either a pair of input FastQ files (BWA), or "
                            "single BAM file (GATK), depending on usage.")
    maingroup.add_argument("-s", "--sample",
                        help="Sample name to prepend to all output files.")
    maingroup.add_argument("-o", "--outdir", default='./', metavar='PATH',
                        help="Output directory. (default: current directory)")
    maingroup.add_argument("--tmpdir", metavar='PATH', help="Tmp directory. "
                                    "(default: system tmp location)")
    maingroup.add_argument("-t", "--threads", default=2, type=int, metavar='T',
                        help="Max number of threads to use. NOTE: not all"
                        "stages use all threads. (default: 2)")
    ####
    bwagroup = parser.add_argument_group('BWA Mem/alignment stage',
        'Additional arguments used when running the BWA stage manually')
    bwagroup.add_argument("--bwa_version",  default='stable',
                            choices=['stable', 'beta'],
                            help="Version of BWA Mem stage to run: Stable "
                            "or more recent, beta version")
    bwagroup.add_argument("--bwa_ind_ref", 
                        default='genomes/Hsapiens/GRCh37/bwa/GRCh37.fa',
                        help="Location of the indexed reference genome file "
                        "for use with BWA")
    ####
    gatkgroup = parser.add_argument_group('GATK stages',
                'Additional arguments needed when running the GATK stage ')
    gatkgroup.add_argument("-j", "--jvm_opts", metavar='J_ARGS', default='',
                help="Arguments passed to Java when running GATK (e.g. max "
                "memory, tmpdirs).")
    gatkgroup.add_argument("-r", "--ref", metavar='FILE', 
                default='genomes/Hsapiens/GRCh37/seq/GRCh37.2bit',
                help="Reference sequence file.")  
    gatkgroup.add_argument("-k", "--knownsites", metavar='FILE', 
                default='genomes/Hsapiens/GRCh37/variation/dbsnp-147.vcf.gz',
                help="One or more databases of known polymorphic sites.")              

    args = parser.parse_args()
    if not args.files and __name__ != "__main__":
        sys.exit("\nusage: bioexcel_align -h for help \n\nbioexcel_align "
                    "error: the following arguments are required: -f/--files")

    if args.tmpdir:
        args.tmpdir = os.path.abspath(args.tmpdir)
    args.outdir = os.path.abspath(args.outdir)
    args.bwadir = os.path.abspath("{0}/BWA_out".format(args.outdir))
    args.gatkdir = os.path.abspath("{0}/GATK_out".format(args.outdir))
    
    args.date = datetime.now().strftime('%Y_%m_%d')
    args.files = get_files(args.files)
    if not args.sample: 
        args.sample=args.date

    return args

def make_paths(dirpath):
    """
    Create paths required for run of SeqQC pipeline
    """
    if not os.path.exists(dirpath):
        os.makedirs(dirpath, exist_ok=True)

    return

def get_files(filelist):
    """
    Return list of files to pass through SeqQC pipeline, throws error
    """

    # make sure files exist
    for checkfile in filelist:
        if not os.path.isfile(checkfile):
            sys.exit("Error: Input file {} does not exist.".format(checkfile))
    # expand paths to files (now we know the all exist)
    infiles = [os.path.abspath(x) for x in filelist]
    return infiles

if __name__ == "__main__":
    description = ("This script contains a series of useful functions for the "
                "Alignment stage of the BioExcel Cancer Genome Workflow. Run "
                "to print out example args.")
    args = parse_command_line(description)
    print(args)
    