#!/usr/bin/env python

"""
This script runs the main stage of the Alignment stage. The script opens a
series of tools (BWA MEM, Samtools, etc) with the correct parameters.
"""

import shlex
import subprocess as sp
import Alignment as al

def run_main(arglist, outdir, infiles):
    """
    Create and run subprocess for fastqc
    """
    command = "fastqc -o {0} -d {1} -t {2} --extract {3}".format(outdir,
                             arglist.tmpdir, arglist.threads,
                             ' '.join(infiles))

    cmdargs = shlex.split(command)
    print(command)
    print(cmdargs)

    p = sp.Popen(cmdargs)
    #p = sp.Popen('date')

    return p

def main(arglist):
    """
    Main function to run standalone FastQC instance
    """
    print("Hello!")
    print(arglist)
    #run_fqc(args)

if __name__ == "__main__":
    description = ("This script runs the FastQC step of SeqQC")
    args = SeqQC.parse_command_line(description)
    main(args)
    