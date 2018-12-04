#!/usr/bin/env python

"""
This script runs GATK BaseRecalibrator and/or ApplyBQSR. The script opens a
GATK process with the correct parameters.
"""

import shlex
import subprocess as sp
import bioexcel_align.alignutils as au

def baserecal(jopts, threads, ref, infile, knownsites, gatkdir, sample):
    '''
    Create and run command for GATK BaseRecalibratorSpark Local mode
    '''

    au.make_paths(gatkdir)

    command = str("gatk BaseRecalibratorSpark \
        --java-options {0} \
        --spark-master local[{1}] \
        -R {2} \
        -I {3} \
        --known-sites {4} \
        -O {5}/{6}.recal.table".format(jopts, threads, ref, infile, knownsites,
                            gatkdir, sample))

    cmdargs = shlex.split(command)
    print(command)
    print(cmdargs)

    p = sp.Popen(command, shell=True, executable='/bin/bash')

    return p

def applybqsr(jopts, threads, infile, gatkdir, sample):
    '''
    Create and run command for GATK ApplyBQSRSpark Local mode
    '''

    au.make_paths(gatkdir)

    command = str("gatk ApplyBQSRSpark \
    --java-options {0} \
    --spark-master local[{1}] \
    -I {2} \
    --bqsr-recal-file {4}/{5}.recal.table \
    -O {4}/{5}.final.bam".format(jopts, threads, infile, gatkdir, sample))

    cmdargs = shlex.split(command)
    print(command)
    print(cmdargs)

    p = sp.Popen(command, shell=True, executable='/bin/bash')

    return p

if __name__ == "__main__":
    description = ("This script runs GATK BaseRecalibrator and/or ApplyBQSR")
    args = au.parse_command_line(description)

    args.files = au.get_files(args)

    pbr = baserecal(args.jvm_opts, args.threads, args.ref, args.files, 
                                        args.knownsites, args.gatkdir, args.sample)
    pbr.wait()

    pab = applybqsr(args.jvm_opts, args.threads, args.files, args.recal, 
                                                    args.gatkdir, args.sample)
    pab.wait()
    