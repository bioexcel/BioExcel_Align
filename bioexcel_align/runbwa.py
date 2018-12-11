#!/usr/bin/env python

"""
This script runs the main stage of the Alignment stage. The script opens a
series of tools (BWA MEM, Samtools, etc) with the correct parameters.
"""

import sys
import subprocess as sp
import bioexcel_align.alignutils as au

def bwamem_stable(bwa_ref, threads, date, sample, fqfiles, bwadir):
    """
    Create and run process for the original implementation of BWA-Mem alignment
    """

    au.make_paths(bwadir)

    # Future: Split this across several processes properly. Shell=True is not
    # safe.
    command = str('bwa mem -R "@RG\\tID:{0}\\tPL:illumina\\tPU:{0}\\tSM:{1}" '
            ' -t {2} -M {3} {4} | samblaster '
            '--splitterFile >(samtools view -hSu /dev/stdin | samtools sort '
            '-@ {2} /dev/stdin > {5}/{1}.sr.bam) --discordantFile >(samtools '
            'view -hSu /dev/stdin | samtools sort -@ {2} /dev/stdin > '
            '{5}/{1}.disc.bam) | samtools view -hSu /dev/stdin | samtools '
            'sort -@ {2} /dev/stdin > {5}/{1}.raw.bam'.format(date, sample,
                                threads, bwa_ref, ' '.join(fqfiles), bwadir))

    print(command)

    p = sp.Popen(command, shell=True, executable='/bin/bash')

    return p

def bwamem_beta(bwa_ref, threads, date, sample, fqfiles, bwadir):
    """
    Create and run process for the new implementation of BWA-Mem alignment
    """

    au.make_paths(bwadir)

    command = str('bwa mem -R "@RG\\tID:{0}\\tPL:illumina\\tPU:{0}\\tSM:{1}" '
    '-c 250 -M -t {2} {3} {4} | bamsormadup inputformat=sam threads={2} '
    'tmpfile={5}/tmp_{1} SO=coordinate '
    'indexfilename={5}/{1}.raw.bam.bai > {5}/{1}.raw.bam'.format(date, sample,
                                threads, bwa_ref, ' '.join(fqfiles), bwadir))

    print(command)

    p = sp.Popen(command, shell=True, executable='/bin/bash')

    return p

def samtools_index(infile):
    """
    Run samtools index on newly created, aligned bam files
    """
    command = str('samtools index {}'.format(infile))

    print(command)

    p = sp.Popen(command, shell=True, executable='/bin/bash')

    return p

if __name__ == "__main__":
    description = ("This script runs the BWA-Mem step of Bioexcel Align")
    args = au.parse_command_line(description)
    if args.bwa_version == 'stable':
        pbwa = bwamem_stable(args.bwa_ind_ref, args.threads, args.date,
                                        args.sample, args.files, args.bwadir)
        pbwa.wait()
    elif args.bwa_version == 'beta':
        pbwa = bwamem_beta(args.bwa_ind_ref, args.threads, args.date,
                                        args.sample, args.files, args.bwadir)
        pbwa.wait()

    else:
        sys.exit('No BWA Mem version selected')

    psamidx = samtools_index('{0}/{1}.raw.bam'.format(args.bwadir, args.sample))
    psamidx.wait()
