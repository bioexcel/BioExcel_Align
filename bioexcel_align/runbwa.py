#!/usr/bin/env python

"""
This script runs the main stage of the Alignment stage. The script opens a
series of tools (BWA MEM, Samtools, etc) with the correct parameters.
"""

import shlex
import subprocess as sp
import bioexcel_align.alignutils as au

def bwamem_stable(bwa_ref, threads, date, sample, fqfiles, bwadir):
    """
    Create and run process for the original implementation of BWA-Mem alignment
    """

    au.make_paths(bwadir)

    command = str('bwa mem -R "@RG\tID:{0}\tPL:illumina\tPU:{0}\tSM:{1}" '
            ' -t {2} -M {3} {4} | samblaster '
            '--splitterFile >(samtools view -hSu /dev/stdin | samtools sort '
            '-@ {2} /dev/stdin > {5}/{1}.sr.bam) --discordantFile >(samtools '
            'view -hSu /dev/stdin | samtools sort -@ {2} /dev/stdin > '
            '{5}/{1}.disc.bam) | samtools view -hSu /dev/stdin | samtools '
            'sort -@ {2} /dev/stdin > {5}/{1}.raw.bam'.format(date, sample, 
                                threads, bwa_ref, ' '.join(fqfiles), bwadir))

    cmdargs = shlex.split(command)
    print(command)
    print(cmdargs)

    p = sp.Popen(cmdargs)

    return p

def bwamem_beta(arglist):
    """
    Create and run process for the new implementation of BWA-Mem alignment
    """
    print("Hello!")
    print(arglist)
    #/anaconda/bin/bwa mem   -c 250 -M -t 15  -R '@RG\tID:SHGSOC001N\tPL:illumina\tPU:1_20180823_SHGSOC\tSM:SHGSOC001N' -v 1 /gpfs/igmmfs01/eddie/bioinfsvice/ameynert/software/bcbio-1.0.7/genomes/Hsapiens/hg38/bwa/hg38.fa align_prep/SHGSOC001N_SHGSOC001N_R1.fastq.gz align_prep/SHGSOC001N_SHGSOC001N_R2.fastq.gz | /anaconda/bin/bamsormadup inputformat=sam threads=15 tmpfile=bcbiotx/tmpMwNtCb/SHGSOC001N-sort-sorttmp-markdup SO=coordinate indexfilename=bcbiotx/tmpMwNtCb/SHGSOC001N-sort.bam.bai > bcbiotx/tmpMwNtCb/SHGSOC001N-sort.bam

def samtools_index(infile):
    """
    Run samtools index on newly created, aligned bam files
    """
    command = str('samtools index {}'.format(infile))
    cmdargs = shlex.split(command)
    print(command)
    print(cmdargs)

    p = sp.Popen(cmdargs)

    return p

if __name__ == "__main__":
    description = ("This script runs the BWA-Mem step of Bioexcel Align")
    args = au.parse_command_line(description)
    if args.bwa_version == 'stable':
        pbwa = bwamem_stable(args.bwa_ind_ref, args.threads, args.date,                 args.sample, args.files, args.bwadir)
        pbwa.wait()
    # elif args.bwa_version == 'beta':
    #     pbwa = bwamem_beta(args.bwa_ind_ref, args.threads, args.date, 
    #               args.sample, args.fqfiles, args.bwadir)
    #     pbwa.wait()

    else: 
        sys.exit('No BWA Mem version selected')

    psamidx = samtools_index('{}.raw.bam'.format(args.sample))
    psamidx.wait()
    
    