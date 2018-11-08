#!/usr/bin/env python

"""
This script runs GATK BaseRecalibrator and/or ApplyBQSR. The script opens a
GATK process with the correct parameters.
"""

import shlex
import subprocess as sp
import bioexcel_seqqc.alignutils as au

if __name__ == "__main__":
    description = ("This script runs GATK BaseRecalibrator and/or ApplyBQSR")
    args = au.parse_command_line(description)

    args.files = au.get_files(args)

    ## DECIDE ON LOGIC FOR RUNNING GATK STAGES