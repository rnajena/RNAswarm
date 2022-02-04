#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

def run_samtools_flagstat(sam_file, flagstat_file):
    """
    Run samtools flagstat.
    """
    cmd = f'samtools flagstat {sam_file} > {flagstat_file}'
    print(cmd)
    os.system(cmd)

def get_samfiles(sam_dir):
    """
    Get sam files in sam_dir.
    """
    sam_files = []
    for root, dirs, files in os.walk(sam_dir):
        for file in files:
            if file.endswith(".sam"):
                sam_files.append(os.path.join(root, file))
    return sam_files

def run_samtools_flagstat_all(sam_dir, flagstat_dir):
    """
    Run samtools flagstat for all sam files in sam_dir.
    """
    sam_files = get_samfiles(sam_dir)
    if not os.path.exists(flagstat_dir):
        os.makedirs(flagstat_dir)
    for sam_file in sam_files:
        flagstat_file = sam_file.replace(sam_dir, flagstat_dir)
        flagstat_file = flagstat_file.replace(".sam", ".flagstat")
        run_samtools_flagstat(sam_file, flagstat_file)

run_samtools_flagstat_all('/home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/02-mappings/segemehl', '/home/ru27wav/Projects/gl_iav-splash_freiburg/results/dadonaite_2019/02-mappings/samtools_flagstat')