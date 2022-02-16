#!/bin/sh

# Look up subdirectories and concatenate reads
#
# Usage:
#   ./scripts/file_ops.sh <path> <out_dir> <out_name>
#
# Example:
#   ./scripts/file_ops.sh /home/user/data/ /home/user/data/out/
#
# Input:
#   <path> - path to directory containing subdirectories
#   <out_dir> - path to output directory
#   <out_name> - name of output file
#
# Output:
#   <out_dir>/<out_name> - concatenated reads
#
# Dependencies:
#   - samtools
#   - awk
#   - cat
#
# Assumptions:
#   - <path> is a directory
#   - <out_dir> is a directory
#   - <out_name> is a string



