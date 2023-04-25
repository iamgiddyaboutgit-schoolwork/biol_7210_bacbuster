#!/usr/bin/env python3
"""
Fix problems with paired-end sequencing reads stored in .fastq format.
"""

import sys
import argparse

import dnaio

def main():
    user_args = get_user_args()
    fixer(
        in_read_1=user_args.i1,
        in_read_2=user_args.i2
    )

def get_user_args() -> argparse.Namespace:
    """Get arguments from the command line and validate them."""
    parser = argparse.ArgumentParser(description = "Perfom basic sanity checks on two .fastq files.")
    # Optional arguments are prefixed by single or double dashes.
    # The remaining arguments are positional.
    parser.add_argument("-i1", required = True, metavar="<Input file 1>",
        help="1st FASTQ file path")
    parser.add_argument("-i2", required = True, metavar="<Input file 2>",
        help="2nd FASTQ file path")
    parser.add_argument("-o1", required=True, metavar="<Output file 1>",
        help="Output path for 1st FASTQ file")
    parser.add_argument("-o2", required=True, metavar="<Output file 2>",
        help="Output path for 2nd FASTQ file")
    args = parser.parse_args()

    return args

def fixer(in_read_1, in_read_2):
    # https://dnaio.readthedocs.io/en/latest/tutorial.html#paired-end-data
    
    with dnaio.open(
        in_read_1, 
        in_read_2, 
        fileformat="fastq",
        mode="r",
        qualities=True,
        open_threads=2
    ) as reader:
        bp = 0
        for r1, r2 in reader:
            bp += len(r1) + len(r2)
        print(f"The paired-end input contains {bp/1E6:.1f} Mbp")

###############################################################################
# if this module is being run directly and not imported
if __name__ == "__main__":
    # Call the function main.
    main()
