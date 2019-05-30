#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:25:44 2019

@author: aokada
"""

import sys
import argparse
#from ecsub import __version__
__version__ = "0.0.1"

import parser
import fetch_bp

def main():

    PROG = "nanomon-sv"
    ap = argparse.ArgumentParser(prog = PROG)
    ap.add_argument("--version", action = "version", version = PROG + "-" + __version__)
    subparsers = ap.add_subparsers()
    
    ##########
    # parse
    parse_parser = subparsers.add_parser("parse", help = "parse junction")
    parse_parser.add_argument("-i", "--input_bam", metavar = "path/to/input.bam", help = "output temporary data", type = str)
    parse_parser.add_argument("-o", "--output_prefix", metavar = "path/to/output-dir/", help = "output-directory or prefix", type = str)

    parse_parser.add_argument("--min_major_clip_size", metavar = 100, help = "(count as clipping) >= min_major_clip_size (bps)", type = int, default = 100)
    parse_parser.add_argument("--insertion_thres", metavar = 20, help = "(count as insertion) >= insertion_thres", type = int, default = 20)
    parse_parser.add_argument("--deletion_thres", metavar = 50, help = "(count as deletion) >= deletion_thres", type = int, default = 50)
    
    # mappQ
    parse_parser.add_argument("--min_mapping_qual", metavar = 0, help = "no use", type = int, default = 0) # TODO: 1
    # Primary Read: large clippoing at both ends, use or not?
    parse_parser.add_argument("--max_minor_clip_size", metavar = 0, help = "no use", type = int, default = 0) # TODO: 100
    # SA Read: small mapping size, use or not?
    parse_parser.add_argument("--min_sa_mapping_size", metavar = 0, help = "no use", type = int, default = 0) # TODO: 500
    
    parse_parser.set_defaults(func = parser.main)
    
    ##########
    # fetch break-point
    fetch_parser = subparsers.add_parser("fetch", help = "parse junction")
    fetch_parser.add_argument("-i", "--input_bp", metavar = "path/to/breakpoints.txt", help = "input file", type = str)
    fetch_parser.add_argument("-t", "--parsed_file_tumor", metavar = "path/to/tumor/xxx.junction.sort.gz", help = "parsed file", type = str)
    fetch_parser.add_argument("-n", "--parsed_file_normal", metavar = "path/to/normal/xxx.junction.sort.gz", help = "parsed file", type = str)
    fetch_parser.add_argument("-o", "--output_prefix", metavar = "path/to/output-dir/xxx.fetch", help = "output-directory or prefix", type = str)
    fetch_parser.add_argument("--margin", metavar = 30, help = "fetch margin size (bps)", type = int, default = 30)
    
    fetch_parser.set_defaults(func = fetch_bp.main)
    
    argv = sys.argv[1:]
    if len(argv) < 1:
        argv = [""]
        
    args = ap.parse_args(argv)
    
    return args.func(args)
    
if __name__ == "__main__":
    print (">>> " + " ".join(sys.argv))
    sys.exit(main())
