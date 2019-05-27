#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:24:22 2019

@author: Okada
"""

import subprocess

if not "pysam" in locals():
    import pysam
    
def sort_file(input_file, output_file, sort_option):
    
    subprocess.check_call("sort %s %s > %s" % (sort_option, input_file, output_file), shell=True)

def tabix_compress_index(input_file, output_file, seq_col, start_col, end_col):
    
    pysam.tabix_compress(input_file, output_file, force = True)
    pysam.tabix_index(output_file, force = True, seq_col = seq_col, start_col = start_col, end_col = end_col)

def bed_from_tabix(input_file, output_file, wide):
    
    def __output(f, chrom, start, end, num):
        f.write("%s\t%d\t%d\t%d\n" % (chrom, start, end, num))
        
    tbx_sa = pysam.TabixFile(input_file)
    f_out = open(output_file, "w")
    
    chrom = ""
    start = 0
    end = 0
    member = 0
    
    for row in tbx_sa.fetch():
        tpl = row.split("\t")
        if len(tpl) < 2:
            continue

        current_chr = tpl[0]
        current_start = int(tpl[1]) - 1
        current_end = int(tpl[1])
        
        if chrom != current_chr:
            if chrom != "":
                __output(f_out, chrom, start, end, member)
            chrom = current_chr
            start = current_start
            end = current_end
            member = 1
            continue
        if current_end - start > wide:
            __output(f_out, chrom, start, end, member)
            start = current_start
            end = current_end
            member = 1
            continue
     
        if current_end > end:
            end = current_end
            member += 1
            continue
        member += 1
        
    if chrom != "":
        __output(f_out, chrom, start, end, member)
    
    f_out.close()
    tbx_sa.close()
    
if __name__ == "__main__":
    pass
