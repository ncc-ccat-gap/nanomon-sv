#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 11:38:37 2019

@author: aokada
"""

import os
import time
import datetime

if not "pysam" in locals():
    import pysam

def fetch_breakpoint(input_bp, src_bp_tbx, output_bp, margin = [0,0], header = False, comment = "#"):

    def __fetch_tabix(tbx, bp1_chr, bp1_pos, bp1_dir, bp2_chr, bp2_pos, bp2_dir, margin):
        
        margin1 = margin[0]
        margin2 = margin[1]
        
        excepts = []
        try:
            itr = tbx.fetch(reference = bp1_chr, start = bp1_pos-1 - margin1, end = bp1_pos + margin1)
        except Exception:
            return excepts
        
        for row in itr:
            tpl = row.split("\t")
            if len(tpl) < 10:
                continue
            
            src1_chr = tpl[0]
            src1_pos = int(tpl[2])
            src1_dir = tpl[8]
            src2_chr = tpl[3]
            src2_pos = int(tpl[5])
            src2_dir = tpl[9]
            qname = tpl[6]
            
            if bp1_dir != src1_dir:
                continue
            if bp2_dir != src2_dir:
                continue
            if bp2_chr != src2_chr:
                continue
            
            d1 = src1_pos - bp1_pos
            d2 = src2_pos - bp2_pos
            
            if d2 < d1 - margin2:
                continue
            if d2 > d1 + margin2:
                continue
            
            excepts.append([src1_chr, str(src1_pos), src1_dir, 
                            src2_chr, str(src2_pos), src2_dir, 
                            str(d1), str(d2), qname]
                          + tpl[10:])
                
        return excepts
    
    tbx = pysam.TabixFile(src_bp_tbx)
    f_bp = open(output_bp, "w")
    for row in open(input_bp).readlines():
        if row.startswith(comment):
            continue
        
        items = row.rstrip("\n").split("\t")
        if len(items) < 8:
            continue
        
        if header:
            f_bp.write(("{header}\t{vtype}\tsrc1_chr\tsrc1_pos\tsrc1_dir\tsrc2_chr\tsrc2_pos\tsrc2_dir\td1\td2"
                        + "\tqname\tvtype\tH1\tM\tD\tI\tH2\tH1'\tM'\tD'\tI'\tH2"
                        + "\tRL\tRL'\tHcl-RL'\tHcl'-RL\n").format(
                header = "\t".join(items[0:6]),
                vtype = items[7]
            ))
            header = False
            continue
        
        bp1_chr = items[0]
        bp1_pos = int(items[1])
        bp1_dir = items[2]
        bp2_chr = items[3]
        bp2_pos = int(items[4])
        bp2_dir = items[5]
        vtype = items[7]
        
        find1 = __fetch_tabix(tbx, bp1_chr, bp1_pos, bp1_dir, bp2_chr, bp2_pos, bp2_dir, margin)
        find2 = __fetch_tabix(tbx, bp2_chr, bp2_pos, bp2_dir, bp1_chr, bp1_pos, bp1_dir, margin)
        find = find1 + find2
        
        if len(find) == 0:
            f_bp.write(("{data}\t{vtype}\t{find}\n").format(
                data = "\t".join(items[0:6]),
                vtype = vtype,
                find = ""
            ))
            continue
        
        for i in find:
            f_bp.write(("{data}\t{vtype}\t{find}\n").format(
                data = "\t".join(items[0:6]),
                vtype = vtype,
                find = "\t".join(i)
            ))
    f_bp.close()
    tbx.close()

def fetch_tabix_debug(junction_text, margin, bp_chr, bp_pos, bp_chr2 = None, bp_pos2 = None):
    def __is_match(src_chr, src_pos, bp_chr, bp_pos, margin):
        if src_chr == bp_chr and src_pos >= (bp_pos - margin) and src_pos <= (bp_pos + margin):
            return True
        return False
    
    def __to_str(items):
        return ",".join([items[0], items[2], items[8], items[3], items[5], items[9], items[6], items[10]])
    
    excepts = []
    f = open(junction_text)
    for row in f.readlines():
        items = row.rstrip("\n").split("\t")
        if len (items) < 10: continue
        if __is_match(items[0], int(items[2]), bp_chr, bp_pos, margin):
            if bp_chr2 != None:
                if __is_match(items[3], int(items[5]), bp_chr2, bp_pos2, margin):
                    excepts.append(__to_str(items))
            else:
                excepts.append(__to_str(items))
                
        elif __is_match(items[3], int(items[5]), bp_chr, bp_pos, margin):
            if bp_chr2 != None:
                if __is_match(items[0], int(items[2]), bp_chr2, bp_pos2, margin):
                    excepts.append(__to_str(items))
            else:
                excepts.append(__to_str(items))
            
    f.close()

    import pprint            
    pprint.pprint(excepts)

def __count_breakpoint(fetched, counted, header = False, comment = "#"):
    
    f_bp = open(counted, "w")
    
    long_reads = 0
    last_key = ""
    data = {}
    for row in open(fetched).readlines():
        if row.startswith(comment):
            continue
        
        items = row.rstrip("\n").split("\t")
        if len(items) < 1:
            continue
        
        if header:
            f_bp.write("Chr_1\tPos_1\tDir_1\tChr_2\tPos_2\tDir_2\tVariant_Type\tReadCount\n")
            header = False
            continue
        
        bp_key = "\t".join(items[0:7])
        if last_key == "":
            last_key = bp_key
            
        elif bp_key != last_key:
            f_bp.write("{key}\t{j}\n".format(key = last_key, j = long_reads))
            data[last_key] = long_reads
            
            long_reads = 0
            last_key = bp_key
            
        if len(items) > 8:
            long_reads += 1
    
    f_bp.write("{key}\t{j}\n".format(key = last_key, j = long_reads))
    data[last_key] = long_reads
    
    f_bp.close()
    
    return data

def count_breakpoint(input_bp, fetched_tumor, fetched_normal, counted, header = False, comment = "#"):
    
    data_tumor = __count_breakpoint(fetched_tumor, counted + ".tumor", header = header, comment = comment)
    data_normal = __count_breakpoint(fetched_normal, counted + ".normal", header = header, comment = comment)
                       
    f_bp = open(counted, "w")
    
    for row in open(input_bp).readlines():
        if row.startswith(comment):
            continue
        
        items = row.rstrip("\n").split("\t")
        if len(items) < 1:
            continue
        
        if header:
            f_bp.write("{header}\ttumor_reads\tnormal_reads\n".format(header = "\t".join(items)))
            header = False
            continue
        
        bp_key = "\t".join(items[0:6] + [items[7]])
        f_bp.write("{data}\t{i}\t{j}\n".format(data = "\t".join(items), i = data_tumor[bp_key], j = data_normal[bp_key]))
        
    f_bp.close()
    
def main(args):
    
    os.makedirs(os.path.dirname(args.output_prefix), exist_ok = True)
    
    if not os.path.exists(args.parsed_file_tumor + ".tbi"):
        print("%s.tbi is not exists." % (args.parsed_file_tumor))
        return 1
    
    if not os.path.exists(args.parsed_file_normal + ".tbi"):
        print("%s.tbi is not exists." % (args.parsed_file_normal))
        return 1
    
    t = time.time()
    print ("[%s] (%.3f) start parse_junction_from_bam." % (datetime.datetime.now(), 0))
    fetch_breakpoint(args.input_bp, 
                     args.parsed_file_tumor, 
                     output_bp = args.output_prefix + ".tumor.txt", 
                     margin = [args.margin, args.margin], 
                     header = True, 
                     comment = "#")
    
    fetch_breakpoint(args.input_bp, 
                     args.parsed_file_normal, 
                     output_bp = args.output_prefix + ".normal.txt", 
                     margin = [args.margin, args.margin], 
                     header = True, 
                     comment = "#")
    
    print ("[%s] (%.3f) end fetch_breakpoint." % (datetime.datetime.now(), (time.time()-t)))
    count_breakpoint(args.input_bp,
                     fetched_tumor = args.output_prefix + ".tumor.txt", 
                     fetched_normal = args.output_prefix + ".normal.txt", 
                     counted = args.output_prefix + ".count.txt",
                     header = True,
                     comment = "#")
    
    print ("[%s] (%.3f) end count_breakpoint." % (datetime.datetime.now(), (time.time()-t)))
    
    return 0

if __name__ == "__main__":
    pass
    
    #__count_breakpoint("./output/genomon.fetch3.normal.txt", "./output/genomon.fetch4.count.txt")
                       