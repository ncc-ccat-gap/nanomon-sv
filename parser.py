#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:24:22 2019

@author: Okada
"""

import re
import datetime
import time
import os
import json

if not "pysam" in locals():
    import pysam

import utils

def _text_to_list(text):
    
    code_list = re.split('[0-9]+', text)
    num_list = re.split('[^0-9]', text)
    if len(num_list) != len(code_list):
        return None
    
    li = []
    for i in range(len(num_list)-1):
        li.append([code_list[i+1], int(num_list[i])])
    return li

def _sum_cigar_label (data, label):
    l = 0    
    for c in data:
        if c[0] in label:
            l += c[1]
    return l

def _sum_M (data):
    return _sum_cigar_label (data, ["M", 0])

def _sum_D (data):
    return _sum_cigar_label (data, ["D", 2])

def _sum_I (data):
    return _sum_cigar_label (data, ["I", 1])

def _sum_cigar(data):

    m = 0
    d = 0
    i = 0
    for c in data:
        if c[0] == 0:   m += c[1]
        elif c[0] == 2: d += c[1]
        elif c[0] == 1: i += c[1]
        
    s1 = (data[0][1] if data[0][0] in [4, 5] else 0)
    s2 = (data[-1][1] if data[-1][0] in [4, 5] else 0)
    return [s1, m, d, i, s2]

def _cigar_to_juction(cigar, min_major_clip_size, max_minor_clip_size):
    # get the clipping size in the both side
    # M 	BAM_CMATCH 	0
    # I 	BAM_CINS 	1
    # D 	BAM_CDEL 	2
    # N 	BAM_CREF_SKIP 	3
    # S 	BAM_CSOFT_CLIP 	4
    # H 	BAM_CHARD_CLIP 	5
    # P 	BAM_CPAD 	6
    # = 	BAM_CEQUAL 	7
    # X 	BAM_CDIFF 	8
    # B 	BAM_CBACK 	9
    
    left_clipping = (cigar[0][1] if cigar[0][0] in [4, 5] else 0)
    right_clipping = (cigar[-1][1] if cigar[-1][0] in [4, 5] else 0)

    if left_clipping < min_major_clip_size and right_clipping < min_major_clip_size:
        return None

    # for comparing with the previous script (this will removed soon)
#    if left_clipping >= max_minor_clip_size and right_clipping >= max_minor_clip_size:
#        return None
    
    return {"left_clipping": left_clipping, "right_clipping": right_clipping}

# _pop_sa_tags(read.get_tags(), read.qname)
    
def _pop_sa_tag(tags, qname, sa_reads):
    
    sa_list = []
    if not qname in sa_reads:
        return sa_list
    
    for tag in tags:
        if tag[0] != "SA": continue
        
        for sa_str in tag[1].split(";"):
            
            sa = sa_str.split(",")
            if len(sa) < 6: continue
        
            sa_chr = sa[0]
            sa_position = int(sa[1])
            #sa_dir_strand = sa[2]
            #sa_cigar = sa[3]
            #sa_mq = int(sa[4])
            #sa_nm = int(sa[5])
            
            if sa_chr == "hs37d5": continue
            
            for r in sa_reads[qname]:
                if r["chr"] == sa_chr and r["reference"] == sa_position:
                    sa_list.append(r)

    return sa_list

def parse_junction_from_bam(input_bam, 
                         output_file_junction, 
                         min_mapping_qual, 
                         min_major_clip_size, 
                         max_minor_clip_size, 
                         min_sa_mapping_size,
                         insertion_thres,
                         deletion_thres):

    """
    This function utilizes the SA tags (SA:Z:rname, pos, strand, CIGAR, mapQ, number of mismatch).
    The strand of the supplementary alignment information in the SA tag is determined by the orignal sequence (before taking complement).
    Therefore, please not that when the primary alignment is in the reverse direction, the sequence shown in the bam file does not match
    to the SA tags..
    """
    junction_template = (
        "{bp1_chr}\t{bp1_pos1}\t{bp1_pos2}"
    + "\t{bp2_chr}\t{bp2_pos1}\t{bp2_pos2}"
    + "\t{rname}.{junction_index}"
    + "\t{score}"
    + "\t{bp1_dir}\t{bp2_dir}"
    + "\t{vtype}"
    + "\t{cigar_s1_a}\t{cigar_m_a}\t{cigar_d_a}\t{cigar_i_a}\t{cigar_s2_a}"
    + "\t{cigar_s1_b}\t{cigar_m_b}\t{cigar_d_b}\t{cigar_i_b}\t{cigar_s2_b}"
    + "\t{cigar_rl_a}\t{cigar_rl_b}\t{scl_a_rl_b}\t{scl_b_rl_a}"
    + "\n"
    )

    t = time.time()
    print ("[%s] (%.3f) start parse_junction_from_bam." % (datetime.datetime.now(), 0))
    
    if not os.path.exists(output_file_junction + ".sa.json"):
        bamfile = pysam.AlignmentFile(input_bam, "rb")
        sa_reads = {}
        for read in bamfile.fetch():
    
            if not read.is_supplementary: continue
            if len(read.cigar) == 1: continue
            if (int(read.mapq) < min_mapping_qual): continue
            if read.reference_name == "hs37d5": continue
        
            [sum_s1, sum_m, sum_d, sum_i, sum_s2] = _sum_cigar (read.cigar)
            if sum_s1 < min_major_clip_size and sum_s2 < min_major_clip_size: continue
        
            #[sum_m, sum_d, sum_i] = _sum_cigar (read.cigar)
            #clip = _cigar_to_juction(read.cigar, min_major_clip_size, max_minor_clip_size)
            #if clip == None: continue
            
            if not read.qname in sa_reads:
                sa_reads[read.qname] = []
    
            if sum_s2 >= min_major_clip_size:
                sa_reads[read.qname].append({
                    "chr": read.reference_name,
                    "reference": int(read.reference_start) + 1,
                    "pos": int(read.reference_end),
                    "dir": "+",
                    "cigar": {
                        "s1": sum_s1,
                        "s2": sum_s2,
                        "m": sum_m,
                        "d": sum_d,
                        "i": sum_i,
                        "Scl": sum_s2,
                        "RL": sum_m + sum_i + sum_s1
                    }
                })
            if sum_s1 >= min_major_clip_size:
                sa_reads[read.qname].append({
                    "chr": read.reference_name,
                    "reference": int(read.reference_start) + 1,
                    "pos": int(read.reference_start) + 1,
                    "dir": "-",
                    "cigar": {
                        "s1": sum_s1,
                        "s2": sum_s2,
                        "m": sum_m,
                        "d": sum_d,
                        "i": sum_i,
                        "Scl": sum_s1,
                        "RL": sum_m + sum_i + sum_s2
                    }
                })
                    
            
        bamfile.close()
        json.dump(sa_reads, open(output_file_junction + ".sa.json", "w"))
    
    else:
        sa_reads = json.load(open(output_file_junction + ".sa.json"))
        
    print ("[%s] (%.3f) loaded supplymentaly reads." % (datetime.datetime.now(), (time.time()-t)))

    bamfile = pysam.AlignmentFile(input_bam, "rb")
    f_junction = open(output_file_junction, "w")
    
    # maybe add the regional extraction of bam files
    for read in bamfile.fetch():

        # skip supplementary alignment
        if read.is_secondary: continue
        if read.is_supplementary: continue
        if len(read.cigar) == 1: continue
        if (int(read.mapq) < min_mapping_qual): continue
        if read.reference_name == "hs37d5": continue
        
        junction_index = 0
        
        # ---------------
        # insertion deletion
        # ---------------
        if True:
            indel_pos = int(read.reference_start) + 1
            cigar_index = 0
            for c in read.cigar:
                if c[0] in [0, 2]:
                    indel_pos += c[1]
                    
                vtype = ""
                pos1 = 0
                pos2 = 0
                if c[0] == 1 and c[1] >= insertion_thres:
                    vtype = "insertion"
                    pos1 = indel_pos
                    pos2 = indel_pos
                    
                elif c[0] == 2 and c[1] >= deletion_thres:
                    vtype = "deletion"
                    pos1 = indel_pos - c[1]
                    pos2 = indel_pos

                if vtype != "":                
                    junction_index += 1
                    [sum_s1_a, sum_m_a, sum_d_a, sum_i_a, sum_s2_a] = _sum_cigar (read.cigar[0:cigar_index])
                    [sum_s1_b, sum_m_b, sum_d_b, sum_i_b, sum_s2_b] = _sum_cigar (read.cigar[cigar_index+1:])
                    
                    f_junction.write(junction_template.format(
                        bp1_chr = read.reference_name,
                        bp1_pos1 = pos1-1,
                        bp1_pos2 = pos1,
                        bp2_chr = read.reference_name,
                        bp2_pos1 = pos2-1,
                        bp2_pos2 = pos2,
                        rname = read.qname,
                        junction_index = junction_index,
                        score = c[1],
                        bp1_dir = "+",
                        bp2_dir = "-",
                        vtype = vtype,
                        cigar_s1_a = sum_s1_a,
                        cigar_m_a = sum_m_a,
                        cigar_d_a = sum_d_a,
                        cigar_i_a = sum_i_a,
                        cigar_s2_a = sum_s2_a,
                        cigar_s1_b = sum_s1_b,
                        cigar_m_b = sum_m_b,
                        cigar_d_b = sum_d_b,
                        cigar_i_b = sum_i_b,
                        cigar_s2_b = sum_s2_b,
                        cigar_rl_a = "",
                        cigar_rl_b = "",
                        scl_a_rl_b = "",
                        scl_b_rl_a = ""
                    ))
                cigar_index += 1
                
        # --------------------
        # clipping
        # --------------------
        if True:
            [sum_s1, sum_m, sum_d, sum_i, sum_s2] = _sum_cigar (read.cigar)
            if sum_s1 < min_major_clip_size and sum_s2 < min_major_clip_size: continue

            #clip = _cigar_to_juction(read.cigar, min_major_clip_size, max_minor_clip_size)
            #if clip == None: continue
        
            sa_list = _pop_sa_tag(read.get_tags(), read.qname, sa_reads)
            
            for r in sa_list:

                if sum_s2 >= min_major_clip_size:
                    junction_index += 1
                    
                    p_pos = int(read.reference_end)
                    p_dir_clip = "+"
                    scl_a = sum_s2
                    sop_a = sum_s1
                    rl_a = sum_m + sum_i + sop_a
                    
                    f_junction.write(junction_template.format(
                        bp1_chr = read.reference_name,
                        bp1_pos1 = p_pos-1,
                        bp1_pos2 = p_pos,
                        bp2_chr = r["chr"],
                        bp2_pos1 = r["pos"]-1,
                        bp2_pos2 = r["pos"],
                        rname = read.qname,
                        junction_index = junction_index,
                        score = 0,
                        bp1_dir = p_dir_clip,
                        bp2_dir = r["dir"],
                        vtype = "clipping",
                        cigar_s1_a = sum_s1,
                        cigar_s2_a = sum_s2,
                        cigar_m_a = sum_m,
                        cigar_d_a = sum_d,
                        cigar_i_a = sum_i,
                        cigar_s1_b = r["cigar"]["s1"],
                        cigar_s2_b = r["cigar"]["s2"],
                        cigar_m_b = r["cigar"]["m"],
                        cigar_d_b = r["cigar"]["d"],
                        cigar_i_b = r["cigar"]["i"],
                        cigar_rl_a = rl_a,
                        cigar_rl_b = r["cigar"]["RL"],
                        scl_a_rl_b = scl_a - r["cigar"]["RL"],
                        scl_b_rl_a = r["cigar"]["Scl"] - rl_a
                    ))
                    
                if sum_s1 >= min_major_clip_size:
                    junction_index += 1
                    
                    p_pos = int(read.reference_start) + 1
                    p_dir_clip = "-"
                    scl_a = sum_s1
                    sop_a = sum_s2
                    rl_a = sum_m + sum_i + sop_a
                    
                    f_junction.write(junction_template.format(
                        bp1_chr = read.reference_name,
                        bp1_pos1 = p_pos-1,
                        bp1_pos2 = p_pos,
                        bp2_chr = r["chr"],
                        bp2_pos1 = r["pos"]-1,
                        bp2_pos2 = r["pos"],
                        rname = read.qname,
                        junction_index = junction_index,
                        score = 0,
                        bp1_dir = p_dir_clip,
                        bp2_dir = r["dir"],
                        vtype = "clipping",
                        cigar_s1_a = sum_s1,
                        cigar_s2_a = sum_s2,
                        cigar_m_a = sum_m,
                        cigar_d_a = sum_d,
                        cigar_i_a = sum_i,
                        cigar_s1_b = r["cigar"]["s1"],
                        cigar_s2_b = r["cigar"]["s2"],
                        cigar_m_b = r["cigar"]["m"],
                        cigar_d_b = r["cigar"]["d"],
                        cigar_i_b = r["cigar"]["i"],
                        cigar_rl_a = rl_a,
                        cigar_rl_b = r["cigar"]["RL"],
                        scl_a_rl_b = scl_a - r["cigar"]["RL"],
                        scl_b_rl_a = r["cigar"]["Scl"] - rl_a
                    ))
        
    bamfile.close()
    f_junction.close()
    print ("[%s] (%.3f) end parse_junction_from_bam." % (datetime.datetime.now(), (time.time()-t)))
   
def main(args):
    
    os.makedirs(os.path.dirname(args.output_prefix), exist_ok = True)
    
    parse_junction_from_bam(args.input_bam, 
                            args.output_prefix + "junction.unsort.txt", 
                            args.min_mapping_qual, 
                            args.min_major_clip_size, 
                            args.max_minor_clip_size, 
                            args.min_sa_mapping_size, 
                            args.insertion_thres, 
                            args.deletion_thres)
    
    utils.sort_file(args.output_prefix + "junction.unsort.txt", 
              args.output_prefix + "junction.sort.txt", 
              "-k1,1 -k2,2n -k4,4 -k5,5n")
    
    utils.tabix_compress_index(args.output_prefix + "junction.sort.txt", 
                         args.output_prefix + "junction.sort.gz", 
                         seq_col = 0, 
                         start_col = 1, 
                         end_col = 1)

    return 0

if __name__ == "__main__":
    pass
