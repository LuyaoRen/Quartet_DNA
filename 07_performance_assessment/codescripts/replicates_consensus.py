from __future__ import division
from glob import glob
import sys, argparse, os
import fileinput
import re
import pandas as pd
from operator import itemgetter
from collections import Counter
from itertools import islice
from numpy import *
import statistics

# input arguments
parser = argparse.ArgumentParser(description="this script is to merge mendelian and vcfinfo, and extract high_confidence_calls")

parser.add_argument('-prefix', '--prefix', type=str, help='prefix of output file',  required=True)
parser.add_argument('-vcf', '--vcf', type=str, help='merged multiple sample vcf',  required=True)


args = parser.parse_args()
prefix = args.prefix
vcf = args.vcf

# input files
vcf_dat = pd.read_table(vcf)

# all info
all_file_name = prefix + "_all_summary.txt"
all_sample_outfile = open(all_file_name,'w')
all_info_col = 'CHROM\tPOS\tREF\tALT\tLCL5_consensus_calls\tLCL5_detect_number\tLCL5_same_diff\tLCL6_consensus_calls\tLCL6_detect_number\tLCl6_same_diff\tLCL7_consensus_calls\tLCL7_detect_number\tLCL7_same_diff\tLCL8_consensus_calls\tLCL8_detect_number\tLCL8_same_diff\n'
all_sample_outfile.write(all_info_col)

# filtered info
vcf_header = '''##fileformat=VCFv4.2
##fileDate=20200501
##source=high_confidence_calls_intergration(choppy app)
##reference=GRCh38.d1.vd1
##INFO=<ID=VOTED,Number=1,Type=Integer,Description="Number mendelian consisitent votes">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
'''
consensus_file_name = prefix + "_consensus.vcf"
consensus_outfile = open(consensus_file_name,'w')
consensus_col = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLCL5_consensus_call\tLCL6_consensus_call\tLCL7_consensus_call\tLCL8_consensus_call\n'
consensus_outfile.write(vcf_header)
consensus_outfile.write(consensus_col)

# function
def decide_by_rep(vcf_list):
	consensus_rep = ''
	gt = [x.split(':')[0] for x in vcf_list]
	gt_num_dict = Counter(gt)
	highest_gt = gt_num_dict.most_common(1)
	candidate_gt = highest_gt[0][0]
	freq_gt = highest_gt[0][1]
	if freq_gt >= 2:
		consensus_rep = candidate_gt
	else:
		consensus_rep = 'inconGT'
	return consensus_rep

def consensus_call(vcf_info_list):
	consensus_call = '.'
	detect_number = '.'
	same_diff = '.'
	# pcr
	SEQ2000 = decide_by_rep(vcf_info_list[0:3])
	XTen_ARD = decide_by_rep(vcf_info_list[18:21])
	XTen_NVG = decide_by_rep(vcf_info_list[21:24])
	XTen_WUX = decide_by_rep(vcf_info_list[24:27])
	Nova_WUX = decide_by_rep(vcf_info_list[15:18])
	pcr_sequence_site = [SEQ2000,XTen_ARD,XTen_NVG,XTen_WUX,Nova_WUX]
	pcr_sequence_dict = Counter(pcr_sequence_site)
	pcr_highest_sequence = pcr_sequence_dict.most_common(1)
	pcr_candidate_sequence = pcr_highest_sequence[0][0]
	pcr_freq_sequence = pcr_highest_sequence[0][1]
	if pcr_freq_sequence > 3:
		pcr_consensus = pcr_candidate_sequence
	else:
		pcr_consensus = 'inconSequenceSite'
	# pcr-free
	T7_WGE = decide_by_rep(vcf_info_list[3:6])
	Nova_ARD_1 = decide_by_rep(vcf_info_list[6:9])
	Nova_ARD_2 = decide_by_rep(vcf_info_list[9:12])
	Nova_BRG = decide_by_rep(vcf_info_list[12:15])
	sequence_site = [T7_WGE,Nova_ARD_1,Nova_ARD_2,Nova_BRG]
	sequence_dict = Counter(sequence_site)
	highest_sequence = sequence_dict.most_common(1)
	candidate_sequence = highest_sequence[0][0]
	freq_sequence = highest_sequence[0][1]
	if freq_sequence > 2:
		pcr_free_consensus = candidate_sequence
	else:
		pcr_free_consensus = 'inconSequenceSite'
	gt = [x.split(':')[0] for x in vcf_info_list]
	gt = [x.replace('./.','.') for x in gt]
	detected_num = 27 - gt.count('.')
	gt_remain = [e for e in gt if e not in {'.'}]
	gt_set = set(gt_remain)
	if len(gt_set) == 1:
		same_diff = 'same'
	else:
		same_diff = 'diff'
	tag = ['inconGT','inconSequenceSite']
	if (pcr_consensus == pcr_free_consensus) and (pcr_consensus not in tag):
		consensus_call = pcr_consensus
	elif (pcr_consensus in tag) and (pcr_free_consensus in tag):
		consensus_call = 'notAgree'				

	else:
		consensus_call = 'notConsensus'
	return consensus_call, detected_num, same_diff

	elif (pcr_consensus in tag) and (pcr_free_consensus in tag):
		consensus_call = 'filtered'				
	elif ((pcr_consensus == './.') or (pcr_consensus in tag)) and ((pcr_free_consensus not in tag) and (pcr_free_consensus != './.')):
		consensus_call = 'pcr-free-speicifc'				
	elif ((pcr_consensus != './.') or (pcr_consensus not in tag)) and ((pcr_free_consensus in tag) and (pcr_free_consensus == './.')):
		consensus_call = 'pcr-speicifc'			
	elif (pcr_consensus == '0/0') and (pcr_free_consensus == '0/0'):
		consensus_call = '0/0'								
	else:
		consensus_call = 'filtered'


for row in vcf_dat.itertuples():
# length
#alt
	if ',' in row.ALT:
		alt = row.ALT.split(',')
		alt_len = [len(i) for i in alt]
		alt_max = max(alt_len)
	else:
		alt_max = len(row.ALT)
#ref
	ref_len = len(row.REF)
	if (alt_max > 50) or (ref_len > 50):
		pass
	else:
# consensus
		lcl5_list = [row.Quartet_DNA_BGI_SEQ2000_BGI_LCL5_1_20180518,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL5_2_20180530,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL5_3_20180530, \
						row.Quartet_DNA_BGI_T7_WGE_LCL5_1_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL5_2_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL5_3_20191105, \
						row.Quartet_DNA_ILM_Nova_ARD_LCL5_1_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL5_2_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL5_3_20181108, \
						row.Quartet_DNA_ILM_Nova_ARD_LCL5_4_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL5_5_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL5_6_20190111, \
						row.Quartet_DNA_ILM_Nova_BRG_LCL5_1_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL5_2_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL5_3_20180930, \
						row.Quartet_DNA_ILM_Nova_WUX_LCL5_1_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL5_2_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL5_3_20190917, \
						row.Quartet_DNA_ILM_XTen_ARD_LCL5_1_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL5_2_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL5_3_20170403, \
						row.Quartet_DNA_ILM_XTen_NVG_LCL5_1_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL5_2_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL5_3_20170329, \
						row.Quartet_DNA_ILM_XTen_WUX_LCL5_1_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL5_2_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL5_3_20170216]
		lcl6_list = [row.Quartet_DNA_BGI_SEQ2000_BGI_LCL6_1_20180518,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL6_2_20180530,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL6_3_20180530, \
						row.Quartet_DNA_BGI_T7_WGE_LCL6_1_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL6_2_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL6_3_20191105, \
						row.Quartet_DNA_ILM_Nova_ARD_LCL6_1_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL6_2_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL6_3_20181108, \
						row.Quartet_DNA_ILM_Nova_ARD_LCL6_4_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL6_5_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL6_6_20190111, \
						row.Quartet_DNA_ILM_Nova_BRG_LCL6_1_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL6_2_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL6_3_20180930, \
						row.Quartet_DNA_ILM_Nova_WUX_LCL6_1_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL6_2_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL6_3_20190917, \
						row.Quartet_DNA_ILM_XTen_ARD_LCL6_1_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL6_2_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL6_3_20170403, \
						row.Quartet_DNA_ILM_XTen_NVG_LCL6_1_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL6_2_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL6_3_20170329, \
						row.Quartet_DNA_ILM_XTen_WUX_LCL6_1_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL6_2_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL6_3_20170216]
		lcl7_list = [row.Quartet_DNA_BGI_SEQ2000_BGI_LCL7_1_20180518,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL7_2_20180530,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL7_3_20180530, \
						row.Quartet_DNA_BGI_T7_WGE_LCL7_1_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL7_2_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL7_3_20191105, \
						row.Quartet_DNA_ILM_Nova_ARD_LCL7_1_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL7_2_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL7_3_20181108, \
						row.Quartet_DNA_ILM_Nova_ARD_LCL7_4_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL7_5_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL7_6_20190111, \
						row.Quartet_DNA_ILM_Nova_BRG_LCL7_1_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL7_2_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL7_3_20180930, \
						row.Quartet_DNA_ILM_Nova_WUX_LCL7_1_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL7_2_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL7_3_20190917, \
						row.Quartet_DNA_ILM_XTen_ARD_LCL7_1_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL7_2_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL7_3_20170403, \
						row.Quartet_DNA_ILM_XTen_NVG_LCL7_1_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL7_2_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL7_3_20170329, \
						row.Quartet_DNA_ILM_XTen_WUX_LCL7_1_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL7_2_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL7_3_20170216]
		lcl8_list = [row.Quartet_DNA_BGI_SEQ2000_BGI_LCL8_1_20180518,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL8_2_20180530,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL8_3_20180530, \
						row.Quartet_DNA_BGI_T7_WGE_LCL8_1_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL8_2_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL8_3_20191105, \
						row.Quartet_DNA_ILM_Nova_ARD_LCL8_1_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL8_2_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL8_3_20181108, \
						row.Quartet_DNA_ILM_Nova_ARD_LCL8_4_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL8_5_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL8_6_20190111, \
						row.Quartet_DNA_ILM_Nova_BRG_LCL8_1_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL8_2_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL8_3_20180930, \
						row.Quartet_DNA_ILM_Nova_WUX_LCL8_1_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL8_2_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL8_3_20190917, \
						row.Quartet_DNA_ILM_XTen_ARD_LCL8_1_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL8_2_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL8_3_20170403, \
						row.Quartet_DNA_ILM_XTen_NVG_LCL8_1_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL8_2_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL8_3_20170329, \
						row.Quartet_DNA_ILM_XTen_WUX_LCL8_1_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL8_2_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL8_3_20170216]
		# LCL5
		LCL5_consensus_call, LCL5_detected_num, LCL5_same_diff = consensus_call(lcl5_list)
		# LCL6
		LCL6_consensus_call, LCL6_detected_num, LCL6_same_diff = consensus_call(lcl6_list)
		# LCL7
		LCL7_consensus_call, LCL7_detected_num, LCL7_same_diff = consensus_call(lcl7_list)
		# LCL8
		LCL8_consensus_call, LCL8_detected_num, LCL8_same_diff = consensus_call(lcl8_list)
		# all data
		all_output = row._1 + '\t' + str(row.POS) + '\t' + row.REF + '\t' + row.ALT + '\t' +  LCL5_consensus_call + '\t' + str(LCL5_detected_num) + '\t' + LCL5_same_diff + '\t' +\
						 LCL6_consensus_call + '\t' + str(LCL6_detected_num) + '\t' + LCL6_same_diff + '\t' +\
						 LCL7_consensus_call + '\t' + str(LCL7_detected_num) + '\t' + LCL7_same_diff + '\t' +\
						 LCL8_consensus_call + '\t' + str(LCL8_detected_num) + '\t' + LCL8_same_diff + '\n'
		all_sample_outfile.write(all_output)
		#consensus vcf
		one_position = [LCL5_consensus_call,LCL6_consensus_call,LCL7_consensus_call,LCL8_consensus_call]
		if ('notConsensus' in one_position) or (((len(set(one_position)) == 1) and ('./.' in set(one_position))) or ((len(set(one_position)) == 1) and ('0/0' in set(one_position))) or ((len(set(one_position)) == 2) and ('0/0' in set(one_position) and ('./.' in set(one_position))))):
			pass
		else:
			consensus_output = row._1 + '\t' + str(row.POS) + '\t' + '.' + '\t' + row.REF + '\t' + row.ALT + '\t' + '.' + '\t' + '.' + '\t' +'.' + '\t' + 'GT' + '\t' + LCL5_consensus_call + '\t' + LCL6_consensus_call + '\t' + LCL7_consensus_call + '\t' + LCL8_consensus_call +'\n'
			consensus_outfile.write(consensus_output)








