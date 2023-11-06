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

parser.add_argument('-folder', '--folder', type=str, help='directory that holds all the mendelian info',  required=True)
parser.add_argument('-vcf', '--vcf', type=str, help='merged multiple sample vcf',  required=True)


args = parser.parse_args()
folder = args.folder
vcf = args.vcf

# input files
folder = folder + '/*.txt'
filenames = glob(folder)
dataframes = []
for filename in filenames:
	dataframes.append(pd.read_table(filename,header=None))

dfs = [df.set_index([0, 1]) for df in dataframes]
merged_mendelian = pd.concat(dfs, axis=1).reset_index()
family_name = [i.split('/')[-1].replace('.txt','') for i in filenames]
columns = ['CHROM','POS'] + family_name
merged_mendelian.columns = columns

vcf_dat = pd.read_table(vcf)

merged_df = pd.merge(merged_mendelian, vcf_dat,  how='outer', left_on=['CHROM','POS'], right_on = ['#CHROM','POS'])
merged_df = merged_df.fillna('nan')

vcf_header = '''##fileformat=VCFv4.2
##fileDate=20200501
##source=high_confidence_calls_intergration(choppy app)
##reference=GRCh38.d1.vd1
##INFO=<ID=VOTED,Number=1,Type=Integer,Description="Number mendelian consisitent votes">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sum depth of all samples">
##FORMAT=<ID=ALT,Number=1,Type=Integer,Description="Sum alternative depth of all samples">
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
# output files
benchmark_LCL5 = open('LCL5_voted.vcf','w')
benchmark_LCL6 = open('LCL6_voted.vcf','w')
benchmark_LCL7 = open('LCL7_voted.vcf','w')
benchmark_LCL8 = open('LCL8_voted.vcf','w')

all_sample_outfile = open('all_sample_information.txt','w')

# write VCF
LCL5_col = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLCL5_benchmark_calls\n'
LCL6_col = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLCL6_benchmark_calls\n'
LCL7_col = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLCL7_benchmark_calls\n'
LCL8_col = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tLCL8_benchmark_calls\n'

benchmark_LCL5.write(vcf_header)
benchmark_LCL5.write(LCL5_col)
benchmark_LCL6.write(vcf_header)
benchmark_LCL6.write(LCL6_col)
benchmark_LCL7.write(vcf_header)
benchmark_LCL7.write(LCL7_col)
benchmark_LCL8.write(vcf_header)
benchmark_LCL8.write(LCL8_col)

# all info
all_info_col = 'CHROM\tPOS\tLCL5_pcr_consensus\tLCL5_pcr_free_consensus\tLCL5_mendelian_num\tLCL5_consensus_call\tLCL5_consensus_alt_seq\tLCL5_alt\tLCL5_dp\tLCL5_detected_num\tLCL6_pcr_consensus\tLCL6_pcr_free_consensus\tLCL6_mendelian_num\tLCL6_consensus_call\tLCL6_consensus_alt_seq\tLCL6_alt\tLCL6_dp\tLCL6_detected_num\tLCL7_pcr_consensus\tLCL7_pcr_free_consensus\tLCL7_mendelian_num\t LCL7_consensus_call\tLCL7_consensus_alt_seq\tLCL7_alt\tLCL7_dp\tLCL7_detected_num\tLCL8_pcr_consensus\tLCL8_pcr_free_consensus\tLCL8_mendelian_num\tLCL8_consensus_call\tLCL8_consensus_alt_seq\tLCL8_alt\tLCL8_dp\tLCL8_detected_num\n'
all_sample_outfile.write(all_info_col)

# function
def decide_by_rep(vcf_list,mendelian_list):
	consensus_rep = ''
	gt = [x.split(':')[0] for x in vcf_list]
	# mendelian consistent?
	mendelian_dict = Counter(mendelian_list)
	highest_mendelian = mendelian_dict.most_common(1)
	candidate_mendelian = highest_mendelian[0][0]
	freq_mendelian = highest_mendelian[0][1]
	if (candidate_mendelian == '1:1:1') and (freq_mendelian >= 2):
		con_loc = [i for i in range(len(mendelian_list)) if mendelian_list[i] == '1:1:1']
		gt_con = itemgetter(*con_loc)(gt)
		gt_num_dict = Counter(gt_con)
		highest_gt = gt_num_dict.most_common(1)
		candidate_gt = highest_gt[0][0]
		freq_gt = highest_gt[0][1]
		if (candidate_gt != './.') and (freq_gt >= 2):
			consensus_rep = candidate_gt
		elif (candidate_gt == './.') and (freq_gt >= 2):
			consensus_rep = 'noGTInfo'
		else:
			consensus_rep = 'inconGT'
	elif (candidate_mendelian == 'nan') and (freq_mendelian >= 2):
		consensus_rep = 'noMenInfo'
	else:
		consensus_rep = 'inconMen'
	return consensus_rep

def consensus_call(vcf_info_list,mendelian_list,alt_seq):
	pcr_consensus = '.'
	pcr_free_consensus = '.'
	mendelian_num = '.'
	consensus_call = '.'
	consensus_alt_seq = '.'
	# pcr
	SEQ2000 = decide_by_rep(vcf_info_list[0:3],mendelian_list[0:3])
	XTen_ARD = decide_by_rep(vcf_info_list[18:21],mendelian_list[18:21])
	XTen_NVG = decide_by_rep(vcf_info_list[21:24],mendelian_list[21:24])
	XTen_WUX = decide_by_rep(vcf_info_list[24:27],mendelian_list[24:27])
	pcr_sequence_site = [SEQ2000,XTen_ARD,XTen_NVG,XTen_WUX]
	pcr_sequence_dict = Counter(pcr_sequence_site)
	pcr_highest_sequence = pcr_sequence_dict.most_common(1)
	pcr_candidate_sequence = pcr_highest_sequence[0][0]
	pcr_freq_sequence = pcr_highest_sequence[0][1]
	if pcr_freq_sequence > 2:
		pcr_consensus = pcr_candidate_sequence
	else:
		pcr_consensus = 'inconSequenceSite'
	# pcr-free
	T7_WGE = decide_by_rep(vcf_info_list[3:6],mendelian_list[3:6])
	Nova_ARD_1 = decide_by_rep(vcf_info_list[6:9],mendelian_list[6:9])
	Nova_ARD_2 = decide_by_rep(vcf_info_list[9:12],mendelian_list[9:12])
	Nova_BRG = decide_by_rep(vcf_info_list[12:15],mendelian_list[12:15])
	Nova_WUX = decide_by_rep(vcf_info_list[15:18],mendelian_list[15:18])
	sequence_site = [T7_WGE,Nova_ARD_1,Nova_ARD_2,Nova_BRG,Nova_WUX]
	sequence_dict = Counter(sequence_site)
	highest_sequence = sequence_dict.most_common(1)
	candidate_sequence = highest_sequence[0][0]
	freq_sequence = highest_sequence[0][1]
	if freq_sequence > 3:
		pcr_free_consensus = candidate_sequence
	else:
		pcr_free_consensus = 'inconSequenceSite'
	# net alt, dp
	# alt
	AD = [x.split(':')[1] for x in vcf_info_list]
	ALT = [x.split(',')[1] for x in AD]
	ALT = [int(x) for x in ALT]
	ALL_ALT = sum(ALT)
	# dp
	DP = [x.split(':')[2] for x in vcf_info_list]
	DP = [int(x) for x in DP]
	ALL_DP = sum(DP)
	# detected number
	gt = [x.split(':')[0] for x in vcf_info_list]
	gt = [x.replace('0/0','.') for x in gt]
	gt = [x.replace('./.','.') for x in gt]
	detected_num = 27 - gt.count('.')
	# decide consensus calls
	tag = ['inconGT','noMenInfo','inconMen','inconSequenceSite','noGTInfo']
	if (pcr_consensus != '0/0') and (pcr_consensus == pcr_free_consensus) and (pcr_consensus not in tag):
		consensus_call = pcr_consensus
		gt = [x.split(':')[0] for x in vcf_info_list]
		indices = [i for i, x in enumerate(gt) if x == consensus_call]
		matched_mendelian = itemgetter(*indices)(mendelian_list)
		mendelian_num = matched_mendelian.count('1:1:1')
		# Delete multiple alternative genotype to necessary expression
		alt_gt = alt_seq.split(',')
		if len(alt_gt) > 1:
			allele1 = consensus_call.split('/')[0]
			allele2 = consensus_call.split('/')[1]
			if allele1 == '0':
				allele2_seq = alt_gt[int(allele2) - 1]
				consensus_alt_seq = allele2_seq
				consensus_call = '0/1'
			else:
				allele1_seq = alt_gt[int(allele1) - 1]
				allele2_seq = alt_gt[int(allele2) - 1]
				if int(allele1) > int(allele2):
					consensus_alt_seq = allele2_seq + ',' + allele1_seq
					consensus_call = '1/2'
				elif int(allele1) < int(allele2):
					consensus_alt_seq = allele1_seq + ',' + allele2_seq
					consensus_call = '1/2'
				else:
					consensus_alt_seq = allele1_seq 
					consensus_call = '1/1'
		else:
			consensus_alt_seq = alt_seq
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
	return pcr_consensus, pcr_free_consensus, mendelian_num, consensus_call, consensus_alt_seq, ALL_ALT, ALL_DP, detected_num


for row in merged_df.itertuples():
	mendelian_list = [row.Quartet_DNA_BGI_SEQ2000_BGI_1_20180518,row.Quartet_DNA_BGI_SEQ2000_BGI_2_20180530,row.Quartet_DNA_BGI_SEQ2000_BGI_3_20180530, \
					row.Quartet_DNA_BGI_T7_WGE_1_20191105,row.Quartet_DNA_BGI_T7_WGE_2_20191105,row.Quartet_DNA_BGI_T7_WGE_3_20191105, \
					row.Quartet_DNA_ILM_Nova_ARD_1_20181108,row.Quartet_DNA_ILM_Nova_ARD_2_20181108,row.Quartet_DNA_ILM_Nova_ARD_3_20181108, \
					row.Quartet_DNA_ILM_Nova_ARD_4_20190111,row.Quartet_DNA_ILM_Nova_ARD_5_20190111,row.Quartet_DNA_ILM_Nova_ARD_6_20190111, \
					row.Quartet_DNA_ILM_Nova_BRG_1_20180930,row.Quartet_DNA_ILM_Nova_BRG_2_20180930,row.Quartet_DNA_ILM_Nova_BRG_3_20180930, \
					row.Quartet_DNA_ILM_Nova_WUX_1_20190917,row.Quartet_DNA_ILM_Nova_WUX_2_20190917,row.Quartet_DNA_ILM_Nova_WUX_3_20190917, \
					row.Quartet_DNA_ILM_XTen_ARD_1_20170403,row.Quartet_DNA_ILM_XTen_ARD_2_20170403,row.Quartet_DNA_ILM_XTen_ARD_3_20170403, \
					row.Quartet_DNA_ILM_XTen_NVG_1_20170329,row.Quartet_DNA_ILM_XTen_NVG_2_20170329,row.Quartet_DNA_ILM_XTen_NVG_3_20170329, \
					row.Quartet_DNA_ILM_XTen_WUX_1_20170216,row.Quartet_DNA_ILM_XTen_WUX_2_20170216,row.Quartet_DNA_ILM_XTen_WUX_3_20170216]
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
	LCL5_pcr_consensus, LCL5_pcr_free_consensus, LCL5_mendelian_num, LCL5_consensus_call, LCL5_consensus_alt_seq, LCL5_alt, LCL5_dp, LCL5_detected_num = consensus_call(lcl5_list,mendelian_list,row.ALT)
	if LCL5_mendelian_num != '.':
		LCL5_output = row.CHROM + '\t' + str(row.POS) + '\t' + '.' + '\t' + row.REF + '\t' + LCL5_consensus_alt_seq + '\t' + '.' + '\t' + '.' + '\t' +'VOTED=' + str(LCL5_mendelian_num) + '\t' + 'GT:ALT:DP' + '\t' + LCL5_consensus_call + ':' + str(LCL5_alt) + ':' + str(LCL5_dp) +  '\n'
		benchmark_LCL5.write(LCL5_output)
	# LCL6
	LCL6_pcr_consensus, LCL6_pcr_free_consensus, LCL6_mendelian_num, LCL6_consensus_call, LCL6_consensus_alt_seq, LCL6_alt, LCL6_dp, LCL6_detected_num = consensus_call(lcl6_list,mendelian_list,row.ALT)
	if LCL6_mendelian_num != '.':
		LCL6_output = row.CHROM + '\t' + str(row.POS) + '\t' + '.' + '\t' + row.REF + '\t' + LCL6_consensus_alt_seq + '\t' + '.' + '\t' + '.' + '\t' +'VOTED=' + str(LCL6_mendelian_num) + '\t' + 'GT:ALT:DP' + '\t' + LCL6_consensus_call + ':' + str(LCL6_alt) + ':' + str(LCL6_dp) +  '\n'
		benchmark_LCL6.write(LCL6_output)
	# LCL7
	LCL7_pcr_consensus, LCL7_pcr_free_consensus, LCL7_mendelian_num, LCL7_consensus_call, LCL7_consensus_alt_seq, LCL7_alt, LCL7_dp, LCL7_detected_num = consensus_call(lcl7_list,mendelian_list,row.ALT)
	if LCL7_mendelian_num != '.':
		LCL7_output = row.CHROM + '\t' + str(row.POS) + '\t' + '.' + '\t' + row.REF + '\t' + LCL7_consensus_alt_seq + '\t' + '.' + '\t' + '.' + '\t' +'VOTED=' + str(LCL7_mendelian_num) + '\t' + 'GT:ALT:DP' + '\t' + LCL7_consensus_call + ':' + str(LCL7_alt) + ':' + str(LCL7_dp) +  '\n'
		benchmark_LCL7.write(LCL7_output)
	# LCL8
	LCL8_pcr_consensus, LCL8_pcr_free_consensus, LCL8_mendelian_num, LCL8_consensus_call, LCL8_consensus_alt_seq, LCL8_alt, LCL8_dp, LCL8_detected_num = consensus_call(lcl8_list,mendelian_list,row.ALT)
	if LCL8_mendelian_num != '.':
		LCL8_output = row.CHROM + '\t' + str(row.POS) + '\t' + '.' + '\t' + row.REF + '\t' + LCL8_consensus_alt_seq + '\t' + '.' + '\t' + '.' + '\t' +'VOTED=' + str(LCL8_mendelian_num) + '\t' + 'GT:ALT:DP' + '\t' + LCL8_consensus_call + ':' + str(LCL8_alt) + ':' + str(LCL8_dp) +  '\n'
		benchmark_LCL8.write(LCL8_output)
	# all data
	all_output = row.CHROM + '\t' + str(row.POS) + '\t' + LCL5_pcr_consensus + '\t' + LCL5_pcr_free_consensus + '\t' + str(LCL5_mendelian_num) + '\t' + LCL5_consensus_call + '\t' + LCL5_consensus_alt_seq + '\t' + str(LCL5_alt) + '\t' + str(LCL5_dp)  + '\t' + str(LCL5_detected_num) + '\t' +\
					LCL6_pcr_consensus + '\t' + LCL6_pcr_free_consensus + '\t' + str(LCL6_mendelian_num) + '\t' + LCL6_consensus_call + '\t' + LCL6_consensus_alt_seq + '\t' + str(LCL6_alt) + '\t' + str(LCL6_dp)  + '\t' + str(LCL6_detected_num) + '\t' +\
					LCL7_pcr_consensus + '\t' + LCL7_pcr_free_consensus + '\t' + str(LCL7_mendelian_num) + '\t' + LCL7_consensus_call + '\t' + LCL7_consensus_alt_seq + '\t' + str(LCL7_alt) + '\t' + str(LCL7_dp) + '\t' + str(LCL7_detected_num) + '\t' +\
					LCL8_pcr_consensus + '\t' + LCL8_pcr_free_consensus + '\t' + str(LCL8_mendelian_num) + '\t' + LCL8_consensus_call + '\t' + LCL8_consensus_alt_seq + '\t' + str(LCL8_alt) + '\t' + str(LCL8_dp) + '\t' + str(LCL8_detected_num) + '\n'
	all_sample_outfile.write(all_output)

