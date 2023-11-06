from __future__ import division 
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
parser = argparse.ArgumentParser(description="this script is to count voting number")

parser.add_argument('-vcf', '--multi_sample_vcf', type=str, help='The VCF file you want to count the voting number',  required=True)
parser.add_argument('-dup', '--dup_list', type=str, help='Duplication list',  required=True)
parser.add_argument('-sample', '--sample_name', type=str, help='which sample of quartet',  required=True)
parser.add_argument('-prefix', '--prefix', type=str, help='Prefix of output file name',  required=True)

args = parser.parse_args()
multi_sample_vcf = args.multi_sample_vcf
dup_list = args.dup_list
prefix = args.prefix
sample_name = args.sample_name

vcf_header = '''##fileformat=VCFv4.2
##fileDate=20200331
##source=high_confidence_calls_intergration(choppy app)
##reference=GRCh38.d1.vd1
##INFO=<ID=location,Number=1,Type=String,Description="Repeat region">
##INFO=<ID=DETECTED,Number=1,Type=Integer,Description="Number of detected votes">
##INFO=<ID=VOTED,Number=1,Type=Integer,Description="Number of consnesus votes">
##INFO=<ID=FAM,Number=1,Type=Integer,Description="Number mendelian consisitent votes">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sum depth of all samples">
##FORMAT=<ID=ALT,Number=1,Type=Integer,Description="Sum alternative depth of all samples">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency, sum alternative depth / sum depth">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Average genotype quality">
##FORMAT=<ID=QD,Number=1,Type=Float,Description="Average Variant Confidence/Quality by Depth">
##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Average mapping quality">
##FORMAT=<ID=FS,Number=1,Type=Float,Description="Average Phred-scaled p-value using Fisher's exact test to detect strand bias">
##FORMAT=<ID=QUALI,Number=1,Type=Float,Description="Average variant quality">
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

vcf_header_all_sample = '''##fileformat=VCFv4.2
##fileDate=20200331
##reference=GRCh38.d1.vd1
##INFO=<ID=location,Number=1,Type=String,Description="Repeat region">
##INFO=<ID=DUP,Number=1,Type=Flag,Description="Duplicated variant records">
##INFO=<ID=DETECTED,Number=1,Type=Integer,Description="Number of detected votes">
##INFO=<ID=VOTED,Number=1,Type=Integer,Description="Number of consnesus votes">
##INFO=<ID=FAM,Number=1,Type=Integer,Description="Number mendelian consisitent votes">
##INFO=<ID=ALL_ALT,Number=1,Type=Float,Description="Sum of alternative reads of all samples">
##INFO=<ID=ALL_DP,Number=1,Type=Float,Description="Sum of depth of all samples">
##INFO=<ID=ALL_AF,Number=1,Type=Float,Description="Allele frequency of net alternatice reads and net depth">
##INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description="Mean of genotype quality of all samples">
##INFO=<ID=QD_MEAN,Number=1,Type=Float,Description="Average Variant Confidence/Quality by Depth">
##INFO=<ID=MQ_MEAN,Number=1,Type=Float,Description="Mean of mapping quality of all samples">
##INFO=<ID=FS_MEAN,Number=1,Type=Float,Description="Average Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=QUAL_MEAN,Number=1,Type=Float,Description="Average variant quality">
##INFO=<ID=PCR,Number=1,Type=String,Description="Consensus of PCR votes">
##INFO=<ID=PCR_FREE,Number=1,Type=String,Description="Consensus of PCR-free votes">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Consensus calls">
##INFO=<ID=CONSENSUS_SEQ,Number=1,Type=String,Description="Consensus sequence">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=String,Description="Depth">
##FORMAT=<ID=ALT,Number=1,Type=Integer,Description="Alternative Depth">
##FORMAT=<ID=AF,Number=1,Type=String,Description="Allele frequency">
##FORMAT=<ID=GQ,Number=1,Type=String,Description="Genotype quality">
##FORMAT=<ID=MQ,Number=1,Type=String,Description="Mapping quality">
##FORMAT=<ID=TWINS,Number=1,Type=String,Description="1 is twins shared, 0 is twins discordant ">
##FORMAT=<ID=TRIO5,Number=1,Type=String,Description="1 is LCL7, LCL8 and LCL5 mendelian consistent, 0 is mendelian vioaltion">
##FORMAT=<ID=TRIO6,Number=1,Type=String,Description="1 is LCL7, LCL8 and LCL6 mendelian consistent, 0 is mendelian vioaltion">
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
# read in duplication list
dup = pd.read_table(dup_list,header=None)
var_dup = dup[0].tolist()

# output file
benchmark_file_name = prefix + '_voted.vcf'
benchmark_outfile = open(benchmark_file_name,'w')

all_sample_file_name = prefix + '_all_sample_information.vcf'
all_sample_outfile = open(all_sample_file_name,'w')

# write VCF
outputcolumn = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample_name + '_benchmark_calls\n'
benchmark_outfile.write(vcf_header)
benchmark_outfile.write(outputcolumn)

outputcolumn_all_sample = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+ \
'Quartet_DNA_BGI_SEQ2000_BGI_1_20180518\tQuartet_DNA_BGI_SEQ2000_BGI_2_20180530\tQuartet_DNA_BGI_SEQ2000_BGI_3_20180530\t' + \
'Quartet_DNA_BGI_T7_WGE_1_20191105\tQuartet_DNA_BGI_T7_WGE_2_20191105\tQuartet_DNA_BGI_T7_WGE_3_20191105\t' + \
'Quartet_DNA_ILM_Nova_ARD_1_20181108\tQuartet_DNA_ILM_Nova_ARD_2_20181108\tQuartet_DNA_ILM_Nova_ARD_3_20181108\t' + \
'Quartet_DNA_ILM_Nova_ARD_4_20190111\tQuartet_DNA_ILM_Nova_ARD_5_20190111\tQuartet_DNA_ILM_Nova_ARD_6_20190111\t' + \
'Quartet_DNA_ILM_Nova_BRG_1_20180930\tQuartet_DNA_ILM_Nova_BRG_2_20180930\tQuartet_DNA_ILM_Nova_BRG_3_20180930\t' + \
'Quartet_DNA_ILM_Nova_WUX_1_20190917\tQuartet_DNA_ILM_Nova_WUX_2_20190917\tQuartet_DNA_ILM_Nova_WUX_3_20190917\t' + \
'Quartet_DNA_ILM_XTen_ARD_1_20170403\tQuartet_DNA_ILM_XTen_ARD_2_20170403\tQuartet_DNA_ILM_XTen_ARD_3_20170403\t' + \
'Quartet_DNA_ILM_XTen_NVG_1_20170329\tQuartet_DNA_ILM_XTen_NVG_2_20170329\tQuartet_DNA_ILM_XTen_NVG_3_20170329\t' + \
'Quartet_DNA_ILM_XTen_WUX_1_20170216\tQuartet_DNA_ILM_XTen_WUX_2_20170216\tQuartet_DNA_ILM_XTen_WUX_3_20170216\n'
all_sample_outfile.write(vcf_header_all_sample)
all_sample_outfile.write(outputcolumn_all_sample)



#function
def replace_nan(strings_list):
	updated_list = []
	for i in strings_list:
		if i == '.':
			updated_list.append('.:.:.:.:.:.:.:.:.:.:.:.')
		else:
			updated_list.append(i)
	return updated_list

def remove_dot(strings_list):
	updated_list = []
	for i in strings_list:
		if i == '.':
			pass
		else:
			updated_list.append(i)
	return updated_list	

def detected_number(strings):
	gt = [x.split(':')[0] for x in strings]
	percentage = 27 - gt.count('.')
	return(str(percentage))

def vote_number(strings,consensus_call):
	gt = [x.split(':')[0] for x in strings]
	gt = [x.replace('.','0/0') for x in gt]
	gt = list(map(gt_uniform,[i for i in gt]))
	vote_num = gt.count(consensus_call)
	return(str(vote_num))

def family_vote(strings,consensus_call):
	gt = [x.split(':')[0] for x in strings]
	gt = [x.replace('.','0/0') for x in gt]
	gt = list(map(gt_uniform,[i for i in gt]))
	mendelian = [':'.join(x.split(':')[1:4]) for x in strings]
	indices = [i for i, x in enumerate(gt) if x == consensus_call]
	matched_mendelian = itemgetter(*indices)(mendelian)
	mendelian_num = matched_mendelian.count('1:1:1')
	return(str(mendelian_num))

def gt_uniform(strings):
	uniformed_gt = ''
	allele1 = strings.split('/')[0]
	allele2 = strings.split('/')[1]
	if int(allele1) > int(allele2):
		uniformed_gt = allele2 + '/' + allele1
	else:
		uniformed_gt = allele1 + '/' + allele2
	return uniformed_gt

def decide_by_rep(strings):
	consensus_rep = ''
	mendelian = [':'.join(x.split(':')[1:4]) for x in strings]
	gt = [x.split(':')[0] for x in strings]
	gt = [x.replace('.','0/0') for x in gt]
	# modified gt turn 2/1 to 1/2
	gt = list(map(gt_uniform,[i for i in gt]))
	# mendelian consistent?
	mendelian_dict = Counter(mendelian)
	highest_mendelian = mendelian_dict.most_common(1)
	candidate_mendelian = highest_mendelian[0][0]
	freq_mendelian = highest_mendelian[0][1]
	if (candidate_mendelian == '1:1:1') and (freq_mendelian >= 2):
		gt_num_dict = Counter(gt)
		highest_gt = gt_num_dict.most_common(1)
		candidate_gt = highest_gt[0][0]
		freq_gt = highest_gt[0][1]
		if (candidate_gt != '0/0') and (freq_gt >= 2):
			consensus_rep = candidate_gt
		elif (candidate_gt == '0/0') and (freq_gt >= 2):
			consensus_rep = '0/0'
		else:
			consensus_rep = 'inconGT'
	elif (candidate_mendelian == '') and (freq_mendelian >= 2):
		consensus_rep = 'noInfo'
	else:
		consensus_rep = 'inconMen'
	return consensus_rep


def main():
	for line in fileinput.input(multi_sample_vcf):
		headline = re.match('^\#',line)
		if headline is not None:
			pass
		else:
			line = line.strip()
			strings = line.split('\t')
			variant_id = '_'.join([strings[0],strings[1]])
			# check if the variants location is duplicated
			if variant_id in var_dup:
				strings[7] = strings[7] + ';DUP'
				outLine = '\t'.join(strings) + '\n'
				all_sample_outfile.write(outLine)
			else:
				# pre-define
				pcr_consensus = '.'
				pcr_free_consensus = '.'
				consensus_call = '.'
				consensus_alt_seq = '.'
				# pcr 
				strings[9:] = replace_nan(strings[9:])
				pcr = itemgetter(*[9,10,11,27,28,29,30,31,32,33,34,35])(strings)
				SEQ2000 = decide_by_rep(pcr[0:3])
				XTen_ARD = decide_by_rep(pcr[3:6])
				XTen_NVG = decide_by_rep(pcr[6:9])
				XTen_WUX = decide_by_rep(pcr[9:12])
				sequence_site = [SEQ2000,XTen_ARD,XTen_NVG,XTen_WUX]
				sequence_dict = Counter(sequence_site)
				highest_sequence = sequence_dict.most_common(1)
				candidate_sequence = highest_sequence[0][0]
				freq_sequence = highest_sequence[0][1]
				if freq_sequence > 2:
					pcr_consensus = candidate_sequence
				else:
					pcr_consensus = 'inconSequenceSite'
				# pcr-free
				pcr_free = itemgetter(*[12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])(strings)
				T7_WGE = decide_by_rep(pcr_free[0:3])
				Nova_ARD_1 = decide_by_rep(pcr_free[3:6])
				Nova_ARD_2 = decide_by_rep(pcr_free[6:9])
				Nova_BRG = decide_by_rep(pcr_free[9:12])
				Nova_WUX = decide_by_rep(pcr_free[12:15])
				sequence_site = [T7_WGE,Nova_ARD_1,Nova_ARD_2,Nova_BRG,Nova_WUX]
				highest_sequence = sequence_dict.most_common(1)
				candidate_sequence = highest_sequence[0][0]
				freq_sequence = highest_sequence[0][1]
				if freq_sequence > 3:
					pcr_free_consensus = candidate_sequence
				else:
					pcr_free_consensus = 'inconSequenceSite'
				# pcr and pcr-free
				tag = ['inconGT','noInfo','inconMen','inconSequenceSite']
				if (pcr_consensus == pcr_free_consensus) and (pcr_consensus not in tag) and (pcr_consensus != '0/0'):
					consensus_call = pcr_consensus
					VOTED = vote_number(strings[9:],consensus_call)
					strings[7] = strings[7] + ';VOTED=' + VOTED
					DETECTED = detected_number(strings[9:])
					strings[7] = strings[7] + ';DETECTED=' + DETECTED
					FAM = family_vote(strings[9:],consensus_call)
					strings[7] = strings[7] + ';FAM=' + FAM
					# Delete multiple alternative genotype to necessary expression
					alt = strings[4]
					alt_gt = alt.split(',')
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
						consensus_alt_seq = alt
					# GT:DP:ALT:AF:GQ:QD:MQ:FS:QUAL
					# GT:TWINS:TRIO5:TRIO6:DP:ALT:AF:GQ:QD:MQ:FS:QUAL:rawGT
					# DP
					DP = [x.split(':')[4] for x in strings[9:]]
					DP = remove_dot(DP)
					DP = [int(x) for x in DP]
					ALL_DP = sum(DP)
					# AF
					ALT = [x.split(':')[5] for x in strings[9:]]
					ALT = remove_dot(ALT)
					ALT = [int(x) for x in ALT]
					ALL_ALT = sum(ALT)
					ALL_AF = round(ALL_ALT/ALL_DP,2)
					# GQ
					GQ = [x.split(':')[7] for x in strings[9:]]
					GQ = remove_dot(GQ)
					GQ = [int(x) for x in GQ]
					GQ_MEAN = round(mean(GQ),2)
					# QD
					QD = [x.split(':')[8] for x in strings[9:]]
					QD = remove_dot(QD)
					QD = [float(x) for x in QD]
					QD_MEAN = round(mean(QD),2)
					# MQ
					MQ = [x.split(':')[9] for x in strings[9:]]
					MQ = remove_dot(MQ)
					MQ = [float(x) for x in MQ]
					MQ_MEAN = round(mean(MQ),2)
					# FS
					FS = [x.split(':')[10] for x in strings[9:]]
					FS = remove_dot(FS)
					FS = [float(x) for x in FS]
					FS_MEAN = round(mean(FS),2)
					# QUAL
					QUAL = [x.split(':')[11] for x in strings[9:]]
					QUAL = remove_dot(QUAL)
					QUAL = [float(x) for x in QUAL]
					QUAL_MEAN = round(mean(QUAL),2)
					# benchmark output 
					output_format = consensus_call + ':' + str(ALL_DP) + ':' + str(ALL_ALT) + ':' + str(ALL_AF) + ':' + str(GQ_MEAN) + ':' + str(QD_MEAN) + ':' + str(MQ_MEAN) + ':' + str(FS_MEAN) + ':' + str(QUAL_MEAN)
					outLine = strings[0] + '\t' + strings[1] + '\t' + strings[2] + '\t' + strings[3] + '\t' + consensus_alt_seq + '\t' + '.' + '\t' + '.' + '\t' + strings[7] + '\t' + 'GT:DP:ALT:AF:GQ:QD:MQ:FS:QUAL' + '\t' + output_format + '\n'
					benchmark_outfile.write(outLine)
					# all sample output
					strings[7] = strings[7] + ';ALL_ALT=' + str(ALL_ALT) + ';ALL_DP=' + str(ALL_DP) + ';ALL_AF=' + str(ALL_AF) \
					+ ';GQ_MEAN=' + str(GQ_MEAN) + ';QD_MEAN=' + str(QD_MEAN) + ';MQ_MEAN=' + str(MQ_MEAN) + ';FS_MEAN=' + str(FS_MEAN) \
					+ ';QUAL_MEAN=' + str(QUAL_MEAN) + ';PCR=' + consensus_call + ';PCR_FREE=' + consensus_call + ';CONSENSUS=' + consensus_call \
					+ ';CONSENSUS_SEQ=' + consensus_alt_seq
					all_sample_outLine = '\t'.join(strings) + '\n'
					all_sample_outfile.write(all_sample_outLine)
				elif (pcr_consensus in tag) and (pcr_free_consensus in tag):
					consensus_call = 'filtered'
					DETECTED = detected_number(strings[9:])
					strings[7] = strings[7] + ';DETECTED=' + DETECTED
					strings[7] = strings[7] + ';CONSENSUS=' + consensus_call
					all_sample_outLine = '\t'.join(strings) + '\n'
					all_sample_outfile.write(all_sample_outLine)					
				elif ((pcr_consensus == '0/0') or (pcr_consensus in tag)) and ((pcr_free_consensus not in tag) and (pcr_free_consensus != '0/0')):
					consensus_call = 'pcr-free-speicifc'
					DETECTED = detected_number(strings[9:])
					strings[7] = strings[7] + ';DETECTED=' + DETECTED
					strings[7] = strings[7] + ';CONSENSUS=' + consensus_call
					all_sample_outLine = '\t'.join(strings) + '\n'
					all_sample_outfile.write(all_sample_outLine)					
				elif ((pcr_consensus != '0/0') or (pcr_consensus not in tag)) and ((pcr_free_consensus in tag) and (pcr_free_consensus == '0/0')):
					consensus_call = 'pcr-speicifc'
					DETECTED = detected_number(strings[9:])
					strings[7] = strings[7] + ';DETECTED=' + DETECTED
					strings[7] = strings[7] + ';CONSENSUS=' + consensus_call + ';PCR=' + pcr_consensus + ';PCR_FREE=' + pcr_free_consensus
					all_sample_outLine = '\t'.join(strings) + '\n'
					all_sample_outfile.write(all_sample_outLine)					
				elif (pcr_consensus == '0/0') and (pcr_free_consensus == '0/0'):
					consensus_call = 'confirm for parents'				
					DETECTED = detected_number(strings[9:])
					strings[7] = strings[7] + ';DETECTED=' + DETECTED
					strings[7] = strings[7] + ';CONSENSUS=' + consensus_call
					all_sample_outLine = '\t'.join(strings) + '\n'
					all_sample_outfile.write(all_sample_outLine)					
				else:
					consensus_call = 'filtered'
					DETECTED = detected_number(strings[9:])
					strings[7] = strings[7] + ';DETECTED=' + DETECTED
					strings[7] = strings[7] + ';CONSENSUS=' + consensus_call
					all_sample_outLine = '\t'.join(strings) + '\n'
					all_sample_outfile.write(all_sample_outLine)					

if __name__ == '__main__':
	main()












