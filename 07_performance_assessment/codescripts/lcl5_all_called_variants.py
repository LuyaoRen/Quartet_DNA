from __future__ import division
import sys, argparse, os
import pandas as pd
from collections import Counter

# input arguments
parser = argparse.ArgumentParser(description="this script is to merge mendelian and vcfinfo, and extract high_confidence_calls")
parser.add_argument('-vcf', '--vcf', type=str, help='merged multiple sample vcf',  required=True)


args = parser.parse_args()
vcf = args.vcf

lcl5_outfile = open('LCL5_all_variants.txt','w')
filtered_outfile = open('LCL5_filtered_variants.txt','w')

vcf_dat = pd.read_table(vcf)


for row in vcf_dat.itertuples():
	lcl5_list = [row.Quartet_DNA_BGI_SEQ2000_BGI_LCL5_1_20180518,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL5_2_20180530,row.Quartet_DNA_BGI_SEQ2000_BGI_LCL5_3_20180530, \
					row.Quartet_DNA_BGI_T7_WGE_LCL5_1_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL5_2_20191105,row.Quartet_DNA_BGI_T7_WGE_LCL5_3_20191105, \
					row.Quartet_DNA_ILM_Nova_ARD_LCL5_1_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL5_2_20181108,row.Quartet_DNA_ILM_Nova_ARD_LCL5_3_20181108, \
					row.Quartet_DNA_ILM_Nova_ARD_LCL5_4_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL5_5_20190111,row.Quartet_DNA_ILM_Nova_ARD_LCL5_6_20190111, \
					row.Quartet_DNA_ILM_Nova_BRG_LCL5_1_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL5_2_20180930,row.Quartet_DNA_ILM_Nova_BRG_LCL5_3_20180930, \
					row.Quartet_DNA_ILM_Nova_WUX_LCL5_1_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL5_2_20190917,row.Quartet_DNA_ILM_Nova_WUX_LCL5_3_20190917, \
					row.Quartet_DNA_ILM_XTen_ARD_LCL5_1_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL5_2_20170403,row.Quartet_DNA_ILM_XTen_ARD_LCL5_3_20170403, \
					row.Quartet_DNA_ILM_XTen_NVG_LCL5_1_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL5_2_20170329,row.Quartet_DNA_ILM_XTen_NVG_LCL5_3_20170329, \
					row.Quartet_DNA_ILM_XTen_WUX_LCL5_1_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL5_2_20170216,row.Quartet_DNA_ILM_XTen_WUX_LCL5_3_20170216]
	lcl5_vcf_gt = [x.split(':')[0] for x in lcl5_list]
	lcl5_gt=[item.replace('./.', '0/0') for item in lcl5_vcf_gt]
	gt_dict = Counter(lcl5_gt)
	highest_gt = gt_dict.most_common(1)
	candidate_gt = highest_gt[0][0]
	freq_gt = highest_gt[0][1]
	output = row._1 + '\t' + str(row.POS) + '\t' + '\t'.join(lcl5_gt) + '\n'
	if (candidate_gt == '0/0') and (freq_gt == 27):
		filtered_outfile.write(output)
	else:
		lcl5_outfile.write(output)


