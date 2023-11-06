from __future__ import division
import pandas as pd
import sys, argparse, os
import fileinput
import re

# input arguments
parser = argparse.ArgumentParser(description="this script is to extract mendelian concordance information")

parser.add_argument('-LCL5', '--LCL5', type=str, help='LCL5 family info',  required=True)
parser.add_argument('-LCL6', '--LCL6', type=str, help='LCL6 family info',  required=True)
parser.add_argument('-genotype', '--genotype', type=str, help='Genotype information of a set of four family members',  required=True)
parser.add_argument('-family', '--family', type=str, help='family name',  required=True)


args = parser.parse_args()
lcl5 = args.LCL5
lcl6 = args.LCL6
genotype = args.genotype
family = args.family


# output file
family_name = family + '.txt'

family_file = open(family_name,'w')

# input files
lcl5_dat = pd.read_table(lcl5)
lcl6_dat = pd.read_table(lcl6)
genotype_dat = pd.read_table(genotype)
merged_df = pd.merge(lcl5_dat, lcl6_dat,  how='outer', left_on=['#CHROM','POS'], right_on = ['#CHROM','POS'])
merged_genotype_df = pd.merge(merged_df, genotype_dat,  how='outer', left_on=['#CHROM','POS'], right_on = ['#CHROM','POS'])

merged_genotype_df_sub = merged_genotype_df.iloc[:,[0,1,23,24,29,30,31,32,7,17]]
merged_genotype_df_sub.columns = ['CHROM', 'POS', 'REF', 'ALT','LCL5','LCL6','LCL7','LCL8', 'TRIO5', 'TRIO6']

for row in merged_genotype_df_sub.itertuples():
	# sister
	if row.LCL5 == row.LCL6:
		if row.LCL5 == './.':
			mendelian = 'noInfo'
			sister_count = "no"
		elif row.LCL5 == '0/0':
			mendelian = 'Ref'
			sister_count = "no"
		else:
			mendelian = '1'
			sister_count = "yes_same"
	else:
		mendelian = '0'
		if (row.LCL5 == './.' or row.LCL5 == '0/0') and (row.LCL6 == './.' or row.LCL6 == '0/0'):
			sister_count = "no"
		else:
			sister_count = "yes_diff"
	# family trio5
	if row.LCL5 == row. LCL7 == row.LCL8 == './.':
		mendelian = mendelian + ':noInfo'
	elif row.LCL5 == row. LCL7 == row.LCL8 == '0/0':
		mendelian = mendelian + ':Ref'
	elif pd.isnull(row.TRIO5) == True:
		mendelian = mendelian + ':unVBT'
	else:
		mendelian = mendelian + ':' + row.TRIO5.split('=')[1]
	# family trio6
	if row.LCL6 == row. LCL7 == row.LCL8 == './.':
		mendelian = mendelian + ':noInfo'
	elif row.LCL6 == row. LCL7 == row.LCL8 == '0/0':
		mendelian = mendelian + ':Ref'
	elif pd.isnull(row.TRIO6) == True:
		mendelian = mendelian + ':unVBT'
	else:
		mendelian =  mendelian + ':' + row.TRIO6.split('=')[1]
	# not count into family
	if (row.LCL5 == './.' or row.LCL5 == '0/0') and (row.LCL6 == './.' or row.LCL6 == '0/0') and (row.LCL7 == './.' or row.LCL7 == '0/0') and (row.LCL8 == './.' or row.LCL8 == '0/0'):
		mendelian_count = "no"
	else:
		mendelian_count = "yes"
	outline = row.CHROM + '\t' + str(row.POS) + '\t' + row.REF + '\t' + row.ALT + '\t' + row.LCL5 + '\t' + row.LCL6 + '\t' + row.LCL7 + '\t' + row.LCL8 + '\t' + str(row.TRIO5) + '\t' + str(row.TRIO6) + '\t' + str(mendelian) + '\t' + str(mendelian_count) + '\t' + str(sister_count) + '\n'
	family_file.write(outline)
