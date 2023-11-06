from __future__ import division
import pandas as pd
import sys, argparse, os
import fileinput
import re

# input arguments
parser = argparse.ArgumentParser(description="this script is to get final high confidence calls and information of all replicates")

parser.add_argument('-vcfInfo', '--vcfInfo', type=str, help='The txt file of variants information, this file is named as prefix__variant_quality_location.txt',  required=True)
parser.add_argument('-mendelianInfo', '--mendelianInfo', type=str, help='The merged mendelian information of all samples',  required=True)
parser.add_argument('-sample', '--sample_name', type=str, help='which sample of quartet',  required=True)


args = parser.parse_args()
vcfInfo = args.vcfInfo
mendelianInfo = args.mendelianInfo
sample_name = args.sample_name


#GT:TWINS:TRIO5:TRIO6:DP:AF:GQ:QD:MQ:FS:QUAL

vcf_header = '''##fileformat=VCFv4.2
##fileDate=20200331
##source=high_confidence_calls_intergration(choppy app)
##reference=GRCh38.d1.vd1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=TWINS,Number=1,Type=Flag,Description="1 for sister consistent, 0 for sister different">
##FORMAT=<ID=TRIO5,Number=1,Type=Flag,Description="1 for LCL7, LCL8 and LCL5 mendelian consistent, 0 for family violation">
##FORMAT=<ID=TRIO6,Number=1,Type=Flag,Description="1 for LCL7, LCL8 and LCL6 mendelian consistent, 0 for family violation">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=ALT,Number=1,Type=Integer,Description="Alternative Depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##FORMAT=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##FORMAT=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="variant quality">
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

# output file
file_name = sample_name + '_mendelian_vcfInfo.vcf'
outfile = open(file_name,'w')

outputcolumn = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +  sample_name + '\n'
outfile.write(vcf_header)
outfile.write(outputcolumn)

# input files
vcf_info = pd.read_table(vcfInfo)
mendelian_info = pd.read_table(mendelianInfo)

merged_df = pd.merge(vcf_info, mendelian_info,  how='outer', left_on=['#CHROM','POS'], right_on = ['#CHROM','POS'])
merged_df = merged_df.fillna('.')

#
def parse_INFO(info):
	strings = info.strip().split(';')
	keys = []
	values = []
	for i in strings:
		kv = i.split('=')
		if kv[0] == 'DB':
			keys.append('DB')
			values.append('1')
		else:
			keys.append(kv[0])
			values.append(kv[1])
	infoDict = dict(zip(keys, values))
	return infoDict
#
for row in merged_df.itertuples():
	if row[18] != '.':
		# format
		# GT:TWINS:TRIO5:TRIO6:DP:AF:GQ:QD:MQ:FS:QUAL
		FORMAT_x = row[10].split(':')
		ALT = int(FORMAT_x[1].split(',')[1])
		if int(FORMAT_x[2]) != 0:
			AF = round(ALT/int(FORMAT_x[2]),2)
		else:
			AF = '.'
		INFO_x = parse_INFO(row.INFO_x)
		if FORMAT_x[2] == '0':
			INFO_x['QD'] = '.'
		else:
			pass
		FORMAT = row[18] + ':' + FORMAT_x[2] + ':' + str(ALT) + ':' + str(AF) + ':' + FORMAT_x[3] + ':' + INFO_x['QD'] + ':' + INFO_x['MQ'] + ':' + INFO_x['FS'] + ':' + str(row.QUAL_x)
		# outline
		outline = row._1 + '\t' + str(row.POS) + '\t' + row.ID_x + '\t' + row.REF_y + '\t' + row.ALT_y + '\t' + '.' + '\t' + '.' + '\t' + '.' + '\t' + 'GT:TWINS:TRIO5:TRIO6:DP:ALT:AF:GQ:QD:MQ:FS:QUAL' + '\t' + FORMAT + '\n'
	else:
		rawGT = row[10].split(':')
		FORMAT_x = row[10].split(':')
		ALT = int(FORMAT_x[1].split(',')[1])
		if int(FORMAT_x[2]) != 0:
			AF = round(ALT/int(FORMAT_x[2]),2)
		else:
			AF = '.'
		INFO_x = parse_INFO(row.INFO_x)
		if FORMAT_x[2] == '0':
			INFO_x['QD'] = '.'
		else:
			pass
		FORMAT = '.:.:.:.' + ':' + FORMAT_x[2] + ':' + str(ALT) + ':' + str(AF) + ':' + FORMAT_x[3] + ':' + INFO_x['QD'] + ':' + INFO_x['MQ'] + ':' + INFO_x['FS'] + ':' + str(row.QUAL_x) + ':' + rawGT[0]
		# outline
		outline = row._1 + '\t' + str(row.POS) + '\t' + row.ID_x + '\t' + row.REF_x + '\t' + row.ALT_x + '\t' + '.' + '\t' + '.' + '\t' + '.' + '\t' + 'GT:TWINS:TRIO5:TRIO6:DP:ALT:AF:GQ:QD:MQ:FS:QUAL:rawGT' + '\t' + FORMAT + '\n'
	outfile.write(outline)
