# import modules
import sys, argparse, os
import fileinput
import re

parser = argparse.ArgumentParser(description="This script is to split samples in VCF files and rewrite to the right style")

parser.add_argument('-vcf', '--familyVCF', type=str, help='VCF with sister and mendelian infomation',  required=True)
parser.add_argument('-name', '--familyName', type=str, help='Family name of the VCF file',  required=True)

args = parser.parse_args()

# Rename input:
inputFile = args.familyVCF
family_name = args.familyName

# output filename
LCL5_name = family_name + '.LCL5.vcf'
LCL5file = open(LCL5_name,'w')
LCL6_name = family_name + '.LCL6.vcf'
LCL6file = open(LCL6_name,'w')
LCL7_name = family_name + '.LCL7.vcf'
LCL7file = open(LCL7_name,'w')
LCL8_name = family_name + '.LCL8.vcf'
LCL8file = open(LCL8_name,'w')
family_filename = family_name + '.vcf'
familyfile = open(family_filename,'w')

# default columns, which will be included in the included in the calssifier
vcfheader = '''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="the same genotype between twin sister and mendelian consistent in 578 and 678">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=TWINS,Number=0,Type=Flag,Description="0 for sister consistent, 1 for sister inconsistent">
##FORMAT=<ID=TRIO5,Number=0,Type=Flag,Description="0 for trio consistent, 1 for trio inconsistent">
##FORMAT=<ID=TRIO6,Number=0,Type=Flag,Description="0 for trio consistent, 1 for trio inconsistent">
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
# write VCF
LCL5colname = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+family_name+'_LCL5'+'\n'
LCL5file.write(vcfheader)
LCL5file.write(LCL5colname)

LCL6colname = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+family_name+'_LCL6'+'\n'
LCL6file.write(vcfheader)
LCL6file.write(LCL6colname)

LCL7colname = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+family_name+'_LCL7'+'\n'
LCL7file.write(vcfheader)
LCL7file.write(LCL7colname)

LCL8colname = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+family_name+'_LCL8'+'\n'
LCL8file.write(vcfheader)
LCL8file.write(LCL8colname)

familycolname = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'LCL5\t'+'LCL6\t'+'LCL7\t'+'LCL8'+'\n'
familyfile.write(vcfheader)
familyfile.write(familycolname)

# reform VCF
def process(oneLine):
	line = oneLine.rstrip()
	strings = line.strip().split('\t')
	# replace .
	# LCL6 uniq
	if strings[11] == '.':
		strings[11] = '0/0'
		strings[9] = strings[12]
		strings[10] = strings[13]
	else:
		pass
	# LCL5 uniq
	if strings[14] == '.':
		strings[14] = '0/0'
		strings[12] = strings[9]
		strings[13] = strings[10]
	else:
		pass
	# sister
	if strings[11] == strings[14]:
		add_format = ":1"
	else:
		add_format = ":0"
	# trioLCL5
	if strings[15] == 'MD=1':
		add_format = add_format + ":1"
	else:
		add_format = add_format + ":0"
	# trioLCL6
	if strings[7] == 'MD=1':
		add_format = add_format + ":1"
	else:
		add_format = add_format + ":0"
	# filter
	if (strings[11] == strings[14]) and (strings[15] == 'MD=1') and (strings[7] == 'MD=1'):
		strings[6] = 'PASS'
	else:
		strings[6] = '.'
	# output LCL5
	LCL5outLine = strings[0]+'\t'+strings[1]+'\t'+strings[2]+'\t'+strings[3]+'\t'+strings[4]+'\t'+'.'+'\t'+strings[6]+'\t'+ '.' +'\t'+ 'GT:TWINS:TRIO5:TRIO6' + '\t' + strings[14] + add_format + '\n'
	LCL5file.write(LCL5outLine)
	# output LCL6
	LCL6outLine = strings[0]+'\t'+strings[1]+'\t'+strings[2]+'\t'+strings[3]+'\t'+strings[4]+'\t'+'.'+'\t'+strings[6]+'\t'+ '.' +'\t'+ 'GT:TWINS:TRIO5:TRIO6' + '\t' + strings[11] + add_format + '\n'
	LCL6file.write(LCL6outLine)
	# output LCL7
	LCL7outLine = strings[0]+'\t'+strings[1]+'\t'+strings[2]+'\t'+strings[3]+'\t'+strings[4]+'\t'+'.'+'\t'+strings[6]+'\t'+ '.' +'\t'+ 'GT:TWINS:TRIO5:TRIO6' + '\t' + strings[10] + add_format + '\n'
	LCL7file.write(LCL7outLine)
	# output LCL8
	LCL8outLine = strings[0]+'\t'+strings[1]+'\t'+strings[2]+'\t'+strings[3]+'\t'+strings[4]+'\t'+'.'+'\t'+strings[6]+'\t'+ '.' +'\t'+ 'GT:TWINS:TRIO5:TRIO6' + '\t' + strings[9] + add_format + '\n'
	LCL8file.write(LCL8outLine)
	# output family
	familyoutLine = strings[0]+'\t'+strings[1]+'\t'+strings[2]+'\t'+strings[3]+'\t'+strings[4]+'\t'+ '.'+'\t'+strings[6]+'\t'+ '.' +'\t'+ 'GT:TWINS:TRIO5:TRIO6' + '\t' + strings[14] + add_format +'\t' + strings[11] + add_format + '\t' + strings[10] + add_format +'\t' + strings[9] + add_format + '\n'
	familyfile.write(familyoutLine)


for line in fileinput.input(inputFile):
	m = re.match('^\#',line)
	if m is not None:
		pass
	else:
		process(line)
		

