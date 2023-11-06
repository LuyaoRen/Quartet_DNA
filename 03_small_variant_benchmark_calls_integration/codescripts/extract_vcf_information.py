import sys,getopt
import os
import re
import fileinput
import pandas as pd

def usage():
	print(
		"""
Usage: python extract_vcf_information.py -i input_merged_vcf_file -o parsed_file

This script will extract SNVs and Indels information from the vcf files and output a tab-delimited files.

Input:
-i the selected vcf file

Output:
-o tab-delimited parsed file
		""")

# select supported small variants
def process(oneLine):
	line = oneLine.rstrip()
	strings = line.strip().split('\t')
	infoParsed = parse_INFO(strings[7])
	formatKeys = strings[8].split(':')
	formatValues = strings[9].split(':')
	for i in range(0,len(formatKeys) -1) :
		if formatKeys[i] == 'AD':
			ra = formatValues[i].split(',')
			infoParsed['RefDP'] = ra[0]
			infoParsed['AltDP'] = ra[1]
			if (int(ra[1]) + int(ra[0])) != 0:
				infoParsed['af'] = float(int(ra[1])/(int(ra[1]) + int(ra[0])))
			else:
				pass
		else:
			infoParsed[formatKeys[i]] = formatValues[i]
	infoParsed['chromo'] = strings[0]
	infoParsed['pos'] = strings[1]
	infoParsed['id'] = strings[2]
	infoParsed['ref'] = strings[3]
	infoParsed['alt'] = strings[4]
	infoParsed['qual'] = strings[5]
	return infoParsed


def parse_INFO(info):
	strings = info.strip().split(';')
	keys = []
	values = []
	for i in strings:
		kv = i.split('=')
		if kv[0] == 'DB':
			keys.append('DB')
			values.append('1')
		elif kv[0] == 'AF':
			pass
		else:
			keys.append(kv[0])
			values.append(kv[1])
	infoDict = dict(zip(keys, values))
	return infoDict
	

opts,args = getopt.getopt(sys.argv[1:],"hi:o:") 
for op,value in opts:
	if op == "-i":
		inputFile=value
	elif op == "-o":
		outputFile=value	
	elif op == "-h":
		usage()
		sys.exit()

if len(sys.argv[1:]) < 3:
	usage()
	sys.exit()

allDict = []
for line in fileinput.input(inputFile):
	m = re.match('^\#',line)
	if m is not None:
		pass
	else:
		oneDict = process(line)
		allDict.append(oneDict)

allTable = pd.DataFrame(allDict)

allTable.to_csv(outputFile,sep='\t',index=False)

