# import modules
import sys, argparse, os
import fileinput
from operator import itemgetter   

parser = argparse.ArgumentParser(description="this script is to vote callable bed region")

parser.add_argument('-bed', '--multiSampleBED', type=str, help='The bed file to get high confidence region',  required=True)
parser.add_argument('-prefix', '--prefix', type=str, help='The output file you want to name',  required=True)

args = parser.parse_args()

# Rename input:
input_file = args.multiSampleBED
prefix = args.prefix

consensus_filename = prefix + '.27consensus.bed'
outCONSENSUS = open(consensus_filename,'w')
filter_filename = prefix + '.filtered.bed'
outFiltered = open(filter_filename,'w')
#initial
#sequence_tech = ['chr','start','end','number','sample','SEQ2000','SEQ2000','SEQ2000','SEQT7','SEQT7','SEQT7','Nova','Nova','Nova','Nova','Nova','Nova','Nova','Nova','Nova','Nova','Nova','Nova','Nova','XTen','XTen','XTen','XTen','XTen','XTen','XTen','XTen','XTen']
#sequence_site = ['chr','start','end','number','sample','BGI','BGI','BGI','WGE','WGE','WGE','ARD','ARD','ARD','ARD','ARD','ARD','BRG','BRG','BRG','WUX','WUX','WUX','ARD','ARD','ARD','NVG','NVG','NVG','WUX','WUX','WUX']

def consensus_bed(oneLine):
	line = oneLine.strip()
	strings = line.split('\t')
	# replicate
	SEQ2000_BGI = 1 if strings[5:8].count('1') > 1 else 0
	T7_WGE = 1 if strings[8:11].count('1') > 1 else 0 
	Nova_ARD_1 = 1 if strings[11:14].count('1') > 1 else 0
	Nova_ARD_2 = 1 if strings[14:17].count('1') > 1 else 0		
	Nova_BRG = 1 if strings[17:20].count('1') > 1 else 0
	Nova_WUX = 1 if strings[20:23].count('1') > 1 else 0
	XTen_ARD = 1 if strings[23:26].count('1') >1 else 0
	XTen_NVG = 1 if strings[26:29].count('1') > 1 else 0
	XTen_WUX = 1 if strings[29:32].count('1') > 1 else 0
	# library
	pcr = 1 if [SEQ2000_BGI,Nova_WUX,XTen_ARD,XTen_WUX,XTen_NVG].count(1) > 3 else 0
	pcr_free = 1 if [T7_WGE,Nova_ARD_1,Nova_ARD_2,Nova_BRG].count(1) > 2 else 0
	voted = 1 if [pcr,pcr_free].count(1) > 1 else 0 
	# get consensus bed and tech specific bed
	if voted == 1:
		outCONSENSUS.write(oneLine)
	else:
		outFiltered.write(oneLine)

for oneLine in fileinput.input(input_file):
	consensus_bed(oneLine)


outCONSENSUS.close()
outFiltered.close()


