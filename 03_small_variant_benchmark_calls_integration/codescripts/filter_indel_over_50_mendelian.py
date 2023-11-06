from itertools import islice
import fileinput
import sys, argparse, os


# input arguments
parser = argparse.ArgumentParser(description="this script is to exclude indel over 50bp")

parser.add_argument('-i', '--mergedGVCF', type=str, help='merged gVCF txt with only chr, pos, ref, alt and genotypes',  required=True)
parser.add_argument('-prefix', '--prefix', type=str, help='prefix of output file',  required=True)


args = parser.parse_args()
input_dat = args.mergedGVCF
prefix = args.prefix


# output file
output_name = prefix + '.indel.lessthan50bp.txt'
outfile = open(output_name,'w')


def process(line):
	strings = line.strip().split('\t')
#d5
	if ',' in strings[3]:
		alt = strings[3].split(',')
		alt_len = [len(i) for i in alt]
		alt_max = max(alt_len)
	else:
		alt_max = len(strings[3])
#ref
	ref_len = len(strings[2])
	if (alt_max > 50) or (ref_len > 50):
		pass
	else:
		outfile.write(line)


for line in fileinput.input(input_dat):
	process(line)


