import pandas as pd
import sys, argparse, os
from operator import itemgetter   

parser = argparse.ArgumentParser(description="this script is to annotate high confidence calls")

parser.add_argument('-info', '--info', type=str, help='The infomation file',  required=True)
parser.add_argument('-vcf', '--vcf', type=str, help='The vcf file',  required=True)
parser.add_argument('-prefix', '--prefix', type=str, help='The outputname',  required=True)


args = parser.parse_args()

# Rename input:
info = args.info
vcf = args.vcf
prefix = args.prefix

info = pd.read_table(info,header=None)
vcf = pd.read_table(vcf,header=None)
merged_df = pd.merge(vcf, info,  how='inner', left_on=[0,1], right_on = [0,1])

filename = prefix + '.annotated.txt'
merged_df.to_csv(filename,header=None,index=None,sep="\t")
