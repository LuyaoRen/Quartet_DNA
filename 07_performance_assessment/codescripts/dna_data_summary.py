import pandas as pd
from functools import reduce
import sys, argparse, os
import glob

parser = argparse.ArgumentParser(description="This script is to get information from multiqc and sentieon, output the raw fastq, bam and variants calling (precision and recall) quality metrics")


parser.add_argument('-pre', '--pre_alignment_path', type=str, help='pre-alignment directory')
parser.add_argument('-post', '--post_alignment_path', type=str, help='post-alignment directory')
parser.add_argument('-vcf', '--benchmark_path', type=str, help='benchmark directory')
parser.add_argument('-men', '--mendelian_path', type=str, help='mendelian directory')


args = parser.parse_args()


# Rename input:
pre_alignment_path = args.pre_alignment_path
post_alignment_path = args.post_alignment_path
benchmark_path = args.benchmark_path
mendelian_path = args.mendelian_path


###pre-alignment
pre_alignment_df = pd.concat(map(pd.read_table,glob.glob(os.path.join(pre_alignment_path,'*.txt'))))

pre_alignment_df.to_csv('pre-alignment-all.txt',sep="\t",index=0)

###post-alignment

post_alignment_df = pd.concat(map(pd.read_table,glob.glob(os.path.join(post_alignment_path,'*.txt'))))

post_alignment_df.to_csv('post-alignment-all.txt',sep="\t",index=0)

###variant-calling qc

## detail
variant_calling_df = pd.concat(map(pd.read_table,glob.glob(os.path.join(benchmark_path,'*.txt'))))

variant_calling_df.to_csv('variant-calling-all.txt',sep="\t",index=0)

## mean + sd

####snv
a = variant_calling_df["SNV precision"].mean().round(2)
b = variant_calling_df["SNV precision"].std().round(2)

c = variant_calling_df["SNV recall"].mean().round(2)
d = variant_calling_df["SNV precision"].std().round(2)

variant_calling_df["SNV F1 score"] = 2 * variant_calling_df["SNV precision"] * variant_calling_df["SNV recall"] / (variant_calling_df["SNV precision"] + variant_calling_df["SNV recall"])

e = variant_calling_df["SNV F1 score"].mean().round(2)
f = variant_calling_df["SNV F1 score"].std().round(2)

#### indel

a2 = variant_calling_df["INDEL precision"].mean().round(2)
b2 = variant_calling_df["INDEL precision"].std().round(2)

c2 = variant_calling_df["INDEL recall"].mean().round(2)
d2 = variant_calling_df["INDEL precision"].std().round(2)

variant_calling_df["INDEL F1 score"] = 2 * variant_calling_df["INDEL precision"] * variant_calling_df["INDEL recall"] / (variant_calling_df["INDEL precision"] + variant_calling_df["INDEL recall"])

e2 = variant_calling_df["INDEL F1 score"].mean().round(2)
f2 = variant_calling_df["INDEL F1 score"].std().round(2)

data = {'precision':[str(a),str(a2)], 'precision_sd':[str(b),str(b2)], 'recall':[str(c),str(c2)], 'recall_sd':[str(d),str(d2)], 'F1-score':[str(e),str(e2)], 'F1-score_sd':[str(f),str(f2)] }  

df = pd.DataFrame(data, index =['SNV', 'INDEL'])  

df.to_csv('benchmark_summary.txt',sep="\t")


### Mendelian

mendelian_df = pd.concat(map(pd.read_table,glob.glob(os.path.join(mendelian_path,'*.txt'))))

mendelian_df.to_csv('mendelian-all.txt',sep="\t",index=0)

### snv
snv = mendelian_df[mendelian_df['Family'].str.contains("SNV")]

indel = mendelian_df[mendelian_df['Family'].str.contains("INDEL")]

g = snv["Mendelian_Concordance_Rate"].mean().round(2)
h = snv["Mendelian_Concordance_Rate"].std().round(2)

g2 = indel["Mendelian_Concordance_Rate"].mean().round(2)
h2 = indel["Mendelian_Concordance_Rate"].std().round(2)


data = {'MCR':[str(g),str(g2)], 'MCR_sd':[str(h),str(h2)]}  

df = pd.DataFrame(data, index =['SNV', 'INDEL'])  

df.to_csv('mendelian_summary.txt',sep="\t")



