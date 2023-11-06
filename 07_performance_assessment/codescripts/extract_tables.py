import json
import pandas as pd
from functools import reduce
import sys, argparse, os

parser = argparse.ArgumentParser(description="This script is to get information from multiqc and sentieon, output the raw fastq, bam and variants calling (precision and recall) quality metrics")


parser.add_argument('-quality', '--quality_yield', type=str, help='*.quality_yield.txt')
parser.add_argument('-depth', '--wgs_metrics', type=str, help='*deduped_WgsMetricsAlgo.txt')
parser.add_argument('-aln', '--aln_metrics', type=str, help='*_deduped_aln_metrics.txt')
parser.add_argument('-is', '--is_metrics', type=str, help='*_deduped_is_metrics.txt')

parser.add_argument('-fastqc', '--fastqc', type=str, help='multiqc_fastqc.txt')
parser.add_argument('-fastqscreen', '--fastqscreen', type=str, help='multiqc_fastq_screen.txt')
parser.add_argument('-hap', '--happy', type=str, help='multiqc_happy_data.json',  required=True)

parser.add_argument('-project', '--project_name', type=str, help='project_name')

args = parser.parse_args()

if args.quality_yield:
	# Rename input:
	quality_yield_file = args.quality_yield
	wgs_metrics_file = args.wgs_metrics
	aln_metrics_file = args.aln_metrics
	is_metrics_file = args.is_metrics
	fastqc_file = args.fastqc
	fastqscreen_file = args.fastqscreen
	hap_file = args.happy
	project_name = args.project_name
	#############################################
	# fastqc
	fastqc = pd.read_table(fastqc_file)
	# fastqscreen
	dat = pd.read_table(fastqscreen_file)
	fastqscreen = dat.loc[:, dat.columns.str.endswith('percentage')]
	dat['Sample'] = [i.replace('_screen','') for i in dat['Sample']]
	fastqscreen.insert(loc=0, column='Sample', value=dat['Sample'])
	# pre-alignment
	pre_alignment_dat = pd.merge(fastqc,fastqscreen,how="outer",left_on=['Sample'],right_on=['Sample'])
	pre_alignment_dat['FastQC_mqc-generalstats-fastqc-total_sequences'] = pre_alignment_dat['FastQC_mqc-generalstats-fastqc-total_sequences']/1000000
	del pre_alignment_dat['FastQC_mqc-generalstats-fastqc-percent_fails']
	del pre_alignment_dat['FastQC_mqc-generalstats-fastqc-avg_sequence_length']
	del pre_alignment_dat['ERCC percentage']
	del pre_alignment_dat['Phix percentage']
	del pre_alignment_dat['Mouse percentage']
	pre_alignment_dat = pre_alignment_dat.round(2)
	pre_alignment_dat.columns = ['Sample','%Dup','%GC','Total Sequences (million)','%Human','%EColi','%Adapter','%Vector','%rRNA','%Virus','%Yeast','%Mitoch','%No hits']
	pre_alignment_dat.to_csv('pre_alignment.txt',sep="\t",index=0)
	############################
	dat = pd.read_table(aln_metrics_file,index_col=False)
	dat['PCT_ALIGNED_READS'] = dat["PF_READS_ALIGNED"]/dat["TOTAL_READS"]
	aln_metrics = dat[["Sample", "PCT_ALIGNED_READS","PF_MISMATCH_RATE"]]
	aln_metrics = aln_metrics * 100
	aln_metrics['Sample'] = [x[-1] for x in aln_metrics['Sample'].str.split('/')]
	dat = pd.read_table(is_metrics_file,index_col=False)
	is_metrics = dat[['Sample', 'MEDIAN_INSERT_SIZE']]
	is_metrics['Sample'] = [x[-1] for x in is_metrics['Sample'].str.split('/')]
	dat = pd.read_table(quality_yield_file,index_col=False)
	dat['%Q20'] = dat['Q20_BASES']/dat['TOTAL_BASES']
	dat['%Q30'] = dat['Q30_BASES']/dat['TOTAL_BASES']
	quality_yield = dat[['Sample','%Q20','%Q30']]
	quality_yield = quality_yield * 100
	quality_yield['Sample'] = [x[-1] for x in quality_yield['Sample'].str.split('/')]
	dat = pd.read_table(wgs_metrics_file,index_col=False)
	wgs_metrics = dat[['Sample','MEDIAN_COVERAGE','PCT_1X', 'PCT_5X', 'PCT_10X','PCT_20X','PCT_30X']]
	wgs_metrics['PCT_1X'] = wgs_metrics['PCT_1X'] * 100
	wgs_metrics['PCT_5X'] = wgs_metrics['PCT_5X'] * 100
	wgs_metrics['PCT_10X'] = wgs_metrics['PCT_10X'] * 100
	wgs_metrics['PCT_20X'] = wgs_metrics['PCT_20X'] * 100
	wgs_metrics['PCT_30X'] = wgs_metrics['PCT_30X'] * 100
	wgs_metrics['Sample'] = [x[-1] for x in wgs_metrics['Sample'].str.split('/')]
	data_frames = [aln_metrics, is_metrics, quality_yield, wgs_metrics]
	post_alignment_dat = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='outer'), data_frames)
	post_alignment_dat.columns = ['Sample', '%Mapping', '%Mismatch Rate', 'Mendelian Insert Size','%Q20', '%Q30', 'Median Coverage', 'PCT_1X', 'PCT_5X', 'PCT_10X','PCT_20X','PCT_30X']
	post_alignment_dat = post_alignment_dat.round(2)
	post_alignment_dat.to_csv('post_alignment.txt',sep="\t",index=0)
	#########################################
	# variants calling
	with open(hap_file) as hap_json:
		happy = json.load(hap_json)
	dat =pd.DataFrame.from_records(happy)
	dat = dat.loc[:, dat.columns.str.endswith('ALL')]
	dat_transposed = dat.T
	dat_transposed = dat_transposed.loc[:,['sample_id','QUERY.TOTAL','METRIC.Precision','METRIC.Recall']]
	indel = dat_transposed[['INDEL' in s for s in dat_transposed.index]]
	snv = dat_transposed[['SNP' in s for s in dat_transposed.index]]
	indel.reset_index(drop=True, inplace=True)
	snv.reset_index(drop=True, inplace=True)
	benchmark = pd.concat([snv, indel], axis=1)
	benchmark = benchmark[["sample_id", 'QUERY.TOTAL', 'METRIC.Precision', 'METRIC.Recall']]
	benchmark.columns = ['Sample','sample_id','SNV number','INDEL number','SNV precision','INDEL precision','SNV recall','INDEL recall']
	benchmark = benchmark[['Sample','SNV number','INDEL number','SNV precision','INDEL precision','SNV recall','INDEL recall']]
	benchmark['SNV precision'] = benchmark['SNV precision'].astype(float)
	benchmark['INDEL precision'] = benchmark['INDEL precision'].astype(float)
	benchmark['SNV recall'] = benchmark['SNV recall'].astype(float)
	benchmark['INDEL recall'] = benchmark['INDEL recall'].astype(float)
	benchmark['SNV precision'] = benchmark['SNV precision'] * 100
	benchmark['INDEL precision'] = benchmark['INDEL precision'] * 100
	benchmark['SNV recall'] = benchmark['SNV recall'] * 100
	benchmark['INDEL recall'] = benchmark['INDEL recall']* 100
	benchmark = benchmark.round(2)
	benchmark.to_csv('variants.calling.qc.txt',sep="\t",index=0)
else:
	hap_file = args.happy
	with open(hap_file) as hap_json:
		happy = json.load(hap_json)
	dat =pd.DataFrame.from_records(happy)
	dat = dat.loc[:, dat.columns.str.endswith('ALL')]
	dat_transposed = dat.T
	dat_transposed = dat_transposed.loc[:,['sample_id','QUERY.TOTAL','METRIC.Precision','METRIC.Recall']]
	indel = dat_transposed[['INDEL' in s for s in dat_transposed.index]]
	snv = dat_transposed[['SNP' in s for s in dat_transposed.index]]
	indel.reset_index(drop=True, inplace=True)
	snv.reset_index(drop=True, inplace=True)
	benchmark = pd.concat([snv, indel], axis=1)
	benchmark = benchmark[["sample_id", 'QUERY.TOTAL', 'METRIC.Precision', 'METRIC.Recall']]
	benchmark.columns = ['Sample','sample_id','SNV number','INDEL number','SNV precision','INDEL precision','SNV recall','INDEL recall']
	benchmark = benchmark[['Sample','SNV number','INDEL number','SNV precision','INDEL precision','SNV recall','INDEL recall']]
	benchmark['SNV precision'] = benchmark['SNV precision'].astype(float)
	benchmark['INDEL precision'] = benchmark['INDEL precision'].astype(float)
	benchmark['SNV recall'] = benchmark['SNV recall'].astype(float)
	benchmark['INDEL recall'] = benchmark['INDEL recall'].astype(float)
	benchmark['SNV precision'] = benchmark['SNV precision'] * 100
	benchmark['INDEL precision'] = benchmark['INDEL precision'] * 100
	benchmark['SNV recall'] = benchmark['SNV recall'] * 100
	benchmark['INDEL recall'] = benchmark['INDEL recall']* 100
	benchmark = benchmark.round(2)
	benchmark.to_csv('variants.calling.qc.txt',sep="\t",index=0)
	



