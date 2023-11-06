import json
import pandas as pd
import sys, argparse, os
import statistics 

parser = argparse.ArgumentParser(description="This script is to summary information for pre-alignment QC")

parser.add_argument('-general', '--general_stat', type=str, help='multiqc_general_stats.txt',  required=True)
parser.add_argument('-is', '--is_metrics', type=str, help='_is_metrics.txt',  required=True)
parser.add_argument('-wgsmetrics', '--WgsMetricsAlgo', type=str, help='deduped_WgsMetricsAlgo',  required=True)
parser.add_argument('-qualityyield', '--QualityYield', type=str, help='deduped_QualityYield',  required=True)
parser.add_argument('-aln', '--aln_metrics', type=str, help='aln_metrics.txt',  required=True)

args = parser.parse_args()

general_file = args.general_stat
is_file = args.is_metrics
wgsmetrics_file = args.wgsmetrics
qualityyield_file = args.qualityyield
aln_file = args.aln_metrics

##### Table
## general stat: % GC
dat = pd.read_table(general_file)
qualimap = dat.loc[:, dat.columns.str.startswith('QualiMap')]
qualimap.insert(loc=0, column='Sample', value=dat['Sample'])
qualimap_stat = qualimap.dropna()
part1 = fastqc_stat.loc[:,['Sample', 'FastQC_mqc-generalstats-fastqc-percent_duplicates','FastQC_mqc-generalstats-fastqc-total_sequences']]

## is_metrics: median insert size
## deduped_WgsMetricsAlgo: 1x, 5x, 10x, 30x, median coverage
with open(html_file) as file:
	origDict = json.load(file)
newdict = {(k1, k2):v2 for k1,v1 in origDict.items() \
                       for k2,v2 in origDict[k1].items()}
df = pd.DataFrame([newdict[i] for i in sorted(newdict)],
                  index=pd.MultiIndex.from_tuples([i for i in sorted(newdict.keys())]))
gc = []
at = []
for i in part1['Sample']:
	sub_df = df.loc[i,:]
	gc.append(statistics.mean(sub_df['g']/sub_df['c']))
	at.append(statistics.mean(sub_df['a']/sub_df['t']))

## fastq_screen
dat = pd.read_table(fastqscreen_file)
fastqscreen = dat.loc[:, dat.columns.str.endswith('percentage')]
del fastqscreen['ERCC percentage']
del fastqscreen['Phix percentage']

### merge all information
part1.insert(loc=3, column='G/C ratio', value=gc)
part1.insert(loc=4, column='A/T ratio', value=at)
part1.reset_index(drop=True, inplace=True)
fastqscreen.reset_index(drop=True, inplace=True)
df = pd.concat([part1, fastqscreen], axis=1)
df = df.append(df.mean(axis=0),ignore_index=True)
df = df.fillna('Batch average value')
df.columns = ['Sample','Total sequences (million)','% Dup','G/C ratio','A/T ratio','% Human','% EColi','% Adapter' , '% Vector','% rRNA' , '% Virus','% Yeast' ,'% Mitoch' ,'% No hits']
df.to_csv('per-alignment_table_summary.txt',sep='\t',index=False)

##### Picture
## cumulative genome coverage
with open(json_file) as file:
	all_dat = json.load(file)
genome_coverage_json = all_dat['report_plot_data']['qualimap_genome_fraction']['datasets'][0]
dat =pd.DataFrame.from_records(genome_coverage_json)
genome_coverage = pd.DataFrame(index=pd.DataFrame(dat.loc[0,'data'])[0])
for i in range(dat.shape[0]):
	one_sample = pd.DataFrame(dat.loc[i,'data'])
	one_sample.index = one_sample[0]
	genome_coverage[dat.loc[i,'name']] = one_sample[1]
genome_coverage = genome_coverage.transpose() 
genome_coverage['Sample'] = genome_coverage.index
genome_coverage.to_csv('post-alignment_genome_coverage.txt',sep='\t',index=False)

## insert size histogram
insert_size_json = all_dat['report_plot_data']['qualimap_insert_size']['datasets'][0]
dat =pd.DataFrame.from_records(insert_size_json)
insert_size = pd.DataFrame(index=pd.DataFrame(dat.loc[0,'data'])[0])
for i in range(dat.shape[0]):
	one_sample = pd.DataFrame(dat.loc[i,'data'])
	one_sample.index = one_sample[0]
	insert_size[dat.loc[i,'name']] = one_sample[1]
insert_size = insert_size.transpose() 
insert_size['Sample'] = insert_size.index
insert_size.to_csv('post-alignment_insert_size.txt',sep='\t',index=False)

## GC content distribution
gc_content_json = all_dat['report_plot_data']['qualimap_gc_content']['datasets'][0]
dat =pd.DataFrame.from_records(gc_content_json)
gc_content = pd.DataFrame(index=pd.DataFrame(dat.loc[0,'data'])[0])
for i in range(dat.shape[0]):
	one_sample = pd.DataFrame(dat.loc[i,'data'])
	one_sample.index = one_sample[0]
	gc_content[dat.loc[i,'name']] = one_sample[1]
gc_content = gc_content.transpose() 
gc_content['Sample'] = gc_content.index
gc_content.to_csv('post-alignment_gc_content.txt',sep='\t',index=False)


