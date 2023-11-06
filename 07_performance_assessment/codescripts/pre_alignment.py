import json
import pandas as pd
import sys, argparse, os
import statistics 

parser = argparse.ArgumentParser(description="This script is to summary information for pre-alignment QC")

parser.add_argument('-general', '--general_stat', type=str, help='multiqc_general_stats.txt',  required=True)
parser.add_argument('-html', '--html', type=str, help='multiqc_report.html',  required=True)
parser.add_argument('-fastqscreen', '--fastqscreen', type=str, help='multiqc_fastq_screen.txt',  required=True)
parser.add_argument('-json', '--json', type=str, help='multiqc_happy_data.json',  required=True)

args = parser.parse_args()

general_file = args.general_stat
html_file = args.html
fastqscreen_file = args.fastqscreen
json_file = args.json

##### Table
## general stat: 1. Total sequences; 2. %Dup
dat = pd.read_table(general_file)

fastqc = dat.loc[:, dat.columns.str.startswith('FastQC')]
fastqc.insert(loc=0, column='Sample', value=dat['Sample'])
fastqc_stat = fastqc.dropna()
part1 = fastqc_stat.loc[:,['Sample', 'FastQC_mqc-generalstats-fastqc-percent_duplicates','FastQC_mqc-generalstats-fastqc-total_sequences']]

## report html: 1. G/C ratio; 2. A/T ratio
## cat multiqc_report.html | grep 'fastqc_seq_content_data = ' | sed s'/fastqc_seq_content_data\ =\ //g' | sed 's/^[ \t]*//g' | sed s'/;//g' > fastqc_sequence_content.json
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
## mean quality scores
with open(json_file) as file:
	all_dat = json.load(file)
mean_quality_json = all_dat['report_plot_data']['fastqc_per_base_sequence_quality_plot']['datasets'][0]
dat =pd.DataFrame.from_records(mean_quality_json)
mean_quality = pd.DataFrame(index=pd.DataFrame(dat.loc[0,'data'])[0])
for i in range(dat.shape[0]):
	one_sample = pd.DataFrame(dat.loc[i,'data'])
	one_sample.index = one_sample[0]
	mean_quality[dat.loc[i,'name']] = one_sample[1]
mean_quality = mean_quality.transpose() 
mean_quality['Sample'] = mean_quality.index
mean_quality.to_csv('pre-alignment_mean_quality.txt',sep='\t',index=False)

## per sequence GC content

gc_content_json = all_dat['report_plot_data']['fastqc_per_sequence_gc_content_plot']['datasets'][0]
dat =pd.DataFrame.from_records(gc_content_json)
gc_content = pd.DataFrame(index=pd.DataFrame(dat.loc[0,'data'])[0])
for i in range(dat.shape[0]):
	one_sample = pd.DataFrame(dat.loc[i,'data'])
	one_sample.index = one_sample[0]
	gc_content[dat.loc[i,'name']] = one_sample[1]
gc_content = gc_content.transpose() 
gc_content['Sample'] = gc_content.index
gc_content.to_csv('pre-alignment_gc_content.txt',sep='\t',index=False)

# fastqc and qualimap
dat = pd.read_table(fastqc_qualimap_file)

fastqc = dat.loc[:, dat.columns.str.startswith('FastQC')]
fastqc.insert(loc=0, column='Sample', value=dat['Sample'])
fastqc_stat = fastqc.dropna()

# qulimap
qualimap = dat.loc[:, dat.columns.str.startswith('QualiMap')]
qualimap.insert(loc=0, column='Sample', value=dat['Sample'])
qualimap_stat = qualimap.dropna()

# fastqc
dat = pd.read_table(fastqc_file)

fastqc_module = dat.loc[:, "per_base_sequence_quality":"kmer_content"]
fastqc_module.insert(loc=0, column='Sample', value=dat['Sample'])
fastqc_all = pd.merge(fastqc_stat,fastqc_module,  how='outer', left_on=['Sample'], right_on = ['Sample'])

# fastqscreen
dat = pd.read_table(fastqscreen_file)
fastqscreen = dat.loc[:, dat.columns.str.endswith('percentage')]
dat['Sample'] = [i.replace('_screen','') for i in dat['Sample']]
fastqscreen.insert(loc=0, column='Sample', value=dat['Sample'])

# benchmark
with open(hap_file) as hap_json:
	happy = json.load(hap_json)
dat =pd.DataFrame.from_records(happy)
dat = dat.loc[:, dat.columns.str.endswith('ALL')]
dat_transposed = dat.T
benchmark = dat_transposed.loc[:,['sample_id','METRIC.Precision','METRIC.Recall']]
benchmark['sample_id'] = benchmark.index
benchmark.columns = ['Sample','Precision','Recall']

#output
fastqc_all.to_csv('fastqc.final.result.txt',sep="\t",index=0)
fastqscreen.to_csv('fastqscreen.final.result.txt',sep="\t",index=0)
qualimap_stat.to_csv('qualimap.final.result.txt',sep="\t",index=0)
benchmark.to_csv('benchmark.final.result.txt',sep="\t",index=0)




