# import modules
import numpy as np
import pandas as pd
from sklearn import svm
from sklearn import preprocessing
import sys, argparse, os
from vcf2bed import position_to_bed,padding_region



parser = argparse.ArgumentParser(description="this script is to preform one calss svm on each chromosome")

parser.add_argument('-train', '--trainDataset', type=str, help='training dataset generated from extracting vcf information part, with mutaitons supported by callsets',  required=True)
parser.add_argument('-test', '--testDataset', type=str, help='testing dataset generated from extracting vcf information part, with mutaitons not called by all callsets',  required=True)
parser.add_argument('-name', '--sampleName', type=str, help='sample name for output file name',  required=True)

args = parser.parse_args()

# Rename input:
train_input = args.trainDataset
test_input = args.testDataset
sample_name = args.sampleName

# default columns, which will be included in the included in the calssifier
chromosome = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15' ,'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
feature_heter_cols = ['AltDP','BaseQRankSum','DB','DP','FS','GQ','MQ','MQRankSum','QD','ReadPosRankSum','RefDP','SOR','af']
feature_homo_cols = ['AltDP','DB','DP','FS','GQ','MQ','QD','RefDP','SOR','af']


# import datasets sepearate the records with or without BaseQRankSum annotation, etc.
def load_dat(dat_file_name):
    dat = pd.read_table(dat_file_name)
    dat['DB'] = dat['DB'].fillna(0)
    dat = dat[dat['DP'] != 0]
    dat['af'] = dat['AltDP']/(dat['AltDP'] + dat['RefDP'])
    homo_rows = dat[dat['BaseQRankSum'].isnull()]
    heter_rows = dat[dat['BaseQRankSum'].notnull()]
    return homo_rows,heter_rows


train_homo,train_heter = load_dat(train_input)
test_homo,test_heter = load_dat(test_input)
clf = svm.OneClassSVM(nu=0.05,kernel='rbf', gamma='auto_deprecated',cache_size=500)

def prepare_dat(train_dat,test_dat,feature_cols,chromo):
	chr_train = train_dat[train_dat['chromo'] == chromo]
	chr_test = test_dat[test_dat['chromo'] == chromo]
	
	train_dat = chr_train.loc[:,feature_cols]
	test_dat = chr_test.loc[:,feature_cols]
	train_dat_scaled = preprocessing.scale(train_dat)
	test_dat_scaled = preprocessing.scale(test_dat)
	return chr_test,train_dat_scaled,test_dat_scaled

def oneclass(X_train,X_test,chr_test):
	clf.fit(X_train)
	y_pred_test = clf.predict(X_test)
	test_true_dat = chr_test[y_pred_test == 1]
	test_false_dat = chr_test[y_pred_test == -1]
	return test_true_dat,test_false_dat

predicted_true = pd.DataFrame(columns=train_homo.columns)
predicted_false = pd.DataFrame(columns=train_homo.columns)

for chromo in chromosome:
	# homo datasets
	chr_test_homo,X_train_homo,X_test_homo = prepare_dat(train_homo,test_homo,feature_homo_cols,chromo)
	test_true_homo,test_false_homo = oneclass(X_train_homo,X_test_homo,chr_test_homo)
	predicted_true = predicted_true.append(test_true_homo)
	predicted_false = predicted_false.append(test_false_homo)
	# heter datasets
	chr_test_heter,X_train_heter,X_test_heter = prepare_dat(train_heter,test_heter,feature_heter_cols,chromo)
	test_true_heter,test_false_heter = oneclass(X_train_heter,X_test_heter,chr_test_heter)
	predicted_true = predicted_true.append(test_true_heter)
	predicted_false = predicted_false.append(test_false_heter)

predicted_true_filename = sample_name + '_predicted_true.txt'
predicted_false_filename = sample_name + '_predicted_false.txt'

predicted_true.to_csv(predicted_true_filename,sep='\t',index=False)
predicted_false.to_csv(predicted_false_filename,sep='\t',index=False)

# output the bed file and padding bed region 50bp

predicted_true_bed_filename = sample_name + '_predicted_true.bed'
predicted_false_bed_filename = sample_name + '_predicted_false.bed'
padding_filename = sample_name + '_padding.bed'

predicted_true_bed = open(predicted_true_bed_filename,'w')
predicted_false_bed = open(predicted_false_bed_filename,'w')
padding = open(padding_filename,'w')

#
for index,row in predicted_false.iterrows():
	chromo,pos1,pos2 = position_to_bed(row['chromo'],row['pos'],row['ref'],row['alt'])
	outline_pos = chromo + '\t' + str(pos1) + '\t' + str(pos2) + '\n'
	predicted_false_bed.write(outline_pos)
	chromo,pad_pos1,pad_pos2,pad_pos3,pad_pos4 = padding_region(chromo,pos1,pos2,50)
	outline_pad_1 = chromo + '\t' + str(pad_pos1) + '\t' + str(pad_pos2) + '\n'
	outline_pad_2 = chromo + '\t' + str(pad_pos3) + '\t' + str(pad_pos4) + '\n'
	padding.write(outline_pad_1)
	padding.write(outline_pad_2)

for index,row in predicted_true.iterrows():
	chromo,pos1,pos2 = position_to_bed(row['chromo'],row['pos'],row['ref'],row['alt'])
	outline_pos = chromo + '\t' + str(pos1) + '\t' + str(pos2) + '\n'
	predicted_true_bed.write(outline_pos)


