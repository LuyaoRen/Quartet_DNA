# import modules
import numpy as np
import pandas as pd
from sklearn import svm
from sklearn import preprocessing
import sys, argparse, os



parser = argparse.ArgumentParser(description="This script is to preform one calss svm on each chromosome")

parser.add_argument('-train', '--trainDataset', type=str, help='Training dataset is Mendelian consistent variants in Quartet samples',  required=True)
parser.add_argument('-test', '--testDataset', type=str, help='Testing dataset is all variants in tested samples',  required=True)
parser.add_argument('-name', '--sampleName', type=str, help='Sample name for output file name',  required=True)

args = parser.parse_args()

# Rename input:
train_input = args.trainDataset
test_input = args.testDataset
sample_name = args.sampleName

# default columns, which will be included in the included in the calssifier
chromosome = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15' ,'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
feature_heter_cols = ['AltDP','BaseQRankSum','DB','DP','FS','GQ','MQ','MQRankSum','QD','ReadPosRankSum','RefDP','SOR']
feature_homo_cols = ['AltDP','DB','DP','FS','GQ','MQ','QD','RefDP','SOR']


# import datasets sepearate the records with or without BaseQRankSum annotation, etc.
# seperate snvs and indels
def load_dat(dat_file_name):
    dat = pd.read_table(dat_file_name)
    dat['DB'] = dat['DB'].fillna(0)
    dat = dat[dat['DP'] != 0]
    dat['af'] = dat['AltDP']/(dat['AltDP'] + dat['RefDP'])
    ref_len = dat['ref'].apply(len)
    alt_len = dat['alt'].apply(len)
    ##snv
    snv_dat = dat.loc[(ref_len==1)&(alt_len==1)]
    snv_homo_rows = snv_dat[snv_dat['BaseQRankSum'].isnull()]
    snv_heter_rows = snv_dat[snv_dat['BaseQRankSum'].notnull()]
    ##indel
    indel_dat = dat.loc[(ref_len!=1)|(alt_len!=1)]
    indel_homo_rows = indel_dat[indel_dat['BaseQRankSum'].isnull()]
    indel_heter_rows = indel_dat[indel_dat['BaseQRankSum'].notnull()]
    return snv_homo_rows,snv_heter_rows,indel_homo_rows,indel_heter_rows



train_snv_homo,train_snv_heter,train_indel_homo,train_indel_heter = load_dat(train_input)
test_snv_homo,test_snv_heter,test_indel_homo,test_indel_heter = load_dat(test_input)
clf = svm.OneClassSVM(nu=0.0001,kernel='rbf')

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
	chr_test['score'] = clf.score_samples(X_test)
	test_true_dat = chr_test[y_pred_test == 1]
	test_false_dat = chr_test[y_pred_test == -1]
	return test_true_dat,test_false_dat

predicted_true = pd.DataFrame(columns=test_snv_homo.columns)
predicted_false = pd.DataFrame(columns=test_snv_homo.columns)

for chromo in chromosome:
	###snv
	# homo datasets
	chr_test_homo,X_train_homo,X_test_homo = prepare_dat(train_snv_homo,test_snv_homo,feature_homo_cols,chromo)
	test_true_homo,test_false_homo = oneclass(X_train_homo,X_test_homo,chr_test_homo)
	predicted_true = predicted_true.append(test_true_homo)
	predicted_false = predicted_false.append(test_false_homo)
	# heter datasets
	chr_test_heter,X_train_heter,X_test_heter = prepare_dat(train_snv_heter,test_snv_heter,feature_heter_cols,chromo)
	test_true_heter,test_false_heter = oneclass(X_train_heter,X_test_heter,chr_test_heter)
	predicted_true = predicted_true.append(test_true_heter)
	predicted_false = predicted_false.append(test_false_heter)
	##indel
	# homo datasets
	chr_test_homo,X_train_homo,X_test_homo = prepare_dat(train_indel_homo,test_indel_homo,feature_homo_cols,chromo)
	test_true_homo,test_false_homo = oneclass(X_train_homo,X_test_homo,chr_test_homo)
	predicted_true = predicted_true.append(test_true_homo)
	predicted_false = predicted_false.append(test_false_homo)
	# heter datasets
	chr_test_heter,X_train_heter,X_test_heter = prepare_dat(train_indel_heter,test_indel_heter,feature_heter_cols,chromo)
	test_true_heter,test_false_heter = oneclass(X_train_heter,X_test_heter,chr_test_heter)
	predicted_true = predicted_true.append(test_true_heter)
	predicted_false = predicted_false.append(test_false_heter)


predicted_true_filename = sample_name + '_predicted_true.txt'
predicted_false_filename = sample_name + '_predicted_false.txt'

predicted_true.to_csv(predicted_true_filename,sep='\t',index=False)
predicted_false.to_csv(predicted_false_filename,sep='\t',index=False)


