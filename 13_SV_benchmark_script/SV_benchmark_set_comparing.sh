#!/bin/bash

helpdoc(){
    cat <<EOF
Description:
    This shell script is used for comparing query SV set with benchmark SV set or comparison within quartet.
Prerequisites:
    Jasmine, a tool to merge structural variants (SVs) across samples, should be in the environment variable.
    bedtools is required for benchmark set comparison. Bcftools is required for quartet comparison.

Usage:
SV discovery:
    $0 -p <Bench> -i <VCF list files> -c <bed file>
SV genotyping:
    $0 -p <Bench> -i <VCF list files> -c <bed file> -g <T/True>
SV reproducibility:
    $0 -p <Quartet> -i <VCF list files>
SV Mendelian concordance rate:
    $0 -p <Quartet> -i <VCF list files> -g <T/True>

Option:
    -p    The performance evaluation (Bench: benchmark set comparison; Quartet: quartet comparison) (Required)
    -i    VCF list files in order: Bench: benchmark set and query set; Quartet: D5, D6, F7 and M8 (Required)
    -t    SV type to compare (choose one from INS and DEL) [default: all] (Optional)
    -g    Consider the genotype concordance of the variants (T/True) [default: False] (Optional)
    -c    The BED file to considers variations in this region using benchmark set comparison (Required)
    -e    The power to use in kd-tree distances (1 is Manhattan, 2 is Euclidean, etc.) [default: 2] (Optional)
    -m    The maximum distance of variants when being merged [default: 1000] (Optional)
    -s    The minimum sequence identity for two insertions to be merged [default: 0] (Optional)	
    -k    The kmer size to use when computing Jaccard similarity of insertions [default: 9] (Optional)	
EOF
}

if [ $# = 0 ]
then
    helpdoc
    exit 1
fi

while getopts "h:p:i:t:g:c:e:m:s:k:" opt
do
    case $opt in
        h)
            helpdoc
            exit 0
            ;;
        p)
            per=$OPTARG
            ;;	
        i)
            vcf_list=$OPTARG
            ;;
        t)
            type=$OPTARG
            ;;
        g)
            geno=$OPTARG
            ;;
        c)
            confidence=$OPTARG
            ;;
        e)
            kd_tree=$OPTARG
            ;;
        m)
            max=$OPTARG
            ;;
        s)
            seq=$OPTARG
            ;;
        k)
            kmer=$OPTARG
            ;;
        ?)
            echo "Unknown option: $opt"
            helpdoc
            exit 1
            ;;
    esac
done


kd_tree_norm=${kd_tree-2}
max_dist=${max-1000}
min_seq_id=${seq-0}
k_jaccard=${kmer-9}

echo "${vcf_list%.*}.txt----running"

#Using SV benchmark
if [ "$per" = "Bench" ];then

#Considerning SV genotype
if [ "$geno" = "T" ] || [ "$geno" = "True" ]
then
jasmine \
file_list=${vcf_list} \
kd_tree_norm=${kd_tree_norm} \
max_dist=${max_dist} \
min_seq_id=${min_seq_id} \
k_jaccard=${k_jaccard} \
--ignore_strand --output_genotypes \
out_file=${vcf_list%.*}.vcf >${vcf_list%.*}.log 2>&1

if [ "$type" = "" ]
then
    echo -e "\nThe precision and recall of SV genotype:\nAll SV type (INS and DEL)"
    cat ${vcf_list%.*}.vcf|bedtools intersect -a - -b $confidence|sort|uniq|grep -E "#|SVTYPE=INS|SVTYPE=DEL"| \
	awk -v OFS="\t" \
    'BEGIN{Benchmark=0;Query=0;shared=0;geno=0}{match($0,/SUPP_VEC=([0-9]+)/,a);split(a[1],vec,"");Benchmark+=vec[1];Query+=vec[2];if(vec[1]+vec[2]==2){num=2+9;if($num !~ /\.\/\./){shared+=1}{split($10,truth,":");split($num,sample,":")}{if(truth[1]==sample[1]){geno+=1}}}}END{print "Bench_num","Query_num","TP_num","Precision","Recall","F1_score","\n",Benchmark,Query,geno,substr(geno/Query,1,6),substr(geno/Benchmark,1,6),substr(2*(geno/Benchmark)*(geno/Query)/((geno/Benchmark)+(geno/Query)),1,6)}'
else
    echo -e "\nThe precision and recall of SV genotype:\n$type type"
    cat ${vcf_list%.*}.vcf|bedtools intersect -a - -b $confidence|sort|uniq|grep -E SVTYPE=$type| \
	awk -v OFS="\t" \
	'BEGIN{Benchmark=0;Query=0;shared=0;geno=0}{match($0,/SUPP_VEC=([0-9]+)/,a);split(a[1],vec,"");Benchmark+=vec[1];Query+=vec[2];if(vec[1]+vec[2]==2){num=2+9;if($num !~ /\.\/\./){shared+=1}{split($10,truth,":");split($num,sample,":")}{if(truth[1]==sample[1]){geno+=1}}}}END{print "Bench_num","Query_num","TP_num","Precision","Recall","F1_score","\n",Benchmark,Query,geno,substr(geno/Query,1,6),substr(geno/Benchmark,1,6),substr(2*(geno/Benchmark)*(geno/Query)/((geno/Benchmark)+(geno/Query)),1,6)}'
fi

#Not considerning SV genotype
else
jasmine \
file_list=${vcf_list} \
kd_tree_norm=${kd_tree_norm} \
max_dist=${max_dist} \
--ignore_strand \
min_seq_id=${min_seq_id} \
k_jaccard=${k_jaccard} \
out_file=${vcf_list%.*}.vcf >${vcf_list%.*}.log 2>&1

if [ "$type" = "" ]
then
    echo -e "\nThe precision and recall of SV discovery:\nAll SV type (INS and DEL)"
    cat ${vcf_list%.*}.vcf|bedtools intersect -a - -b $confidence|sort|uniq|grep -E "#|SVTYPE=INS|SVTYPE=DEL"| \
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' -| \
	awk -v OFS="\t" -v FS="" \
	'BEGIN{Benchmark=0;Query=0;True=0}{Benchmark+=$1;Query+=$2;if($1+$2==2){True+=1}}END{print "Bench_num","Query_num","TP_num","Precision","Recall","F1_score","\n",Benchmark,Query,True,substr(True/Query,1,6),substr(True/Benchmark,1,6),substr(2*(True/Benchmark)*(True/Query)/((True/Benchmark)+(True/Query)),1,6)}'
else
    echo -e "\nThe precision and recall of SV discovery:\n$type type"
    cat ${vcf_list%.*}.vcf|bedtools intersect -a - -b $confidence|sort|uniq|grep -E SVTYPE=$type| \
	perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' -| \
	awk -v OFS="\t" -v FS="" \
	'BEGIN{Benchmark=0;Query=0;True=0}{Benchmark+=$1;Query+=$2;if($1+$2==2){True+=1}}END{print "Bench_num","Query_num","TP_num","Precision","Recall","F1_score","\n",Benchmark,Query,True,substr(True/Query,1,6),substr(True/Benchmark,1,6),substr(2*(True/Benchmark)*(True/Query)/((True/Benchmark)+(True/Query)),1,6)}'
fi

fi

#Using quartet family
elif [ "$per" = "Quartet" ];then

#Considerning SV genotype
if [ "$geno" = "T" ] || [ "$geno" = "True" ]
then
jasmine \
file_list=${vcf_list} \
kd_tree_norm=${kd_tree_norm} \
max_dist=${max_dist} \
min_seq_id=${min_seq_id} \
k_jaccard=${k_jaccard} \
--ignore_strand --output_genotypes \
out_file=${vcf_list%.*}.vcf >${vcf_list%.*}.log 2>&1

if [ "$type" = "" ]
then
    echo -e "\nThe Mendelian concordance rate in quartet:\nAll SV type"
    input=${vcf_list%.*}.vcf
    Mendelian=`bcftools +mendelian <(cat <(grep "#" ${vcf_list%.*}.vcf|sed 's/\tINFO.*/\tINFO\tFORMAT\tmother\tfather\tchild/g') <(grep -E "SUPP_VEC=1111" $input|grep -E -v "\.\/\."|grep -E -v "0\/0.*0\/0.*0\/0.*0\/0.*"| \
    grep -v -E "chrY|chrM"|\
    awk -v OFS="\t" '/#/{print}{split($10,LCL5,":");split($11,LCL6,":")}{if(LCL5[1]==LCL6[1]){print $0}}' -| \
    perl -lane 'if(/SVTYPE=(.*?);.*?GT.*?(.\/.):.*?(.\/.):.*?(.\/.):.*?(.\/.):.*?/){$type="SVTYPE=".$1;$m=$5;$f=$4;$c=$2;print join("\t",@F[0..5],"PASS",$type,"GT",$m,$f,$c)}')) \
    -t mother,father,child -l+|grep -v "#"|wc -l`
    Quartet_union=`grep -v "#" $input|grep -v -E "chrY|chrM"|wc -l`
    awk -v OFS="\t" -v Me=$Mendelian -v Qu=$Quartet_union 'BEGIN{print "Quartet_union","Mendelian","MCR","\n",Qu,Me,substr(Me/Qu,1,6)}'
else
    echo -e "\nThe Mendelian concordance rate in quartet:\n$type type"
    SV=$type
    input=${vcf_list%.*}.vcf
    Mendelian=`bcftools +mendelian <(cat <(grep "#" ${vcf_list%.*}.vcf|sed 's/\tINFO.*/\tINFO\tFORMAT\tmother\tfather\tchild/g') <(grep -E "SUPP_VEC=1111" $input|grep -E -v "\.\/\."|grep -E -v "0\/0.*0\/0.*0\/0.*0\/0.*"| \
    grep -v -E "chrY|chrM"|grep -E SVTYPE=$SV|\
    awk -v OFS="\t" '/#/{print}{split($10,LCL5,":");split($11,LCL6,":")}{if(LCL5[1]==LCL6[1]){print $0}}' -| \
    perl -lane 'if(/SVTYPE=(.*?);.*?GT.*?(.\/.):.*?(.\/.):.*?(.\/.):.*?(.\/.):.*?/){$type="SVTYPE=".$1;$m=$5;$f=$4;$c=$2;print join("\t",@F[0..5],"PASS",$type,"GT",$m,$f,$c)}')) \
    -t mother,father,child -l+|grep -v "#"|wc -l`
    Quartet_union=`grep -v "#" $input|grep -v -E "chrY|chrM"|grep -E SVTYPE=$SV|wc -l`
    awk -v OFS="\t" -v Me=$Mendelian -v Qu=$Quartet_union -v Ty=$SV 'BEGIN{print "Quartet_union","Mendelian","MCR","\n",Qu,Me,substr(Me/Qu,1,6)}'
fi

#Not considerning SV genotype
else
jasmine \
file_list=${vcf_list} \
kd_tree_norm=${kd_tree_norm} \
max_dist=${max_dist} \
--ignore_strand \
min_seq_id=${min_seq_id} \
k_jaccard=${k_jaccard} \
out_file=${vcf_list%.*}.vcf >${vcf_list%.*}.log 2>&1

if [ "$type" = "" ]
then
    echo -e "\nThe reproducibility between monozygotic twins:\nAll SV type"
    cat ${vcf_list%.*}.vcf| \
    perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' -| \
    awk -v OFS="\t" -v FS="" \
    'BEGIN{D5=0;D6=0;share=0}{D5+=$1;D6+=$2;if($1+$2==2){share+=1}}END{print "D5_num","D6_num","Shared_num","Union_num","Reproducibility","\n",D5,D6,share,D5+D6-share,substr(share/(D5+D6-share),1,6)}'
else
    echo -e "\nThe reproducibility between monozygotic twins:\n$type type"
    cat ${vcf_list%.*}.vcf|grep -E SVTYPE=$type| \
    perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' -| \
    awk -v OFS="\t" -v FS="" \
    'BEGIN{D5=0;D6=0;share=0}{D5+=$1;D6+=$2;if($1+$2==2){share+=1}}END{print "D5_num","D6_num","Shared_num","Union_num","Reproducibility","\n",D5,D6,share,D5+D6-share,substr(share/(D5+D6-share),1,6)}'
fi

fi

else

echo -e "\nPlease choose one kind of performance evaluation. See parameter -p"

fi

echo "${vcf_list%.*}.vcf----done"
