#高置信突变位点的整合

> Author： Run Luyao
>
> E-mail：18110700050@fudan.edu.cn
>
> Git：http://choppy.3steps.cn/renluyao/high_confidence_calls_manuscript.git
>
> Last Updates: 18/03/2020

## 安装指南

```bash
# 激活choppy环境
source activate choppy
# 安装app
choppy install renluyao/high_confidence_calls_manuscript
```

## App概述

中华家系1号全基因组高置信small variants（SNVs和Indels）的整合流程。



## 流程与参数

####1. variantsNorm

保留chr1-22，X上的突变

用bcftools norm进行突变格式的统一

####2. mendelian

LCL5、LCL7和LCL8为三口之家，进行trio-analysis，分析符合孟德尔和不符合孟德尔遗传规律的突变位点

LCL6、LCL7和LCL8为三口之家，进行trio-analysis，分析符合孟德尔和不符合孟德尔遗传规律的突变位点

得到LCL5和LCL6两个家系合并的vcf文件

####3. zipIndex

对LCL5和LCL6两个家系合并的文件压缩和检索引

####4. VCFrename

将VBT的输出结果，VCF文件中的MOTHER FATHER CHILD改成对应的样本名

####5. mergeSister

将LCL5和LCL6修改过名字后的家系VCF文件合并

####6. reformVCF

根据两个三口之家和姐妹的孟德尔遗传的信息，将之前和合并的VCF分解成4个人的vcf，并且包含了家系遗传的信息

####7. familyzipIndex

将上一步输出的4个文件进行压缩和检索引

####8. merge

将所有注释后的LCL5 vcf文件合并

将所有注释后的LCL6 vcf文件合并

将所有注释后的LCL7 vcf文件合并

将所有注释后的LCL8 vcf文件合并

####9. vote

投票选择高置信的突变位点

t

## App输入变量与输入文件

inputSamplesFile的格式如下

```bash
#LCL5_VCF	#LCL6_VCF	#LCL7_VCF	#LCL8_VCF	#LCL5_sampleName	#LCL6_sampleName	#LCL7_sampleName	#LCL8_sampleName	#familyName
```

最终版的整合文件包括：



## App输出文件



## 结果展示与解读



## CHANGELOG



## FAQ

