
#canu
cd canu-1.8/
wget -c https://github.com/marbl/canu/releases/download/v1.8/canu-1.8.Linux-amd64.tar.xz
xz -d canu-1.8.Linux-amd64.tar.xz
tar xvf canu-1.8.Linux-amd64.tar

canu -p prefix \
-d prefix genomeSize=3.1g -pacbio-raw pacbio.fasta.gz \
-haplotypeF7 F7.NGS.fastq.gz -haplotypeM8 M8.NGS.fastq.gz \
gridEngineMemoryOption="-l mem_total=MEMORY" -s spec.txt

#spect.txt
corFilter="quick"
ovbMemory=8g
maxMemory=500g
maxThreads=18
ovsMemory=8-500g

corThreads=4
cormhapThreads=4
obtmhapThreads=4
utgmhapThreads=4

ovsThreads=2
oeaThreads=2
ovlThreads=3
merylThreads=5
redThreads=4
utgmhapThreads=4
utgovlThreads=4
cnsThreads=4
gfaThreads=4
corovlThreads=4

#svanalyzer
conda install svanalyzer
v0.36

#mummer-4.0.0beta2
wget -c https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
tar -xvzf mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2
./configure --prefix=/software/mummer-4.0.0beta2
make
make install

dir=/software/mummer-4.0.0beta2
$dir/nucmer -t 20 --maxmatch -l 100 -c 500 –p case_ nucmer ref.fa case_assembly.fa


#Assemblytics
wget -c https://github.com/MariaNattestad/Assemblytics/archive/1.2.1.tar.gz
tar -xvzf 1.2.1.tar.gz
cd /software/Assemblytics-1.2.1/scripts
cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region
dir=/software/Assemblytics-1.2.1/scripts

SURVIVOR convertAssemblytics D5_haplotypeF7.Assemblytics_structural_variants.bed 19 D5_haplotypeF7_Assemblytics_SV_20bp.vcf
SURVIVOR convertAssemblytics D5_haplotypeF7.Assemblytics_structural_variants.bed 49 D5_haplotypeF7_Assemblytics_SV_50bp.vcf
SURVIVOR convertAssemblytics D5_haplotypeM8.Assemblytics_structural_variants.bed 19 D5_haplotypeM8_Assemblytics_SV_20bp.vcf
SURVIVOR convertAssemblytics D5_haplotypeM8.Assemblytics_structural_variants.bed 49 D5_haplotypeM8_Assemblytics_SV_50bp.vcf

cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/D6_haplotypeF7
SURVIVOR convertAssemblytics D6_haplotypeF7.Assemblytics_structural_variants.bed 49 D6_haplotypeF7_Assemblytics_SV_50bp.vcf
cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/D6_haplotypeM8
SURVIVOR convertAssemblytics D6_haplotypeM8.Assemblytics_structural_variants.bed 49 D6_haplotypeM8_Assemblytics_SV_50bp.vcf
cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/Assemblytics/M8
SURVIVOR convertAssemblytics M8.Assemblytics_structural_variants.bed 49 M8_Assemblytics_SV_50bp.vcf


#svmu-0.4
git clone https://github.com/mahulchak/svmu.git
cd svmu/
make

dir=/software/mummer-4.0.0beta2
#$dir/nucmer --threads 4 --prefix D5_M8 $ref $LCL5
$dir/nucmer -t 12 -l 100 -c 500 --prefix=D5_M8 $ref $LCL5
/software/svmu/svmu D5_M8.delta $ref $LCL5 l last_out.txt D5-haplotypeM8

cd /Quartet_data/SV_combine/pacbio_ont_pacbio2_combine_3/Assembly_SV_Region/svmu/D5_haplotypeM8
i=sv.M8.txt
awk -v OFS='\t' '{if($9>=50)print}' $i| \
perl -lane 'if(/INS/){$num+=1;print "$F[0]\t$F[1]\tsvmu_M8_$num\tN\t<INS>\t.\tPASS\tSVMETHOD=svmu;CHR2=$F[0];END=$F[1];SVTYPE=INS;SVLEN=$F[8]\tGT\t./."}elsif(/DEL/){$num+=1;print "$F[0]\t$F[1]\tsvmu_M8_$num\tN\t<DEL>\t.\tPASS\tSVMETHOD=svmu;CHR2=$F[0];END=$F[2];SVTYPE=DEL;SVLEN=$F[8]\tGT\t./."}' - \
> M8_svmu_50bp.vcf

#paftools
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cd k8-0.2.4

/k8-0.2.4/k8-Linux
minimap2 -cx asm5 -t8 --cs ref.fa asm.fa > asm.paf  # keeping this file is recommended; --cs required!
sort -k6,6 -k8,8n asm.paf > asm.srt.paf             # sort by reference start coordinate
k8 paftools.js call asm.srt.paf > asm.var.txt

/miniconda3/bin/paftools.js version
2.17-r941

i=LCL8_asm.var.txt
perl -lane 'if(length($F[6])<length($F[7])){if(length($F[7])>=50){$num+=1;$LEN=length($F[7]);@ref=split(//,$F[7]);print "$F[1]\t$F[2]\tpaftools_M8_$num\t$ref[0]\t$F[7]\t$F[4]\t$F[5]\tSVMETHOD=paftools;CHR2=$F[1];END=$F[2];SVTYPE=INS;SVLEN=$LEN\tGT\t./."}}elsif(length($F[6])>length($F[7])){if(length($F[6])>=50){$num+=1;$LEN=length($F[6]);@ref=split(//,$F[6]);print "$F[1]\t$F[2]\tpaftools_M8_$num\t$F[6]\t$ref[0]\t$F[4]\t$F[5]\tSVMETHOD=paftools;CHR2=$F[1];END=$F[3];SVTYPE=DEL;SVLEN=$LEN\tGT\t./."}}else{next}' \
$i > M8_paftools_50bp.vcf