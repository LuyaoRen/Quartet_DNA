# GATK CallableLoci

[GATK (v3.8-1) CallableLoci](<https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_coverage_CallableLoci.php>) is used to identify callable regions to create high confidence region  for Quartet reference small variants datasets.

Command line used in this app:

```bash
java -jar /usr/GenomeAnalysisTK.jar \
		-T CallableLoci \
		-R ${ref_dir}/${fasta} \
		-I ${bam} \
		--maxDepth 300 \
		--maxFractionOfReadsWithLowMAPQ 0.1 \
		--maxLowMAPQ 1 \
		--minBaseQuality 20 \
		--minMappingQuality 20 \
		--minDepth 5 \
		--minDepthForLowMAPQ 10 \
		-summary ${sample}_table.txt \
		-o ${sample}_callable_status.bed
```

#####Parameters explanation:

`--maxDepth` 300, If the QC+ depth exceed this value, the site is considered to have EXCESSIVE_DEPTH.

`--maxFractionOfReadsWithLowMAPQ` 0.1, If the number of reads at this site is greater than minDepthForLowMAPQ and the fraction of reads with low mapping quality exceeds this fraction then the site has POOR_MAPPING_QUALITY.

`--maxLowMAPQ` 1, **Maximum value for ([mapping quality](<https://genome.sph.umich.edu/wiki/Mapping_Quality_Scores>)) MAPQ to be considered a problematic mapped read.**
The gap between this value and mmq are reads that are not sufficiently well mapped for calling but aren't indicative of mapping problems. For example, if maxLowMAPQ = 1 and mmq = 20, then reads with MAPQ == 0 are poorly mapped, MAPQ >= 20 are considered as contributing to calling, where reads with MAPQ >= 1 and < 20 are not bad in and of themselves but aren't sufficiently good to contribute to calling. In effect this reads are invisible, driving the base to the NO_ or LOW_COVERAGE states.

`--minBaseQuality` 20, **Minimum quality of bases to count towards depth.**
Bases with less than minBaseQuality are viewed as not sufficiently high quality to contribute to the PASS state.

`--minMappingQuality` 20, **Minimum mapping quality of reads to count towards depth.**
Reads with MAPQ > minMappingQuality are treated as usable for variation detection, contributing to the PASS state.

`--minDepth` 5, **Minimum QC+ read depth before a locus is considered callable**
If the number of QC+ bases (on reads with MAPQ > minMappingQuality and with base quality > minBaseQuality) exceeds this value and is less than maxDepth the site is considered PASS.

`--minDepthForLowMAPQ` 10, **Minimum read depth before a locus is considered a potential candidate for poorly mapped**
We don't want to consider a site as POOR_MAPPING_QUALITY just because it has two reads, and one is MAPQ. We won't assign a site to the POOR_MAPPING_QUALITY state unless there are at least minDepthForLowMAPQ reads covering the site.

##### Per sample running time 

Settings:

Disk size: 400

Cluster: OnDemand ecs.sn1ne.8xlarge img-ubuntu-vpc

sample: 30x WGS 

2h30min