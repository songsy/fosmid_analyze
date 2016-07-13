Analyze fosmid pool sequencing data
==============
## Find clones
Here are some scripts to discover clones, calculate clone statistics and do filtering
* discover clones 

>python find_clone.py #.mpileup #.bam index > #.log

* generate useful command lines

>python write_poststep.py sample_name

* calculate clone statistics and do filtering

>python fosmid_clone_stats.py sample_name

## Perform SNP calling

* Generate GVCF files using GATK HaplotypeCaller, then generate variant vcf and run VQSR

> process-SNP-calling-HC-by-chr.py --bam # --sample # --chr #

> process-SNP-calling-HC-VQSR.py --bam # --sample # --log # --level

Then generate callable sites (mask.bed.gz) and the vcf files with variants that either pass VQSR or present in reference panel, by each chromosome.
The vcf.gz and mask.bed.gz serves as input for MSMC
 
> gVCFcaller.py --depth 26 --gvcf #.chr#.all.sites.vcf.gz --chr chr1 --mask #_chr#.mask.bed.gz --legend_file ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz --log #.log | gzip -c > #.chr#.vcf.gz

*  For each pool, genotype the heterzygous SNPs discovered in whole genome sequence

> python process-pool-SNP-calling-HC.py —sample # —pool #


