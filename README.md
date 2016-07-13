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

# Perform SNP calling

* Generate GVCF files using GATK HaplotypeCaller, then generate variant vcf and run VQSR
 
> process-SNP-calling-HC-by-chr.py --bam # --sample # --chr #

> process-SNP-calling-HC-VQSR.py --bam # --sample # --log # --level



*  For each pool, genotype the heterzygous SNPs discovered in whole genome sequence

> python process-pool-SNP-calling-HC.py —sample # —pool #

