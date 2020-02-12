Analyze fosmid pool sequencing data, get physically phased haplotypes, run PSMC/MSMC
==============

Scripts are for `Song S, Sliwerska E, Emery S, Kidd JM. Modeling human population separation history using physically phased genomes. Genetics. 2017 Jan 1;205(1):385-95.`

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
 
> gVCFcaller.py --depth 26 --gvcf #.chr#.all.sites.vcf.gz --chr chr1 --mask #_chr#.mask.bed.gz --legend_file #.legend.gz --log #.log | gzip -c > #.chr#.vcf.gz

* Then generate a #_het_snp.vcf that include all heterzygous SNPs that pass the filter. For each pool, genotype those SNPs.

> python process-pool-SNP-calling-HC.py —sample # —pool #

## Perform physically phasing

> python assign_parental_in_block_gvcf.py

Assign parental allele within each blocks, calculate switch error, correct some switch errors, make matrix2png plots, assign clone into parental alleles

> python calc_concordance_mec_cor.py

For each pair of SNP, calc their physical distance, bin physical distance into categories. Calculate pairwise LD between each pair of SNP, bin into categories.
In each category of LD, calc the genotype concordance, also calc the MEC value; draw LD vs genotype concordance, MEC value, and the correlation of genotype concordance and MEC value

> python phasing_compare.py

## Run PRISM

> python make_prism_input.py

> run_prism.sh

> python prism_posterior.py

> python prism_merge.py

> python prism_interleave_prestep.py

> python prism_interleave_check.py

## Run PSMC

* make psmc input from simulated sequences using Macs

> make-psmc-from-macs-unphase.py

## Run MSMC

> python fa_to_msmc.py

> python msmc_plot.py
