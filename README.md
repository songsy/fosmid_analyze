Analyze fosmid pool sequencing data
==============
# Find clones
Here are some scripts to discover clones, calculate clone statistics and do filtering
* discover clones 
>python find_clone.py #.mpileup #.bam index > #.log
* generate useful command lines
>python write_poststep.py sample_name
* calculate clone statistics and do filtering
>python fosmid_clone_stats.py sample_name
# Perform SNP calling
* call whole genome snp: process-SNP-calling-HC.py
* python process-pool-SNP-calling-HC.py —sample —pool

