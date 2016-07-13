#!/usr/bin/env python
# Shiya Song
# 19rd August 2014
# process-SNP-calling-HC-VQSR.py


import subprocess
import argparse
import os

# Helper function to run commands, handle return values and print to log file
def runCMD(cmd,logFile):
	val = subprocess.Popen(cmd, shell=True).wait()
	if val == 0:
		dt = subprocess.Popen('date', shell=True,stdout=subprocess.PIPE).communicate()[0]
		dt = dt.rstrip()
		s = 'echo \'%s CMD %s\' >> %s' % (dt, cmd, logFile)
		val = subprocess.Popen(s, shell=True).wait()
	else:
		print 'command failed'
		print cmd
		ns = 'echo \'FAILED %s\' >> %s' % (cmd, logFile)
		val = subprocess.Popen(ns, shell=True).wait()
		exit(1)

def runCMDNoFail(cmd):
	val = subprocess.Popen(cmd, shell=True).wait()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='assign parental allele in block')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument('--bam',dest='bamFile',help='bam file after markduplicates')
	parser.add_argument('--level',dest='level',type = float, default = 1, help='level to start at, default = 0 ')
	parser.add_argument('--log',dest='log',help='log file')
	parser.add_argument('--ref_version',dest='ref_version',default='hg19',help='hg19 or b37')
	args = parser.parse_args()

	inputDir = '/'.join(args.bamFile.split('/')[:-2])
	print inputDir

	outDir = inputDir+'/'+args.sample + '/gVCF_calls/'

	gatkDir = '/home/jmkidd/kidd-lab/progs/'

	dt = subprocess.Popen('date', shell=True,stdout=subprocess.PIPE).communicate()[0]
	dt = dt.rstrip()
	cmd = 'echo \'%s STARTING %s level %i\' >> %s' % (args.bamFile, dt, args.level, args.log)
	#print cmd + '\n\n'
	runCMD(cmd,args.log)

	#Setup tmpdir to use
	tmpDir = '/home/jmkidd/kidd-lab-scratch/shiya-projects/temp'

	#Step 0 gVCF to vcf
	if args.ref_version=='b37':
		known_indel1 = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/Mills_and_1000G_gold_standard.indels.b37.vcf'
		known_indel2 = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/1000G_phase1.indels.b37.vcf'
		reference = '/home/jmkidd/kidd-lab/genomes/hg19/human_g1k_v37.fasta'
		dbSNP = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/dbsnp_138.b37.vcf'
		HapmapSNP = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/hapmap_3.3.b37.sites.vcf'
		KGSNP = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/1000G_omni2.5.b37.sites.vcf'
		KGSNP2 = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/1000G_phase1.snps.high_confidence.b37.vcf'
	elif args.ref_version=='hg19':
		known_indel1 = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/Mills_and_1000G_gold_standard.indels.hg19.vcf'
		known_indel2 = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/1000G_phase1.indels.hg19.vcf'
		reference = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa'
		dbSNP = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/dbsnp_138.b37.refmt.vcf'
		HapmapSNP = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/hapmap_3.3.b37.sites.refmt.vcf'
		KGSNP = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/1000G_omni2.5.hg19.vcf'
		KGSNP2 = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/1000G_phase1.snps.high_confidence.hg19.sites.vcf'
	if args.level <=0:
		merge_cmd = 'java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g ' + gatkDir + 'GenomeAnalysisTK.jar -T GenotypeGVCFs '
		merge_cmd += '-R %s ' %(reference)
		for chr in range(1,23)+['X']:
			if args.ref_version=='hg19':
				chr = 'chr'+str(chr)
			out_vcf = outDir + args.sample + '.%s.all.sites.vcf' %(chr)
			if os.path.isfile(out_vcf) is False:
				runCMD('gunzip '+out_vcf+'.gz',args.log)
			merge_cmd += '--variant %s ' %(out_vcf)
		all_vcf = outDir + args.sample + ".raw.SNPs.vcf"
		merge_cmd += '-o ' + all_vcf
		runCMD(merge_cmd,args.log)
		
		for chr in range(1,23)+['X']:
			if args.ref_version=='hg19':
				chr = 'chr'+str(chr)
			out_vcf = outDir + args.sample + '.%s.all.sites.vcf' %(chr)
			runCMD('gzip '+out_vcf,args.log)
		
	if args.level <=1:
		out_vcf = outDir + args.sample + ".raw.SNPs.vcf"
		recal_file = outDir + args.sample + ".raw.SNPs.recal"
		tranch_file = outDir + args.sample + ".raw.SNPs.tranches"
		R_file = outDir + args.sample + ".snp.plots.R"
		cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T VariantRecalibrator '
		cmd += '-R %s -input ' %(reference) + out_vcf 
		cmd += ' -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %s ' %(HapmapSNP)
		cmd += '-resource:omni,known=false,training=true,truth=true,prior=12.0 %s ' %(KGSNP)
		cmd += '-resource:1000G,known=false,training=true,truth=false,prior=10.0 %s ' %(KGSNP2)
		cmd += '-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %s ' %(dbSNP)		
		cmd += '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an DP -mode SNP '
		cmd += '-recalFile '+recal_file+' -tranchesFile '+tranch_file +' -rscriptFile '+ R_file

		runCMD(cmd,args.log)
#		print cmd

		recal_vcf = outDir + args.sample + ".recal.SNPs.vcf"
		cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T ApplyRecalibration '
		cmd += '-R %s -input ' %(reference) + out_vcf  + ' -recalFile '+recal_file+' -tranchesFile '+tranch_file + ' -mode SNP --ts_filter_level 99.0 -o ' + recal_vcf
		runCMD(cmd,args.log)
#		print cmd
	
	if args.level <=2:
		recal_file = outDir + args.sample + ".raw.INDELs.recal"
		tranch_file = outDir + args.sample + ".raw.INDELs.tranches"
		R_file = outDir + args.sample + ".indel.plots.R"
		recal_vcf = outDir + args.sample + ".recal.SNPs.vcf"
		cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T VariantRecalibrator '
		cmd += '-R %s -input ' %(reference) + recal_vcf 
		cmd += ' --maxGaussians 4 -resource:mills,known=false,training=true,truth=true,prior=12.0 %s -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 %s -an QD -an DP -an FS -an ReadPosRankSum -an MQRankSum -mode INDEL ' %(known_indel1,dbSNP)
		cmd += '-recalFile '+recal_file+' -tranchesFile '+tranch_file +' -rscriptFile '+ R_file
		runCMD(cmd,args.log)
#		print cmd

		filtered_vcf = outDir + args.sample + ".recal.ALL.vcf"
		cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T ApplyRecalibration '
		cmd += '-R %s -input ' %(reference) + recal_vcf  + ' -recalFile '+recal_file+' -tranchesFile '+tranch_file + ' -mode INDEL -ts_filter_level 99.0 -o ' + filtered_vcf
		runCMD(cmd,args.log)
#		print cmd


		
		
