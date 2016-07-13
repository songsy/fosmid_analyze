#!/usr/bin/env python
# python write_poststep.py sample_name
# Shiya Song
# 26th June 2014
# This script generate cmds for different purpose

import sys
import commands
import os
import argparse
import re
def chr_list():
	list = []
	for i in range(1,23):
		list.append('chr'+str(i))
	list.append('chrX')
	return list

def get_pool_name(file,pool_name):
	inFile = open(file,'r')
	for line in inFile:
		line = line.rstrip()
		line = line.split()
		sample = line[0]
		if sample == 'runName':
			continue
		pool_name[sample]=1
	return pool_name

def find_clone_cmds(pool_name,sample_name):
	outDirBase = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools/%s/' %(sample_name)
	outFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/process/script/%s_%s_findclone.cmds' %(sample_name,commands.getoutput('date +"%m%d"'))
	f_out = open(outFile,'w')
	for pool in pool_name.keys():
		baseDir = outDirBase+pool + '/'
		mpileupFile = baseDir + pool + '.markdup.mpileup'
		findClonesLog = mpileupFile + '.findclones.log'
		markDupBam = baseDir + pool + '.markdup.bam'
		cmd = 'python /home/songsy/script/fosmid_analyze/find_clone.py '
		cmd += '%s %s %s >> %s\n ' % (mpileupFile,markDupBam,pool,findClonesLog)
		f_out.write(cmd)

def SNP_calling_cmds(pool_name,sample_name):
	outFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/process/script/%s_%s_SNPcalling.cmds' %(sample_name,commands.getoutput('date +"%m%d"'))
	f_out = open(outFile,'w')
	for pool in pool_name.keys():
		cmd = 'python /home/songsy/script/fosmid_analyze/process-pool-SNP-calling-HC.py --sample %s --pool %s' %(sample_name,pool)
		f_out.write(cmd+'\n')

def checkbam(pool_name,sample_name):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools'
	for pool in pool_name.keys():
		outDir = inputDir+'/'+sample_name + '/' + pool + '/'
		recalBAM = outDir + pool + ".markdup.realign.recal.bam"
		bamFile = outDir + pool + ".markdup.bam"
		realignInterval = outDir + pool + ".realign.intervals"
		realignBAM = outDir + pool + ".markdup.realign.bam"
		cmd = "bam validate --in %s --maxErrors 1" %(recalBAM)
		line=commands.getoutput(cmd)
		if line.strip()[-8:-1]=='SUCCESS':
			cmd = 'rm '+realignBAM
			os.popen(cmd)
			print cmd
			cmd = 'rm '+realignInterval
			os.popen(cmd)
			print cmd
		else:
			print "Wrong", recalBAM

def checkvcf(pool_name,sample_name):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools'
#	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools'
	for pool in pool_name.keys():
		outDir = inputDir+'/'+sample_name + '/' + pool + '/'
#		outDir = inputDir+'/'+ pool + '/'
		vcf = outDir + pool + "_gvcf_all.vcf"
		if os.path.isfile(vcf) is False:
			print vcf
#			print 'python /home/songsy/script/fosmid_analyze/process-pool-SNP-calling-HC.py --sample HG03428 --pool %s' %(pool)
		else:
			cmd = "tail -n 1 %s" %(vcf)
			line=commands.getoutput(cmd)
			if line[:4]!='chrX':
				print 'python /home/songsy/script/fosmid_analyze/process-pool-SNP-calling-HC.py --sample HG03428 --pool %s' %(pool)
			print line.strip()

def remove(pool_name,sample_name):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools'
	for pool in pool_name.keys():
		outDir = inputDir+'/'+sample_name + '/' + pool + '/'
		file = outDir + pool + ".*.grp"
		cmd = 'rm -rf '+ file
		os.popen(cmd)
		print cmd
		file = outDir + pool + "_wgs_all.vcf*"
		cmd = 'rm -rf '+ file
#		os.popen(cmd)
#		print cmd
		file = outDir + pool + "*.bai"
		cmd = 'rm -rf '+ file
#		os.popen(cmd)
#		print cmd
		file = outDir + pool + ".markdup.1k_coverage"
		cmd = 'gzip '+ file
#		os.popen(cmd)
#		print cmd
		file = outDir + pool + ".markdup.mpileup.gz"
		cmd = 'rm '+ file
#		os.popen(cmd)
#		print cmd

def wgs_all_sites(sample,coverage):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align'
	outFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/process/script/%s_%s_wgs_all_sites.cmds' %(sample_name,commands.getoutput('date +"%m%d"'))
	f = open(outFile,'w')
	for chr in chr_list():
		legend = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz' %(chr)
		cmds = 'samtools mpileup -q 20 -Q 20 -C 50 -u -r %s -f /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa ' %(chr)
		cmds +='%s/%s/%s.realign.recal.bam | bcftools view -cgI - | /home/songsy/script/msmc_tools/bamCaller.py ' %(inputDir,sample,sample)
		cmds +='%i %s/%s/all_sites/%s_%s.mask.bed.gz --legend_file %s| gzip -c > %s/%s/wgs_all_sites/%s.%s.vcf.gz' %(coverage,inputDir,sample,sample,chr,legend,inputDir,sample,sample,chr)
		print >>f,cmds

def wgs_all_sites_v2(sample,coverage):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align'
	outFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/process/script/%s_%s_wgs_all_sites.cmds' %(sample_name,commands.getoutput('date +"%m%d"'))
	f = open(outFile,'w')
	for chr in chr_list():
		legend = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz' %(chr)
		cmds = 'samtools mpileup -q 20 -Q 20 -C 50 -u -r %s -f /home/jmkidd/kidd-lab/genomes/hg19/human_g1k_v37.fasta ' %(chr[3:])
		cmds +='%s/%s/high_coverage/%s.chrom%s.ILLUMINA.bwa.YRI.high_coverage.20100311.bam | bcftools view -cgI - | /home/songsy/script/msmc_tools/bamCaller_v2.py ' %(inputDir,sample,sample,chr[3:])
		cmds +='%i %s/%s/1KG_all_sites/%s_%s.mask.bed.gz --legend_file %s| gzip -c > %s/%s/1KG_all_sites/%s.%s.vcf.gz' %(coverage,inputDir,sample,sample,chr,legend,inputDir,sample,sample,chr)
		print >>f,cmds

def gVCFcaller(sample,coverage):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align'
	outFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/process/script/%s_%s_gVCFcaller.cmds' %(sample,commands.getoutput('date +"%m%d"'))
	f = open(outFile,'w')
	vcf = '%s/%s/%s.recal.ALL.vcf.gz' %(inputDir,sample,sample)
	for chr in chr_list():
		legend = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz' %(chr)
		gvcf = '%s/%s/gVCF_calls/%s.%s.all.sites.vcf.gz' %(inputDir,sample,sample,chr)
		mask_file = '%s/%s/gVCF_calls/%s_%s.mask.bed.gz' %(inputDir,sample,sample,chr)
		out = '%s/%s/gVCF_calls/%s.%s.vcf.gz' %(inputDir,sample,sample,chr)
		log = '%s/%s/gVCF_calls/%s_gVCFcaller.log' %(inputDir,sample,sample)
		cmds = '/home/songsy/script/msmc_tools/gVCFcaller.py --depth %s --gvcf %s --vcf %s --chr %s --mask %s --log %s --legend_file %s| gzip -c > %s' %(coverage,gvcf,vcf,chr,mask_file,log,legend,out)
		print >>f,cmds

def gVCFcaller_v2(sample,coverage):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align'
	outFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/process/script/2015/%s_%s_gVCFcaller.cmds' %(sample,commands.getoutput('date +"%m%d"'))
#	f = open(outFile,'w')
	vcf = '%s/%s/%s.recal.ALL.vcf.gz' %(inputDir,sample,sample)
	for chr in chr_list():
		legend = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz' %(chr)
		gvcf = '%s/%s/gVCF_calls/%s.%s.all.sites.vcf.gz' %(inputDir,sample,sample,chr)
		mask_file = '%s/%s/gVCF_calls/%s_%s.mask.v2.bed.gz' %(inputDir,sample,sample,chr)
		out = '%s/%s/gVCF_calls/%s.%s.vcf.v2.gz' %(inputDir,sample,sample,chr)
		log = '%s/%s/gVCF_calls/%s_gVCFcaller_v2.log' %(inputDir,sample,sample)
		if chr!='chrX':
			cmds = '/home/songsy/script/msmc_tools/gVCFcaller_v2.py --depth %s --gvcf %s --vcf %s --chr %s --mask %s --log %s --legend_file %s| gzip -c > %s' %(coverage,gvcf,vcf,chr,mask_file,log,legend,out)
		else:
			cmds = '/home/songsy/script/msmc_tools/gVCFcaller_v2.py --depth %s --gvcf %s --vcf %s --chr %s --mask %s --log %s | gzip -c > %s' %(coverage,gvcf,vcf,chr,mask_file,log,out)
			print cmds

def wgs_all_sites_v3(sample):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align'
	outFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/process/script/%s_%s_wgs_all_sites.cmds' %(sample,commands.getoutput('date +"%m%d"'))
	f = open(outFile,'w')
	for chr in range(1,23)+['X']:
		cmds = 'python /home/songsy/script/snp-calling/process-SNP-calling-HC-all-by-chr.py --bam %s/%s/high_coverage/NA19240.chrom%s.ILLUMINA.bwa.YRI.high_coverage.20100311.bam' %(inputDir,sample,chr)
		cmds += ' --chr %s --sample %s --log %s/%s/high_coverage/%s_wgs_chr%s_snpcalling.log' %(chr,sample,inputDir,sample,sample,chr)
		print >>f,cmds

def wgs_all_sites_v4(sample):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align'
	outFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/process/script/%s_additional_%s_wgs_all_sites.cmds' %(sample,commands.getoutput('date +"%m%d"'))
	f = open(outFile,'w')
	for chr in range(1,23)+['X']:
		cmds = 'python /home/songsy/script/snp-calling/process-SNP-calling-HC-all-by-chr.py --bam %s/%s-additional/NA19240-additional.markdup.merge.bam' %(inputDir,sample)
		cmds += ' --chr chr%s --sample %s --log %s/%s-additional/%s_wgs_chr%s.log' %(chr,sample,inputDir,sample,sample,chr)
		print >>f,cmds

def wgs_all_sites_v5(sample):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align'
	outFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/process/script/%s_%s_wgs_all_sites.cmds' %(sample,commands.getoutput('date +"%m%d"'))
	f = open(outFile,'w')
	for chr in range(1,23)+['X']:
		cmds = 'python /home/songsy/script/snp-calling/process-SNP-calling-HC-all-by-chr.py --bam %s/%s/NA19240.realign.recal.bam' %(inputDir,sample)
		cmds += ' --chr chr%s --sample %s --log %s/%s/wgs/%s_wgs_chr%s_snpcalling.log' %(chr,sample,inputDir,sample,sample,chr)
		print >>f,cmds

def wgs_all_sites_v6(sample):
	inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align'
	outFile = '%s_wgs_all_sites.cmds' %(sample)
	f = open(outFile,'w')
	for chr in range(1,23)+['X']:
		cmds = 'python /home/songsy/script/snp-calling/process-SNP-calling-HC-all-by-chr.py --bam SRR1598471.markdup.bam'
		cmds += ' --chr chr%s --sample %s --log %s_wgs_chr%s.log' %(chr,sample,sample,chr)
		print >>f,cmds

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='assign parental allele in block')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--coverage", dest='coverage',help="sequence coverage")
	args = parser.parse_args()
	wk = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/'
	if args.sample =='HG02799':
		pool_name = {}
		file = wk+'HG02799-fosmids.plate4.txt.info'
		pool_name=get_pool_name(file,pool_name)
		file = wk+'HG02799-fosmids.plates1_2.txt.info'
		pool_name=get_pool_name(file,pool_name)
	if args.sample =='HG03108':
		pool_name = {}
		file = wk+'HG03108-fosmids.plate2.txt.info'
		pool_name=get_pool_name(file,pool_name)
		file = wk+'HG03108-fosmids.plates3_4.txt.info'
		pool_name=get_pool_name(file,pool_name)
	if args.sample =='NA21302':
		pool_name = {}
		file = wk+'NA21302-fosmids.plate1.txt.info'
		pool_name=get_pool_name(file,pool_name)
		file = wk+'NA21302-fosmids.plate2.txt.info'
		pool_name=get_pool_name(file,pool_name)
		file = wk+'HG21302_p3.txt.info'
		pool_name=get_pool_name(file,pool_name)
	if args.sample == 'HG03428':
		pool_name = {}
		a = os.listdir('/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools/HG03428/')
		for i in a:
			if i[-10:]=='clonesizes':
				continue
			m = re.search('.+run759_match.+',i)
			if m is not None:
				continue
			if i[:9]=='HG03428L2' or i[:10]=='HG03428_p2' or i[:10]=='HG03428_p1' or i[:14]=='HG03428_plate1':
				pool_name[i]=1
			if i[:14]=='HG03428_plate2':
				if i[-8:]=='combined':
					pool_name[i]=1
				else:
					if i+'_combined' not in a:
						pool_name[i]=1
#		print len(pool_name)
#		print pool_name
	if args.sample == 'NA19240':
		pool_name = {}
		interval=range(1,289)
		interval.remove(120)
		interval.remove(155)
		for i in interval:
			if i <=9:
				index="00"+str(i)
			elif i<=99:
				index="0"+str(i)
			else:
				index=str(i)
			pool_name['NA19240_pool_'+index]=1
	if args.sample == 'NA20847':
		pool_name={}
		for i in range(1,116):
			pool_name['fraction'+str(i)]=1
#	find_clone_cmds(pool_name,args.sample)
#	SNP_calling_cmds(pool_name,args.sample)
#	checkbam(pool_name,sample_name)
#	checkvcf(pool_name,args.sample)
#	wgs_all_sites_v6(args.sample)
#	remove(pool_name,args.sample)
#	wgs_all_sites(sample_name,coverage)
	gVCFcaller_v2(args.sample,args.coverage)
#	wgs_all_sites_v3(args.sample)
#	wgs_all_sites_v5(args.sample)