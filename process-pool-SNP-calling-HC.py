#!/usr/bin/env python
# Shiya Song
# 26rd June 2014
# process-pool-SNP-calling-HC.py
# Step Zero: Indel Realignment
# Step One: Indel Recalibration
# Step Two: Haplotype Caller
# Step Three: VQSR
# Step Four: Indel calling

# Takes results of bam file after markduplicates and runs through indel realignment and recalibration pipeline.

import subprocess
from optparse import  OptionParser
import os
import re
import commands

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
###############################################################################

# Helper function to run commands,
# doesn't check return or print to log.	 Use for 'grep' so that doesn't
# fail if no items found
def runCMDNoFail(cmd):
	val = subprocess.Popen(cmd, shell=True).wait()


###############################################################################

USAGE = """
process-pool-SNP-calling-HC.py  --sample <sample name> --pool <pool name>
						--level <step to begin at, default = 0> 

										 
Takes results of bam file after markduplicates and runs through indel realignment and recalibration pipeline.

Use --redoCalib to make second qual reclaibration plot

"""

parser = OptionParser(USAGE)
parser.add_option('--pool',dest='poolName',help='bam file after markduplicates')
parser.add_option('--sample',dest='sample',help='sample name')
parser.add_option('--level',dest='level',type = float, default = 2, help='level to start at, default = 0 ')



(options,args)=parser.parse_args()

if options.poolName is None:
	parser.error('pool name not given')
if options.sample is None:
	parser.error('sample not given')

"""
if options.outDir[-1] != '/':
	options.outDir = options.outDir + '/'
	
if options.tmpDir[-1] != '/':
	options.tmpDir = options.tmpDir + '/'
"""	   

###############################################################################
#val = subprocess.Popen('ls -lh', shell=True,stdout=subprocess.PIPE).communicate()[0]
#print 'val is',val

inputDir = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools'

outDir = inputDir+'/'+options.sample + '/' + options.poolName + '/'

gatkDir = '/home/jmkidd/kidd-lab/progs/GenomeAnalysisTK-2.8-1-g932cd3a/'

logFile = outDir+options.poolName+".SNPcalling.log"
#logFile = '/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/'+options.poolName + '/' +options.poolName+".SNPcalling.log"

dt = subprocess.Popen('date', shell=True,stdout=subprocess.PIPE).communicate()[0]
dt = dt.rstrip()
cmd = 'echo \'%s STARTING %s level %i\' >> %s' % (options.poolName, dt, options.level, logFile)
#print cmd + '\n\n'
runCMD(cmd,logFile)

#Setup tmpdir to use
tmpDir = '/home/jmkidd/kidd-lab-scratch/shiya-projects/temp'


###############################################################
#Step 0 Indel Realignment
bamFile = outDir + options.poolName + ".markdup.bam"
realignInterval = outDir + options.poolName + ".realign.intervals"
realignBAM = outDir + options.poolName + ".markdup.realign.bam"
known_indel1 = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/Mills_and_1000G_gold_standard.indels.hg19.vcf'
known_indel2 = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/1000G_phase1.indels.hg19.vcf'

if options.level <=0:
	cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T RealignerTargetCreator '
	cmd += '-R /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa -I ' + bamFile + ' -o ' + realignInterval
	runCMD(cmd,logFile)
#	print cmd

	cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T IndelRealigner '
	cmd += '-R /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa -I ' + bamFile + ' -targetIntervals '+ realignInterval + ' -known ' + known_indel1 + ' -known ' + known_indel2 +' -o ' + realignBAM
	runCMD(cmd,logFile)
#	print cmd

###############################################################
#Step 1 Indel Recalibration
recalInterval = outDir + options.poolName + ".recal_data.grp"
recalBAM = outDir + options.poolName + ".markdup.realign.recal.bam"
dbSNP_file = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/dbsnp_138.b37.refmt.vcf'
after_recal = outDir + options.poolName + ".recal_after.grp"
recal_plots = outDir + options.poolName + ".recal_plots.pdf"
checkbam = outDir + options.poolName + ".bam.check"
if options.level <=1:
	cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T BaseRecalibrator '
	cmd += '-R /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa -I ' + realignBAM + ' -knownSites '+ dbSNP_file + ' -knownSites ' + known_indel1 + ' -knownSites ' + known_indel2 + ' -o ' + recalInterval
	runCMD(cmd,logFile)
#	print cmd

	cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T PrintReads '
	cmd += '-R /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa -I ' + realignBAM + ' -BQSR '+ recalInterval + ' -o ' + recalBAM
	runCMD(cmd,logFile)
#	print cmd
	
	cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T BaseRecalibrator '
	cmd += '-R /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa -I ' + realignBAM + ' -knownSites '+ dbSNP_file + ' -knownSites ' + known_indel1 + ' -knownSites ' + known_indel2 
	cmd += ' -BQSR ' + recalInterval + ' -o ' + after_recal
	runCMD(cmd,logFile)
#	print cmd
	
	cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T AnalyzeCovariates '
	cmd += '-R /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa -before ' + recalInterval + ' -after ' + after_recal + ' -plots ' + recal_plots
#	runCMD(cmd,logFile)
#	print cmd
	
	cmd = "bam validate --in %s --maxErrors 1" %(recalBAM)
	line=commands.getoutput(cmd)
	if line.strip()[-8:-1]=='SUCCESS':
		cmd = 'rm '+realignBAM
		runCMD(cmd,logFile)
		print cmd
		cmd = 'rm '+realignInterval
		runCMD(cmd,logFile)
		print cmd
		cmd = 'rm '+realignBAM.replace('bam','bai')
		runCMD(cmd,logFile)
		print cmd

###############################################################

wgs_vcfFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/' + options.sample+ '/'+ options.sample + '.recal.PASS.vcf'
if options.level <=2:
    out_vcf = outDir + options.poolName + "_wgs_all.vcf"
    cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T HaplotypeCaller '
    cmd += '-R /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa -I ' + recalBAM + ' --genotyping_mode GENOTYPE_GIVEN_ALLELES --alleles ' + wgs_vcfFile
    cmd += ' --dbsnp '+ dbSNP_file + ' -o ' + out_vcf
    runCMD(cmd,logFile)
    #print cmd

###############################################################
# This is for NA19240
wgs_vcfFile = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/' + options.sample+ '/'+ options.sample + '.recal.PASS.vcf'
outDir = '/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/'+options.poolName + '/'
recalBAM = outDir + options.poolName + ".markdup.realign.recal.bam"
if options.level ==3:
    out_vcf = outDir + options.poolName + "_wgs_all.vcf"
    cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T HaplotypeCaller '
    cmd += '-R /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa -I ' + recalBAM + ' --genotyping_mode GENOTYPE_GIVEN_ALLELES --alleles ' + wgs_vcfFile
    cmd += ' --dbsnp '+ dbSNP_file + ' -o ' + out_vcf
    runCMD(cmd,logFile)
    #print cmd