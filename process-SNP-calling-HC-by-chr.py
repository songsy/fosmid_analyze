#!/usr/bin/env python
# Shiya Song
# 23rd June 2014
# process-SNP-calling-HC-by-chr.py
# Step Zero: Indel Realignment
# Step One: Indel Recalibration
# Step Two: Haplotype Caller
# Step Three: VQSR
# Step Four: Indel calling

# Takes results of bam file after markduplicates and runs through indel realignment and recalibration pipeline.

import subprocess
from optparse import  OptionParser
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
###############################################################################

# Helper function to run commands,
# doesn't check return or print to log.	 Use for 'grep' so that doesn't
# fail if no items found
def runCMDNoFail(cmd):
	val = subprocess.Popen(cmd, shell=True).wait()


###############################################################################

USAGE = """
process-SNP-VQSR.py --bam <bam file> --sample <sample name> --chr <chromosome>
						--log <file for log output>
						--level <step to begin at, default = 0> 

										 
Takes results of bam file after markduplicates and runs through indel realignment and recalibration pipeline.

Use --redoCalib to make second qual reclaibration plot

"""

parser = OptionParser(USAGE)
parser.add_option('--bam',dest='bamFile',help='bam file after markduplicates')
parser.add_option('--sample',dest='sample',help='sample name')
parser.add_option('--chr',dest='chr',help='chromosome')
parser.add_option('--level',dest='level',type = float, default = 0, help='level to start at, default = 0 ')
parser.add_option('--log',dest='log',help='log file')


(options,args)=parser.parse_args()

if options.bamFile is None:
	parser.error('bam file not given')
if options.sample is None:
	parser.error('sample not given')
if options.log is None:
	parser.error('log file not given')

"""
if options.outDir[-1] != '/':
	options.outDir = options.outDir + '/'
	
if options.tmpDir[-1] != '/':
	options.tmpDir = options.tmpDir + '/'
"""	   

###############################################################################
#val = subprocess.Popen('ls -lh', shell=True,stdout=subprocess.PIPE).communicate()[0]
#print 'val is',val

inputDir = '/'.join(options.bamFile.split('/')[:-2])
print inputDir

outDir = inputDir+'/'+options.sample + '/gVCF_calls'

gatkDir = '/home/jmkidd/kidd-lab/progs/'


dt = subprocess.Popen('date', shell=True,stdout=subprocess.PIPE).communicate()[0]
dt = dt.rstrip()
cmd = 'echo \'%s STARTING %s level %i\' >> %s' % (options.bamFile, dt, options.level, options.log)
#print cmd + '\n\n'
runCMD(cmd,options.log)

#Setup tmpdir to use
tmpDir = '/home/jmkidd/kidd-lab-scratch/shiya-projects/temp'

dbSNP_file = '/home/jmkidd/kidd-lab/genomes/hg19/vqsr-rods/dbsnp_138.b37.refmt.vcf'
#Step 3 Haplotypecaller
if options.level <=0:
	out_vcf = outDir+'/' + options.sample + ".%s.all.sites.vcf" %(options.chr)
	cmd = "java -jar -Djava.io.tmpdir=/home/jmkidd/kidd-lab-scratch/shiya-projects/temp -Xmx4g " + gatkDir + 'GenomeAnalysisTK.jar -T HaplotypeCaller '
	cmd += '-L ' +"/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/process/interval_list/" + options.chr + '.interval_list '
	cmd += '-R /home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 -I ' + options.bamFile + ' --dbsnp '+ dbSNP_file + ' -o ' + out_vcf + ' -stand_call_conf 30 -stand_emit_conf 10'
	runCMD(cmd,options.log)
#	print cmd
	runCMD('gzip '+out_vcf,options.log)
