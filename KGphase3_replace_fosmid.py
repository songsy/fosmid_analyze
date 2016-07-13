#!/usr/bin/env python
# python KGphase3_replace_fosmid.py
# Shiya Song
# 2 Feb 2015
# 1000 Genomes phase3 vcf file, fill in the missing data in fosmid haplotypes

import sys
from NGS_utils import *
import argparse
import pandas as pd
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pickle
import os

class VcfIteratorV2:
	def __init__(self, filename):
		if filename[-2:]=='gz':
			self.file = gzip.open(filename, "r")
		else:
			self.file = open(filename, "r")
	def __iter__(self):
		return self  
	def __next__(self):
		return self.next()  	
	def next(self):
		line = next(self.file)
		while line[0] == "#":
			line = next(self.file)
		fields = line.strip().split()
		chrom = fields[0]
		pos = int(fields[1])
		alleles = [fields[3]]
		for alt_a in fields[4].split(","):
			alleles.append(alt_a)
		geno = fields[9][:3]
		phased = geno[1]=='|'
		return (chrom, pos, tuple(alleles), geno)

def create_snp_list(file1):
	snp_pos = []
	snp_geno = []
	snp_phase = []
	snp_phase2 = []
	for i in VcfIteratorV2(file1):
		if i[3][0]!=i[3][2]: 
			snp_pos.append(i[1])				# store the position of the snp, used for fast search 
			snp_geno.append(i[2])		# store the ref and alt allele of the snp
			snp_phase.append(i[3])			# store the phasing information, like ["0011","TTGG","TG"]
			snp_phase2.append(())
	snp = {}
	snp["geno"] = pd.Series(snp_geno,index=snp_pos)
	snp["phase1"] = pd.Series(snp_phase,index=snp_pos)
	snp["phase2"] = pd.Series(snp_phase2,index=snp_pos)
	snp = pd.DataFrame(snp)
	return snp

def read_snp(SNP,file2):		# read in hapmap phasing after liftOver
	same = 0
	opposite = 0
	for i in VcfIteratorV2(file2):
		if i[3][1]=='|' and i[3][0]!=i[3][2]: 
			try:
				info = SNP["geno"][i[1]]
			except KeyError:
				continue
			if SNP["geno"][i[1]]!=i[2]:
				phase = []
				allele = i[2][int(i[3][0])]+i[2][int(i[3][2])]
				if info[0]+info[1]==allele[0]+allele[1] or info[0]+info[1]==allele[1]+allele[0]:
					if int(i[3][0])<2:
						phase.append(i[3][0])
					elif i[3][0]==2:
						phase.append('1')
					if int(i[3][2])<2:
						phase.append(i[3][2])
					elif i[3][0]==2:
						phase.append('1')
					SNP["phase2"][i[1]]='|'.join(phase)
					if info[0]+info[1]==allele[0]+allele[1]:
						same +=1
					elif info[0]+info[1]==allele[1]+allele[0]:
						opposite+=1
			else:
				SNP["phase2"][i[1]]=i[3]
				if SNP["phase1"][i[1]][0]==SNP["phase2"][i[1]][0] and SNP["phase1"][i[1]][2]==SNP["phase2"][i[1]][2] and SNP["phase1"][i[1]][1]=='|' and SNP["phase2"][i[1]][1]=='|':
					same +=1
				elif SNP["phase1"][i[1]][0]==SNP["phase2"][i[1]][2] and SNP["phase1"][i[1]][2]==SNP["phase2"][i[1]][0] and SNP["phase1"][i[1]][1]=='|' and SNP["phase2"][i[1]][1]=='|':
					opposite +=1
	print same,opposite,same+opposite
	flip = 0 if same >= opposite else 1
	return SNP,flip

def write_vcf(wgs_vcf,SNP,flip):
	f = gzip.open(wgs_vcf,'r')
	f_out=gzip.open("%s%s/gVCF_calls/%s.%s.fosmid.1KGphase3.phased.vcf.gz" %(args.wgs_dir,SAMPLE,SAMPLE,CHROM),'w')
	extra = 0
	for line in f:
		line = line.strip()
		if line[0] == "#":
			f_out.write(line+'\n')
			continue
		line1 = line.split('\t')
		chr = line1[0]
		if chr==CHROM:
			ref = line1[3]
			alt = line1[4].split(',')
			pos = int(line1[1])
			info = line1[9]
			info = info.split(':')
			if info[0][1] == "/" and info[0][0]!=info[0][2]:
				try:
					phase = SNP["phase2"][pos]
				except KeyError:
					print >>f_out,line
					continue
				if len(phase)<3:
					print >>f_out,line
					continue
				if phase[1]=='|':
					if flip:
						info[0]=phase[2]+phase[1]+phase[0]
					else:
						info[0]=phase
					extra +=1
					print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
			else:
				print >>f_out,line		
		else:
			print line
	print 'extra:',extra

if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--name",help="Input VCF files")
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	parser.add_argument("--phase_panel", dest='phase_panel',default='1KGphase1',help="reference panel")
	parser.add_argument("--prism", dest='prism',default=0,help="prism or not")
	parser.add_argument("--phase", dest='phase',default=0,help="prism or not")
	args = parser.parse_args()
	
#	for chr in range(1,23):
	for chr in ['X']:
		chr = 'chr'+str(chr)
		SAMPLE = args.name
		CHROM = chr
		if args.prism=='1':
			file1 = '%s%s/gVCF_calls/%s.%s.prism.v2.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
		elif args.name=='NA12878':
			file1 = '%s%s/all_sites/%s.%s.fosmid.phase.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
		else:
			file1 = '%s%s/gVCF_calls/%s.%s.fosmid.v2.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
		if args.phase=='shapeit':
			if args.phase_panel=='1KGphase1':
				file2 = '%s%s/gVCF_calls/%s.%s.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
				if args.name=='NA12878':
					file2 = '%s%s/all_sites/%s.%s.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
			elif args.phase_panel=='1KGphase3':
				file2 = '%s%s/gVCF_calls/%s.%s.1KGphase3.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
				if args.name=='NA12878':
					file2 = '%s%s/all_sites/%s.%s.1KGphase3.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
		else:
			file2 = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/%s/%s.%s.1KGphase3.snp.vcf.gz' %(args.name,args.name,chr)
		print file1,file2
		SNP = create_snp_list(file1)
		SNP,flip = read_snp(SNP,file2)
		print chr,flip
		write_vcf(file1,SNP,flip)