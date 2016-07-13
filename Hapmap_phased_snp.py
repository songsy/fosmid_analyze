#!/usr/bin/env python
# python Hapmap_phased_snp.py
# Shiya Song
# 30 July 2014
# Deal with Hapmap3 phased data, either from unrelated individuals or trios, liftover to hg19
# python /home/songsy/script/fosmid_analyze/Hapmap_phased_snp.py --pop YRI > YRI_phased_snp.txt

import sys
import argparse
import gzip
import os
import subprocess

class BedIterator:
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
	while line[:4] == "rsID":
		line = next(self.file)
	fields = line.strip().split()
	return tuple(fields)

def chrom_list():
	chr_list=[]
	for i in range(1,24):
		chr_list.append("chr"+str(i))
	return chr_list
chr_list = chrom_list()

def phase_by_transmission(m,outfile):
	f_out=open(outfile,'w')
	for i in m:
		if i[3] is False:
			continue
		index=get_phase(i[3])
		if index is False:
			print 'Wrong:','\t'.join(map(str,i))
			continue
		elif index == 'undetermine':
			print 'Undetermine:','\t'.join(map(str,i))
			continue
		try:
			assert "".join(map(str,i[3][0]))==str(index[0])+str(index[2]) or "".join(map(str,i[3][0]))==str(index[2])+str(index[0])
		except AssertionError:
			print 'AssertionError:%s\t%i\t%i\t%s\t%s' %(i[0],i[1]-1,i[1],"".join(map(str,index)),i[3])
		genotype = [i[2][j] for j in index]
		print >>f_out,'%s\t%i\t%i\t%s\t%s' %(i[0],i[1]-1,i[1],"".join(map(str,index)),"".join(genotype))

def get_phase(geno):
	child = geno[0]
	dad = geno[1]
	mom = geno[2]
	paternal = None
	maternal = None
	if child[0]!=child[1]:
		if child[0] in dad and child[0] not in mom:
			m=dad.index(child[0])
			paternal = (dad[m],dad[1-m])
			maternal = mom
		elif child[0] in mom and child[0] not in dad:
			m=mom.index(child[0])
			maternal = (mom[m],mom[1-m])
			paternal = dad
		if child[1] in dad and child[1] not in mom:
			m=dad.index(child[1])
			paternal = (dad[m],dad[1-m])
			maternal = mom
		elif child[1] in mom and child[1] not in dad:
			m=mom.index(child[1])
			maternal = (mom[m],mom[1-m])
			paternal = dad
		if child==dad and child==mom:
			return 'undetermine'
		if paternal is None:
			return False
	else:
		try:
			m=dad.index(child[0])
			paternal = (dad[m],dad[1-m])
			m=mom.index(child[0])
			maternal = (mom[m],mom[1-m])
		except ValueError:
			return False
	return paternal+maternal

def read_ped(file):
	family = []
	f=open(file,'r')
	for line in f:
		line = line.strip().split()
		family.append((line[0],line[1],line[2],line[3]))
	return family

def runCMD(cmd):
	val = subprocess.Popen(cmd, shell=True).wait()
	if val == 0:
		cmd = cmd
	else:
		print 'command failed:',cmd

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='get Hapmap phased data')
	parser.add_argument("--pop", dest='population',help="population")
	parser.add_argument("--dir", dest='dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/',help="directory")
	parser.add_argument("--ped", dest='ped_file',help="ped file")
	
	args = parser.parse_args()
	f_out = open('%s_phased_sites.bed' %(args.population),'w')
	dict = {}
	for chr in chr_list:
		file = 'hapmap3_r2_b36_fwd.consensus.qc.poly.CHROM_POP.unr.phased.gz'
		file = file.replace("CHROM",chr)
		file = file.replace("POP",args.population.lower())
		if chr=='chr23':
			chr = 'chrX'
		dict[chr]={}
		for info in BedIterator(file):
			print >>f_out,"%s\t%i\t%i\t%s" %(chr,int(info[1])-1,int(info[1]),info[0])
			dict[chr][info[0]]=" ".join(info[2:])
#	cmds = 'liftOver %s /home/jmkidd/kidd-lab/genomes/lift-over/hg18ToHg19.over.chain %s unmapped.bed' %('%s_phased_sites.bed' %(args.population),'%s_phased_sites_hg19.bed' %(args.population))
#	runCMD(cmds)

	for info in BedIterator('%s_phased_sites_hg19.bed' %(args.population)):
		print '%s %s %s %s' %(info[0],info[3],info[2],dict[info[0]][info[3]])
	
	
	