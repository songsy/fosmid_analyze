#!/usr/bin/env python
# python Affy_phased_snp.py
# Shiya Song
# 28 July 2014
# Deal with phase3_affy/
# format 0000 is p_transmit+p_untransmit+m_transmit+m_untransmit

import sys
import argparse
import gzip
class VcfIterator:
  def __init__(self, filename,name):
	self.file = gzip.open(filename, "r")
  	self.name = name
	self.index = [0,0,0]

  def __iter__(self):
	return self
  
  def __next__(self):
	return self.next()
	
  def next(self):
	line = next(self.file)
	while line[:2] == "##":
		line = next(self.file)
	while line[0] == "#":
		col = line.strip().split('\t')
		for i in range(len(self.name)):
			try:
				self.index[i]=col.index(self.name[i])	
			except ValueError:
				print self.name[i],'not found'	
		line = next(self.file)
	fields = line.strip().split()
	chrom = fields[0]
	while chrom =='MT' or chrom == 'Y':
		line = next(self.file)
		fields = line.strip().split()
		chrom = fields[0]
	if chrom =='MT' or chrom == 'Y':
		pass
	pos = int(fields[1])
	alleles = [fields[3]]
	for alt_a in fields[4].split(","):
		alleles.append(alt_a)
	geno=[(0,0),(0,0),(0,0)]
	for i in range(len(self.name)):
		if chrom == 'X' and len(fields[self.index[i]])==1:
			if fields[self.index[i]][0]=='.':
				return (chrom, pos, tuple(alleles), False)
			else:
				geno[i] = (int(fields[self.index[i]][0]),int(fields[self.index[i]][0]))
		elif len(fields[self.index[i]])==1:
			print line
			if fields[self.index[i]][0]=='.':
				return (chrom, pos, tuple(alleles), False)
			else:
				geno[i] = (int(fields[self.index[i]][0]),int(fields[self.index[i]][0]))
		elif fields[self.index[i]][0]=='.' or fields[self.index[i]][2]=='.':
			return (chrom, pos, tuple(alleles), False)
		else:
			geno[i] = (int(fields[self.index[i]][0]), int(fields[self.index[i]][2]))
	return (chrom, pos, tuple(alleles), geno)

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
	
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='get trio phasing info from Affy')
	parser.add_argument("--vcf", dest='vcf',help="vcf file")
	parser.add_argument("--ped", dest='ped_file',help="ped file")

	args = parser.parse_args()
	family = read_ped(args.ped_file)
	for trio in family:
		print trio
		outfile = trio[0]+'_Affy_snps.bed'
		phase_by_transmission(VcfIterator(args.vcf,trio[1:]),outfile)
	