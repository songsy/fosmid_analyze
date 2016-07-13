#!/usr/bin/env python
# python SNP_stats_from_vcf.py
# Shiya Song
# 23 July 2014
# 

import sys
import argparse
import gzip

class VcfIterator:
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
	filter = fields[6]
	if filter=='PASS':
		filter=True
	else:
		filter=False
	return (chrom, pos, tuple(alleles), (int(geno[0]), int(geno[2])), filter)

def calc_vcf_stats(m):
	snp=[0,0,0,0,0,0]	  #tot,tot_filter,homo,homo_filter,het,het_filter
	indel=[0,0,0,0,0,0]
	for i in m:
		chrom = i[0]
		if chrom=='chrY' or chrom=='EBV' or chrom=='pCC1FOS' or chrom=='eColiK12':
			continue
		a=i[2][i[3][0]]
		b=i[2][i[3][1]]
		if len(a)==1 and len(b)==1:
			snp[0]+=1
			if i[3][0]==i[3][1]:
				snp[2]+=1
				c = 3
			else:
				snp[4]+=1
				c = 5
			if i[4] is True:
				snp[1]+=1
				snp[c]+=1				
		else:
			indel[0]+=1
			if i[3][0]==i[3][1]:
				indel[2]+=1
				c = 3
			else:
				indel[4]+=1
				c = 5
			if i[4] is True:
				indel[1]+=1
				indel[c]+=1
	return snp,indel	
			
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='summarize snp stats from vcf')
	parser.add_argument("--vcf", dest='vcf',help="vcf file")
	parser.add_argument("--sample", dest='sample',help="sample name")
	args = parser.parse_args()
	snp,indel=calc_vcf_stats(VcfIterator(args.vcf))
	f_out=open('/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/wgs_snp_stats.txt','a')
	print >>f_out,args.sample
	print >>f_out,'\t'.join(map(str,snp))
	print >>f_out,'\t'.join(map(str,indel))