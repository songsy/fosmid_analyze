#!/usr/bin/env python
# python read_1000G_hapmap.py
# Shiya Song
# 12th July 2014
# This script is used to access switch error comparing with hapmap trio phased data

import gzip
import sys
import os
from NGS_utils import *

def find_pos_v2(pos,list):
	left = 0
	right = len(list)
	find = -1
	while right - left >0:
		midpoint = (right + left)/2
		if pos < list[midpoint][0]:
			right = midpoint
		elif pos > list[midpoint][0]:
			left = midpoint+1
		elif pos == list[midpoint][0]:
			left = midpoint
			find = midpoint
			break
	return find

def compare_allele(nuc1,nuc2):
	if nuc1[0]==nuc2[0] and nuc1[1]==nuc2[1]:
		return True

def compare_allele_v2(nuc1,nuc2):
	if nuc1[0]==nuc2[1] and nuc1[0]==nuc2[1]:
		return True
		
def fosmid_replace_phase3(phase3_file,locus):
	if phase3_file[-2:]=='gz':
		f=gzip.open(phase3_file, "r")
	else:
		f=open(phase3_file,"r")
	f_out1=gzip.open(phase3_file.replace('_phased_snps.vcf.gz','_first_phased_snps.vcf.gz'),"w")
	not_share = 0
	same = 0
	reverse = 0
	unphased = 0
	different = 0
	for line in f:
		if line[0]=='#':
			f_out1.write(line)
			continue
		line = line.strip().split('\t')
		if line[0][0]!='c':
			chr = 'chr'+line[0]
		else:
			chr = line[0]
		pos = int(line[1])
		ref=line[3]
		alt=line[4].split(',')
		GT = line[9].split(':')[0]
		line[9] = GT
		line[8]='GT'
		if GT[0]!=GT[2]:
			find=find_pos_v2(pos,locus)
			if find<0:
				not_share +=1
				f_out1.write('\t'.join(line)+'\n')
				continue
			else:
			 	allele = ''
			 	if GT[0]=='0':
			 		allele += ref
			 	else:
			 		allele += alt[int(GT[0])-1]
				if GT[2]=='0':
			 		allele += ref
			 	else:
			 		allele += alt[int(GT[2])-1]
			 	if GT[1]=='|':
					if compare_allele(allele,locus[find][2]):
						same +=1
					elif compare_allele_v2(allele,locus[find][2]):
						print 'reverse',pos,allele,GT,locus[find]
						reverse +=1
						line[9]= GT[2]+'|'+GT[0]
					else:
						print 'different',pos,allele,GT,locus[find]
						different +=1
						f_out1.write('\t'.join(line)+'\n')
						if locus[find][1][0]=='1' and locus[find][1][1]=='0':
							line[4]=locus[find][2][0]
							line[9]='|'.join(locus[find][1])
						elif locus[find][1][1]=='1' and locus[find][1][0]=='0':
							line[4]=locus[find][2][1]
							line[9]='|'.join(locus[find][1])
						elif locus[find][1]=='12' or locus[find][1]=='21':
							if locus[find][1][1]=='1':
								assert ord(locus[find][2][0])<ord(locus[find][2][1])
								line[4]=locus[find][2][0]+','+locus[find][2][1]
							elif locus[find][1][2]=='1':
								assert ord(locus[find][2][1])<ord(locus[find][2][0])
								line[4]=locus[find][2][1]+','+locus[find][2][0]
							line[9]='|'.join(locus[find][1])
					f_out1.write('\t'.join(line)+'\n')
				elif GT[1]=='/':
					unphased +=1
					if locus[find][1][0]=='1' and locus[find][1][1]=='0':
						line[4]=locus[find][2][0]
						line[9]='|'.join(locus[find][1])
					elif locus[find][1][1]=='1' and locus[find][1][0]=='0':
						line[4]=locus[find][2][1]
						line[9]='|'.join(locus[find][1])
					elif locus[find][1]=='12' or locus[find][1]=='21':
						if locus[find][1][1]=='1':
							assert ord(locus[find][2][0])<ord(locus[find][2][1])
							line[4]=locus[find][2][0]+','+locus[find][2][1]
						elif locus[find][1][2]=='1':
							assert ord(locus[find][2][1])<ord(locus[find][2][0])
							line[4]=locus[find][2][1]+','+locus[find][2][0]
						line[9]='|'.join(locus[find][1])
					f_out1.write('\t'.join(line)+'\n')
		else:
			f_out1.write('\t'.join(line)+'\n')
	return not_share,same,reverse,different,unphased

if __name__=="__main__":

	chr = sys.argv[1]
	fosmid_vcf = '/home/jmkidd/kidd-lab-scratch/shiya-projects/phased-genomes/NA19240.vcf.gz'
	locus=vcf_to_bed_by_chr(fosmid_vcf,chr,0)

	phase3_file = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/BLOCK/NA19240/NA19240_%s_phased_snps.vcf.gz' %(chr)
	not_share,same,reverse,different,unphased=fosmid_replace_phase3(phase3_file,locus)
	print chr,not_share,same,reverse,different,unphased
