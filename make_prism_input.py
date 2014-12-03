#!/usr/bin/env python
# python make_prism_input.py
# Shiya Song
# 26th August 2014
# This script generate input file for prism. Especially work on 1000G phase1 reference panel

import sys
import commands
import os
import argparse
from NGS_utils import *
import numpy as np
import pandas as pd
import gzip
from scipy.interpolate import interp1d

def genetic_map(chr,file):
	POS = []
	MAP = []
	interpolate=[]
	f = open(file,"r")
	for line in f:
		if line[:8]=="position":
			continue
		line=line.strip().split(" ")
		pos = int(line[0])
		map = float(line[2])
		POS.append(pos)
		MAP.append(map)
	POS=np.array(POS)
	MAP=np.array(MAP)
	map=pd.Series(MAP,index=POS)
	interpolate=interp1d(POS,MAP)
	return interpolate,map	

def parse_reference_panel(vcf_file,legend_file,haplotype_file,interpolate,map):
	f_out = open(args.out_dir + args.sample+'.chr'+args.chr+'.positions','w')
	pos_list = []
	geno = []
	allele = []
	for i in VcfIterator(vcf_file):
		pos_list.append(i[1])
		geno.append(str(i[3][0])+'/'+str(i[3][1]))
		allele.append(i[2][0]+i[2][1])
	index_to_pos = {}
	f = gzip.open(legend_file,'r')
	for i, line in enumerate(f):
		line = line.strip().split()
		if line[0]=='id':
			continue
		if line[4]=='SNP':
			pos = int(line[1])
			assert line[2] in "ACTG" and line[3] in "ACTG"
			find = find_pos(pos,pos_list)
			if find>=0:
				index_to_pos[i-1]=int(pos)
				if not compare_allele(allele[find],line[2]+line[3]):
					print line,allele[find]
	f = gzip.open(haplotype_file,'r')
	for i, line in enumerate(f):
		try:
			pos = index_to_pos[i]
		except KeyError:
			continue
		find = find_pos(pos,pos_list)
		genotype = geno[find]
		if genotype=='0/0':
			tag = "ref_only"
		else:
			tag = "common"
		try:
			cm = interpolate(pos)
		except ValueError:
			if pos<=map.index[0]:
				cm = map[map.index[0]]
			else:
				cm = map[map.index[-1]]
		print >>f_out,"%s,%s,(%s),%s,%s" %(pos,tag,genotype,cm,line.strip())
			
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='make input file for Prism')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--chr", dest='chr',help="chromosome")
	parser.add_argument("--out_dir", dest='out_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/Prism/',help="directory for output")
	parser.add_argument("--ref_dir", dest='ref_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/',help="directory for reference panel")
	args = parser.parse_args()
	
	wk = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/'
	vcf_file = '%s%s/gVCF_calls/%s.chr%s.vcf.gz' %(wk,args.sample,args.sample,args.chr)
	map_file = args.ref_dir+'genetic_map_chr%s_combined_b37.txt' %(args.chr)
	interpolate,map=genetic_map(args.chr,map_file)
	
	legend_file = args.ref_dir+'ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz' %(args.chr)
	haplotype_file = args.ref_dir+'ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz' %(args.chr)
	parse_reference_panel(vcf_file,legend_file,haplotype_file,interpolate,map)
	