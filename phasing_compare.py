#!/usr/bin/env python
# python phasing_compare.py
# Shiya Song
# 23 Sept 2014
# Comparing phasing statistics in terms of switch error and mean distance between switch errors

import sys
from NGS_utils import *
import argparse
import pandas as pd
import numpy as np
import math

def create_snp_list(file1):
	snp_pos = []
	snp_geno = []
	snp_phase = []
	snp_phase2 = []
	for i in VcfIterator(file1):
		if i[4] is True and i[3][0]!=i[3][1]: 
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
	for i in VcfIterator(file2):
		if i[4] is True and i[3][0]!=i[3][1]: 
			try:
				info = SNP["geno"][i[1]]
			except KeyError:
				continue
			if SNP["geno"][i[1]]!=i[2]:
				phase = []
				allele = i[2][int(i[3][0])]+i[2][int(i[3][1])]
				if info[0]+info[1]==allele[0]+allele[1] or info[0]+info[1]==allele[1]+allele[0]:
					if i[3][0]<2:
						phase.append(i[3][0])
					elif i[3][0]==2:
						phase.append(1)
					if i[3][1]<2:
						phase.append(i[3][1])
					elif i[3][0]==2:
						phase.append(1)
					SNP["phase2"][i[1]]=tuple(phase)
				else:
					print i,SNP["geno"][i[1]]
			else:
				SNP["phase2"][i[1]]=i[3]
	return SNP

def compare_phasing(SNP,switch_distance,switch_length):
	list = SNP.index
	tot = 0
	opposite = 0
	is_switch = []
	switch = []
	switch_position = []
	for j in range(len(list)):
		if SNP["phase1"][list[j]]!=() and SNP["phase2"][list[j]]!=():
			tot +=1
			if SNP["phase1"][list[j]]!=SNP["phase2"][list[j]]:
#				print list[j],SNP["phase1"][list[j]],SNP["phase2"][list[j]]
				opposite +=1
				is_switch.append(1)
				switch.append(0)
				switch_position.append(list[j])
			else:
				is_switch.append(0)
				switch.append(0)
				switch_position.append(list[j])
	k = 0
	s = 0
	e = 0
	switch_error = 0
	last_error_pos = None
	if is_switch.count(0)>=is_switch.count(1):
		for k in range(len(is_switch)):
			if k==0:
				last = is_switch[k]
				s = k
			elif is_switch[k]==last:
				continue
			else:
				e = k
				if is_switch[s]==1:
					switch_error +=1
					switch_length.append(int(switch_position[e-1])-int(switch_position[s]))
					if last_error_pos is None:
						last_error_pos = switch_position[s]
					else:
						switch_distance.append(switch_position[s]-last_error_pos)
						last_error_pos = switch_position[s]
				s = k
				last = is_switch[k]
		e = k
		if is_switch[s]==1:
			switch_error +=1
			switch_length.append(int(switch_position[e-1])-int(switch_position[s]))
			if last_error_pos is None:
				last_error_pos = switch_position[s]
			else:
				switch_distance.append(switch_position[s]-last_error_pos)
				last_error_pos = switch_position[s]
	else:
		for k in range(len(is_switch)):
			if k==0:
				last = is_switch[k]
				s = k
			elif is_switch[k]==last:
				continue
			else:
				e = k
				if is_switch[s]==0:
					switch_error +=1
					switch_length.append(int(switch_position[e-1])-int(switch_position[s]))
					if last_error_pos is None:
						last_error_pos = switch_position[s]
					else:
						switch_distance.append(switch_position[s]-last_error_pos)
						last_error_pos = switch_position[s]
				s = k
				last = is_switch[k]
		e = k
		if is_switch[s]==0:
			switch_error +=1
			switch_length.append(int(switch_position[e-1])-int(switch_position[s]))
			if last_error_pos is None:
				last_error_pos = switch_position[s]
			else:
				switch_distance.append(switch_position[s]-last_error_pos)
				last_error_pos = switch_position[s]
	print chr,tot,opposite,switch_error,len(switch_distance),len(switch_length),np.mean(switch_distance),np.mean(switch_length)
	return switch_distance,switch_length,tot,opposite,switch_error

if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--name",help="Input VCF files")
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	args = parser.parse_args()
	
	
	switch_distance = []
	switch_length = []
	tot = 0
	different = 0
	switch_error = 0
	for chr in range(1,23):
		chr = 'chr'+str(chr)
#		print chr
		file1 = '%s%s/gVCF_calls/%s.%s.prism.v2.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
		file2 = '%s%s/gVCF_calls/%s.%s.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
#		file2 = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/%s/%s.%s.1KGphase3.snp.vcf.gz' %(args.name,args.name,chr)
#		file2 = '/share/jmkidd/songsy/complete-genomics/YRI_trio/%s.%s.phased.vcf.gz' %(args.name,chr)
		print file1,file2
		SNP = create_snp_list(file1)
		SNP = read_snp(SNP,file2)
		switch_distance,switch_length,a,b,c=compare_phasing(SNP,switch_distance,switch_length)
		if b<a-b:
			different +=b
		else:
			different += a-b
		tot +=a
		switch_error += c
	print tot,different,switch_error
			
			
