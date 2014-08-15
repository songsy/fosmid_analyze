#!/usr/bin/env python
# Shiya Song
# 24rd June 2014
# indel_phase.py
# from hap1_hap2_unknown_wgs_indel.vcf file, together with wgs indel calls, get phased indels, try to infer

import sys
import re
import os
import signal
import glob
import math
import pickle
from NGS_utils import *

chromOrder = get_chromOrder("human")
def chrom_list():
	chr_list=[]
	for i in range(1,23):
		chr_list.append("chr"+str(i))
	chr_list.append("chrX")
	return chr_list
chr_list = chrom_list()

def read_whole_genome_vcf(whole_genome_vcf):
	f = open(whole_genome_vcf,"r")
	whole_indel = {}
	for i in chr_list:
		whole_indel[i]=[]
	for line in f:
		if line[0]=="#":
			continue
		line = line.strip().split("\t")
		chr = line[0]
		pos = int(line[1])
		ref = line[3]
		alt = line[4].split(",")
		geno = line[9][:3]
		if chr not in chr_list:
			break
		whole_indel[chr].append([pos,ref,alt,geno])
	return whole_indel

def find_pos(pos,list):
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

def compare_vcf_wgs_INDEL(vcf_file,whole_indel,INDEL_same,INDEL_whole):
	gc = 'gunzip -c ' + vcf_file
	f = os.popen(gc, 'r')
	tot = 0
	not_shown= 0
	hap1_indel = 0
	hap2_indel = 0
	homo_indel = 0
	het = 0
	not_indel =0
	biallele_indel = 0
	unknown_indel = 0
	for line in f:
		if line[0]=="#":
			continue
		tot +=1
		line = line.strip().split("\t")
		chr = line[0]
		pos = int(line[1])
		ref = line[3]
		alt = line[4].split(",")
		find = find_pos(pos,whole_indel[chr])
		hap1=line[9]
		hap2=line[10]
		unknown=line[11]
		find2 = find_pos(pos,INDEL_same[chr])
		if find <0:
			print >>f_log,"not_found",line
		else:	
			if hap1[:3]=="./." and hap2[:3]=="./." and unknown[:3]=="./.":
				not_shown +=1
			elif hap1[:3]=="1/1" and hap2[:3]=="0/0" and (unknown[:3]=="./." or unknown[:3]=="0/0"):
				hap1_indel +=1
				if find2>=0:
					INDEL_whole.append([pos,ref,alt,INDEL_same[chr][find2][3],whole_indel[chr][find][3],"1|0"])			#KG,whole_genome,fosmid
				else:
					INDEL_whole.append([pos,ref,alt,"",whole_indel[chr][find][3],"1|0"])			#KG,whole_genome,fosmid
			elif hap2[:3]=="1/1" and hap1[:3]=="0/0" and (unknown[:3]=="./." or unknown[:3]=="0/0"):
				hap2_indel +=1
				if find2>=0:
					INDEL_whole.append([pos,ref,alt,INDEL_same[chr][find2][3],whole_indel[chr][find][3],"0|1"])			#KG,whole_genome,fosmid
				else:
					INDEL_whole.append([pos,ref,alt,"",whole_indel[chr][find][3],"0|1"])			#KG,whole_genome,fosmid						
			elif hap1[:3]=="1/1" and hap2[:3]=="./." and (unknown[:3]=="./." or unknown[:3]=="0/0"):
				geno = ''
				if whole_indel[chr][find][3]=="0/1":
					hap1_indel +=1
					geno = '1|0'										
				elif whole_indel[chr][find][3]=="1/1":
					homo_indel +=1
					geno = '1|1'
				else:
					geno = ''
				if geno!='':
					if find2>=0:
						INDEL_whole.append([pos,ref,alt,INDEL_same[chr][find2][3],whole_indel[chr][find][3],geno])			#KG,whole_genome,fosmid
					else:
						INDEL_whole.append([pos,ref,alt,"",whole_indel[chr][find][3],geno])			#KG,whole_genome,fosmid			
			elif hap2[:3]=="1/1" and hap1[:3]=="./." and (unknown[:3]=="./." or unknown[:3]=="0/0"):
				geno = ''
				if whole_indel[chr][find][3]=="1/0":
					hap2_indel +=1
					geno = '0|1'										
				elif whole_indel[chr][find][3]=="1/1":
					homo_indel +=1
					geno = '1|1'
				else:
					geno = ''
				if geno!='':
					if find2>=0:
						INDEL_whole.append([pos,ref,alt,INDEL_same[chr][find2][3],whole_indel[chr][find][3],geno])			#KG,whole_genome,fosmid
					else:
						INDEL_whole.append([pos,ref,alt,"",whole_indel[chr][find][3],geno])			#KG,whole_genome,fosmid			
			elif hap1[:3]=="1/1" and hap2[:3]=="1/1":
				homo_indel +=1
				if find2>=0:
					INDEL_whole.append([pos,ref,alt,INDEL_same[chr][find2][3],whole_indel[chr][find][3],"1|1"])			#KG,whole_genome,fosmid
				else:
					INDEL_whole.append([pos,ref,alt,"",whole_indel[chr][find][3],"1|1"])			#KG,whole_genome,fosmid
			elif (hap1[:3]=="1/1" and hap2[:3]=="2/2") or (hap2[:3]=="1/1" and hap1[:3]=="2/2") and unknown[:3]=="./.":
				biallele_indel +=1
				geno = ''
				if hap1[:3]=="1/1":
					geno = '1|2'
				if hap1[:3]=="2/2":
					geno = '2|1'
				if find2>=0:
					INDEL_whole.append([pos,ref,alt,INDEL_same[chr][find2][3],whole_indel[chr][find][3],geno])			#KG,whole_genome,fosmid
				else:
					INDEL_whole.append([pos,ref,alt,"",whole_indel[chr][find][3],geno])			#KG,whole_genome,fosmid
			elif hap1[:3]=="0/1" or hap2[:3]=="0/1" or hap1[:3]=="1/2" or hap2[:3]=="1/2" or hap1[:3]=="0/2" or hap2[:3]=="0/2":
				geno = ''
				if whole_indel[chr][find][3]=="0/1":
					if hap1[:3]=="0/1" and hap2[:3]=='0/0':
						geno = '1|0'
						hap1_indel +=1
					elif hap1[:3]=="0/1" and hap2[:3]=='1/1':
						geno = '0|1'
						hap2_indel +=1
					elif hap2[:3]=="0/1" and hap1[:3]=='0/0':
						geno = '0|1'
						hap2_indel +=1
					elif hap2[:3]=="0/1" and hap1[:3]=='1/1':
						geno = '1|0'
						hap1_indel +=1
				elif whole_indel[chr][find][3]=="1/2":
					if hap1[:3]=="1/2" and hap2[:3]=='1/1':
						geno = '2|1'
						biallele_indel +=1
					elif hap1[:3]=="1/2" and hap2[:3]=='2/2':
						geno = '1|2'
						biallele_indel +=1
					elif hap2[:3]=="1/2" and hap1[:3]=='1/1':
						geno = '1|2'
						biallele_indel +=1
					elif hap2[:3]=="1/2" and hap1[:3]=='2/2':
						geno = '2|1'
						biallele_indel +=1
				elif whole_indel[chr][find][3]=="1/1":
					het +=1
					geno = ''
				else:
					geno = ''
				if geno!='':
					if find2>=0:
						INDEL_whole.append([pos,ref,alt,INDEL_same[chr][find2][3],whole_indel[chr][find][3],geno])			#KG,whole_genome,fosmid
					else:
						INDEL_whole.append([pos,ref,alt,"",whole_indel[chr][find][3],geno])			#KG,whole_genome,fosmid						
			elif hap1[:3]=="0/0" or hap2[:3]=="0/0" or unknown[:3]=="0/0":
				not_indel +=1
			elif unknown[:3]=="1/1" or unknown[:3]=="1/2" or unknown[:3]=="2/2" or unknown[:3]=="0/2" or unknown[:3]=="0/1":
				unknown_indel +=1
	m = [tot,hap1_indel,hap2_indel,homo_indel,biallele_indel,not_shown,not_indel,unknown_indel,het]
	return m,INDEL_whole

if __name__=="__main__":	
	whole_genome_vcf = "../../../merge/NA19240_combined_indel.vcf"
	whole_indel = read_whole_genome_vcf(whole_genome_vcf)
	info = {}
	dbfile=open('NA19240_indel_same_wgs_pickle','rb')    # This is from 1000G phased indels
	INDEL_same=pickle.load(dbfile)
	f=open("NA19240_indel_stats_1KG_v5.txt","a")
	INDEL_whole = {}
	for chr in chr_list:
		INDEL_whole[chr]=[]
		vcf_file="NA19240_hap1_hap2_unknown_wgs_indel_"+chr+".vcf.gz"
		print >>f,chr,"\tsame with 1kG\t",len(INDEL_same[chr])
		m,INDEL_whole[chr]=compare_vcf_wgs_INDEL(vcf_file,whole_indel,INDEL_same,INDEL_whole[chr])
		print >>f,chr,"\t","tot,hap1_indel,hap2_indel,homo_indel,biallele_indel,not_shown,not_indel,unknown_indel,het"
		print >>f,chr,"\t","\t".join(map(str,m))
		info[chr]=m
	stats = [0,0,0,0,0,0,0,0,0]  #	tot,hap1_indel,hap2_indel,homo_indel,biallele_indel,not_shown,not_indel,unknown_indel,het
	for i in info.keys():
		for j in range(len(stats)):
			stats[j] += info[i][j]
	print >>f,"tot,hap1_indel,hap2_indel,homo_indel,biallele_indel,not_shown,not_indel,unknown_indel,het"
	print >>f,"\t".join(map(str,stats))
	dbfile=open('NA19240_fosmid_of_wgs_indel_infer_pickle','wb')    # pos,ref,alt,KG,whole_genome,fosmid
	pickle.dump(INDEL_whole,dbfile)

	