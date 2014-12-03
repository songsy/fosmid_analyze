#!/usr/bin/env python
# python assign_parental_in_block_gvcf_v2.py
# Shiya Song
# 29 July 2014
# Assign parental allele within each blocks, calculate switch error, correct some switch errors, make matrix2png plots, assign clone into parental alleles
# Used for NA19240

import argparse
import os
import signal
import glob
import pandas as pd
import numpy as np
import math
from NGS_utils import *
import pickle
import subprocess
import genutils

chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
chromLens = genutils.read_chrom_len(chromLenFile)
chromOrder = get_chromOrder("human")

def create_snp_list(all_var_file):
	f = open(all_var_file,"r")
	snp_pos = []
	snp_geno = []
	snp_phase = []
	for line in f:
		line = line.strip().split("\t")
		chr = line[0]
#		if chr != "chr1" and chr!="chr2":
#			break
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			snp_pos.append(line[1])				# store the position of the snp, used for fast search 
			snp_geno.append(line[3]+line[4])		# store the ref and alt allele of the snp
			snp_phase.append("")			# store the phasing information, like ["0011","TTGG","TG"]
		else:
			break
	snp = {}
	SNP={}	
	snp["geno"] = pd.Series(snp_geno,index=snp_pos)
	snp["phase"] = pd.Series(snp_phase,index=snp_pos)
	snp["phase_1000G"] = pd.Series(snp_phase,index=snp_pos)
	SNP[CHROM]= pd.DataFrame(snp)
	return SNP

def SNPs_not_shown(file):
	f = open(file,"r")
	SNP_not_shown = {}
	SNP_not_shown[CHROM]=[]
	for line in f:
		line = line.strip().split("\t")
		chr = line[1]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			SNP_not_shown[line[1]].append(int(line[2]))
		else:
			break
	return SNP_not_shown

def read_haplotype(hap_file,SNP_not_shown):
	f = open(hap_file,"r")
	BLOCK = {}
	block_list=[]
	decode = {}
	mec = {}
	start = True
	for i in range(1,2800000):
		decode[i] = 0
	for line in f:
		line = line.strip()
		if line[0]=="B":
			start = True
		elif line[0]=="*":
			if start == False:
				block_info = {}
				block_info["hap"] = pd.Series(hap,index=block)
				BLOCK[chr+" "+pos] = pd.DataFrame(block_info)
		else:
			col = line.split("\t")
			if start == True:
				pos = col[4]
				chr = col[3]
				if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
					continue
				elif chr==CHROM:
					find = find_pos(int(pos),SNP_not_shown[chr])
					if find <0:			# SNPs are shown in the matrix
						start = False			
						pos = col[4]
						block = []				# store the position in each block
						block_list.append(chr+" "+pos)		# store the first position of each block
						hap = []				# store the 1/0 in each block
						mec[chr+" "+pos] = [0,0,0,0,0,0,0,0,0,0]		# number of snp, number of snp from Affy, number of support from Affy, MEC value before, MEC value after correction, switch_error, switch_error_corrected, source type
					else:		# SNPs not shown in the matrix
						continue
				else:
					break
			find = find_pos(int(col[4]),SNP_not_shown[chr])
			if find >=0:
				continue
			if col[1] == "-":
				block.append(col[4])
				hap.append(0.5)
				if decode[int(col[0])]==0:
					decode[int(col[0])] = chr +" "+ pos
			else:
				block.append(col[4])	  # just the position of the snp
				hap.append(int(col[1]))			# "1" or "0" for the first haplotype
				decode[int(col[0])] = chr +" "+ pos		# link the refhap result (# of the snp) to the first snp of that block	
	f.close()
	return BLOCK,block_list,decode,mec
	
def decode_number_pos(file):
	f = open(file,"r")
	decode_pos = {}
	decode_ref_alt = {}
	for i,l in enumerate(f):
		line = l.strip().split("\t")
		decode_pos[i+1] = line[1]
#		decode_ref_alt[line[0]+" "+line[1]]=line[3]+line[4]
#	return decode_pos,decode_ref_alt
	return decode_pos

def read_matrix(matrix_file,BLOCK,decode,decode_pos,mec):
	f_mat = open(matrix_file,"r")
	MEC = 0
	done = 0
	for line in f_mat:
		done += 1
		line = line.strip()
		col = line.split(" ")
		name = col[1]
		if SAMPLE=='NA19240':
			chr = ("_").join(col[1].split("_")[3:-1])
			pool = ("_").join(col[1].split("_")[:3])
		else:
			chr = ("_").join(col[1].split("_")[-2:-1])
			pool = ("_").join(col[1].split("_")[:-2])
		pos = col[1].split("_")[-1]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr!=CHROM:
			break
		position = []
		new_position = []
		quality = col[-1]
		qual = [ord(m)-33 for m in quality]
		qual1 = 0
		qual2 = 0
		allele = []
		match = 0
		not_match = 0
		not_found = 0
		hap1 = False
		hap = ""
		assign = False
		for i in range(int(col[0])):
			num = int(col[2+2*i])
			string = col[3+2*i]
			for j in range(len(string)):
				position.append(num+j)
				allele.append(string[j])
		first = True
		for i in range(len(position)):
			if decode[position[i]] == 0:
				not_found +=1
				print "not_found",position[i]
				continue
			assign = True
			real_position = decode_pos[position[i]]
			new_position.append(real_position)			# the actual position of each clone
			if first == True:
				code = decode[position[i]]
				if BLOCK[code]["hap"][real_position]!=0.5:
					first = False
				else:
					not_found+=1
				if BLOCK[code]["hap"][real_position] == int(allele[i]):
					match += 1
					qual1+= qual[i]
				else:
					not_match +=1
					qual2+= qual[i]
			else: 
				if BLOCK[code]["hap"][real_position] == int(allele[i]):
					match += 1
					qual1+= qual[i]
					assert code == decode[position[i]]
				elif BLOCK[code]["hap"][real_position]==0.5:
					not_found+=1
				else:
					not_match += 1
					qual2+= qual[i]
					assert code == decode[position[i]]
#		print done,name,len(position),not_found,match
		if done % 10000 == 0:
			print done,"done"
		if assign == False:
			continue
		if match > not_match:
			hap1 = True
			mismatch = not_match
		elif match == not_match:
			if qual1 >= qual2:
				hap1 = True
				mismatch = not_match
			else:
				hap1 = False
				mismatch = match
		else:
			hap1 = False
			mismatch = match
		if hap1 is True:
			name += "_1"
		else:
			name += "_2"
		try:
			BLOCK[code][name]= pd.Series(allele,index=new_position)
		except AssertionError:
			print line,allele,new_position,position
		BLOCK[code][name+"_qual"]= pd.Series([m for m in quality],index=new_position)
		MEC += mismatch
		mec[code][3] += mismatch
	return MEC,mec,BLOCK

def make_matrix2png(BLOCK,block_list,mec):
	list = ["hap","transmit","untransmit"]
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		pos = block_list[i].split(" ")[1]
		f_matrix = open("BLOCK/%s_gvcf_matrix2png/%s_" %(SAMPLE,SAMPLE)+chr+"_"+pos+".matrix","w")
		match = 0
		switch_match = 0
		switch = False
		print >>f_matrix,"locus\t","\t".join(position for position in BLOCK[block_list[i]].index)
		print >>f_matrix,"refhap_hap1\t","\t".join(map(str,(BLOCK[block_list[i]]["hap"][position] for position in BLOCK[block_list[i]].index)))
		print >>f_matrix,"refhap_hap2\t","\t".join(map(str,((1-BLOCK[block_list[i]]["hap"][position]) for position in BLOCK[block_list[i]].index)))
		string1,string2= [],[]
		for j in BLOCK[block_list[i]].index:
			string1.append("-")
			string2.append("-")
		mec[block_list[i]][0]=len(BLOCK[block_list[i]].index)
		print >>f_matrix,"transmit\t","\t".join(map(str,(item for item in string1)))
		print >>f_matrix,"untransmit\t","\t".join(map(str,(item for item in string2)))
		for each_clone in BLOCK[block_list[i]]:
			if each_clone in list or each_clone[-5:]=="_qual":
				continue
			string = []
			for j in BLOCK[block_list[i]].index:
				if math.isnan(float(BLOCK[block_list[i]][each_clone][j])) is True:
					string.append("-")
				else:
					string.append(BLOCK[block_list[i]][each_clone][j])
			print >>f_matrix,each_clone,"\t","\t".join(map(str,(item for item in string)))
		print chr+"_"+pos+".matrix done"
		f_matrix.close()					
	return mec

def run_matrix2png():
	wk = "BLOCK/%s_gvcf_matrix2png/" %(SAMPLE)
	for i in glob.glob(wk + "*.matrix"):
		j = i.replace(wk,"")
		j = j.replace(".matrix",".png")
		cmd = "matrix2png -data %s -size 10:10 -map 18 -r -c -s -d -range 0:1 > %s" %(i,wk+j)
		os.popen(cmd)
		
def make_track(hap_file,block_list,mec):
	f = open(hap_file,"r")
	f_out = open("BLOCK/%s_gvcf/%s_refhap_track_info_" %(SAMPLE,SAMPLE)+CHROM+".txt","w")
	for line in f:
		line = line.strip()
		if line[0]=="B":
			start = 0
			end = 0
			chr = ""
			continue
		elif line[0] == "*":
			print >>f_out,"%s\t%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i" %(chr,start,end,MEC[7],MEC[0],MEC[1],MEC[2],MEC[3],MEC[4],MEC[5],MEC[6],MEC[8],MEC[9])
		elif start == 0:
			chr = line.split("\t")[3]
			if chr!=CHROM:
				break
			start = line.split("\t")[4]
			find = find_pos(int(start),SNP_not_shown[chr])
			if find <0:  # SNP shows in the matrix
				end = line.split("\t")[4]
				MEC = mec[chr + " " + start]
				mec[chr + " " + start][7]="None"
			else:
				start = 0
		elif start != "":
			end = line.split("\t")[4]		

def make_prism_input(BLOCK,mec):
	f = open("Prism/%s/%s.%s.positions" %(SAMPLE,SAMPLE,CHROM),'r')
	f_track = "BLOCK/%s_gvcf/%s_refhap_track_info_" %(SAMPLE,SAMPLE)+CHROM+".txt"
	position = []
	for line in f:
		position.append(int(line.strip().split(',')[0]))
	snp_decode={}
	f_out1 = open("Prism/%s/%s.%s.local.blocks.all" %(SAMPLE,SAMPLE,CHROM),"w")
#	f_out2 = open("Prism/%s/%s.%s.local.blocks" %(SAMPLE,SAMPLE,CHROM),"w")
	f_out3 = open("Prism/%s/%s.%s.local.interleaving.blocks" %(SAMPLE,SAMPLE,CHROM),"w")
	f_out4 = open("Prism/%s/%s.%s.info" %(SAMPLE,SAMPLE,CHROM),"w")
	id = 0
	last_item = ""
	prev_none = False
	prev_range=None
	current_range=None
	first = True
	for i in BedIterator(f_track):
		interleaving = False
		chr = i[0]
		pos = i[1]
		pos2 = i[2]
		block_name = chr+" "+pos
		if first is True:
			prev_range=(pos,pos2)
			first = False
		else:
			current_range=(pos,pos2)
			if int(current_range[0])<int(prev_range[1]):
				interleaving = True
			prev_range=(pos,pos2)
		last_pos = pos
		string = []
		for j in BLOCK[block_name].index:
			if BLOCK[block_name]["hap"][j] !=0.5:
				find = find_pos(int(j),position)
				if find <0:
					allele1 = int(BLOCK[block_name]["hap"][j])
					allele2 = 1-int(BLOCK[block_name]["hap"][j])
					snp_decode[j]=str(allele1)+'|'+str(allele2)
					continue
				last_pos = j
				string.append(str(j)+":"+str(int(BLOCK[block_name]["hap"][j]))+":"+str(1-int(BLOCK[block_name]["hap"][j])))
				allele1 = int(BLOCK[block_name]["hap"][j])
				allele2 = 1-int(BLOCK[block_name]["hap"][j])
				snp_decode[j]=str(allele1)+'|'+str(allele2)
			else:
				find = find_pos(int(j),position)
				if find <0:
					snp_decode[j]='0/1'
					continue
				id +=1
				tag = "ID=%i:LEFT_CLOUD_IDS=:RIGHT_CLOUD_IDS="  %(id)
				print 'interleave',block_name,id,j
				print >>f_out3,"%s\t%s\t%s\t%s\t%s" %(chr,str(j),str(j),tag,str(j)+":0:1")
				print >>f_out4,"%s\t%s\t%s\t%s\tinterleaving" %(chr,str(j),str(j),id)
				snp_decode[j]='0/1'
		if len(string)==0:
			id +=1
			print 'miss block',block_name
			print >>f_out4,"%s\t%s\t%s\t%s\t%s" %(chr,pos,pos2,id,mec[block_name][7]+",drop")
			continue
		if interleaving is True:
			id +=1
			tag = "ID=%i:LEFT_CLOUD_IDS=:RIGHT_CLOUD_IDS="  %(id)
			print >>f_out3,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,tag,"\t".join(string))
			print >>f_out4,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,id,mec[block_name][7]+",interleave")
			continue
		id +=1
		tag = "ID=%i:LEFT_CLOUD_IDS=:RIGHT_CLOUD_IDS="  %(id)
		print >>f_out1,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,tag,"\t".join(string))
		print >>f_out4,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,id,mec[block_name][7])	
	return snp_decode

def make_prism_interleave(SNP_not_shown):
#	f = open("Prism/%s/%s.%s.local.interleaving.blocks" %(SAMPLE,SAMPLE,CHROM),"r")
	f_out = open("Prism/%s/%s.%s.local.interleaving.v2.blocks" %(SAMPLE,SAMPLE,CHROM),"w")
	f_info = open("Prism/%s/%s.%s.info" %(SAMPLE,SAMPLE,CHROM),"r")
	index = 0
	string = []
	id = 0
	id_index = 0
	old_id = 0
	line = f_info.readline().rstrip()
	col = line.split('\t')
	pos1 = int(col[1])
	pos2 = int(col[2])
	while True:
		if index>=len(SNP_not_shown[CHROM]):
			break
		snp_pos=SNP_not_shown[CHROM][index]
		if int(snp_pos)<pos1:
			if old_id==id:
				id_index +=1
			else:
				id_index=0
			tag = "ID=%s:LEFT_CLOUD_IDS=:RIGHT_CLOUD_IDS="  %(str(id)+'_'+str(id_index))
			string=[str(snp_pos)+":0:1"]
			print >>f_out,"%s\t%s\t%s\t%s\t%s" %(CHROM,snp_pos,snp_pos,tag,"\t".join(string))
		else:
			if int(snp_pos)<pos2:
				id = int(col[3])
				if old_id==id:
					id_index +=1
				else:
					id_index=0
				tag = "ID=%s:LEFT_CLOUD_IDS=:RIGHT_CLOUD_IDS="  %(str(id)+'_'+str(id_index))
				string=[str(snp_pos)+":0:1"]
				print >>f_out,"%s\t%s\t%s\t%s\t%s" %(CHROM,snp_pos,snp_pos,tag,"\t".join(string))
			else:
				while True:
					line = f_info.readline().rstrip()
					if line=='':
						break
					col = line.split('\t')
					pos1 = int(col[1])
					pos2 = int(col[2])
					id = int(col[3])
					id_index = 0
					if int(snp_pos)<pos1:
						if old_id==id:
							id_index +=1
						else:
							id_index=0
						tag = "ID=%s:LEFT_CLOUD_IDS=:RIGHT_CLOUD_IDS="  %(str(id)+'_'+str(id_index))
						string=[str(snp_pos)+":0:1"]
						print >>f_out,"%s\t%s\t%s\t%s\t%s" %(CHROM,snp_pos,snp_pos,tag,"\t".join(string))
						break
		old_id = id
		index +=1

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='assign parental allele in block')
	parser.add_argument("--chr", dest='chr',help="chromosome")
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--source", dest='source',help="source:1KG,Affy or HapMap")
	parser.add_argument("--pop", dest='pop',default='MKK',help="population")
	parser.add_argument("--affy_dir",dest='affy_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/phase3_affy/',help='Affy directory')
	parser.add_argument("--hapmap_dir",dest='hapmap_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/',help='HapMap directory')
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab-scratch/shiya-projects/Reich_fosmid/wgs/',help="directory for whole genome sequencing file")
	parser.add_argument("--pool_dir", dest='pool_dir',default='/home/jmkidd/kidd-lab-scratch/shiya-projects/Reich_fosmid/pools/',help="directory for each pool data file")

	args = parser.parse_args()
	
	SAMPLE = args.sample
	CHROM = args.chr
	SOURCE = args.source
	POP = args.pop
	all_var_file= args.sample+"_combined_genome_allvars_gvcf_SNP"

	SNP = create_snp_list(all_var_file)
	print "SNP list set up"
	
	snp_not_shown_file = args.sample + "_gvcf_SNPs_not_shown"
	SNP_not_shown=SNPs_not_shown(snp_not_shown_file)

#	read refhap result and get information
	haplotype_file = args.sample+'_gvcf_snp_haplotype_detail'
	BLOCK,block_list,decode,mec=read_haplotype(haplotype_file,SNP_not_shown)
	print "haplotype info stored"
	print "Number of blocks",len(block_list)
	print block_list[0],block_list[1]
	decode_pos = decode_number_pos(all_var_file)
	matrix_file = args.sample + '_all_clone_all_fragment_gvcf_snp_sorted.matrix'
	MEC,mec,BLOCK=read_matrix(matrix_file,BLOCK,decode,decode_pos,mec)
	print BLOCK[block_list[0]]

	mec=make_matrix2png(BLOCK,block_list,mec)
	run_matrix2png()

	haplotype_file = 'Refhap/%s_%s.gvcf.snp.haplotype_detail' %(SAMPLE,CHROM)
	make_track(haplotype_file,block_list,mec)

	snp_decode = make_prism_input(BLOCK,mec)
	'''
	SAMPLE = args.sample
	CHROM = args.chr
	snp_not_shown_file = args.sample + "_gvcf_SNPs_not_shown"
	SNP_not_shown=SNPs_not_shown(snp_not_shown_file)
	make_prism_interleave(SNP_not_shown)
	'''	