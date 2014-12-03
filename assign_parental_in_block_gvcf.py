#!/usr/bin/env python
# python assign_parental_in_block.py
# Shiya Song
# 29 July 2014
# Assign parental allele within each blocks, calculate switch error, correct some switch errors, make matrix2png plots, assign clone into parental alleles
# Used for HG02799, HG03108, NA21302
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

def read_heter_snp_affy(SNP,locus,index):
	extra = 0
	different = 0
	same = 0
	het = 0
	for i in locus:
		chr = 'chr'+i[0]
		pos = i[2]
		info1 = i[3]
		info2 = i[4]
		if chr==CHROM:
			if index == 0:
				info = info1[0]+info1[2]
				geno = info2[0]+info2[2]
			elif index == 1:
				info = info1[0]+info1[1]
				geno = info2[0]+info2[1]
			elif index == 2:
				info = info1[2]+info1[3]
				geno = info2[2]+info2[3]
			if info[0]!=info[1]:
				het +=1
			else:
				continue
			if info[0]=="0" and info[1]=="1":
				ref_alt = geno[0]+geno[1]
			elif info[0]=="1" and info[1]=="0":
				ref_alt = geno[1]+geno[0]
			else:
				print "wrong:","\t".join(i),"\t"
			try:
				snp_info = SNP[chr]["geno"][pos]
			except KeyError:
				extra +=1
				continue
			if compare_allele(snp_info,ref_alt) is True:		# compare the child genotype with the ref/alt genotype
				SNP[chr]["phase"][pos]=info
				same +=1
			elif compare_allele_v2(snp_info,ref_alt) is True:	# same genotype but miss up the ref/alt, change it
				info = info.replace("0","B")
				info = info.replace("1","D")
				info = info.replace("B","1")
				info = info.replace("D","0")
				SNP[chr]["phase"][pos]=info
				different +=1
			else:
				print "\t".join(i),"\t",SNP[chr]["geno"][pos]
		else:
			continue
	print 'Affy',CHROM,het,extra,same,different
	return SNP,same+different

def read_heter_snp_hapmap(SNP,locus,index1,index2):		
	extra = 0
	ABAB = 0
	same = 0
	het = 0
	for i in locus:
		chr = i[0]
		pos = i[2]
		paternal = i[index1]+i[index1+1]
		maternal = i[index2]+i[index2+1]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			if paternal[0]!=maternal[0]:
				het +=1
			else:
				continue
			try:
				snp_info = SNP[chr]["geno"][pos]
			except KeyError:
				extra +=1
				continue
			if paternal[0]!=maternal[0] and paternal[1]==maternal[0] and paternal[0]==maternal[1]:
				ABAB +=1
				continue
			map = {SNP[chr]["geno"][pos][0]:'0',SNP[chr]["geno"][pos][1]:'1'}
			try:
				info1 = map[paternal[0]]+map[maternal[0]]
				SNP[chr]["phase"][pos]=info1
				same +=1
			except KeyError:
				continue
		else:
			break
	print 'HapMap',CHROM,het,extra,ABAB,same
	return SNP,ABAB,same

def type_heter_snp(string):
	if string.count("1")==3 or string.count("0")==3:
		return("AAAB")
	elif string[0]==string[1] and string[2]==string[3]:
		return("AABB")
	elif string[0]!=string[1] and string[2]!=string[3]:
		return("ABAB")

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
						mec[chr+" "+pos] = [0,0,0,0,0,0,0,0]		# number of snp, number of snp from Affy, number of support from Affy, MEC value before, MEC value after correction, switch_error, switch_error_corrected, source type
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
			chr = ("_").join(col[1].split("_")[4:-1])
			pool = ("_").join(col[1].split("_")[:4])
		pos = col[1].split("_")[-1]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr!=CHROM:
			break
		position = []
		new_position = []
		new_allele = []
		new_quality = []
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
#				print "not_found",position[i]
				continue
			assign = True
			real_position = decode_pos[position[i]]
			new_position.append(real_position)			# the actual position of each clone
#			new_allele.append(allele[i])
#			new_quality.append(quality[i])
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
		BLOCK[code][name]= pd.Series(allele,index=new_position)
		BLOCK[code][name+"_qual"]= pd.Series([m for m in quality],index=new_position)
		MEC += mismatch
		mec[code][3] += mismatch
	return MEC,mec,BLOCK

def calc_entropy(block,clone_name,pos):
	entropy = 0
	tot = 0
	for i in block.index:
		if math.isnan(float(block[clone_name][i])) is False:
			tot +=1
			if i>=pos:
				entropy+=1
	if entropy>tot-entropy:
		entropy=tot-entropy
	return entropy

def find_break(block,switch_position,start,end):
	list = ["hap","transmit","untransmit","switch"]
	Index = block.index
	k = [m for m in Index if int(m)>int(switch_position[start]) and int(m)<int(switch_position[end])]
	if len(k)==0:
		entropy = 0
		for j in block:
			if j not in list and j[-5:]!="_qual":
				entropy+= calc_entropy(block,j,switch_position[end])
		return switch_position[end],entropy
	else:
		entropy = pd.Series(0,index=k)
	for i in range(len(Index)):
		if int(Index[i])>int(switch_position[start]) and int(Index[i])<int(switch_position[end]):
			for j in block:
				if j not in list and j[-5:]!="_qual":
					entropy[Index[i]]+= calc_entropy(block,j,Index[i])
	print "entropy",entropy
	min=entropy.argmin()
	print "minimum entropy at",entropy.index[min],entropy[min]
	return entropy.index[min],entropy[min]
								
def correct_switch_error(block,start,end,switch_position,resolve,corrected_error):
	list = ["hap","transmit","untransmit","switch"]
	if start!=0:
		start_pos,start_entropy=find_break(block,switch_position,start-1,start)
	else:
		start_pos=0
		start_entropy = 0
	if end!=len(switch_position)-1:
		end_pos,end_entropy=find_break(block,switch_position,end,end+1)
	else:
		end_pos=5000000000
		end_entropy =0
	if start_entropy+end_entropy>=50:
		print 'Too much entropy, do not switch',start_entropy,end_entropy
		return block,corrected_error
	else:
		print 'entropy:',start_entropy,end_entropy
		print "Correct from start",start_pos,"End",end_pos
		for i in block.index:
			if int(i)>=int(start_pos) and int(i)<int(end_pos) and i not in resolve:
				block["hap"][i] = 1 - block["hap"][i]
		return block,corrected_error+1

def calc_switch_error(SNP,BLOCK,block_list,mec):
	f_hap = open("BLOCK/%s_gvcf/%s_" %(SAMPLE,SAMPLE)+CHROM,"w")
	Switch_error = 0
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		assert chr==CHROM
		pos = block_list[i].split(" ")[1]
		match = 0
		switch_match = 0
		is_switch=[]
		switch = []
		print >>f_hap,"BLOCK: %s %s" %(chr,pos)
		list = BLOCK[block_list[i]].index
		resolve = []
		switch_position=[]
		hapmap_1KG_position=[]
		hapmap_phase = 0
		for j in range(len(list)):
			if SNP[chr]["phase"][list[j]]!="":
				hapmap_phase+=1
		if hapmap_phase > 0:
			phase_tag = "phase"
			mec[block_list[i]][7]=SOURCE
		else:
			phase_tag = "phase"
			mec[block_list[i]][7]="None"
		for j in range(len(list)):
			if SNP[chr][phase_tag][list[j]]!="":
				if BLOCK[block_list[i]]["hap"][list[j]]==int(SNP[chr][phase_tag][list[j]][0]):
					match += 1
					is_switch.append(0)
					switch.append(0)
					switch_position.append(list[j])
					hapmap_1KG_position.append(list[j])
					print >>f_hap,"%s\t%s\tnon_switch\t%s\t%i" %(chr,list[j],SNP[chr][phase_tag][list[j]],BLOCK[block_list[i]]["hap"][list[j]])
				elif BLOCK[block_list[i]]["hap"][list[j]]==int(SNP[chr][phase_tag][list[j]][1]):
					switch_match += 1
					is_switch.append(1)
					switch.append(0)
					switch_position.append(list[j])
					hapmap_1KG_position.append(list[j])
					print >>f_hap,"%s\t%s\tswitch\t%s\t%i" %(chr,list[j],SNP[chr][phase_tag][list[j]],BLOCK[block_list[i]]["hap"][list[j]])
				elif BLOCK[block_list[i]]["hap"][list[j]]==0.5:
					BLOCK[block_list[i]]["hap"][list[j]]=int(SNP[chr][phase_tag][list[j]][0])
					resolve.append(list[j])
					hapmap_1KG_position.append(list[j])
					print >>f_hap,"%s\t%s\tresolve\t%s\t%i" %(chr,list[j],SNP[chr][phase_tag][list[j]],BLOCK[block_list[i]]["hap"][list[j]])
			else:
				print >>f_hap,"%s\t%s\tnone\t\t\t" %(chr,list[j])
		BLOCK[block_list[i]]["transmit"]=pd.Series([SNP[chr][phase_tag][n][0] for n in hapmap_1KG_position],index=hapmap_1KG_position)
		BLOCK[block_list[i]]["untransmit"]=pd.Series([SNP[chr][phase_tag][n][1] for n in hapmap_1KG_position],index=hapmap_1KG_position)
		BLOCK[block_list[i]]["switch"]=pd.Series(is_switch,index=switch_position)
		mec[block_list[i]][0]=len(list)
		mec[block_list[i]][1]=match+switch_match
		mec[block_list[i]][2]=match
		if match == 0 and switch_match==0:
			print "non-phased: block %s %s" %(chr,pos)
		elif match >= switch_match:
			print "phased by %s: block %s %s non-switch %i %i %i" %(mec[block_list[i]][7],chr,pos,mec[block_list[i]][0],mec[block_list[i]][1],mec[block_list[i]][2])
			if switch_match != 0 :
				switch_error = 0
				corrected_error = 0
				print is_switch
				last = ""
				k = 0
				s = 0
				e = 0
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
							section_pos = switch_position[s:e]
							section_switch = is_switch[s:e]
							print "start",s,"end",e
							print "position",section_pos
							print "switch info",section_switch
							BLOCK[block_list[i]],corrected_error = correct_switch_error(BLOCK[block_list[i]],s,e-1,switch_position,resolve,corrected_error)
						s = k
						last = is_switch[k]
				e = k
				if is_switch[s]==1:
					switch_error +=1
					section_pos = switch_position[s:e+1]
					section_switch = is_switch[s:e+1]
					print "start",s,"end",e+1
					print "position",section_pos
					print "switch info",section_switch
					BLOCK[block_list[i]],corrected_error = correct_switch_error(BLOCK[block_list[i]],s,e,switch_position,resolve,corrected_error)
				print "switch final:","\t".join(map(str,[BLOCK[block_list[i]]["switch"][n] for n in switch_position]))
				print "switch_error",switch_error
				print "corrected_error",corrected_error
				mec[block_list[i]][5]=switch_error
				mec[block_list[i]][6]=corrected_error
				Switch_error += switch_error
		else:
			print "phased by %s: block %s %s switch %i %i %i" %(mec[block_list[i]][7],chr,pos,mec[block_list[i]][0],mec[block_list[i]][1],mec[block_list[i]][2])
			if match != 0 :
				switch_error = 0
				corrected_error = 0
				print is_switch
				last = ""
				k = 0
				s = 0
				e = 0
				section = []
				last_shift = False
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
							section_pos = switch_position[s:e]
							section_switch = is_switch[s:e]
							print "start",s,"end",e
							print "position",section_pos
							print "switch info",section_switch
							BLOCK[block_list[i]],corrected_error = correct_switch_error(BLOCK[block_list[i]],s,e-1,switch_position,resolve,corrected_error)
						s = k
						last = is_switch[k]
				e = k
				if is_switch[s]==0:
					switch_error +=1
					section_pos = switch_position[s:e+1]
					section_switch = is_switch[s:e+1]
					print "start",s,"end",e+1
					print "position",section_pos
					print "switch info",section_switch
					BLOCK[block_list[i]],corrected_error = correct_switch_error(BLOCK[block_list[i]],s,e,switch_position,resolve,corrected_error)
				print "switch_error",switch_error
				print "corrected_error",corrected_error
				mec[block_list[i]][5]=switch_error
				mec[block_list[i]][6]=corrected_error
				Switch_error += switch_error
			for k in BLOCK[block_list[i]].index:
				BLOCK[block_list[i]]["hap"][k] = 1- BLOCK[block_list[i]]["hap"][k]
	return BLOCK,mec,Switch_error

def make_matrix2png(BLOCK,block_list):
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
			if math.isnan(float(BLOCK[block_list[i]]["transmit"][j])) is True:
				string1.append("-")
				string2.append("-")
			else:
				string1.append(BLOCK[block_list[i]]["transmit"][j])
				string2.append(BLOCK[block_list[i]]["untransmit"][j])
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

def run_matrix2png():
	wk = "BLOCK/%s_gvcf_matrix2png/" %(SAMPLE)
	for i in glob.glob(wk + "*.matrix"):
		j = i.replace(wk,"")
		j = j.replace(".matrix",".png")
		cmd = "matrix2png -data %s -size 10:10 -map 18 -r -c -s -d -range 0:1 > %s" %(i,wk+j)
		os.popen(cmd)

def determine_clone_hap(block,clone,decode_pos):
	match = 0
	mismatch = 0
	unknown = 0
	hap1 = False
	hap = ""
	qual1 = 0
	qual2 = 0
	for i in block.index:
		if math.isnan(float(block[clone][i])) is False:
			if int(block[clone][i])==block["hap"][i]:
				match +=1
				qual1 += ord(block[clone+"_qual"][i])-33
			elif int(block[clone][i])==1-block["hap"][i]:
				mismatch +=1
				qual2 += ord(block[clone+"_qual"][i])-33
			elif block["hap"][i]==0.5:
				unknown +=1
	if match > mismatch:
		hap1 = True
		hap = "hap1"
		mec = mismatch
	elif match == mismatch:
		if qual1>=qual2:
			hap = "hap1"
			mec = mismatch
		else:
			hap = "hap2"
			mec = match
	else:
		hap = "hap2"
		mec = match
	return hap,mec,match+mismatch

def separate_clone_into_hap(BLOCK,block_list,mec,decode_pos):
	MEC_fraction = []
	f_hap1=open("BLOCK/%s_gvcf/%s_hap1_" %(SAMPLE,SAMPLE)+CHROM+"_clone.bed","w")
	f_hap2=open("BLOCK/%s_gvcf/%s_hap2_" %(SAMPLE,SAMPLE)+CHROM+"_clone.bed","w")
	f_unknown_hap1=open("BLOCK/%s_gvcf/%s_unknown_hap1_" %(SAMPLE,SAMPLE)+CHROM+"_clone.bed","w")
	f_unknown_hap2=open("BLOCK/%s_gvcf/%s_unknown_hap2_" %(SAMPLE,SAMPLE)+CHROM+"_clone.bed","w")
	for i in range(len(block_list)):
		for clone in BLOCK[block_list[i]]:
			if clone not in ["hap","transmit","untransmit","switch"] and clone[-5:]!="_qual":
				hap,clone_mec,clone_length=determine_clone_hap(BLOCK[block_list[i]],clone,decode_pos)
				mec[block_list[i]][4]+=clone_mec
				'''
				if SAMPLE=='NA19240':
					chr = ("_").join(clone.split("_")[3:-2])
					pool = ("_").join(clone.split("_")[:3])
				else:
					chr = ("_").join(clone.split("_")[4:-2])
					pool = ("_").join(clone.split("_")[:4])
				pos = clone.split("_")[-2]
				find = False
				clone_name = args.pool_dir+'%s/' %(SAMPLE) +pool + '/' + pool + '.markdup.clone.sel.10000.0.25_SNP_snp_all' 
				f_clone = open(clone_name,"r")
				for l in f_clone:
					l = l.strip()
					chrom = l.split("\t")[0]
					pos1 = l.split("\t")[1]
					pos2 = l.split("\t")[2]
					if chrom == chr and int(pos) == int(pos1):
						find = True
						break
				f_clone.close()
				if os.path.isdir(args.pool_dir+'%s/' %(SAMPLE) +pool + '/hap_bam') is False:
					os.popen('mkdir '+args.pool_dir+'%s/' %(SAMPLE) +pool + '/hap_bam')
				if mec[block_list[i]][1]!=0:
					file = args.pool_dir+'%s/' %(SAMPLE) +pool + '/hap_bam/' + pool + "." + hap+"."+CHROM + ".gvcf.bed"
					f_out = open(file,"a")
					if hap == "hap1":
						print >>f_out, "%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t255,0,0" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
						print >>f_hap1,"%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t255,0,0" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
					else:
						assert hap == "hap2"
						print >>f_out, "%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t0,0,255" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
						print >>f_hap2,"%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t0,0,255" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
				else:
					file = args.pool_dir+'%s/' %(SAMPLE) +pool + '/hap_bam/' + pool + ".unknown_" + hap+"."+CHROM + ".gvcf.bed"
					f_out = open(file,"a")
					if hap == "hap1":
						print >>f_out, "%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t255,0,0" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
						print >>f_unknown_hap1,"%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t255,0,0" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
					else:
						assert hap == "hap2"
						print >>f_out, "%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t0,0,255" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
						print >>f_unknown_hap2,"%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t0,0,255" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
				f_out.close()
				'''
	return mec
					
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
			print >>f_out,"%s\t%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i" %(chr,start,end,MEC[7],MEC[0],MEC[1],MEC[2],MEC[3],MEC[4],MEC[5],MEC[6])
		elif start == 0:
			chr = line.split("\t")[3]
			if chr!=CHROM:
				break
			start = line.split("\t")[4]
			find = find_pos(int(start),SNP_not_shown[chr])
			if find <0:  # SNP shows in the matrix
				end = line.split("\t")[4]
				MEC = mec[chr + " " + start]
			else:
				start = 0
		elif start != "":
			end = line.split("\t")[4]		

def detect_hapmap_wrong(BLOCK,block_list,mec):
	f2 = open("BLOCK/%s_gvcf/%s_%s_wrong_detail_" %(SAMPLE,SAMPLE,SOURCE)+CHROM,"w")
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		pos = block_list[i].split(" ")[1]
		if mec[block_list[i]][7]=="None":
			continue
		for j in BLOCK[block_list[i]]["hap"].index:
			if math.isnan(float(BLOCK[block_list[i]]["transmit"][j])) is False:
				if BLOCK[block_list[i]]["hap"][j]!=int(BLOCK[block_list[i]]["transmit"][j]):
					print >>f2,"%s\t%i\t%s\t%i\t%s\t" %(chr,int(pos),j,BLOCK[block_list[i]]["hap"][j],BLOCK[block_list[i]]["transmit"][j])

def read_hap(BLOCK,block_list,mec):
	snp_decode={}
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		pos = block_list[i].split(" ")[1]
		if mec[block_list[i]][7]=='None':
			for j in BLOCK[block_list[i]].index:
				snp_decode[j]='0/1'
			continue
		for j in BLOCK[block_list[i]].index:
			if BLOCK[block_list[i]]["hap"][j] !=0.5:
				allele1 = int(BLOCK[block_list[i]]["hap"][j])
				allele2 = 1-int(BLOCK[block_list[i]]["hap"][j])
				snp_decode[j]=str(allele1)+'|'+str(allele2)
			else:
				allele1 = int(BLOCK[block_list[i]]["hap"][j])
				allele2 = 1-int(BLOCK[block_list[i]]["hap"][j])
				snp_decode[j]='0/1'
	return snp_decode

def write_vcf(wgs_vcf,snp_decode):
	f = open(wgs_vcf,'r')
	f_out=gzip.open("BLOCK/%s_gvcf/%s_%s_phased_snps.vcf.gz" %(SAMPLE,SAMPLE,CHROM),'w')
	for line in f:
		line = line.strip()
		if line[0] == "#":
			f_out.write(line+'\n')
			continue
		line1 = line.split('\t')
		chr = line1[0]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			ref = line1[3]
			alt = line1[4].split(',')
			pos = line1[1]
			info = line1[9]
			info = info.split(':')
			if info[0] == "1/2":
				if len(ref)!=1 or len(alt[0])!=1 or len(alt[1])!=1:
					continue
				else:
					try:
						info[0]=str(int(snp_decode[pos][0])+1)+snp_decode[pos][1]+str(int(snp_decode[pos][2])+1)
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
					except KeyError:
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
			elif info[0] == "0/1":
				if len(ref)!=1 or len(alt[0])!=1:
					continue
				else:
					try: 
						info[0]=snp_decode[pos]
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
					except KeyError:
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
		else:
			break

def write_vcf_v2(wgs_vcf,snp_decode):
	f = gzip.open(wgs_vcf,'r')
	f_out=gzip.open("%s%s/gVCF_calls/%s.%s.fosmid.v2.phased.vcf.gz" %(args.wgs_dir,SAMPLE,SAMPLE,CHROM),'w')
	for line in f:
		line = line.strip()
		if line[0] == "#":
			f_out.write(line+'\n')
			continue
		line1 = line.split('\t')
		chr = line1[0]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			ref = line1[3]
			alt = line1[4].split(',')
			pos = line1[1]
			info = line1[9]
			info = info.split(':')
			if info[0] == "1/2":
				if len(ref)!=1 or len(alt[0])!=1 or len(alt[1])!=1:
					continue
				else:
					try:
						info[0]=str(int(snp_decode[pos][0])+1)+snp_decode[pos][1]+str(int(snp_decode[pos][2])+1)
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
					except KeyError:
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
			elif info[0] == "0/1":
				if len(ref)!=1 or len(alt[0])!=1:
					continue
				else:
					try: 
						info[0]=snp_decode[pos]
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
					except KeyError:
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
			elif info[0]=="1/1":
				if len(ref)!=1 or len(alt[0])!=1:
					continue
				info[0] = "1|1"
				print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
		else:
			break

def read_ped(file):
	f=open(file,'r')
	line = f.readline().strip().split()
	family=(line[0],line[1],line[2],line[3])
	index = family.index(SAMPLE)
	return family,family[2],family[3],index-1

def make_prism_input(BLOCK,mec):
	f = open("Prism/%s_new/%s.%s.positions" %(SAMPLE,SAMPLE,CHROM),'r')
	f_track = "BLOCK/%s_gvcf/%s_refhap_track_info_" %(SAMPLE,SAMPLE)+CHROM+".txt"
	position = []
	for line in f:
		position.append(int(line.strip().split(',')[0]))
	snp_decode={}
	f_out1 = open("Prism/%s_new/%s.%s.local.blocks.all" %(SAMPLE,SAMPLE,CHROM),"w")
	f_out2 = open("Prism/%s_new/%s.%s.local.blocks" %(SAMPLE,SAMPLE,CHROM),"w")
	f_out3 = open("Prism/%s_new/%s.%s.local.interleaving.blocks" %(SAMPLE,SAMPLE,CHROM),"w")
	f_out4 = open("Prism/%s_new/%s.%s.info" %(SAMPLE,SAMPLE,CHROM),"w")
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
			else:
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
			print 'miss block',block_name
			print >>f_out4,"%s\t%s\t%s\t%s\t%s" %(chr,pos,pos2,id,mec[block_name][7]+",drop")
			continue
		if interleaving is True:
			id +=1
			tag = "ID=%i:LEFT_CLOUD_IDS=:RIGHT_CLOUD_IDS="  %(id)
			print >>f_out3,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,tag,"\t".join(string))
			print >>f_out4,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,id,mec[block_name][7]+",interleave")
			continue
		if mec[block_name][7]=="None":
			if last_item!="":
				print >>f_out2,last_item
				last_item = ""
			id +=1
			tag = "ID=%i:LEFT_CLOUD_IDS=:RIGHT_CLOUD_IDS="  %(id)
			print >>f_out1,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,tag,"\t".join(string))
			print >>f_out2,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,tag,"\t".join(string))
			prev_none = True
		else:
			id +=1
			tag = "ID=%i:LEFT_CLOUD_IDS=:RIGHT_CLOUD_IDS="  %(id)
			print >>f_out1,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,tag,"\t".join(string))
			if prev_none is True:
				print >>f_out2,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,tag,"\t".join(string))
				prev_none = False
			else:
				last_item = "%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,tag,"\t".join(string))
		print >>f_out4,"%s\t%s\t%s\t%s\t%s" %(chr,pos,last_pos,id,mec[block_name][7])	
	return snp_decode

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='assign parental allele in block')
	parser.add_argument("--chr", dest='chr',help="chromosome")
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--source", dest='source',help="source:1KG,Affy or HapMap")
	parser.add_argument("--pop", dest='pop',default='MKK',help="population")
	parser.add_argument("--affy_dir",dest='affy_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/phase3_affy/',help='Affy directory')
	parser.add_argument("--hapmap_dir",dest='hapmap_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/',help='HapMap directory')
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	parser.add_argument("--pool_dir", dest='pool_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools/',help="directory for each pool data file")

	args = parser.parse_args()
	
	SAMPLE = args.sample
	CHROM = args.chr
	SOURCE = args.source
	POP = args.pop
	all_var_file= args.sample+"_combined_genome_allvars_gvcf_SNP"
	
	'''
	SNP = create_snp_list(all_var_file)
	print "SNP list set up"
	
	ped_file = args.sample + '.ped'
	family,dad,mom,index = read_ped(ped_file)
	
	if SOURCE == 'Affy':
#	read the locus from Affy	
		Affy_file = args.affy_dir + family[0]+'_Affy_snps.bed'
		locus_affy = read_tab_file(Affy_file)
		SNP,others =read_heter_snp_affy(SNP,locus_affy,index)
		print "SNP info stored"
		print "Number of SNPs used in %s:" %(SOURCE),others
	
	if SOURCE == 'HapMap':
#	read the locus from HapMap	
		Affy_file = args.hapmap_dir + '%s/%s_phased_snp.txt.gz' %(POP,POP)
		header = gzip.open(args.hapmap_dir + '%s/hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_%s.phased.gz' %(POP,POP.lower()),'r').readline().strip().split()
		index_dad = header.index(dad+'_A')+1
		index_mom = header.index(mom+'_A')+1
		print dad,mom,index_dad,index_mom
		SNP,ABAB,same =read_heter_snp_hapmap(SNP,BedIterator(Affy_file),index_dad,index_mom)
		print "SNP info stored"
		print "Number of SNPs used in %s:" %(SOURCE),same
	

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
	

	BLOCK,mec,switch_error=calc_switch_error(SNP,BLOCK,block_list,mec)
	print "Switch_error",switch_error
	detect_hapmap_wrong(BLOCK,block_list,mec)
	
	make_matrix2png(BLOCK,block_list)
	run_matrix2png()

	dbfile=open('BLOCK/%s_gvcf/%s_BLOCK_' %(SAMPLE,SAMPLE) +CHROM+'_pickle','wb')
	pickle.dump(BLOCK,dbfile) 
	pickle.dump(block_list,dbfile)
	pickle.dump(mec,dbfile)
	mec = separate_clone_into_hap(BLOCK,block_list,mec,decode_pos)
#	mec = compare_with_1000G(BLOCK,block_list,SNP,mec)
	
	haplotype_file = 'Refhap/%s_%s.gvcf.snp.haplotype_detail' %(SAMPLE,CHROM)
	make_track(haplotype_file,block_list,mec)
	'''
	dbfile=open('BLOCK/%s_gvcf/%s_BLOCK_' %(SAMPLE,SAMPLE) +CHROM+'_pickle','rb')
	BLOCK=pickle.load(dbfile) 
	block_list=pickle.load(dbfile)
	mec=pickle.load(dbfile)

	snp_decode = read_hap(BLOCK,block_list,mec)
#	snp_decode = make_prism_input(BLOCK,mec)
	wgs_vcf = args.wgs_dir + SAMPLE+ '/gVCF_calls/%s.%s.vcf.gz' %(SAMPLE,CHROM)
	write_vcf_v2(wgs_vcf,snp_decode)
	
#	print "Summary info as below:"
#	print "Total number of blocks in ",chr,":",len(block_list)
#	print "Total Switch_error",switch_error

	
