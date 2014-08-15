#!/usr/bin/env python
# python assign_parental_in_block.py
# Shiya Song
# 29 July 2014
# 

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
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
			break
	print het,extra,same,different
	return SNP,same+different
	
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
	for i in range(1,2405814):
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
						mec[chr+" "+pos] = [0,0,0,0,0,0,0,0,0]		# number of snp, number of snp from HARSH, number of support from HARSH, error value from clone support, error value from clone support after correction, switch_error
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

def calc_switch_interval(block,pos1,pos2,resolve):
	interval = []
	find = False
	for i in block.index:
		if int(i)>=int(pos1) and int(i)<=int(pos2):
			interval.append(i)
	if len(interval)>=40:
		j = len(interval)/2
		print pos1,pos2,len(interval),interval[j]
		while j<len(interval):
			if interval[j] not in resolve:
				find = True
				break
			else:
				j +=1
		if find == True:	
			return interval[j]	
		else:
			return False			
	else:
		return False

def check_nearby(block,clone_name,index):
	left,right = [0,0]
	for i in range(len(block.index)):
		if i<index and math.isnan(float(block[clone_name][block.index[i]])) is False:
			left +=1
		elif i>index and math.isnan(float(block[clone_name][block.index[i]])) is False:
			right +=1
	print clone_name,block.index[index],left,right
	if (left>=3 and right >=3) or (left>=2 and right>=10) or (left>=10 and right>=2):
		return 1
	else:
		return 0
	
def check_break(block,pos1,pos2):
	list = ["hap","transmit","untransmit","switch"]
	Index = block.index
	unresolve = 0
	tot = 0
	Q = True
	for i in range(len(block.index)):
		if int(Index[i])>=int(pos1) and int(Index[i])<=int(pos2):
			qualify = 0
			for j in block:
				if j not in list and j[-5:]!="_qual" :
					if int(find_first_pos(block,j))<int(Index[i]) and int(find_last_pos(block,j))>int(Index[i]):
						qualify += check_nearby(block,j,i)
			print block.index[i],qualify
			if qualify <1:
				Q = False
				return Q
	for i in range(len(block.index)):
		if int(Index[i])>int(pos1) and int(Index[i])<int(pos2):
			tot +=1
			if block["hap"][Index[i]]==0.5:
				unresolve +=1
	if unresolve >=3:
		Q = False
	return Q

def check_unresolve(block,pos1,pos2):
	unresolve = 0
	Index = block.index
	for i in range(len(block.index)):
		if int(Index[i])>int(pos1) and int(Index[i])<int(pos2):
			if block["hap"][Index[i]]==0.5:
				unresolve +=1
	if unresolve >=1:
		Q = False
	else:
		Q = True
	return Q

def check_connection(block,is_switch,switch_island,resolve):							# This step check the connection between those positions with known phasing information
	list = ["hap","transmit","untransmit","switch"]
	is_connect = []
	correct_pos = []
	for index in range(len(switch_island)-1):
		pos1 = switch_island[index]
		pos2 = switch_island[index+1]
		total = 0
		consistent = 0
		for i in block:
			if i not in list and i[-5:]!="_qual":
				if math.isnan(float(block[i][pos1])) is False and math.isnan(float(block[i][pos2])) is False:
					total +=1
					if int(block[i][pos2])==block["hap"][pos2] and int(block[i][pos1])==block["hap"][pos1]:
						consistent +=1
					elif int(block[i][pos2])==(1-block["hap"][pos2]) and int(block[i][pos1])==(1-block["hap"][pos1]):
						consistent +=1
#		print "check switch_island",pos1,pos2,total,consistent
		if index == 0 :
			is_connect.append(1)
		if total < 1 or consistent<1:
			inter_pos = calc_switch_interval(block,pos1,pos2,resolve)
			if inter_pos is False:
				if check_break(block,pos1,pos2) is False:
					a = is_connect[len(is_connect)-1]
					is_connect.append(1-a)
				else:
					a = is_connect[len(is_connect)-1]
					is_connect.append(a)
			else:
				print "Use intermediate position to decide"
				total_1 = 0
				total_2 = 0
				consistent_1 = 0
				consistent_2 = 0
				for i in block:
					if i not in list and i[-5:]!="_qual":
						if math.isnan(float(block[i][pos1])) is False and math.isnan(float(block[i][inter_pos])) is False:
							total_1 +=1
							if int(block[i][inter_pos])==block["hap"][inter_pos] and int(block[i][pos1])==block["hap"][pos1]:
								consistent_1 +=1
							elif int(block[i][inter_pos])==(1-block["hap"][inter_pos]) and int(block[i][pos1])==(1-block["hap"][pos1]):
								consistent_1 +=1
						if math.isnan(float(block[i][pos2])) is False and math.isnan(float(block[i][inter_pos])) is False:
							total_2 +=1
							if int(block[i][pos2])==block["hap"][pos2] and int(block[i][inter_pos])==block["hap"][inter_pos]:
								consistent_2 +=1
							elif int(block[i][pos2])==(1-block["hap"][pos2]) and int(block[i][inter_pos])==(1-block["hap"][inter_pos]):
								consistent_2 +=1
#				print "check switch_island",pos1,inter_pos,total_1,consistent_1
#				print "check switch_island",inter_pos,pos2,total_2,consistent_2
				if total_1 >=1 and total_2 >=1 and consistent_1>=1 and consistent_2>=1:
					a = is_connect[len(is_connect)-1]
					is_connect.append(a)
					print "intermediate position helped",pos1,pos2
				else:
					if check_break(block,pos1,pos2) is False:
						a = is_connect[len(is_connect)-1]
						is_connect.append(1-a)
					else:
						a = is_connect[len(is_connect)-1]
						is_connect.append(a)
		else:
			a = is_connect[len(is_connect)-1]
			is_connect.append(a)	
	assert len(is_connect)==len(switch_island)
	print is_connect
	print "\t".join(switch_island)
#	block["connect"]=pd.Series(is_connect,index=switch_island)
	return is_connect

def switch_whole(block,switch_position,switch_error,resolve,list):
	section = []
	first = True
	for i in range(len(switch_position)):
		pos = switch_position[i]
		if block["switch"][pos]==1:
			section.append(i)
		else:
			if len(section)==0:
				continue
			else:
				if len(section)==1:
					snp = 0
					if section[0]>=1 and section[0]<len(switch_position)-1:
						for j in list:
							if int(j)>int(switch_position[section[0]-1]) and int(j)<int(switch_position[section[0]+1]):
								snp +=1
					elif section[0]==0:
						for j in list:
							if int(j)>=int(list[0]) and int(j)<int(switch_position[section[0]+1]):
								snp +=1
					elif section[0]==len(switch_position)-1:
						for j in list:
							if int(j)>int(switch_position[section[0]-1]) and int(j)<=int(list[-1]):
								snp +=1
					if snp < 100:
						print "position needs to be correct",section
						block=correct_switch_error(block,section[0],section[-1],switch_position,resolve)
						switch_error += 1
						section = []
					else:
						print 'Only single snp',switch_position[section[0]],'over such long region ',snp
						section = []
				else:
					print "position needs to be correct",section
					block=correct_switch_error(block,section[0],section[-1],switch_position,resolve)
					switch_error += 1
					section = []
	if len(section)!=0:
		if len(section)==1:
			snp = 0
			if section[0]>=1 and section[0]<len(switch_position)-1:
				for j in list:
					if int(j)>int(switch_position[section[0]-1]) and int(j)<int(switch_position[section[0]+1]):
						snp +=1
			elif section[0]==0:
				for j in list:
					if int(j)>=int(list[0]) and int(j)<int(switch_position[section[0]+1]):
						snp +=1
			elif section[0]==len(switch_position)-1:
				for j in list:
					if int(j)>int(switch_position[section[0]-1]) and int(j)<=int(list[-1]):
						snp +=1
			if snp < 100:
				print "position needs to be correct",section
				block=correct_switch_error(block,section[0],section[-1],switch_position,resolve)
				switch_error += 1
				section = []
			else:
				print 'Only single snp',switch_position[section[0]],'over such long region ',snp
				section = []
		else:
			print "position needs to be correct",section
			block=correct_switch_error(block,section[0],section[-1],switch_position,resolve)
			switch_error += 1
	return block,switch_error
		
def find_first_pos(block,clone_name):
	for i in block.index:
		if math.isnan(float(block[clone_name][i])) is False:
			return i

def find_last_pos(block,clone_name):
	m=range(len(block.index))
	m.reverse()
	for i in m:
		if math.isnan(float(block[clone_name][block.index[i]])) is False:
			return block.index[i]

def find_break_front(block,switch_position,start,end):
	list = ["hap","transmit","untransmit","switch"]
	Index = block.index
	k = [m for m in Index if int(m)>int(switch_position[start]) and int(m)<int(switch_position[end])]
	if len(k)==0:
		return switch_position[end]
	else:
		occurrence = pd.Series(0,index=k)
	min_position = k[0]
	min_occurrence = 1000000
	for i in range(len(k)):
		for j in block:
			if j not in list and j[-5:]!="_qual":
				if math.isnan(float(block[j][k[i]])) is False and find_first_pos(block,j)!=k[i]:
					occurrence[k[i]]+=1
		if occurrence[k[i]]<=min_occurrence:
			min_occurrence=occurrence[k[i]]
			min_position = k[i]
	print "occurrence",occurrence
	return k[i]

def find_break_back(block,switch_position,start,end):
	list = ["hap","transmit","untransmit","switch"]
	Index = block.index
	k = [m for m in Index if int(m)>int(switch_position[start]) and int(m)<int(switch_position[end])]
	if len(k)==0:
		return switch_position[end]
	else:
		occurrence = pd.Series(0,index=k)
	min_position = k[0]
	min_occurrence = 1000000
	for i in range(len(k)):
		for j in block:
			if j not in list and j[-5:]!="_qual":
				if math.isnan(float(block[j][k[i]])) is False and find_first_pos(block,j)!=k[i]:
					occurrence[k[i]]+=1
		if occurrence[k[i]]<min_occurrence:
			min_occurrence=occurrence[k[i]]
			min_position = k[i]
	print "occurrence",occurrence
	return k[i]

def find_break(block,switch_position,start,end):
	list = ["hap","transmit","untransmit","switch"]
	Index = block.index
	k = [m for m in Index if int(m)>int(switch_position[start]) and int(m)<int(switch_position[end])]
	if len(k)==0:
		return switch_position[end]
	else:
		occurrence = pd.Series(0,index=k)
	for i in range(len(Index)):
		if int(Index[i])>int(switch_position[start]) and int(Index[i])<int(switch_position[end]):
			for j in block:
				if j not in list and j[-5:]!="_qual":
					if math.isnan(float(block[j][Index[i]])) is False and find_first_pos(block,j)!=Index[i]:
						occurrence[Index[i]]+=1
	print "occurrence",occurrence
	min=occurrence.argmin()
	return occurrence.index[min]
								
def correct_switch_error(block,start,end,switch_position,resolve):
	list = ["hap","transmit","untransmit","switch"]
	if start!=0:
		start_pos=find_break(block,switch_position,start-1,start)
	else:
		start_pos=0
	if end!=len(switch_position)-1:
		end_pos=find_break(block,switch_position,end,end+1)
	else:
		end_pos=5000000000					
	print "Correct from start",start_pos,"End",end_pos
	for i in block.index:
		if int(i)>=int(start_pos) and int(i)<int(end_pos) and i not in resolve:
			block["hap"][i] = 1 - block["hap"][i]
	return block

def check_phase_error(block,section,Phase_error,zero,one):
	list = ["hap","transmit","untransmit","switch"]
	if zero-one>=2:
		print "more zeros"
		for j in section.index:
			if section[j]==1:
				a,b,c,d=[0,0,0,0]
				for i in block:
					if i not in list and i[-5:]!="_qual":
						if i[-1] == "1":
							hap1 = True
						else:
							i[-1]=="2"
							hap1 = False
						if math.isnan(float(block[i][j])) is False :
							if hap1 is True and int(block[i][j])==block["hap"][j]:
								a +=1
							elif hap1 is True and int(block[i][j])==(1-block["hap"][j]):
								b +=1
							elif hap1 is False and int(block[i][j])==(1-block["hap"][j]):
								d +=1
							elif hap1 is False and int(block[i][j])==block["hap"][j]:
								c +=1
							else:
								print "other possibility",j
				cmds ="TwoSidedFastFishersExactTest %i %i %i %i" %(a,b,c,d)
				print cmds
				if (d==0 and a>0) or (a==0 and d>0) or (a==1 and b==0 and c==0 and d==1):
					if a+b+c+d<=1:
						section[j]=1-section[j]
						block["hap"][j] = 1-block["hap"][j]				# possible phase error
						Phase_error +=1
						print "possible phase error",j,a,b,c,d
					else:
						continue
				else:
					m=subprocess.Popen(["TwoSidedFastFishersExactTest",str(a),str(b),str(c),str(d)],stdout=subprocess.PIPE)
					p=m.communicate()[0].strip()
					m.stdout.close()
					if float(p)>=0.5:
						section[j]=1-section[j]
						block["hap"][j] = 1-block["hap"][j]				# possible phase error
						Phase_error +=1
						print "possible phase error",j,a,b,c,d,p
		is_switch = [section[m] for m in section.index]
		print "switch info within island after correcting phasing error",is_switch
	elif one-zero >=2:
		print "more ones"
		for j in section.index:
			if section[j]==0:
				a,b,c,d=[0,0,0,0]
				for i in block:
					if i not in list and i[-5:]!="_qual":
						if i[-1] == "1":
							hap1 = True
						else:
							i[-1]=="2"
							hap1 = False
						if math.isnan(float(block[i][j])) is False :
							if hap1 is True and int(block[i][j])==block["hap"][j]:
								a +=1
							elif hap1 is True and int(block[i][j])==(1-block["hap"][j]):
								b +=1
							elif hap1 is False and int(block[i][j])==(1-block["hap"][j]):
								d +=1
							elif hap1 is False and int(block[i][j])==block["hap"][j]:
								c +=1
							else:
								print "other possibility",j
				cmds ="TwoSidedFastFishersExactTest %i %i %i %i" %(a,b,c,d)
				print cmds
				if (d==0 and a>0) or (a==0 and d>0) or (a==1 and b==0 and c==0 and d==1):
					if a+b+c+d<=1:
						section[j]=1-section[j]
						block["hap"][j] = 1-block["hap"][j]				# possible phase error
						Phase_error +=1
						print "possible phase error",j,a,b,c,d
					else:
						continue
				else:
					m=subprocess.Popen(["TwoSidedFastFishersExactTest",str(a),str(b),str(c),str(d)],stdout=subprocess.PIPE)
					p=m.communicate()[0].strip()
					m.stdout.close()
					if float(p)>=0.5:
						section[j]=1-section[j]
						block["hap"][j] = 1-block["hap"][j]				# possible phase error
						Phase_error +=1
						print "possible phase error",j,a,b,c,d,p
		is_switch = [section[m] for m in section.index]
		print "switch info within island after correcting phasing error",is_switch
		for j in section.index:
			block["switch"][j]=1
	else:
		print "Cannot decide for now"
		for j in section.index:
			a,b,c,d=[0,0,0,0]
			for i in block:
				if i not in list and i[-5:]!="_qual":
					if i[-1] == "1":
						hap1 = True
					else:
						i[-1]=="2"
						hap1 = False
					if math.isnan(float(block[i][j])) is False :
						if hap1 is True and int(block[i][j])==block["hap"][j]:
							a +=1
						elif hap1 is True and int(block[i][j])==(1-block["hap"][j]):
							b +=1
						elif hap1 is False and int(block[i][j])==(1-block["hap"][j]):
							d +=1
						elif hap1 is False and int(block[i][j])==block["hap"][j]:
							c +=1
						else:
							print "other possibility",j
			cmds ="TwoSidedFastFishersExactTest %i %i %i %i" %(a,b,c,d)
			print cmds
			if (d==0 and a>0) or (a==0 and d>0) or (a==1 and b==0 and c==0 and d==1):
				if a+b+c+d<=1:
					section[j]=1-section[j]
					block["hap"][j] = 1-block["hap"][j]				# possible phase error
					Phase_error +=1
					print "possible phase error",j,a,b,c,d
				else:
					continue
			else:
				m=subprocess.Popen(["TwoSidedFastFishersExactTest",str(a),str(b),str(c),str(d)],stdout=subprocess.PIPE)
				p=m.communicate()[0].strip()
				m.stdout.close()
				if float(p)>=0.5:
					if section[j]==1:
						one -=1
					else:
						zero -=1
		print "zero:",zero,"one:",one
		if one > zero:
			for j in section.index:
				if section[j]==0:
					a,b,c,d=[0,0,0,0]
					for i in block:
						if i not in list and i[-5:]!="_qual":
							if i[-1] == "1":
								hap1 = True
							else:
								i[-1]=="2"
								hap1 = False
							if math.isnan(float(block[i][j])) is False :
								if hap1 is True and int(block[i][j])==block["hap"][j]:
									a +=1
								elif hap1 is True and int(block[i][j])==(1-block["hap"][j]):
									b +=1
								elif hap1 is False and int(block[i][j])==(1-block["hap"][j]):
									d +=1
								elif hap1 is False and int(block[i][j])==block["hap"][j]:
									c +=1
								else:
									print "other possibility",j
					cmds ="TwoSidedFastFishersExactTest %i %i %i %i" %(a,b,c,d)
					print cmds
					if (d==0 and a>0) or (a==0 and d>0) or (a==1 and b==0 and c==0 and d==1):
						if a+b+c+d<=1:
							section[j]=1-section[j]
							block["hap"][j] = 1-block["hap"][j]				# possible phase error
							Phase_error +=1
							print "possible phase error",j,a,b,c,d
						else:
							continue
					else:
						m=subprocess.Popen(["TwoSidedFastFishersExactTest",str(a),str(b),str(c),str(d)],stdout=subprocess.PIPE)
						p=m.communicate()[0].strip()
						m.stdout.close()
						if float(p)>=0.5:
							section[j]=1-section[j]
							block["hap"][j] = 1-block["hap"][j]				# possible phase error
							Phase_error +=1
							print "possible phase error",j,a,b,c,d,p
			is_switch = [section[m] for m in section.index]
			print "switch info within island after correcting phasing error",is_switch
			for j in section.index:
				block["switch"][j]=1
		else:
			for j in section.index:
				if section[j]==1:
					a,b,c,d=[0,0,0,0]
					for i in block:
						if i not in list and i[-5:]!="_qual":
							if i[-1] == "1":
								hap1 = True
							else:
								i[-1]=="2"
								hap1 = False
							if math.isnan(float(block[i][j])) is False :
								if hap1 is True and int(block[i][j])==block["hap"][j]:
									a +=1
								elif hap1 is True and int(block[i][j])==(1-block["hap"][j]):
									b +=1
								elif hap1 is False and int(block[i][j])==(1-block["hap"][j]):
									d +=1
								elif hap1 is False and int(block[i][j])==block["hap"][j]:
									c +=1
								else:
									print "other possibility",j
					cmds ="TwoSidedFastFishersExactTest %i %i %i %i" %(a,b,c,d)
					print cmds
					if (d==0 and a>0) or (a==0 and d>0) or (a==1 and b==0 and c==0 and d==1):
						if a+b+c+d<=1:
							section[j]=1-section[j]
							block["hap"][j] = 1-block["hap"][j]				# possible phase error
							Phase_error +=1
							print "possible phase error",j,a,b,c,d
					else:
						m=subprocess.Popen(["TwoSidedFastFishersExactTest",str(a),str(b),str(c),str(d)],stdout=subprocess.PIPE)
						p=m.communicate()[0].strip()
						m.stdout.close()
						if float(p)>=0.5:
							section[j]=1-section[j]
							block["hap"][j] = 1-block["hap"][j]				# possible phase error
							Phase_error +=1
							print "possible phase error",j,a,b,c,d,p
			is_switch = [section[m] for m in section.index]
			print "switch info within island after correcting phasing error",is_switch
	return block,Phase_error

def count_consecutive(is_switch,switch_position,indicator):
	list = []
	pos = []
	for i in range(len(is_switch)):
		if is_switch[i]==indicator:
			if len(list)==0:
				list.append(i)
				pos.append(switch_position[i])
			elif i==list[len(list)-1]+1:
				list.append(i)
				pos.append(switch_position[i])
			elif len(list)<10:
				list = []
				pos =[]
				list.append(i)
				pos.append(switch_position[i])
	if len(list)<10:
		pos = []
		return pos
	else:
		return pos
	
def determine_switch(block,is_switch,switch_position):
	section = pd.Series(is_switch,index=switch_position)
	zero = is_switch.count(0)
	one = is_switch.count(1)
	if one > zero:
		if zero<10:
			is_one = True
			print "is_one"
			for j in switch_position:
				block["switch"][j]=1
		else:
			pos_list=count_consecutive(is_switch,switch_position,0)
			if len(pos_list)==0:
				is_one = True
				print "is_one"
				for j in switch_position:
					block["switch"][j]=1
			else:
				print "partial is_one,find consecutive inconsistent window",pos_list
				for j in switch_position:
					if j in pos_list:
						block["switch"][j]=0
					else:
						block["switch"][j]=1
	else:
		if one<10:
			is_one = False
			print "is_zero"
		else:
			pos_list=count_consecutive(is_switch,switch_position,1)
			if len(pos_list)==0:
				is_one = False
				print "is_zero"
			else:
				print "partial is_zero,find consecutive inconsistent window",pos_list
				for j in switch_position:
					if j in pos_list:
						block["switch"][j]=1
	return block

def compare_haplotype(SNP,BLOCK,block_list,mec,chrom,sample):
	f_hap = open("BLOCK/%s_" %(sample)+chrom,"w")
	Switch_error = 0
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		assert chr==chrom
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
			mec[block_list[i]][7]="Affy"
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
		BLOCK[block_list[i]]["switch"]=pd.Series(switch,index=switch_position)
		mec[block_list[i]][0]=len(list)
		mec[block_list[i]][1]=match+switch_match
		mec[block_list[i]][2]=match
		if match == 0 and switch_match==0:
			print "non-phased: block %s %s" %(chr,pos)
		elif match >= switch_match:
			print "phased by %s: block %s %s non-switch %i %i %i" %(mec[block_list[i]][7],chr,pos,mec[block_list[i]][0],mec[block_list[i]][1],mec[block_list[i]][2])
			if switch_match != 0 :
				switch_error = 0
				is_connect =check_connection(BLOCK[block_list[i]],is_switch,switch_position,resolve)
				print is_connect
				last = ""
				k = 0
				s = 0
				e = 0
				section = []
				last_shift = False
				for k in range(len(is_connect)):
					if k==0:
						last = is_connect[k]
						s = k
					elif is_connect[k]==last:
						continue
					else:
						e = k
						section_pos = switch_position[s:e]
						section_switch = is_switch[s:e]
						print "start",s,"end",e
						print "position",section_pos
						print "switch info",section_switch
						block = determine_switch(BLOCK[block_list[i]],section_switch,section_pos)
						s = k
						last = is_connect[k]
				e = k
				section_pos = switch_position[s:e+1]
				section_switch = is_switch[s:e+1]
				print "start",s,"end",e+1
				print "position",section_pos
				print "switch info",section_switch
				block = determine_switch(BLOCK[block_list[i]],section_switch,section_pos)
				print "switch final:","\t".join(map(str,[BLOCK[block_list[i]]["switch"][n] for n in switch_position]))
				BLOCK[block_list[i]],switch_error = switch_whole(BLOCK[block_list[i]],switch_position,switch_error,resolve,list)
				print "switch_error",switch_error
				Switch_error += switch_error
		else:
			print "phased by %s: block %s %s switch %i %i %i" %(mec[block_list[i]][7],chr,pos,mec[block_list[i]][0],mec[block_list[i]][1],mec[block_list[i]][2])
			if match != 0 :
				switch_error = 0
				is_connect =check_connection(BLOCK[block_list[i]],is_switch,switch_position,resolve)
				print is_connect
				last = ""
				k = 0
				s = 0
				e = 0
				section = []
				last_shift = False
				for k in range(len(is_connect)):
					if k==0:
						last = is_connect[k]
						s = k
					elif is_connect[k]==last:
						continue
					else:
						e = k
						section_pos = switch_position[s:e]
						section_switch = is_switch[s:e]
						print "start",s,"end",e
						print "position",section_pos
						print "switch info",section_switch
						block = determine_switch(BLOCK[block_list[i]],section_switch,section_pos)
						s = k
						last = is_connect[k]
				e = k
				section_pos = switch_position[s:e+1]
				section_switch = is_switch[s:e+1]
				print "start",s,"end",e+1
				print "position",section_pos
				print "switch info",section_switch
				block = determine_switch(BLOCK[block_list[i]],section_switch,section_pos)
				print "switch final:","\t".join(map(str,[BLOCK[block_list[i]]["switch"][n] for n in switch_position]))	
				BLOCK[block_list[i]],switch_error = switch_whole(BLOCK[block_list[i]],switch_position,switch_error,resolve,list)
				print "switch_error",switch_error
				Switch_error += switch_error
			else:
				for k in BLOCK[block_list[i]].index:
					if k not in resolve:
						BLOCK[block_list[i]]["hap"][k] = 1- BLOCK[block_list[i]]["hap"][k]
	return BLOCK,mec,Switch_error

def calc_switch_error(SNP,BLOCK,block_list,mec):
	f_hap = open("BLOCK/%s_" %(SAMPLE)+CHROM,"w")
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
			mec[block_list[i]][7]="Affy"
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
						section_pos = switch_position[s:e]
						section_switch = is_switch[s:e]
						print "start",s,"end",e
						print "position",section_pos
						print "switch info",section_switch
						if is_switch[s]==1:
							switch_error +=1
						s = k
						last = is_switch[k]
				e = k
				section_pos = switch_position[s:e+1]
				section_switch = is_switch[s:e+1]
				print "start",s,"end",e+1
				print "position",section_pos
				print "switch info",section_switch
				if is_switch[s]==1:
					switch_error +=1
				print "switch_error",switch_error
				mec[block_list[i]][8]=switch_error
				Switch_error += switch_error
		else:
			print "phased by %s: block %s %s switch %i %i %i" %(mec[block_list[i]][7],chr,pos,mec[block_list[i]][0],mec[block_list[i]][1],mec[block_list[i]][2])
			if match != 0 :
				switch_error = 0
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
						section_pos = switch_position[s:e]
						section_switch = is_switch[s:e]
						print "start",s,"end",e
						print "position",section_pos
						print "switch info",section_switch
						if is_switch[s]==0:
							switch_error +=1
						s = k
						last = is_switch[k]
				e = k
				section_pos = switch_position[s:e+1]
				section_switch = is_switch[s:e+1]
				print "start",s,"end",e+1
				print "position",section_pos
				print "switch info",section_switch
				if is_switch[s]==0:
					switch_error +=1
				print "switch_error",switch_error
				mec[block_list[i]][8]=switch_error
				Switch_error += switch_error
			for k in BLOCK[block_list[i]].index:
				BLOCK[block_list[i]]["hap"][k] = 1- BLOCK[block_list[i]]["hap"][k]
	return BLOCK,mec,Switch_error

def make_matrix2png(BLOCK,block_list):
	list = ["hap","transmit","untransmit"]
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		pos = block_list[i].split(" ")[1]
		f_matrix = open("BLOCK/%s_matrix2png/%s_" %(SAMPLE,SAMPLE)+chr+"_"+pos+".matrix","w")
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
	wk = "BLOCK/%s_matrix2png/" %(SAMPLE)
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

def separate_clone_into_hap(BLOCK,block_list,mec,decode_pos,chrom):
	wk = "/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/NA19240_pool_"
	MEC_fraction = []
#	f_hap1=open("/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/BLOCK/NA19240_hap1_"+chrom+"_clone.bed","w")
#	f_hap2=open("/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/BLOCK/NA19240_hap2_"+chrom+"_clone.bed","w")
#	f_unknown_hap1=open("/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/BLOCK/NA19240_unknown_hap1_"+chrom+"_clone.bed","w")
#	f_unknown_hap2=open("/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/BLOCK/NA19240_unknown_hap2_"+chrom+"_clone.bed","w")
	for i in range(len(block_list)):
		for clone in BLOCK[block_list[i]]:
			if clone not in ["hap","transmit","untransmit"] and clone[-5:]!="_qual":
				hap,clone_mec,clone_length=determine_clone_hap(BLOCK[block_list[i]],clone,decode_pos)
				mec[block_list[i]][4]+=clone_mec
				MEC_fraction.append(round(float(clone_mec)/clone_length,2))
				if len(clone.split("_"))==4:
					chr = clone.split("_")[1]
					pool = clone.split("_")[0]
					pos = clone.split("_")[2]
				elif len(clone.split("_"))==5:
					continue
					chr = clone.split("_")[1]+"_"+clone.split("_")[2]
					pool = clone.split("_")[0]
					pos = clone.split("_")[3]
				else:
					continue
					chr = clone.split("_")[1]+"_"+clone.split("_")[2]+"_"+clone.split("_")[3]
					pool = clone.split("_")[0]
					pos = clone.split("_")[4]
				find = False
				"""
				clone_name = wk + pool + "/NA19240_pool_" + pool + ".markdup.clone_v2_filter_SNP_all"
				f_clone = open(clone_name,"r")
				for l in f_clone:
					l = l.strip()
					chrom = l.split("\t")[0]
					pos1 = l.split("\t")[1]
					pos2 = l.split("\t")[2]
#					if chrom == chr and int(pos) >= int(pos1) and int(pos) <= int(pos2):
					if chrom == chr and int(pos) == int(pos1):
						find = True
						break
				f_clone.close()
				if mec[block_list[i]][1]!=0:
					file = wk + pool + "/NA19240_pool_" + pool + "." + hap+"."+chrom + ".bed"
					f_out = open(file,"a")
					if hap == "hap1":
						print >>f_out, "%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t255,0,0" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
						print >>f_hap1,"%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t255,0,0" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
					else:
						assert hap == "hap2"
						print >>f_out, "%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t0,0,255" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
						print >>f_hap2,"%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t0,0,255" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
				else:
					file = wk + pool + "/NA19240_pool_" + pool + ".unknown_" + hap+"."+chrom + ".bed"
					f_out = open(file,"a")
					if hap == "hap1":
						print >>f_out, "%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t255,0,0" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
						print >>f_unknown_hap1,"%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t255,0,0" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
					else:
						assert hap == "hap2"
						print >>f_out, "%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t0,0,255" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
						print >>f_unknown_hap2,"%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\t0,0,255" %(chr,pos1,pos2,clone,clone_mec,clone_length,pos1,pos2)
				f_out.close()
				"""
	"""
	plt.figure()
	plt.hist(MEC_fraction,bins=20,facecolor='blue',alpha=0.75)
	plt.xlim(0,0.6)
	plt.xlabel("mec fration in each clone")
	plt.ylabel("freqeuncy")
	plt.title("MEC fraction distribution in clone data")
	plt.savefig("BLOCK/mec_fraction_in_clone_"+chr+".png",format='png')
	"""
	return mec
					
def make_track(hap_file,block_list,mec):
	f = open(hap_file,"r")
	f_out = open("BLOCK/%s_refhap_track_info_" %(SAMPLE)+CHROM+".txt","w")
	for line in f:
		line = line.strip()
		if line[0]=="B":
			start = 0
			end = 0
			chr = ""
			continue
		elif line[0] == "*":
			print >>f_out,"%s\t%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i" %(chr,start,end,MEC[7],MEC[0],MEC[1],MEC[2],MEC[3],MEC[4],MEC[5],MEC[6],MEC[8])
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

def detect_hapmap_wrong(BLOCK,block_list,mec,chrom,sample):
	ABAB=0
#	f = open("BLOCK/hapmap_wrong_v2v2_"+chrom,"w")
	f2 = open("BLOCK/%s_affy_wrong_detail_" %(sample)+chrom,"w")
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		pos = block_list[i].split(" ")[1]
		if mec[block_list[i]][7]!="Hapmap":
			continue
		for j in BLOCK[block_list[i]]["hap"].index:
			if math.isnan(float(BLOCK[block_list[i]]["transmit"][j])) is False:
				if BLOCK[block_list[i]]["hap"][j]!=int(BLOCK[block_list[i]]["transmit"][j]):
#					print >>f,"%s\t%i\t%s\t%i\t%s\t" %(chr,int(j)-1,j,BLOCK[block_list[i]]["hap"][j],BLOCK[block_list[i]]["p_transmit"][j])
					print >>f2,"%s\t%i\t%s\t%i\t%s\t" %(chr,int(pos),j,BLOCK[block_list[i]]["hap"][j],BLOCK[block_list[i]]["transmit"][j])


def compare_with_1000G(BLOCK,block_list,SNP,mec):
	unknown = 0			# our data is unknown
	hapmap_same = 0		# our data is contradictory with 1000G, but hapmap and 1000G are the same
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		pos = block_list[i].split(" ")[1]
		match = 0			# same with 1000G
		appear = 0			# SNP exists in 1000G
		list = BLOCK[block_list[i]].index
		for j in range(len(list)):
			if SNP[chr]["phase_1000G"][list[j]]!="":
				appear +=1
				if int(SNP[chr]["phase_1000G"][list[j]][0])==BLOCK[block_list[i]]["hap"][list[j]]:
					match +=1
				elif BLOCK[block_list[i]]["hap"][list[j]]==0.5:
					unknown +=1
				elif math.isnan(float(BLOCK[block_list[i]]["p_transmit"][j])) is False:
					if int(SNP[chr]["phase_1000G"][list[j]][0])==int(BLOCK[block_list[i]]["p_transmit"][j]) and mec[block_list[i]][7]=="Hapmap":
						hapmap_same +=1
		mec[block_list[i]][5]+=appear
		mec[block_list[i]][6]+=match
	print "Unknown in our phasing",unknown,"Contradict with our phasing but same with Hapmap",hapmap_same
	return mec

def decode_pos_geno(file):
	f = open(file,"r")
	snp = {}
	for i,l in enumerate(f):
		line = l.strip().split("\t")
		snp[line[0]+" "+line[1]]=line[3]+line[4]
	return snp

def read_hap(BLOCK,block_list,snp_decode,chrom):
	f = open("BLOCK/phased_snps_v3_"+chrom,"w")
	f_unphased = open("BLOCK/phased_block_unresolved_snps_v3_"+chrom,"w")
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		pos = block_list[i].split(" ")[1]
		for j in BLOCK[block_list[i]].index:
			if BLOCK[block_list[i]]["hap"][j] !=0.5:
				allele1 = snp_decode[chr+" "+str(j)][int(BLOCK[block_list[i]]["hap"][j])]
				allele2 = snp_decode[chr+" "+str(j)][1-int(BLOCK[block_list[i]]["hap"][j])]
				print >>f,"%s\t%s\t%s\t%s" %(chr,j,allele1,allele2)
			else:
				allele1 = snp_decode[chr+" "+str(j)][int(BLOCK[block_list[i]]["hap"][j])]
				allele2 = snp_decode[chr+" "+str(j)][1-int(BLOCK[block_list[i]]["hap"][j])]
				print >>f_unphased,"%s\t%s\t%s\t%s" %(chr,j,allele1,allele2)

def read_ped(file):
	f=open(file,'r')
	line = f.readline().strip().split()
	family=(line[0],line[1],line[2],line[3])
	index = family.index(SAMPLE)
	return family,index-1
							
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='assign parental allele in block')
	parser.add_argument("--chr", dest='chr',help="chromosome")
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--affy_dir",dest='affy_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/phase3_affy/',help='Affy directory')
	args = parser.parse_args()
	
	SAMPLE = args.sample
	CHROM = args.chr
	all_var_file= args.sample+"_combined_genome_allvars_SNP"
	SNP = create_snp_list(all_var_file)
	print "SNP list set up"
	
	ped_file = args.sample + '.ped'
	family,index = read_ped(ped_file)
	
#	read the locus from Affy	
	Affy_file = args.affy_dir + family[0]+'_Affy_snps.bed'
	locus_affy = read_tab_file(Affy_file)
	SNP,others =read_heter_snp_affy(SNP,locus_affy,index)
	print "SNP info stored"
	print "Number of SNPs used in hapmap:",others
	

	snp_not_shown_file = args.sample + "_SNPs_not_shown"
	SNP_not_shown=SNPs_not_shown(snp_not_shown_file)
	
#	read refhap result and get information
	haplotype_file = args.sample+'_snp_haplotype_detail'
	BLOCK,block_list,decode,mec=read_haplotype(haplotype_file,SNP_not_shown)
	print "haplotype info stored"
	print "Number of blocks",len(block_list)
	print block_list[0],block_list[1]
	decode_pos = decode_number_pos(all_var_file)
	matrix_file = args.sample + '_all_clone_all_fragment_snp_sorted.matrix'
	MEC,mec,BLOCK=read_matrix(matrix_file,BLOCK,decode,decode_pos,mec)
	print BLOCK[block_list[0]]
	
#	BLOCK,mec,switch_error=compare_haplotype(SNP,BLOCK,block_list,mec,args.chr,args.sample)
	BLOCK,mec,switch_error=calc_switch_error(SNP,BLOCK,block_list,mec)
	print "Switch_error",switch_error
#	detect_hapmap_wrong(BLOCK,block_list,mec,args.chr,args.sample)
	
	make_matrix2png(BLOCK,block_list)
	run_matrix2png()
	

#	mec = separate_clone_into_hap(BLOCK,block_list,mec,decode_pos,chr)
#	mec = compare_with_1000G(BLOCK,block_list,SNP,mec)
	
	haplotype_file = 'Refhap/%s_%s.snp.haplotype_detail' %(SAMPLE,CHROM)
	make_track(haplotype_file,block_list,mec)

	'''	
	snp_decode = decode_pos_geno(all_var_file)
	read_hap(BLOCK,block_list,snp_decode,chr)
	print "Summary info as below:"
	print "Total number of blocks in ",chr,":",len(block_list)
	print "Total Switch_error",switch_error
	'''
	
	

	
