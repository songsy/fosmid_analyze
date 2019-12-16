#!/usr/bin/env python
# python calc_concordance_mec_cor_v2.py
# Shiya Song
# 15 Jan 2015
# For each pair of SNP, calc their physical distance, bin physical distance into categories. 
# Calculate pairwise LD between each pair of SNP, bin into categories.
# In each category of LD, calc the genotype concordance, also calc the MEC value
# draw LD vs genotype concordance, MEC value, and the correlation of genotype concordance and MEC value
# For all possible pair of snp, either within same block or within a certian distance(1MB)
# For HapMap data

import sys
from NGS_utils import *
import argparse
import pandas as pd
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import gzip
import pickle

chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
chromLens = genutils.read_chrom_len(chromLenFile)
chromOrder = get_chromOrder("human")

def create_snp_list(file1):
	snp_pos = []
	snp_geno = []
	snp_phase = []
	snp_phase2 = []
	tot = 0
	for i in VcfIterator(file1):
		if i[4] is True and i[3][0]!=i[3][1]:
			tot +=1
			snp_pos.append(i[1])				# store the position of the snp, used for fast search 
			snp_geno.append(i[2])		# store the ref and alt allele of the snp
			snp_phase.append(i[3][0])			# store the phasing information,just one haplotype, 0 or 1
			snp_phase2.append(-1)
	snp = {}
	snp["geno"] = pd.Series(snp_geno,index=snp_pos)
	snp["phase1"] = pd.Series(snp_phase,index=snp_pos)
	snp["phase2"] = pd.Series(snp_phase2,index=snp_pos)
	snp = pd.DataFrame(snp)
	print tot, len(snp_pos)
	return snp

def read_snp(SNP,file2):		# read in hapmap phasing after liftOver
	tot = 0
	same = 0
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
					SNP["phase2"][i[1]]=phase[0]   # store the phasing information,just one haplotype, 0 or 1
					tot +=1
					if SNP["phase2"][i[1]]==SNP["phase1"][i[1]]:
						same +=1
				else:
					print i,SNP["geno"][i[1]]
			else:
				assert i[3][0]!=i[3][1]
				SNP["phase2"][i[1]]=i[3][0]		# store the phasing information,just one haplotype, 0 or 1
				tot +=1
				if SNP["phase2"][i[1]]==SNP["phase1"][i[1]]:
					same +=1
	flip = 1 if float(same)/tot<0.5 else 0 # need to flip it
	print tot,same,flip
	return SNP,flip

def read_heter_snp_hapmap(SNP,hap_file,index):		
	het = 0
	tot = 0
	same = 0
	f=gzip.open(hap_file,'r')
	for line in f:
		line = line.strip().split(' ')
		chr = line[0]
		pos = int(line[2])
		parental = line[index]+line[index+1]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			if parental[0]!=parental[1]:
				het +=1
			else:
				continue
			try:
				snp_info = SNP["geno"][pos]
			except KeyError:
				continue
			map = {SNP["geno"][pos][0]:'0',SNP["geno"][pos][1]:'1'}
			try:
				info1 = map[parental[0]]
				SNP["phase2"][pos]=info1
			except KeyError:
				continue
			if SNP["phase2"][pos]==SNP["phase1"][pos]:
				same +=1
			tot +=1
		else:
			break
	flip = 1 if float(same)/tot<0.5 else 0 # need to flip it
	print 'HapMap',CHROM,tot,same,flip
	return SNP,flip


def merge_BLOCK_SNP(block_file,SNP):
	block_list_left = []
	block_list_right = []
	for i in BedIterator(block_file):
		if args.prism ==1:
			if i[4]!='interleaving' and i[4]!='None,interleave' and i[1]!=i[2]:
				block_list_left.append(int(i[1]))
				block_list_right.append(int(i[2]))
		else:
			if len(block_list_left)==0:
				block_list_left.append(int(i[1]))
				block_list_right.append(int(i[2]))
			elif int(i[1])>block_list_right[-1]:
				block_list_left.append(int(i[1]))
				block_list_right.append(int(i[2]))
	for i in range(1,len(block_list_left)):
		assert block_list_left[i]>block_list_left[i-1]
	list = SNP.index
	block_index = []
	for j in range(len(list)):
		find = find_interval(list[j],block_list_left)
		try:
			assert list[j]>=block_list_left[find] and list[j]<=block_list_right[find]
			block_index.append(find)   # the index of the block that the SNP belongs to 
		except AssertionError:
#			print 'Assertion Error',list[j],block_list_left[find],block_list_right[find]
			block_index.append(-1)
	print len(list),'-1:',block_index.count(-1)
	SNP["block_index"] = pd.Series(block_index,index=list)   # store the index of the block that this snp belongs to 
	return SNP

def determine_flip(SNP,block):
	tot = 0
	same = 0
	for pos in block:
		tot +=1
		if SNP["phase1"][pos]==SNP["phase2"][pos]:
			same +=1
	flip = 1 if float(same)/tot<0.5 else 0 # need to flip it
	return flip

def read_haplotype(hap_file,SNP):
	list = SNP.index
	SNP["block_index"] = pd.Series(map(int,np.zeros(len(list))),index=list)   # store the index of the block that this snp belongs to
	f = open(hap_file,"r")
	flip_array = {}
	start = True
	for line in f:
		line = line.strip()
		if line[0]=="B":
			start = True
		elif line[0]=="*":
			if start == False:
				if len(block)!=0:
					flip_array[int(pos)]= determine_flip(SNP,block)
				else:
					flip_array[int(pos)] = 0
		else:
			col = line.split("\t")
			if start == True:
				pos = col[4]
				chr = col[3]
				if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
					continue
				elif chr==CHROM:
					start = False			
					pos = col[4]
					block = []				# store the position in each block
				else:
					break
			if col[1] == "-":
				try:
					m = SNP["block_index"][int(col[4])]
				except KeyError:
					continue
				if m ==0:
					block.append(int(col[4]))
					SNP["block_index"][int(col[4])]=int(pos)
			else:
				try:
					m = SNP["block_index"][int(col[4])]
				except KeyError:
					continue
				block.append(int(col[4]))
				SNP["block_index"][int(col[4])]=int(pos)		# link the refhap result (# of the snp) to the first snp of that block	
	f.close()
	return SNP,flip_array

def read_sample(file,popname):
	pop_list=[]
	f = open(file,'r')
	index = 0
	for line in f:
		line = line.rstrip().split(' ')
		if line[0]=='sample' or line[0]=='ID':
			continue
		if popname is None:
			pop_list.append(index)
		else:
			if line[1]==popname:
				pop_list.append(index)
		index+=1
	return pop_list

def read_1000G_legend(legend_file,snp_pos_list):
	snp_pos = []
	snp_geno = []
	snp_LD = []
	f=gzip.open(legend_file)
	pos_index =[]
	index_to_pos=[]
	snp = pd.Series(map(int,np.zeros(len(snp_pos_list))),index=map(int,snp_pos_list))
	for line in f:
		line = line.rstrip().split(' ')
		if line[0]=='id':
			continue
		pos=int(line[1])
		if line[4]!='SNP' and line[4]!='Biallelic_SNP':
			pos_index.append(0)
			index_to_pos.append(pos)
			continue
		try:
			m = snp[pos]
		except KeyError:
			pos_index.append(0)
			index_to_pos.append(pos)
			continue
		pos_index.append(1)
		index_to_pos.append(pos)
	print pos_index.count(0),pos_index.count(1)
	return pos_index,index_to_pos

def read_reference_SNP(SNP,hap_file,pop_list,pos_index,index_to_pos):
	index = 0
	list = SNP.index
	snp_pos = []
	snp_LD = []
	prev_reference_hap=[]
	reference = []
	first = True
	for i in BedIterator(hap_file):
		if pos_index[index]==1:
			pos = index_to_pos[index]
			reference_hap = []
			for j in pop_list:
				reference_hap.append(i[j*2])
				reference_hap.append(i[j*2+1])
			snp_pos.append(pos)
			reference_hap = map(int,reference_hap)
			reference.append(reference_hap)
		index+=1
		if index%10000==0:
			print index,'done'
	try:
		SNP["reference"] = pd.Series([1 for m in snp_pos],index=snp_pos)
	except TypeError:
		new_column = pd.Series([1 for m in snp_pos],index=snp_pos)
		SNP["reference"] = new_column.reset_index(level=0, drop=True)
	return SNP,snp_pos,reference

def read_reference_SNP_v2(SNP,hap_file):
	index = 0
	list = SNP.index
	snp_pos = []
	snp_LD = []
	prev_reference_hap=[]
	reference = []
	first = True
	f=gzip.open(hap_file,'r')
	for line in f:
		line = line.strip().split(' ')
		pos = int(line[2])
		chr = line[0]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			reference_hap = line[3:]
			try:
				code_map = {SNP["geno"][pos][0]:0,SNP["geno"][pos][1]:1}
				reference_hap = [code_map[i] for i in reference_hap]
				reference.append(reference_hap)
				snp_pos.append(pos)
			except KeyError:
				continue
		index+=1
		if index%10000==0:
			print index,'done'
	try:
		SNP["reference"] = pd.Series([1 for m in snp_pos],index=snp_pos)
	except TypeError:
		new_column = pd.Series([1 for m in snp_pos],index=snp_pos)
		SNP["reference"] = new_column.reset_index(level=0, drop=True)
		print 'TypeError'
	return SNP,snp_pos,reference

def calc_LD(hA,hB):
	pA=float(hA.count(0))/len(hA)
	pB=float(hB.count(0))/len(hB)
	AB = [1 if i==0 and j==0 else 2 for (i,j) in zip(hA,hB)]
	pAB = float(AB.count(1))/len(AB)
	D=pAB-pA*pB
	if pA==0 or pB == 0 or pA==1 or pB==1:
		Dprime =1
	else:
#		print D,pA,pB
		Dprime = D/min((1-pA)*pB,pA*(1-pB)) if D>=0 else -D/min((1-pA)*pA,pB*(1-pB))
	return Dprime

def determine_clone_hap(block,clone,pos_list,mec_vec,count_vec):
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
	for i in range(len(pos_list)):
		if math.isnan(float(block[clone][pos_list[i]])) is False:
			if hap == "hap1":
				if int(block[clone][pos_list[i]])!=block["hap"][pos_list[i]]:
					mec_vec[i]+=1
			elif hap == "hap2":
				if int(block[clone][pos_list[i]])==block["hap"][pos_list[i]]:
					mec_vec[i]+=1
			count_vec[i]+=1
	return hap,mec_vec,count_vec

def calc_BLOCK_mec(BLOCK,block_list_new,SNP,flip,flip_array):
	list = SNP.index
	MEC_vec = []
	COUNT_vec = []
	POS_list = []
	index = 0
	TOT = 0
	SAME = 0
	for i in block_list_new:
		same = 0
		tot = 0
		for j in BLOCK[i].index:
			tot +=1
			if args.prism=='1':
				flip = flip_array[index]
			if flip==1:
#				if BLOCK[i]["phase1"][j]==1-BLOCK[i]["phase2"][j]:
				if SNP["phase1"][int(j)]!=SNP["phase2"][int(j)]:
					assert SNP["phase1"][int(j)]==1-SNP["phase2"][int(j)]
					SNP["phase_compare"][int(j)]=1    # same haplotype
					same +=1
			else:
#				if BLOCK[i]["phase1"][j]==BLOCK[i]["phase2"][j]:
				if SNP["phase1"][int(j)]==SNP["phase2"][int(j)]:
					SNP["phase_compare"][int(j)]=1    # same haplotype
					same +=1
		TOT += tot
		SAME += same
		mec_vec=map(int,np.zeros(len(BLOCK[i].index)))
		count_vec=map(int,np.zeros(len(BLOCK[i].index)))
		pos_list = BLOCK[i].index
		for clone in BLOCK[i]:
			if clone not in ["hap","transmit","untransmit","switch","phase1","phase2"] and clone[-5:]!="_qual":
				hap,mec_vec,count_vec=determine_clone_hap(BLOCK[i],clone,pos_list,mec_vec,count_vec)
		MEC_vec += mec_vec
		COUNT_vec += count_vec
		POS_list +=  map(int,pos_list)
		index +=1
	SNP["MEC"] = pd.Series(MEC_vec,index=POS_list)
	SNP["clone_num"] = pd.Series(COUNT_vec,index=POS_list)
	return SNP

def make_range():
	dist_range = [0,5,10,100,1e06]
	LD_range = np.linspace(0,1,num=51)

def find_interval(pos,list):
	left = 0
	right = len(list)
	find = -1
	while right - left >0:
		midpoint = (right + left)/2
		if pos < list[midpoint]:
			right = midpoint
		elif pos > list[midpoint]:
			left = midpoint+1
		elif pos == list[midpoint]:
			left = midpoint
			find = midpoint
			break
		if pos<list[midpoint] and midpoint!=0:
			midpoint=midpoint-1
	return midpoint

def make_curve(SNP,dist_range,LD_range):
	list = SNP[SNP["LD"]<=1].index
	abundance = [[0 for x in range(len(LD_range))] for x in range(len(dist_range))] 
	genotype_concordance = [[[] for x in range(len(LD_range))] for x in range(len(dist_range))] 
	MEC = [[[] for x in range(len(LD_range))] for x in range(len(dist_range))] 
	MEC_per_clone = [[[] for x in range(len(LD_range))] for x in range(len(dist_range))] 
	for i in range(len(list)-1):
		pos_dist = float(list[i+1]-list[i])/1000
		index_i = find_interval(pos_dist,dist_range)
		index_j = find_interval(SNP["LD"][list[i]],map(float,LD_range))
#		print index_i,index_j,pos_dist,SNP["LD"][list[i]],SNP["MEC"][list[i]],SNP["clone_num"][list[i]]
		abundance[index_i][index_j]+=1
		genotype_concordance[index_i][index_j].append(SNP["phase_compare"][list[i]])
		MEC[index_i][index_j].append(SNP["MEC"][list[i]])
		MEC_per_clone[index_i][index_j].append(SNP["MEC"][list[i]]/SNP["clone_num"][list[i]])
#		except IndexError:
#			print index_i,index_j,pos_dist,SNP["LD"][list[i]],SNP["MEC"][list[i]],SNP["clone_num"][list[i]]
	print 'done'
	col = ['r','b','g','y','k']
#	fig = plt.figure()
	for i in range(len(dist_range)):
		Genotype_concordance = []
		for j in range(len(LD_range)):
			Genotype_concordance.append(float(np.sum(genotype_concordance[i][j]))/abundance[i][j])
			print i,j,abundance[i][j],float(np.sum(genotype_concordance[i][j]))/abundance[i][j],np.mean(MEC[i][j]),np.std(MEC[i][j]),np.mean(MEC_per_clone[i][j]),np.std(MEC_per_clone[i][j])
		print Genotype_concordance
		plt.plot(LD_range,Genotype_concordance,col[i])
	plt.show
#	fig.savefig('%s_%s_LD_curve.pdf' %(SAMPLE,CHROM))

def calc_pairwise_LD_concordance(SNP,reference,snp_pos,flip,flip_array):
	dist1 = []
	LD1 = []
	concord1 = []
	dist2 = []
	LD2 = []
	concord2 = []
	list = SNP.index
	for i in range(len(list)):
		for j in range(i+1,len(list)):
			pos_dist=list[j]-list[i]
			if abs(pos_dist)>1000000:
				break
			pos_LD = calc_LD(reference[i],reference[j])
			if SNP["phase1"][list[i]]==SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==SNP["phase2"][list[j]]:
				pos_concord1 = 1 if flip ==0 else 0
			elif SNP["phase1"][list[i]]==1-SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==1-SNP["phase2"][list[j]]:
				pos_concord1 = 1 if flip ==1 else 0
			else:
				pos_concord1 = 0
			dist1.append(pos_dist)
			LD1.append(pos_LD)
			concord1.append(pos_concord1)
			a=SNP["block_index"][list[i]]
			b=SNP["block_index"][list[j]]
			if a!=b or a==0 or b==0:
				continue
			if SNP["phase1"][list[i]]==SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==SNP["phase2"][list[j]]:
				pos_concord2 = 1 if flip_array[a] ==0 else 0
			elif SNP["phase1"][list[i]]==1-SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==1-SNP["phase2"][list[j]]:
				try:
					pos_concord2 = 1 if flip_array[a] ==1 else 0
				except KeyError:
					print list[i],list[j],a,b
			else:
				pos_concord2 = 0
			dist2.append(pos_dist)
			LD2.append(pos_LD)
			concord2.append(pos_concord2)
	Compare1 = pd.DataFrame({'dist':dist1,'LD':LD1,'concord':concord1})
	Compare2 = pd.DataFrame({'dist':dist2,'LD':LD2,'concord':concord2})
	print 'length',len(dist1),len(dist2)
	return Compare1,Compare2

def calc_pairwise_LD_concordance_within_block(SNP,reference,snp_pos,flip_array):
	SNP = SNP[SNP["block_index"]>0]
	dist = []
	LD = []
	concord = []
	list = SNP.index
	for i in range(len(list)):
		for j in range(i,len(list)):
			pos_dist=list[j]-list[i]
			if abs(pos_dist)>1000000:
				break
			a=SNP["block_index"][list[i]]
			b=SNP["block_index"][list[j]]
			if a!=b or a==0 or b==0:
				continue
			pos_LD = calc_LD(reference[i],reference[j])
			if SNP["phase1"][list[i]]==SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==SNP["phase2"][list[j]]:
				pos_concord1 = 1 if flip_array[a] ==0 else 0
			elif SNP["phase1"][list[i]]==1-SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==1-SNP["phase2"][list[j]]:
				pos_concord1 = 1 if flip_array[a] ==1 else 0
			else:
				pos_concord1 = 0
			dist.append(pos_dist)
			LD.append(pos_LD)
			concord.append(pos_concord1)
	Compare = pd.DataFrame({'dist':dist,'LD':LD,'concord':concord})
	print 'length',len(dist)
	return Compare

def calc_pairwise_LD_concordance_v2(SNP,reference,snp_pos,flip):   # use flip instead of flip_array
	dist_range = [0,5,10,100,1e06]
	LD_range = np.linspace(0,1,num=21)
	abundance1 = [[0 for x in range(len(LD_range))] for x in range(len(dist_range))] 
	genotype_concordance1 = [[0 for x in range(len(LD_range))] for x in range(len(dist_range))] 
	SNP = SNP[SNP["block_index"]>0]
	dist = []
	LD = []
	concord = []
	list = SNP.index
	tot = 0
	for i in range(len(list)):
		for j in range(i,len(list)):
			pos_dist=float(list[j]-list[i])/1000
			if abs(pos_dist)>1000:
				break
			a=SNP["block_index"][list[i]]
			b=SNP["block_index"][list[j]]
			if a!=b or a==0 or b==0:
				continue
			pos_LD = calc_LD(reference[i],reference[j])
			if SNP["phase1"][list[i]]==SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==SNP["phase2"][list[j]]:
				pos_concord1 = 1 if flip ==0 else 0
			elif SNP["phase1"][list[i]]==1-SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==1-SNP["phase2"][list[j]]:
				pos_concord1 = 1 if flip ==1 else 0
			else:
				pos_concord1 = 0
			index_i = find_interval(pos_dist,dist_range)
			index_j = find_interval(pos_LD,map(float,LD_range))
			abundance1[index_i][index_j]+=1
			genotype_concordance1[index_i][index_j]+=pos_concord1
			tot +=1
	print 'tot',tot
	return abundance1,genotype_concordance1

def calc_pairwise_LD_concordance_within_block_v2(SNP,reference,snp_pos,flip_array):
	dist_range = [0,5,10,100,1e06]
	LD_range = np.linspace(0,1,num=51)
	abundance1 = [[0 for x in range(len(LD_range))] for x in range(len(dist_range))] 
	genotype_concordance1 = [[0 for x in range(len(LD_range))] for x in range(len(dist_range))] 
	SNP = SNP[SNP["block_index"]>0]
	dist = []
	LD = []
	concord = []
	list = SNP.index
	tot = 0
	for i in range(len(list)):
		for j in range(i,len(list)):
			pos_dist=float(list[j]-list[i])/1000
			if abs(pos_dist)>1000:
				break
			a=SNP["block_index"][list[i]]
			b=SNP["block_index"][list[j]]
			if a!=b or a==0 or b==0:
				continue
			pos_LD = calc_LD(reference[i],reference[j])
			if SNP["phase1"][list[i]]==SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==SNP["phase2"][list[j]]:
				pos_concord1 = 1 if flip_array[a] ==0 else 0
			elif SNP["phase1"][list[i]]==1-SNP["phase2"][list[i]] and SNP["phase1"][list[j]]==1-SNP["phase2"][list[j]]:
				pos_concord1 = 1 if flip_array[a] ==1 else 0
			else:
				pos_concord1 = 0
			index_i = find_interval(pos_dist,dist_range)
			index_j = find_interval(pos_LD,map(float,LD_range))
			abundance1[index_i][index_j]+=1
			genotype_concordance1[index_i][index_j]+=pos_concord1
			tot +=1
	print 'tot',tot
	return abundance1,genotype_concordance1


if __name__=="__main__":
	parser = argparse.ArgumentParser(description='assign parental allele in block')
	parser.add_argument("--chr", dest='chr',help="chromosome")
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	parser.add_argument("--prism", dest='prism',default=0,help="prism or not")   # for NA20847, HG03428, compare within block
	parser.add_argument("--phase", dest='phase',default='1KGphase3',help="phase method")
	parser.add_argument("--phase_panel", dest='phase_panel',default='1KGphase1',help="reference panel")
	parser.add_argument("--ref_dir", dest='ref_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono/',help="directory for 1KG reference file")
	parser.add_argument("--pop", dest='pop',help="population name")
	parser.add_argument("--analyze_dir", dest='analyze_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/',help="directory for BLOCK file")
	parser.add_argument("--within_block",dest='within_block',default=0)
	parser.add_argument("--hapmap_dir",dest='hapmap_dir',default='/home/jmkidd/kidd-lab/genomes/snp-sets/hapmap3/',help='HapMap directory')

	args = parser.parse_args()
	SAMPLE = args.sample
	POP = args.pop
	if True:
		CHROM = 'chr'+str(args.chr)
		if args.prism=='1':
			file1 = '%s%s/gVCF_calls/%s.%s.prism.v2.phased.vcf.gz' %(args.wgs_dir,args.sample,args.sample,CHROM)
		elif args.sample=='NA12878':
			file1 = '%s%s/all_sites/%s.%s.fosmid.phase.vcf.gz' %(args.wgs_dir,args.sample,args.sample,CHROM)
		else:
			file1 = '%s%s/gVCF_calls/%s.%s.fosmid.v2.phased.vcf.gz' %(args.wgs_dir,args.sample,args.sample,CHROM)
		if args.phase=='shapeit':
			if args.phase_panel=='1KGphase1':
				file2 = '%s%s/gVCF_calls/%s.%s.phased.vcf.gz' %(args.wgs_dir,args.sample,args.sample,CHROM)
				if args.sample=='NA12878':
					file2 = '%s%s/all_sites/%s.%s.phased.vcf.gz' %(args.wgs_dir,args.sample,args.sample,CHROM)
			elif args.phase_panel=='1KGphase3':
				file2 = '%s%s/gVCF_calls/%s.%s.1KGphase3.vcf.gz' %(args.wgs_dir,args.sample,args.sample,CHROM)
				if args.sample=='NA12878':
					file2 = '%s%s/all_sites/%s.%s.1KGphase3.vcf.gz' %(args.wgs_dir,args.sample,args.sample,CHROM)
		elif args.phase=='hapmap':
			file2= args.hapmap_dir + '%s/%s_phased_snp.txt.gz' %(POP,POP)
		else:
			file2 = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/%s/%s.%s.1KGphase3.snp.vcf.gz' %(args.sample,args.sample,CHROM)
		
		SNP = create_snp_list(file1)
		
		if args.phase=='hapmap':    #	read the locus from HapMap	
			header = gzip.open(args.hapmap_dir + '%s/hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_%s.unr.phased.gz' %(POP,POP.lower()),'r').readline().strip().split()
			index_sample = header.index(args.sample+'_A')+1
			print 'index',index_sample
			SNP,flip =read_heter_snp_hapmap(SNP,file2,index_sample)
		else:
			SNP,flip = read_snp(SNP,file2)

		print 'original snp set',len(SNP)
		SNP = SNP[SNP["phase1"]>=0]
		print 'phase1 snp set',len(SNP)
		SNP = SNP[SNP["phase2"]>=0]
		print 'phase2 snp set',len(SNP)

		print 'finish loading SNP'
		print 'flip',flip

#		block_file = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/BLOCK/%s_gvcf/%s_refhap_track_info_%s.txt' %(args.sample,args.sample,CHROM)

#		SNP = merge_BLOCK_SNP(block_file,SNP)

		hap_file = args.analyze_dir+args.sample+'_gvcf_snp_haplotype_detail'
		SNP,flip_array = read_haplotype(hap_file,SNP)

		# read reference haplotypes and calculate LD

		if args.phase_panel=='1KGphase3':
			sample_file = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/1000GP_Phase3/1000GP_Phase3.sample'
			haplotype_file = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/1000GP_Phase3/1000GP_Phase3_%s.hap.gz' %(CHROM)
			legend_file = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/1000GP_Phase3/1000GP_Phase3_%s.legend.gz' %(CHROM)
		else:
			sample_file = args.ref_dir+'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample'
			haplotype_file = args.ref_dir+'ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz' %(CHROM)
			legend_file = args.ref_dir + 'ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz' %(CHROM)
		
		if args.phase=='hapmap':
			haplotype_file = file2
			SNP,snp_pos,reference=read_reference_SNP_v2(SNP,haplotype_file)
		else:
			pop_list=read_sample(sample_file,args.pop)
			print 'POP:',args.pop,len(pop_list)
			snp_pos_list = list(SNP.index)
			pos_index,index_to_pos = read_1000G_legend(legend_file,snp_pos_list)
			SNP,snp_pos,reference=read_reference_SNP(SNP,haplotype_file,pop_list,pos_index,index_to_pos)

		SNP = SNP[SNP["reference"]==1]
		print 'reference snp set',len(list(SNP.index)),len(snp_pos)

#		assert list(SNP.index)==snp_pos

		print 'finish load reference haplotypes, LD calulation'
		# calculate genotype comparison, local MEC, local clone coverage for each SNP
		if args.within_block=='1':
#			Compare = calc_pairwise_LD_concordance_within_block(SNP,reference,snp_pos,flip_array)
			if args.prism=='1':
				abundance1,genotype_concordance1=calc_pairwise_LD_concordance_within_block_v2(SNP,reference,snp_pos,flip_array)
			else:
				abundance1,genotype_concordance1=calc_pairwise_LD_concordance_v2(SNP,reference,snp_pos,flip)
			dbfile = open('%s_%s_%s_%s_SNP_pair_LD_within_block_pickle' %(SAMPLE,CHROM,args.phase,args.phase_panel),'wb')
			print abundance1
			print genotype_concordance1
			pickle.dump(abundance1,dbfile)
			pickle.dump(genotype_concordance1,dbfile)
		else:
			Compare1,Compare2 = calc_pairwise_LD_concordance(SNP,reference,snp_pos,flip,flip_array)
			dbfile = open('%s_%s_%s_%s_SNP_pair_pickle' %(SAMPLE,CHROM,args.phase,args.phase_panel),'wb')
			pickle.dump(Compare1,dbfile)
			pickle.dump(Compare2,dbfile)
		print CHROM, 'finished'
	'''
	dbfile = open('%s_%s_SNP_pickle' %(SAMPLE,CHROM),'rb')
	SNP = pickle.load(dbfile)

	dist_range = [0,5,10,50,1e06]
	LD_range = np.linspace(0,1,num=21)
	make_curve(SNP,dist_range,LD_range)
	'''

