#!/usr/bin/env python
# python calc_concordance_mec_cor.py
# Shiya Song
# 7 Jan 2015
# For each pair of SNP, calc their physical distance, bin physical distance into categories. 
# Calculate pairwise LD between each pair of SNP, bin into categories.
# In each category of LD, calc the genotype concordance, also calc the MEC value
# draw LD vs genotype concordance, MEC value, and teh correlation of genotype concordance and MEC value

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
				else:
					print i,SNP["geno"][i[1]]
			else:
				assert i[3][0]!=i[3][1]
				SNP["phase2"][i[1]]=i[3][0]		# store the phasing information,just one haplotype, 0 or 1
	return SNP


def merge_BLOCK_SNP(block_file,BLOCK,SNP):
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
	list = SNP.index
	tot = 0
	same = 0
	is_switch = []
	label = []
	switch = []
	switch_position = []
	block_list_new = []
	flip_array = []    # just within this block, if flip is necessary
	start = True
	position = []
	for j in range(len(list)):
		if block_index[j]==-1:
			continue
		if SNP["phase1"][list[j]]!=-1 and SNP["phase2"][list[j]]!=-1:
			if SNP["phase1"][list[j]] == SNP["phase2"][list[j]]:
				same +=1
			tot +=1
			if start is True:
				cur_label=block_index[j]
				start = False
				position=[list[j]]
				phase1=[SNP["phase1"][list[j]]]
				phase2=[SNP["phase2"][list[j]]]
				continue
			if block_index[j]==cur_label:
				position.append(list[j])
				phase1.append(SNP["phase1"][list[j]])
				phase2.append(SNP["phase2"][list[j]])
			else:
				cur_block=CHROM+" "+str(block_list_left[cur_label])
				block_list_new.append(cur_block)
				BLOCK[cur_block]["phase1"]=pd.Series(phase1,index=map(str,position))
				BLOCK[cur_block]["phase2"]=pd.Series(phase2,index=map(str,position))
				m = [abs(a-b) for a,b in zip(phase1,phase2)]
				if float(sum(m))/len(m)>0.5:
					flip_array.append(1)
				else:
					flip_array.append(0)
				cur_label=block_index[j]
				position=[list[j]]
				phase1=[SNP["phase1"][list[j]]]
				phase2=[SNP["phase2"][list[j]]]
	if len(position)!=0:
		cur_block=CHROM+" "+str(block_list_left[cur_label])
		block_list_new.append(cur_block)
		BLOCK[cur_block]["phase1"]=pd.Series(phase1,index=map(str,position))
		BLOCK[cur_block]["phase2"]=pd.Series(phase2,index=map(str,position))
		m = [abs(a-b) for a,b in zip(phase1,phase2)]
		if float(sum(m))/len(m)>0.5:
			flip_array.append(1)
		else:
			flip_array.append(0)
	print 'tot',tot,'same',same,float(same)/tot
	flip = 1 if float(same)/tot<0.5 else 0 # need to flip it
	return BLOCK,block_list_new,flip,flip_array

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
#		find = find_pos(pos,map(int,snp_pos_list))
#		if find<=0:
#			pos_index.append(0)
#			index_to_pos.append(pos)
#			continue
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
	first = True
	for i in BedIterator(hap_file):
		if pos_index[index]==1:
			pos = index_to_pos[index]
			try:
				m = SNP["phase1"][int(pos)]
			except KeyError:
				continue
			reference_hap = []
			for j in pop_list:
				reference_hap.append(i[j*2])
				reference_hap.append(i[j*2+1])
			snp_pos.append(pos)
			reference_hap = map(int,reference_hap)
			if first is True:
				prev_reference_hap = reference_hap
				first = False
			else:
				LD = calc_LD(prev_reference_hap,reference_hap)
				snp_LD.append(LD)
				prev_reference_hap = reference_hap
		index+=1
		if index%10000==0:
			print index,'done'
	print np.mean(snp_LD),np.std(snp_LD),snp_LD.count(0),snp_LD.count(1),min(snp_LD),max(snp_LD)
	snp_LD.append(snp_LD[-1])
	try:
		SNP["LD"] = pd.Series(snp_LD,index=snp_pos)
	except TypeError:
		print 'TypeError,LD'
		new_column = pd.Series(snp_LD,index=snp_pos)
		SNP["LD"] = new_column.reset_index(level=0, drop=True)
	try:
		SNP["phase_compare"] = pd.Series(map(int,np.zeros(len(snp_pos))),index=snp_pos)
	except TypeError:
		print 'TypeError,phase_compare'
		new_column = pd.Series(map(int,np.zeros(len(snp_pos))),index=snp_pos)
		SNP["phase_compare"] = new_column.reset_index(level=0, drop=True)
#	SNP["MEC"] = pd.Series(map(int,list(np.zeros(len(snp_pos)))),index=snp_pos)
#	SNP["clone_num"] = pd.Series(map(int,list(np.zeros(len(snp_pos)))),index=snp_pos)
	return SNP

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

def get_snp_pos(BLOCK,block_list_new):
	snp_pos_list = []
	for i in block_list_new:
		BLOCK[i]=BLOCK[i][BLOCK[i]["phase1"]>=0]
		snp_pos_list += list(BLOCK[i].index)
	return snp_pos_list

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
	LD_range = np.linspace(0,1,num=21)

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

	args = parser.parse_args()
	SAMPLE = args.sample
#	f_out1 = open(args.sample+"_shapeit_switch_error_mec.txt","w")
#	f_out2 = open(args.sample+"_1KGphase3_switch_error_mec.txt","w")
#	switch_error_file1=args.sample+"_shapeit_switch_error.txt"
#	switch_error_file2=args.sample+"_1KGphase3_switch_error.txt"

#	for chr in range(11,22):
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
		else:
			file2 = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/%s/%s.%s.1KGphase3.snp.vcf.gz' %(args.sample,args.sample,CHROM)
		SNP = create_snp_list(file1)
		SNP = read_snp(SNP,file2)
		dbfile=open('%sBLOCK/%s_gvcf/%s_BLOCK_' %(args.analyze_dir,SAMPLE,SAMPLE) +CHROM+'_pickle','rb')
		BLOCK=pickle.load(dbfile) 
		block_list=pickle.load(dbfile)
#		mec=pickle.load(dbfile)

		print 'finish loading BLOCK'

		block_file = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/BLOCK/%s_gvcf/%s_refhap_track_info_%s.txt' %(args.sample,args.sample,CHROM)

		BLOCK,block_list_new,flip,flip_array = merge_BLOCK_SNP(block_file,BLOCK,SNP)
		print 'flip:',flip,flip_array.count(0),flip_array.count(1)

		snp_pos_list = get_snp_pos(BLOCK,block_list_new)

		print 'finish update BLOCK;',len(snp_pos_list)

		# read reference haplotypes and calculate LD
		if args.phase_panel=='1KGphase3':
			sample_file = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/1000GP_Phase3/1000GP_Phase3.sample'
			haplotype_file = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/1000GP_Phase3/1000GP_Phase3_%s.hap.gz' %(CHROM)
			legend_file = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/1000GP_Phase3/1000GP_Phase3_%s.legend.gz' %(CHROM)
		else:
			sample_file = args.ref_dir+'ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample'
			haplotype_file = args.ref_dir+'ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz' %(CHROM)
			legend_file = args.ref_dir + 'ALL.%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz' %(CHROM)
		pop_list=read_sample(sample_file,args.pop)
		print 'POP:',args.pop,len(pop_list)
		pos_index,index_to_pos = read_1000G_legend(legend_file,snp_pos_list)
		SNP=read_reference_SNP(SNP,haplotype_file,pop_list,pos_index,index_to_pos)

		print 'finish load reference haplotypes, LD calulation'
		# calculate genotype comparison, local MEC, local clone coverage for each SNP
		SNP = calc_BLOCK_mec(BLOCK,block_list_new,SNP,flip,flip_array)

		print 'finish MEC calculation'
		dbfile = open('%s_%s_%s_%s_SNP_pickle' %(SAMPLE,CHROM,args.phase,args.phase_panel),'wb')
		pickle.dump(SNP,dbfile)

		print CHROM, 'finished'
	'''
	dbfile = open('%s_%s_SNP_pickle' %(SAMPLE,CHROM),'rb')
	SNP = pickle.load(dbfile)

	dist_range = [0,5,10,50,1e06]
	LD_range = np.linspace(0,1,num=21)
	make_curve(SNP,dist_range,LD_range)
	'''

