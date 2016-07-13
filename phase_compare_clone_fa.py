#!/usr/bin/env python
# python phase_compare_clone_fa.py
# Shiya Song
# 23 Sept 2014
# Comparing phasing statistics in terms of switch error and mean distance between switch errors

import os
import signal
import glob
import pandas as pd
import numpy as np
import math
from NGS_utils import *
import pickle
import subprocess
import argparse
def compare_seq(seq_clone,seq_hap,position,clone_seq,hap_seq):
	j = 0
	length = 0
#	print "in",len(seq_clone),len(seq_hap)
	for i in range(len(seq_hap)):
		if seq_hap[i]!="-":
			j +=1
			if seq_clone[i]!=seq_hap[i] and seq_clone[i]!="-":
				hap_seq.append(seq_hap[i])
				clone_seq.append(seq_clone[i])
				position.append(j)
				if seq_hap[i]!='N':
					different +=1
		if seq_hap[i]!='N' and seq_hap[i]!="-":
			length +=1
	return position,clone_seq,hap_seq,different,length
	
def compare_align(name,Clone):
	clone_hg19_seq=[]
	clone_hap1_seq=[]
	clone_hap2_seq=[]
	seq_position=[]
	hap1_seq=[]
	hap1_position=[]
	hap2_seq=[]
	hap2_position=[]
	hg19_seq=[]
	hg19_position=[]
	compare={}
	seq_clone = ""
	seq_hg19 = ""
	seq_hap1=""
	seq_hap2=""
	f = open(name+"_hg19.aln","r")
	line_number=0
	while True:
		line = f.readline()
		if line=="":
			break
		if line[0]=="#":
			continue
		line = line.strip()
		if line=="":
			line_number +=1
		else:
			if line[:2]=="AC":
				seq_clone +=line.split(" ")[1]
				line = f.readline()
				line = f.readline().strip()
				seq_hg19 +=line.split(" ")[1]
#				print len(seq_clone),len(seq_hg19),line_number
	hg19_position,clone_hg19_seq,hg19_seq,diff1,length1=compare_seq(seq_clone,seq_hg19,hg19_position,clone_hg19_seq,hg19_seq)
	print name,"hg19",len(seq_clone),len(seq_hg19)
	f = open(name+"_hap1.v2.aln","r")
	line_number=0
	seq_clone = ""
	while True:
		line = f.readline()
		if line=="":
			break
		if line[0]=="#":
			continue
		line = line.strip()
		if line=="":
			line_number +=1
		else:
			if line[:2]=="AC":
				seq_clone +=line.split(" ")[1]
				line = f.readline()
				line = f.readline().strip()
				seq_hap1 +=line.split(" ")[1]
	hap1_position,clone_hap1_seq,hap1_seq,diff2,length2=compare_seq(seq_clone,seq_hap1,hap1_position,clone_hap1_seq,hap1_seq)
	print name,"hap1",len(seq_clone),len(seq_hap1)
	f = open(name+"_hap2.v2.aln","r")
	seq_clone = ""
	while True:
		line = f.readline()
		if line=="":
			break
		if line[0]=="#":
			continue
		line = line.strip()
		if line=="":
			line_number +=1
		else:
			if line[:2]=="AC":
				seq_clone +=line.split(" ")[1]
				line = f.readline()
				line = f.readline().strip()
				seq_hap2 +=line.split(" ")[1]
	hap2_position,clone_hap2_seq,hap2_seq=compare_seq(seq_clone,seq_hap2,hap2_position,clone_hap2_seq,hap2_seq)
	print name,"hap2",len(seq_clone),len(seq_hap2)
	compare["hap1"]=pd.Series(hap1_seq,index=hap1_position)
	compare["hap2"]=pd.Series(hap2_seq,index=hap2_position)
	compare["hg19"]=pd.Series(hg19_seq,index=hg19_position)
	compare["clone_hap1"]=pd.Series(clone_hap1_seq,index=hap1_position)
	compare["clone_hap2"]=pd.Series(clone_hap2_seq,index=hap2_position)
	compare["clone_hg19"]=pd.Series(clone_hg19_seq,index=hg19_position)
	Clone[name] = pd.DataFrame(compare)
	hap1_mismatch = 0
	hap2_mismatch = 0
	both_mismatch = 0
	for i in Clone[name].index:
		hap1_wrong = False
		hap1_N = False
		hap2_wrong = False
		hap2_N = False
		hg19_N = False
		hg19_wrong = False
#		print "%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
		if Clone[name]["hap1"][i]== Clone[name]["hap1"][i]:
			if	Clone[name]["hap1"][i]!="N":
				hap1_wrong = True
			else:
				hap1_N = True
		if Clone[name]["hap2"][i]== Clone[name]["hap2"][i]:
			if Clone[name]["hap2"][i]!="N":
				hap2_wrong = True
			else:
				hap2_N = True
		if Clone[name]["hg19"][i]== Clone[name]["hg19"][i]:
			if Clone[name]["hg19"][i]!="N":
				hg19_wrong = True
			else:
				hg19_N = True
		if hap1_wrong is True :
			if hap2_N is True:
				print "hap2 missing,hap1 right","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
			elif hap2_wrong is True:
				both_mismatch +=1
				print "both wrong","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
			else:
				hap1_mismatch +=1
				print "hap1 wrong","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
		else:
			if hap1_N is True:
				if hap2_wrong is True:
					print "hap1 missing,hap2 wrong","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
				elif hap2_N is True:
					if hg19_N is True:
						print "hg19 also uncallable","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
					else:
						print "uncallable","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
				else:
					print "hap1 missing,hap2 right","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
			else:
				if hap2_wrong is True:
					hap2_mismatch +=1
					print "hap2 wrong","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
				elif hap2_N is True:
					print "hap2 missing,hap1 right","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
				else:
					print "fine","%s\t%s\t%s\t%s\t%s\t%s\t%s" %(i,Clone[name]["clone_hap1"][i],Clone[name]["hap1"][i],Clone[name]["clone_hap2"][i],Clone[name]["hap2"][i],Clone[name]["clone_hg19"][i],Clone[name]["hg19"][i])
	print hap1_mismatch,hap2_mismatch,both_mismatch
	return hap1_mismatch,hap2_mismatch,both_mismatch

def compare_seq_v2(seq_clone,seq_hap1,seq_hap2,seq_hg19,position,clone_seq):
	j = 0
	diff1=0
	diff2=0
	length = 0
#	print "in",len(seq_clone),len(seq_hap)
	for i in range(len(seq_hg19)):
		if seq_hg19[i]!="-":
			j +=1
			if seq_hap1[j-1]!=0 and seq_clone[i]!="-":
				length +=1
			else:
				continue
			if seq_hap1[j-1]=="R" and seq_clone[i]!=seq_hg19[i]:
				clone_seq.append(seq_clone[i]+seq_hg19[i]+seq_hg19[i]+seq_hg19[i])
				position.append(j)
			if seq_hap1[j-1]!=0 and seq_hap1[j-1]!="R" and (seq_clone[i]!=seq_hap1[j-1] or seq_clone[i]!=seq_hap2[j-1]):
				clone_seq.append(seq_clone[i]+seq_hap1[j-1]+seq_hap2[j-1]+seq_hg19[i])
				position.append(j)
	return position,clone_seq,length

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
		if pos<list[midpoint]:
			midpoint=midpoint-1
	return midpoint

def get_seq(vcf_file,mask_file,pos1,pos2,sign):
	seq_hap1=[0 for i in range(pos2-pos1)]
	seq_hap2=[0 for i in range(pos2-pos1)]
	for i in VcfIterator(vcf_file):
		if i[1]>pos1 and i[1]<=pos2:
			if i[4] or i[3][0]==i[3][1]:
				seq_hap1[i[1]-pos1-1]=i[2][i[3][0]]
				seq_hap2[i[1]-pos1-1]=i[2][i[3][1]]
		elif i[1]>pos2:
			break
	left_pos = []
	right_pos = []
	for i in BedIterator(mask_file):
		left_pos.append(int(i[1])+1)
		right_pos.append(int(i[2]))
	find=find_interval(pos1+1,left_pos)
	try:
		assert pos1+1>=left_pos[find]
	except AssertionError:
		print pos1+1,find,left_pos[find]
	for i in range(pos1+1,pos2+1):
		while find<len(left_pos):
			if i<left_pos[find]:
				break
			if i<=right_pos[find]:
				if seq_hap1[i-pos1-1]==0:
					seq_hap1[i-pos1-1]="R"			# means reference allele
					seq_hap2[i-pos1-1]="R"
				break
			find=find+1
	print 'count:',seq_hap1.count(0),seq_hap1.count("R"),len(seq_hap1)
	reverse = {'A':'T','C':'G','T':'A','G':'C','R':'R',0:0}
	pos_list = range(pos2-pos1)
	pos_list.reverse()
	if sign =='-':
		seq_hap1=[reverse[seq_hap1[j]] for j in pos_list]
		seq_hap2=[reverse[seq_hap2[j]] for j in pos_list]
	return seq_hap1,seq_hap2

def compare_align_v2(name,start,end,sign,seq_hap1,seq_hap2,chr):
	position=[]
	clone_seq=[]
	compare={}
	seq_clone = ""
	seq_hg19 = ""
	f = open(name+".aln","r")
	line_number = 0
	while True:
		line = f.readline()
		if line=="":
			break
		if line[0]=="#":
			continue
		line = line.strip()
		if line=="":
			line_number +=1
		else:
			if line[:2]=="AC":
				seq_clone +=line.split(" ")[1]
				line = f.readline()
				line = f.readline().strip()
				seq_hg19 +=line.split(" ")[1]
	print name,"hg19",len(seq_clone),len(seq_hg19),len(seq_hap1),len(seq_hap2)
	position,clone_seq,length=compare_seq_v2(seq_clone,seq_hap1,seq_hap2,seq_hg19,position,clone_seq)
	hap1_mismatch = 0
	hap2_mismatch = 0
	both_mismatch = 0
	snp = 0
	for i in range(len(position)):
		hap1_wrong = False
		hap2_wrong = False
		if clone_seq[i][0]!=clone_seq[i][1]:		#hap1 is different
			hap1_wrong = True
		if clone_seq[i][0]!=clone_seq[i][2]:		#hap2 is different
			hap2_wrong = True
		if clone_seq[i][1]!=clone_seq[i][2]:
			snp +=1
		if hap1_wrong is True :
			if hap2_wrong is True:
				both_mismatch +=1
				print "both wrong","%s\t%s\t%s\t%s\t%s" %(position[i],clone_seq[i][0],clone_seq[i][1],clone_seq[i][2],clone_seq[i][3])
			else:
				hap1_mismatch +=1
				print "hap1 wrong","%s\t%s\t%s\t%s\t%s" %(position[i],clone_seq[i][0],clone_seq[i][1],clone_seq[i][2],clone_seq[i][3])
		else:
			if hap2_wrong is True:
				hap2_mismatch +=1
				print "hap2 wrong","%s\t%s\t%s\t%s\t%s" %(position[i],clone_seq[i][0],clone_seq[i][1],clone_seq[i][2],clone_seq[i][3])	
	print snp,hap1_mismatch,hap2_mismatch,both_mismatch,length
	return snp,hap1_mismatch,hap2_mismatch,both_mismatch,length

def read_phased_snps(file):
	f = open(file,"r")
	SNP = {}
	chr_list=[]
	for i in range(1,23):
		SNP["chr"+str(i)]=[]
		chr_list.append("chr"+str(i))
	for line in f:
		line = line.strip().split("\t")
		chr = line[0]
		if chr in chr_list:
			SNP[line[0]].append(int(line[1]))
	return SNP

def find_pos_v2(pos,list):
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
	return find

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='write msHOT simulation commands')
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	args = parser.parse_args()

	f=open("clone_hg19.bed","r")
	f_out = open("clone_assignment_final_new","w")
	tot_error = 0
	het_error =0.
	het_snp = 0.
	total_error = 0.
	total_length = 0.
	clone = []
	for line in f:
		line = line.strip().split("\t")
		name = line[3]
		start = int(line[1])
		end = int(line[2])
		sign = line[5]
		chr = line[0]
		sample ='NA19240'
		vcf_file = '%s%s/gVCF_calls/%s.%s.fosmid.v2.phased.vcf.gz' %(args.wgs_dir,sample,sample,chr)
		mask_file = '%s%s/gVCF_calls/%s_%s.mask.bed.gz' %(args.wgs_dir,sample,sample,chr)
		print name,"begin",sign
		seq_hap1,seq_hap2=get_seq(vcf_file,mask_file,start,end,sign)
		print "".join(map(str,seq_hap1))
		print "".join(map(str,seq_hap2))
		snp,hap1_mismatch,hap2_mismatch,both_mismatch,length=compare_align_v2(name,start,end,sign,seq_hap1,seq_hap2,chr)
		total_length+=length
		if hap1_mismatch > hap2_mismatch:
			hap="hap2"
			error = float(hap2_mismatch)/snp
			het_snp +=snp
			het_error +=hap2_mismatch
			total_error +=both_mismatch+hap2_mismatch
		elif hap1_mismatch<hap2_mismatch:
			hap="hap1"
			error = float(hap1_mismatch)/snp
			het_error +=hap1_mismatch
			het_snp +=snp
			total_error +=both_mismatch+hap1_mismatch
		else:
			hap="unknown"
			error = float(hap1_mismatch)/snp
			het_error +=hap1_mismatch
			het_snp +=snp
			total_error +=both_mismatch+hap1_mismatch
		print >>f_out,"%s\t%s\t%s\t%s\t%s\t%i\t%i\t%i\t%i\t%f\t%s" %(name.replace("_hg19",""),line[0][3:],line[1],line[2],line[5],snp,min(hap1_mismatch,hap2_mismatch),both_mismatch,length,error,hap)
	print >>f_out,"Total het error,", het_error,"Total het snp",het_snp,"average,",het_error/het_snp
	print >>f_out,"Total error,", total_error,"Total length",total_length,"average,",total_error/total_length

	
				
				
