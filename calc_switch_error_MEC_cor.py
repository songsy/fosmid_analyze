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
import sys


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

def calc_switch_error_MEC(switch_error_file,BLOCK,block_list,f_out):
	for i in BedIterator(switch_error_file):
		if i[0]!=CHROM:
			continue
		block_name=i[0]+" "+i[1]
		pos_list=i[3].split(",")
		mec_vec=np.zeros(len(pos_list))
		count_vec=np.zeros(len(pos_list))
		for clone in BLOCK[block_name]:
			if clone not in ["hap","transmit","untransmit","switch"] and clone[-5:]!="_qual":
				hap,mec_vec,count_vec=determine_clone_hap(BLOCK[block_name],clone,pos_list,mec_vec,count_vec)
		print >>f_out,"%s\t\t%s\t%s\t%s\t%s" %("\t".join(i),np.sum(mec_vec),np.mean(mec_vec),np.mean(count_vec),np.mean([float(a)/b for a,b in zip(mec_vec,count_vec)]))
	
if __name__=="__main__":
	parser = argparse.ArgumentParser(description='assign parental allele in block')
#	parser.add_argument("--chr", dest='chr',help="chromosome")
	parser.add_argument("--sample", dest='sample',help="sample name")
#	parser.add_argument("--file",dest='file',help="switch error file name")

	args = parser.parse_args()
	SAMPLE = args.sample
	f_out1 = open(args.sample+"_shapeit_switch_error_mec.txt","w")
	f_out2 = open(args.sample+"_1KGphase3_switch_error_mec.txt","w")
	switch_error_file1=args.sample+"_shapeit_switch_error.txt"
	switch_error_file2=args.sample+"_1KGphase3_switch_error.txt"
	for i in range(1,23):
		CHROM = 'chr'+str(i)
		dbfile=open('../BLOCK/%s_gvcf/%s_BLOCK_' %(SAMPLE,SAMPLE) +CHROM+'_pickle','rb')
		BLOCK=pickle.load(dbfile) 
		block_list=pickle.load(dbfile)
		mec=pickle.load(dbfile)
		if os.path.isfile(switch_error_file1):
			calc_switch_error_MEC(switch_error_file1,BLOCK,block_list,f_out1)
		if os.path.isfile(switch_error_file2):
			calc_switch_error_MEC(switch_error_file2,BLOCK,block_list,f_out2)
	
			

