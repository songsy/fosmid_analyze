#!/usr/bin/env python
# Shiya Song
# 3rd July 2014
# fa_to_msmc.py

import signal
import os
import sys

def compare_seq(lines,line_number,hom_num):
	String = []
	for i in range(len(lines)-1):
		assert len(lines[i])==len(lines[i+1])
	for i in range(len(lines[0])):
		callable = True
		het = False
		for j in range(len(lines)):
			if lines[j][i]=="N":
				callable = False
		if callable is True:
			for j in range(len(lines)):
				for k in range(j+1,len(lines)):
					if lines[j][i]!=lines[k][i]:
						het = True
						break
				if het is True:
					break
			if het is True:
				pos = (line_number-1)*50+i+1
				geno = "".join([m[i] for m in lines])
				String.append(str(pos)+"\t"+str(hom_num+1)+"\t"+geno)
				hom_num = 0
			else:
				hom_num +=1
	return hom_num,String

def read_hap_fa_snp(file_list,name):
	f_handle=[]
	for filename in file_list:
		signal.signal(signal.SIGPIPE, signal.SIG_DFL)
		gc = 'gunzip -c ' + filename
		f_handle.append(os.popen(gc, 'r'))
	while True:
		lines = []
		for f in f_handle:
			lines.append(f.readline().strip())
		if lines[0]=="":
			break
		if lines[0][0]==">":
			chr = lines[0].replace(">chr","")
			f_out = open("msmc/YRI/"+chr+"_"+name+".txt","w")
			for i in range(len(lines)-1):
				assert lines[i]==lines[i+1]
			line_number = 0
			hom_num = 0
		else:
			line_number +=1
			hom_num,String= compare_seq(lines,line_number,hom_num)
			if len(String)!=0:
				for i in String:
					print >>f_out,"%s\t%s\t" %(chr,i)

if __name__=="__main__":	
	file_list = ["phased-genomes/NA19240_chrX_hap1.fa.gz","phased-genomes/NA19240_chrX_hap2.fa.gz"]
	read_hap_fa_snp(file_list,"NA19240_old_fosmid_phase")

