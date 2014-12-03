#!/usr/bin/env python
# python pool_SNP_info_gvcf.py
# Shiya Song
# 21 July 2014
# Pull information of SNP from each pool, make input file for ReFhap, run it

import sys
import genutils
import glob
import re
import argparse
from write_poststep import *
import subprocess
import commands

chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
chromLens = genutils.read_chrom_len(chromLenFile)
chromOrder = {}
inFile = open(chromLenFile,'r')
for i,l in enumerate(inFile):
	l = l.rstrip()
	c = l.split()[0]
	chromOrder[c] = i
inFile.close()

# Helper function to run commands, handle return values and print to log file
def runCMD(cmd):
	val = subprocess.Popen(cmd, shell=True).wait()
	if val == 0:
		cmd = cmd
	else:
		print 'command failed:',cmd

def read_whole_genome_vcf(vcf_file,sample):
	m=['A','C','G','T']
	f_vcf = open(vcf_file,'r')
	f_var = open(sample+"_combined_genome_allvars_gvcf_SNP",'w')
#	f_var2 = open(sample+"_combined_genome_allvars_True",'w')
	whole_SNP = {}
	i = 0
	GQ = []
	for line in f_vcf:
		line = line.strip()
		if line[0] == "#":
			continue
		line1 = line.split('\t')
		chr = line1[0]
		if chr=='chrY' or chr=='EBV' or chr=='pCC1FOS' or chr=='eColiK12':
			continue
		ref = line1[3]
		alt = line1[4].split(',')
		pos = line1[1]
		info = line1[9]
		info = info.split(':')
# only consider heterzygous SNPs
		tag = 'SNP'
		if info[0] == "1/2":
			if len(ref)!=1 or len(alt[0])!=1 or len(alt[1])!=1:
				continue
				tag = 'INDEL'
				if len(alt[0])!=1 and len(alt[1])!=1:
					alt_allele1 = m.index(alt[0][0])
					if alt_allele1 <3:
						alt_allele2 = m[alt_allele1+1]
					else:
						alt_allele2 = m[0]
				elif len(alt[0])!=1:
					alt_allele2 = alt[1]
					alt_allele1 = m.index(alt_allele2)
					if alt_allele1 <3:
						alt_allele1 = m[alt_allele1+1]
					else:
						alt_allele1 = m[0]
				elif len(alt[1])!=1:
					alt_allele1 = alt[0]
					alt_allele2 = m.index(alt_allele1)
					if alt_allele2 <3:
						alt_allele2 = m[alt_allele2+1]
					else:
						alt_allele2 = m[0]
				print >>f_var,"%s\t%s\t%s\t%s\t%s\tSNP" %(chr,pos,pos,alt_allele1,alt_allele2)
#				print >>f_var2,"%s\t%s\t%s\t%s\t%s\t%s" %(chr,pos,pos,alt[0],alt[1],tag)
			else:
				print >>f_var,"%s\t%s\t%s\t%s\t%s\tSNP" %(chr,pos,pos,alt[0],alt[1])
#				print >>f_var2,"%s\t%s\t%s\t%s\t%s\t%s" %(chr,pos,pos,alt[0],alt[1],tag)
			i += 1
			whole_SNP[chr + " " + pos] = [i,chr,int(pos),info[0],alt[0],alt[1]]
			if len(info)==5:
				GQ.append(int(info[3]))
			elif len(info)==3:
				GQ.append(int(info[1]))
			else:
				print line
		elif info[0] == "0/1":
			if len(ref)!=1 or len(alt[0])!=1:
				continue
				tag = 'INDEL'
				alt_allele = m.index(ref[0])
				if alt_allele <3:
					alt_allele = m[alt_allele+1]
				else:
					alt_allele = m[0]
				print >>f_var,"%s\t%s\t%s\t%s\t%s\tSNP" %(chr,pos,pos,ref[0],alt_allele)
#				print >>f_var2,"%s\t%s\t%s\t%s\t%s\t%s" %(chr,pos,pos,ref,alt[0],tag)
			else:
				print >>f_var,"%s\t%s\t%s\t%s\t%s\tSNP" %(chr,pos,pos,ref[0],alt[0])
#				print >>f_var2,"%s\t%s\t%s\t%s\t%s\t%s" %(chr,pos,pos,ref,alt[0],tag)
			i += 1
			whole_SNP[chr + " " + pos] = [i,chr,int(pos),info[0],ref,alt]
			if len(info)==6:
				GQ.append(int(info[3]))
			else:
				print line			
		else:
			whole_SNP[chr + " " + pos] = 0
	return whole_SNP,GQ,i

def filter_vcf(vcf_file,clone_file,chromOrder,whole_SNP,index):
	f_vcf = open(vcf_file,'r')
	f_clone = open(clone_file,'r')
	tot_SNP = 0
	clone_SNP = 0
	heter_SNP = 0
	homo_SNP = 0
	heter_false_SNP = 0 
	clone_invisible_SNP = 0
	first = True
	clone_out = open(clone_file+"_gvcf_SNP_all",'w')
	matrix_out = open(clone_file + "_gvcf_SNP_matrix",'w')
	print >>clone_out,"chrom\tpos1\tpos2\tcov\tsize\tleft_support\tright_support\theter_SNP\thomo_SNP"
	for line in f_vcf:
		line = line.strip()
		if line[0] == "#":
			continue
		line1 = line.split('\t')
		chr = line1[0]
		if chr=='chrY' or chr=='EBV' or chr=='pCC1FOS':
			continue
		ref = line1[3]
		alt = line1[4]
		pos = int(line1[1])
		info = line1[9]
		lowqual = line1[6]
		if info == "./.":
			continue
		tot_SNP +=1				   # number of SNPs in the list found in vcf file of a pool
		info = info.split(':')
		GT = info[0]
		if len(info)==5:
			GQ=int(info[3])
		elif len(info)==3:
			GQ=int(info[1])
		if 33+GQ >= 126:
			GQ = 126
		else:
			GQ = 33 + GQ
		found = False
		try:
			abc = whole_SNP[chr + " " + str(pos)]
		except KeyError:
			continue
		if first is True:
			clone = f_clone.readline().strip()
			clone_first = True
			bad = False
			clone1 = clone.split('\t')
			chrom = clone1[0]
			pos1 = int(clone1[1])
			pos2 = int(clone1[2])
			last_pos = pos1
			cov = float(clone1[3])
			first = False
			clone_heter_SNP = 0
			clone_homo_SNP = 0
			item = []
		if chr == chrom :
			if pos < pos1:
				continue
			elif pos <= pos2:
				found = True
				clone_SNP +=1						# number of SNPs found in a clone that we identified
				if GT == "0/1" or GT == "1/2" or GT == "0/2":
					heter_SNP +=1
					clone_heter_SNP +=1
#					Heter_SNP[chr + " " + str(pos)].append([index,chrom,pos1,pos2])
					bad = True
				elif whole_SNP[chr + " " + str(pos)]!= 0:
					homo_SNP +=1
					clone_homo_SNP +=1
					if clone_first == True:
						item = []
						clone_name = str(index) + "_" + chr + "_" + str(pos1)
						line_number = whole_SNP[chr + " " + str(pos)][0]
						clone_first = False
						seqstring = ""
						qualstring = ""
						item.append(clone_name)
						item.append(line_number)
						qualstring += str(unichr(GQ))
						if GT == "1/1":
							if whole_SNP[chr + " " + str(pos)][3] == "1/2":
								seqstring += "0"
							elif whole_SNP[chr + " " + str(pos)][3] == "0/1":
								seqstring += "1"
							else:
								print line
						elif GT =="0/0":
							seqstring += "0"
							clone_invisible_SNP +=1
						elif GT == "2/2":
							if whole_SNP[chr + " " + str(pos)][3] == "1/2":
								seqstring += "1"
							else:
								print line
						else:
							print line
					else:
						if whole_SNP[chr + " " + str(pos)][0] == line_number + 1 :
							qualstring += str(unichr(GQ))
							line_number +=1
							if GT == "1/1":
								if whole_SNP[chr + " " + str(pos)][3] == "1/2":
									seqstring += "0"
								elif whole_SNP[chr + " " + str(pos)][3] == "0/1":
									seqstring += "1"
								else:
									print line
							elif GT =="0/0":
								seqstring += "0"
								clone_invisible_SNP +=1
							elif GT == "2/2":
								if whole_SNP[chr + " " + str(pos)][3] == "1/2":
									seqstring += "1"
								else:
									print line
							else:
								print line
						else:
							item.append(seqstring)
							seqstring = ""
							line_number = whole_SNP[chr + " " + str(pos)][0]
							item.append(line_number)
							qualstring += str(unichr(GQ))
							if GT == "1/1":
								if whole_SNP[chr + " " + str(pos)][3] == "1/2":
									seqstring += "0"
								elif whole_SNP[chr + " " + str(pos)][3] == "0/1":
									seqstring += "1"
								else:
									print line
							elif GT =="0/0":
								seqstring += "0"
								clone_invisible_SNP +=1
							elif GT == "2/2":
								if whole_SNP[chr + " " + str(pos)][3] == "1/2":
									seqstring += "1"
								else:
									print line
							else:
								print line
				continue
			else:
				while True:
					if clone!='':
						if clone_homo_SNP != 0:
							item.append(seqstring)
							item.append(qualstring)
							string = ""
							string += str(len(item)/2 - 1)
							for i in range(len(item)):
								string += " " + str(item[i])
							item = []
							if bad == False:
								print >>matrix_out, string
						print >>clone_out, clone,'\t',clone_heter_SNP,'\t',clone_homo_SNP
					clone = f_clone.readline().strip()
					clone_first = True
					bad = False
					clone_heter_SNP = 0
					clone_homo_SNP = 0
					clone_split = 0
					if clone == '':
						break
					clone1 = clone.split('\t')
					chrom = clone1[0]
					pos1 = int(clone1[1])
					pos2 = int(clone1[2])
					last_pos = pos1
					cov = float(clone1[3])
					if chrom != chr:
#						print >>f,line
						break
					else:
						if pos < pos1:
#							print >>f,line
							break
						elif pos <= pos2:
							found = True
							clone_SNP +=1
							if GT == "0/1" or GT == "1/2" or GT == "0/2":
								heter_SNP +=1
								clone_heter_SNP +=1
#								Heter_SNP[chr + " " + str(pos)].append([index,chrom,pos1,pos2])
								bad = True
							elif whole_SNP[chr + " " + str(pos)]!= 0:
								homo_SNP +=1
								clone_homo_SNP +=1
								if clone_first == True:
									item = []
									clone_name = str(index) + "_" + chr + "_" + str(pos1)
									line_number = whole_SNP[chr + " " + str(pos)][0]
									clone_first = False
									seqstring = ""
									qualstring = ""
									item.append(clone_name)
									item.append(line_number)
									qualstring += str(unichr(GQ))
									if GT == "1/1":
										if whole_SNP[chr + " " + str(pos)][3] == "1/2":
											seqstring += "0"
										elif whole_SNP[chr + " " + str(pos)][3] == "0/1":
											seqstring += "1"
										else:
											print line
									elif GT =="0/0":
										seqstring += "0"
										clone_invisible_SNP += 1
									elif GT == "2/2":
										if whole_SNP[chr + " " + str(pos)][3] == "1/2":
											seqstring += "1"
										else:
											print line
									else:
										print line
								else:
									if whole_SNP[chr + " " + str(pos)][0] == line_number + 1 :
										qualstring += str(unichr(GQ))
										line_number +=1
										if GT == "1/1":
											if whole_SNP[chr + " " + str(pos)][3] == "1/2":
												seqstring += "0"
											elif whole_SNP[chr + " " + str(pos)][3] == "0/1":
												seqstring += "1"
											else:
												print line
										elif GT =="0/0":
											seqstring += "0"
											clone_invisible_SNP +=1
										elif GT == "2/2":
											if whole_SNP[chr + " " + str(pos)][3] == "1/2":
												seqstring += "1"
											else:
												print line
										else:
											print line
									else:
										item.append(seqstring)
										seqstring = ""
										line_number = whole_SNP[chr + " " + str(pos)][0]
										item.append(line_number)
										qualstring += str(unichr(GQ))
										if GT == "1/1":
											if whole_SNP[chr + " " + str(pos)][3] == "1/2":
												seqstring += "0"
											elif whole_SNP[chr + " " + str(pos)][3] == "0/1":
												seqstring += "1"
											else:
												print line
										elif GT =="0/0":
											seqstring += "0"
											clone_invisible_SNP +=1
										elif GT == "2/2":
											if whole_SNP[chr + " " + str(pos)][3] == "1/2":
												seqstring += "1"
											else:
												print line
										else:
											print line						
							break
		elif chromOrder[chr] < chromOrder[chrom]:
			continue
		else:
			while True:
					if clone!='':
						if clone_homo_SNP != 0:
							item.append(seqstring)
							item.append(qualstring)
							string = ""
							string += str(len(item)/2 - 1)
							for i in range(len(item)):
								string += " " + str(item[i])
							item = []
							if bad == False:
								print >>matrix_out, string
						print >>clone_out, clone,'\t',clone_heter_SNP,'\t',clone_homo_SNP,'\t',clone_split
					clone = f_clone.readline().strip()
					clone_first = True
					bad = False
					clone_heter_SNP = 0
					clone_homo_SNP = 0
					clone_split = 0
					if clone == '':
#						print >>f,line
						break
					clone1 = clone.split('\t')
					chrom = clone1[0]
					pos1 = int(clone1[1])
					pos2 = int(clone1[2])
					last_pos = pos1
					cov = float(clone1[3])
					if chrom != chr:
						continue
					else:
						if pos < pos1:
#							print >>f,line
							break
						elif pos <= pos2:
							found = True
							clone_SNP +=1
							if GT == "0/1" or GT == "1/2" or GT == "0/2":
								heter_SNP +=1
								clone_heter_SNP +=1
#								Heter_SNP[chr + " " + str(pos)].append([index,chrom,pos1,pos2])
								bad = True
							elif whole_SNP[chr + " " + str(pos)]!= 0:
								homo_SNP +=1
								clone_homo_SNP +=1
								if clone_first == True:
									item = []
									clone_name = str(index) + "_" + chr + "_" + str(pos1)
									line_number = whole_SNP[chr + " " + str(pos)][0]
									clone_first = False
									seqstring = ""
									qualstring = ""
									item.append(clone_name)
									item.append(line_number)
									qualstring += str(unichr(GQ))
									if GT == "1/1":
										if whole_SNP[chr + " " + str(pos)][3] == "1/2":
											seqstring += "0"
										elif whole_SNP[chr + " " + str(pos)][3] == "0/1":
											seqstring += "1"
										else:
											print line
									elif GT =="0/0":
										seqstring += "0"
										clone_invisible_SNP +=1
									elif GT == "2/2":
										if whole_SNP[chr + " " + str(pos)][3] == "1/2":
											seqstring += "1"
										else:
											print line
									else:
										print line
								else:
									if whole_SNP[chr + " " + str(pos)][0] == line_number + 1 :
										qualstring += str(unichr(GQ))
										line_number +=1
										if GT == "1/1":
											if whole_SNP[chr + " " + str(pos)][3] == "1/2":
												seqstring += "0"
											elif whole_SNP[chr + " " + str(pos)][3] == "0/1":
												seqstring += "1"
											else:
												print line
										elif GT =="0/0":
											seqstring += "0"
											clone_invisible_SNP +=1
										elif GT == "2/2":
											if whole_SNP[chr + " " + str(pos)][3] == "1/2":
												seqstring += "1"
											else:
												print line
										else:
											print line
									else:
										item.append(seqstring)
										seqstring = ""
										line_number = whole_SNP[chr + " " + str(pos)][0]
										item.append(line_number)
										qualstring += str(unichr(GQ))
										if GT == "1/1":
											if whole_SNP[chr + " " + str(pos)][3] == "1/2":
												seqstring += "0"
											elif whole_SNP[chr + " " + str(pos)][3] == "0/1":
												seqstring += "1"
											else:
												print line
										elif GT =="0/0":
											seqstring += "0"
											clone_invisible_SNP +=1
										elif GT == "2/2":
											if whole_SNP[chr + " " + str(pos)][3] == "1/2":
												seqstring += "1"
											else:
												print line
										else:
											print line
							break
	return tot_SNP,clone_SNP,heter_SNP,homo_SNP

def refhap_decode(refhap_file,whole_SNP,chr):
	f_refhap = open(refhap_file,"r")
	f_out = open(refhap_file+"_detail","w")
	for l in f_refhap:
		l = l.strip()
		if l[0] == "B" or l[0] =="*":
			print >>f_out, l
		else:
			l = l.split("\t")
			pos = l[0]
			m = whole_SNP[chr + " "+pos]
			print >>f_out,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(m[0],l[1],l[2],m[1],pos,m[3],m[4],m[3]+"/"+m[4])
	f_out.close()

def read_allvars(allvars_file):
	f = open(allvars_file,'r')
	whole_SNP = {}
	i = 0
	GQ = []
	for line in f:
		line = line.strip()
		line1 = line.split('\t')
		chr = line1[0]
		ref = line1[3]
		alt = line1[4]
		pos = line1[1]
		i +=1
		whole_SNP[chr + " " + pos] = [i,chr,pos,ref,alt]
	return whole_SNP

def read_haplotype(hap_file,SNP_num):
	f = open(hap_file,"r")
	block = {}
	hap = {}
	decode = {}
	mec = {}
	resolve = []
	for i in range(1,SNP_num+2):
		resolve.append(0)				
	for i in range(1,SNP_num+2):
		decode[i] = 0
	start = True
	for line in f:
		line = line.strip()
		if line[0]=="B" or line[0]=="*":
			start = True
		else:
			col = line.split("\t")
			if col[1] == "-":
				continue
			if col[1] == "1" or col[1] == "0":
				resolve[int(col[0])] = 1
			if start == True:
				start = False
				chr = col[3]
				pos = col[4]
				block[chr+" "+pos] = []
				hap[chr+" "+pos] = ""
				mec[chr+" "+pos] = 0
			block[chr+" "+pos].append(int(col[0]))
			hap[chr+" "+pos] += col[1]
			decode[int(col[0])] = chr +" "+ pos
	f.close()
	return block,hap,decode,mec,resolve

def read_matrix(matrix_file,block,hap,decode,mec,pool_dir,sample,SNP_num):
	f_mat = open(matrix_file,"r")
	f1_out = open(sample+"_refhap_gvcf_hap1.bed","w")
	f2_out = open(sample+"_refhap_gvcf_hap2.bed","w")
	print >>f1_out,"track name='%s_hap1' description='Hap1' visibility=2 itemRgb='On'" %(sample)
	print >>f2_out,"track name='%s_hap2' description='Hap2' visibility=2 itemRgb='On'" %(sample)
	MEC = 0
	done = 0
	show = []
	for i in range(1,SNP_num+2):
		show.append(0)
	for line in f_mat:
		done += 1
		line = line.rstrip()
		col = line.split(" ")
		name = col[1]
		if sample=='NA19240':
			chr = ("_").join(col[1].split("_")[3:-1])
			pool = ("_").join(col[1].split("_")[:3])
		else:
			chr = col[1].split("_")[-2]
			pool = ("_").join(col[1].split("_")[:-2])
		pos = col[1].split("_")[-1]
		position = []
		allele = ""
		match = 0
		not_found = 0
		hap1 = False
		assign = False
		for i in range(int(col[0])):
			try:
				num = int(col[2+2*i])
				string = col[3+2*i]
				for j in range(len(string)):
					position.append(num+j)
					allele += string[j]
					show[num+j] = 1
			except IndexError:
				print 'Wrong:',line,i,col,len(col)
		first = True
		for i in range(len(position)):
			if decode[position[i]] == 0:
				not_found +=1
				continue
			assign = True
			if first == True:
				code = decode[position[i]]
				first = False
				index = block[code].index(position[i])
				if hap[code][index] == allele[i]:
					match += 1
			elif code == decode[position[i]]:
				index = block[code].index(position[i])
				if hap[code][index] == allele[i]:
					match += 1
			else:
				print "wrong",line
		if done % 10000 == 0:
			print done,"done"
		if assign == False:
			continue
		if match >= (len(position)-not_found)/2.:
			hap1 = True
			mismatch = len(position) - match - not_found
		else:
			mismatch = match
		MEC += mismatch
		mec[code] += mismatch
		find = False
		clone_name = pool_dir+'%s/' %(sample) +pool + '/' + pool + '.markdup.clone.sel.10000.0.25'
		f_clone = open(clone_name,"r")
		for l in f_clone:
			l = l.strip()
			chrom = l.split("\t")[0]
			pos1 = l.split("\t")[1]
			pos2 = l.split("\t")[2]
#			if chrom == chr and int(pos) >= int(pos1) and int(pos) <= int(pos2):
			if chrom == chr and int(pos) == int(pos1):
				find = True
				break
		f_clone.close()
		if hap1 == True:
			print >>f1_out, "%s\t%s\t%s\t%s\t%i\t+\t%s\t%s\t255,0,0" %(chr,pos1,pos2,name,mismatch,pos1,pos2)
		else:
			print >>f2_out, "%s\t%s\t%s\t%s\t%i\t+\t%s\t%s\t0,0,255" %(chr,pos1,pos2,name,mismatch,pos1,pos2)
	return MEC,mec,show

def make_track(hap_file,mec,sample):
	block_size = []
	f = open(hap_file,"r")
	f_out = open(sample+"_refhap_gvcf_track_final.txt","w")
	print >>f_out,"track type=bedGraph name='refhap_hap_block' description='refhap' visibility=2 priority=priority"
	for line in f:
		line = line.strip()
		if line[0]=="B":
			start = 0
			end = 0
			chr = ""
			continue
		elif line[0] == "*":
			print >>f_out,"%s\t%s\t%s\t%s" %(chr,start,end,MEC)
			if int(end)-int(start)!=0:
				block_size.append(int(end)-int(start))
		elif start == 0:
			if line.split("\t")[1] == "-":
				continue
			chr = line.split("\t")[3]
			start = line.split("\t")[4]
			end = line.split("\t")[4]
			MEC = mec[chr + " " + start]
		elif start != "":
			end = line.split("\t")[4]
	return block_size

def calc_N50(num,sample):
	f_out = open(sample+"_gvcf_block_size.txt","w")
	num.sort()
	total = sum(num)
	print "tot:",total
	m = range(1,11)
	a = 0
	j = 0
	for i in num:
		a += i
		print >>f_out,"%i\t%i" %(i,a)
		if float(a) >= float(m[j]*0.1*total) and j<9:
			j +=1
			print "N%i0:%i" %(10-j,i)

def SNP_not_shown(show,resolve,file,sample,SNP_num):
	not_shown = 0
	not_resolved = 0
	for i in range(1,SNP_num+1):
		if show[i] == 0:
			not_shown +=1
		if resolve[i]==0 and show[i]==1:
			not_resolved +=1
	print "not_shown",not_shown
	print "not_resolved",not_resolved
	f= open(file,"r")
	f_out = open(sample+"_gvcf_SNPs_not_shown","w")
	for i,l in enumerate(f):
		if show[i+1]==0:
			print >>f_out,i+1,"\t",l.strip()
	print "phased %f snps" %(1-float(not_shown+not_resolved)/SNP_num)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='make input file for ReFhap')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	parser.add_argument("--pool_dir", dest='pool_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools/',help="directory for each pool data file")
	parser.add_argument("--out_dir", dest='out_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/',help="directory for output")
	args = parser.parse_args()

	wk = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/'
	if args.sample =='HG02799':
		pool_name = {}
		file = wk+'HG02799-fosmids.plate4.txt.info'
		pool_name=get_pool_name(file,pool_name)
		file = wk+'HG02799-fosmids.plates1_2.txt.info'
		pool_name=get_pool_name(file,pool_name)
	if args.sample =='HG03108':
		pool_name = {}
		file = wk+'HG03108-fosmids.plate2.txt.info'
		pool_name=get_pool_name(file,pool_name)
		file = wk+'HG03108-fosmids.plates3_4.txt.info'
		pool_name=get_pool_name(file,pool_name)
	if args.sample =='NA21302':
		pool_name = {}
		file = wk+'NA21302-fosmids.plate1.txt.info'
		pool_name=get_pool_name(file,pool_name)
		file = wk+'NA21302-fosmids.plate2.txt.info'
		pool_name=get_pool_name(file,pool_name)
		file = wk+'HG21302_p3.txt.info'
		pool_name=get_pool_name(file,pool_name)
	if args.sample =='NA19240':
		pool_name = {}
		interval=range(1,289)
		interval.remove(120)
		interval.remove(155)
		for i in interval:
			if i <=9:
				index="00"+str(i)
			elif i<=99:
				index="0"+str(i)
			else:
				index=str(i)
			pool_name['NA19240_pool_'+index]=1
	if args.sample == 'HG03428':
		pool_name = {}
		a = os.listdir('/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools/HG03428/')
		for i in a:
			if i[-10:]=='clonesizes':
				continue
			m = re.search('.+run759_match.+',i)
			if m is not None:
				continue
			if i[:9]=='HG03428L2' or i[:10]=='HG03428_p2' or i[:10]=='HG03428_p1' or i[:14]=='HG03428_plate1':
				pool_name[i]=1
			if i[:14]=='HG03428_plate2':
				if i[-8:]=='combined':
					pool_name[i]=1
				else:
					if i+'_combined' not in a:
						pool_name[i]=1
	if args.sample == 'NA20847':
		pool_name={}
		for i in range(1,116):
			pool_name['fraction'+str(i)]=1

	wgs_vcf = args.wgs_dir + args.sample+ '/'+ args.sample + '_het_snp.vcf'
	whole_SNP,GQ,SNP_num = read_whole_genome_vcf(wgs_vcf,args.sample)
	print "number of heterozygous snp to be phased:",SNP_num
	'''
	f = open(args.out_dir+args.sample+'_pool_gvcf_SNP_stats.txt','w')
	print >>f, "pool\ttot_SNP\tclone_SNP\theter_SNP\thomo_SNP"
	outDirBase= args.pool_dir+'%s/' %(args.sample)
	for pool in pool_name.keys():
		cloneVcf = args.pool_dir+'%s/' %(args.sample) +pool + '/' + pool + '_gvcf_all.vcf'
		cloneBed = args.pool_dir+'%s/' %(args.sample) +pool + '/' + pool + '.markdup.clone.sel.10000.0.25'
		
		tot_SNP,clone_SNP,heter_SNP,homo_SNP = filter_vcf(cloneVcf,cloneBed,chromOrder,whole_SNP,pool)
		print '%s done' %(pool)
		print >>f, "%s\t%i\t%i\t%i\t%i" %(pool,tot_SNP,clone_SNP,heter_SNP,homo_SNP)
	f.close()

	cmd= "cat "
	for pool in pool_name.keys():
		cloneBed = args.pool_dir+'%s/' %(args.sample) +pool + '/' + pool + '.markdup.clone.sel.10000.0.25' 
		cloneMatrix = cloneBed + "_gvcf_SNP_matrix "
		cmd += cloneMatrix
	cmd += "> %s_all_clone_all_fragment_gvcf_snp.matrix" %(args.sample)
	os.popen(cmd)


	f = open("%s_all_clone_all_fragment_gvcf_snp.matrix" %(args.sample),"r")
	window = []
	for line in f:
		line = line.strip()
		col = line.split(" ")
		if len(col)<2:
			print line
		if args.sample=='NA19240':
			window.append([("_").join(col[1].split("_")[3:-1]),int(col[2]),line])
		else:
			window.append([("_").join(col[1].split("_")[-2:-1]),int(col[2]),line])		
	
	# by window start
	window.sort( key=lambda w: w[1] )
	# by chrom name
	window.sort( key=lambda w: chromOrder[w[0]])
	print 'Matrix finished'

	f_out = open("%s_all_clone_all_fragment_gvcf_snp_sorted.matrix" %(args.sample),"w")
	last_chr = window[0][0]
	f_chr = open("Refhap/%s_" %(args.sample)+last_chr + ".gvcf.snp.matrix","w")
	f_cmds = open("Refhap/%s_gvcf_snp_refhap.cmds" %(args.sample),"w")
	cmds = 'java -cp /home/jmkidd/kidd-lab-scratch/shiya-projects/indian_Jacob/softwares/SingleIndividualHaplotyper/SIH.jar mpg.molgen.sih.main.SIH -a Refhap -v %s_combined_genome_allvars_gvcf_SNP -c 2 Refhap/%s_%s.gvcf.snp.matrix Refhap/%s_%s.gvcf.snp.haplotype' %(args.sample,args.sample,last_chr,args.sample,last_chr)
	print >>f_cmds,cmds
	for i in window:
		chr = i[0]
		print >>f_out,i[2]
		if chr==last_chr:
			print >>f_chr,i[2]
		else:
			f_chr.close()
			f_chr = open("Refhap/%s_" %(args.sample)+chr + ".gvcf.snp.matrix","w")
			runCMD(cmds)
			cmds = 'java -cp /home/jmkidd/kidd-lab-scratch/shiya-projects/indian_Jacob/softwares/SingleIndividualHaplotyper/SIH.jar mpg.molgen.sih.main.SIH -a Refhap -v %s_combined_genome_allvars_gvcf_SNP -c 2 Refhap/%s_%s.gvcf.snp.matrix Refhap/%s_%s.gvcf.snp.haplotype' %(args.sample,args.sample,chr,args.sample,chr)
			print >>f_cmds,cmds
			print >>f_chr,i[2]
			last_chr=chr
	runCMD(cmds)
	f_chr.close()
	print 'Refhap finished'
	'''
	whole_SNP = read_allvars(args.sample+"_combined_genome_allvars_gvcf_SNP")
	cmds = 'cat '
	for file in glob.glob("Refhap/%s*.gvcf.snp.haplotype" %(args.sample)):
		m=re.search('\w+(chr\w+).gvcf.snp.haplotype',file)
		chr = m.group(1)
		print chr
		if chr=='chrY':
			continue
		refhap_decode(file,whole_SNP,chr)
		cmds += file +"_detail "
	cmds += '> %s_gvcf_snp_haplotype_detail' %(args.sample)
	os.popen(cmds)


	haplotype_file = args.sample+'_gvcf_snp_haplotype_detail'
	block,hap,decode,mec,resolve=read_haplotype(haplotype_file,SNP_num)
	
	matrix_file = args.sample + '_all_clone_all_fragment_gvcf_snp_sorted.matrix'
	MEC,mec,show = read_matrix(matrix_file,block,hap,decode,mec,args.pool_dir,args.sample,SNP_num)
	print 'MEC:',MEC
	block_size=make_track(haplotype_file,mec,args.sample)
	calc_N50(block_size,args.sample)

	'''
	block_size=[]
	f=open(args.sample+"_refhap_gvcf_track_final.txt","r")
	for l in f:
		if l[0]=='t':
			continue
		l = l.strip().split()
		m = int(l[2])-int(l[1])
		if m!=0:
			block_size.append(m)
	calc_N50(block_size,args.sample)

	matrix_file = args.sample + '_all_clone_all_fragment_gvcf_snp_sorted.matrix'
	haplotype_file = args.sample+'_gvcf_snp_haplotype_detail'
	SNP_num = int(commands.getoutput('less %s_combined_genome_allvars_gvcf_SNP|wc -l' %(args.sample)))
	print SNP_num
	f_mat = open(matrix_file,"r")
	show = []
	for i in range(1,SNP_num+2):
		show.append(0)
	for line in f_mat:
		line = line.strip()
		col = line.split(" ")
		for i in range(int(col[0])):
			num = int(col[2+2*i])
			string = col[3+2*i]
			for j in range(len(string)):
				show[num+j] = 1
	f_mat.close()
	'''
	f = open(haplotype_file,"r")
	resolve = []
	for i in range(1,SNP_num+2):
		resolve.append(0)
	resolve.append(0)
	for line in f:
		line = line.strip()
		if line[0]=="B" or line[0]=="*":
			continue
		else:
			col = line.split("\t")
			if col[1] == "1" or col[1] == "0":
				resolve[int(col[0])] = 1
	f.close()
	
	SNP_not_shown(show,resolve,args.sample+"_combined_genome_allvars_gvcf_SNP",args.sample,SNP_num)