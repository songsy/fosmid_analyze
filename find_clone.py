#!/usr/bin/env python
# python find_clone.py #.mpileup #.bam index
# Shiya Song
# 27 March 2013
#This script takes mpileup result as input. It generates 1k-window file, possible clones(named fraction#_clone)and possible overlapping clones(named fraction#_clone_overlap) 
#in two different files. It better assign clone breakpoint by supporting read nearby and trim some part of the clone. 
# fraction#_clone_overlap file include clones that has more than 2 supporting read in the middle of the clone, which indicate there is some SV or overlapping clone.

import sys
import genutils
import math
import os
import signal
import numpy as np


chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'


def get_chromOrder():
	chromOrder = {}
	inFile = open(chromLenFile,'r')
	for i,l in enumerate(inFile):
		l = l.rstrip()
		c = l.split()[0]
		chromOrder[c] = i
	return chromOrder

def set_1k_windows(windowFile):
	windows = []
	inFile = open(windowFile,'r')
	for line in inFile:
		line = line.rstrip()
		line = line.split("\t")
#		if line[0]!= "chr1":
#			break
	# i'll do everything in normal 1-based coordinates 
		windows.append([line[0],int(line[1])+1,int(line[2]),0,0,0])	  #chr,start,end,cov,number of neg strand, number of pos strand
	inFile.close()
	return windows

def set_1bp_list():
	chromLens = genutils.read_chrom_len(chromLenFile)
	counts = {}
	for chromName in chromLens:
		counts[chromName] = []
#		if chromName != "chr1":
#			continue
		for i in range(0,chromLens[chromName] + 1):
			counts[chromName].append("")
	return counts
	
def find_trans_read(BamFileName,ChrName,index,counts):
	gc = 'samtools view ' + BamFileName
	myFile = os.popen(gc,'r')
	file_name1 = BamFileName.replace("bam","trans_pos")
	file_out1 = open(file_name1,'w')
	file_name2 = BamFileName.replace("bam","trans_neg")
	file_out2 = open(file_name2,'w')
	file_name = BamFileName.replace("bam","trans_read.txt")
	file_out = open(file_name,'w')
	print >>file_out1, "track type=bedGraph name='%s_pos' description='+' visibility=2 priority=priority color=255,0,0" %(index)
	print >>file_out2, "track type=bedGraph name='%s_neg' description='-' visibility=2 priority=priority color=0,0,255" %(index)
	for i in myFile:
		i = i.strip()
		col = i.split('\t')
		chr1 = col[2]
		chr2 = col[6]
		flag = int(col[1])
		if chr2 == ChrName and flag & 0x0004 == 0 and flag & 0x0008 == 0:
			m = 1
			if chr1 != 'EBV' and chr1 != 'eColiK12':
#				if chr1 != "chr1":
#					continue
				if flag & 0x0010 != 0: # this strand is reversed
					print >>file_out2, "%s\t%d\t%d\t%s" %(chr1,int(col[3]),int(col[3])+76,m)
					print >>file_out, "%s\t%d\t%d\t-" %(chr1,int(col[3]),int(col[3])+76)
					counts[chr1][int(col[3])] += "-"
				else: #this strand is positive
					print >>file_out1, "%s\t%d\t%d\t%s" %(chr1,int(col[3]),int(col[3])+76,m)
					print >>file_out, "%s\t%d\t%d\t+" %(chr1,int(col[3]),int(col[3])+76)
					counts[chr1][int(col[3])] += "+"
	return counts
					
def read_trans_read(read_file,counts):
	f_read = open(read_file,'r')
	for line in f_read:
		line = line.strip()
		col = line.split('\t')
#		if col[0]!= "chr1":
#			break
		counts[col[0]][int(col[1])] += col[3]
	return counts

def mpileup_to_window(mpileupFile,windows,windowSize,counts,windowOut,index):
	signal.signal(signal.SIGPIPE, signal.SIG_DFL)
	gc = 'gunzip -c ' + mpileupFile
	inFile = os.popen(gc, 'r')
	outFile = open(windowOut,'w')
	print >>outFile,"track type=bedGraph name='%s_cov' description='1k cov' visibility=2 priority=priority" %(index)
	windowIndex = 0
	current_c = windows[windowIndex][0]
	current_b = windows[windowIndex][1]
	current_e = windows[windowIndex][2]
	current_total= 0
	neg_support = 0
	pos_support = 0
	for line in inFile:
		line = line.rstrip()
		line = line.split()
		c = line[0]
		pos = int(line[1])
		depth = int(line[3])
#		if c != "chr1":
#			break
		if c == current_c and pos >= current_b and pos <= current_e:
			current_total += depth
		else:
			avgDepth = float(current_total) / windowSize
			windows[windowIndex][3] = avgDepth
			for i in range(current_b,current_e + 1):
				neg_support += counts[current_c][i].count('-')
			for i in range(current_b,current_e + 1):
				pos_support += counts[current_c][i].count('+')
			windows[windowIndex][4] = neg_support
			windows[windowIndex][5] = pos_support
			if windows[windowIndex][0]!='EBV' and windows[windowIndex][0] != 'eColiK12' and windows[windowIndex][0] != 'pCC1FOS':
				outFile.write('%s\t%i\t%i\t%f\t%i\t%i\n' % (windows[windowIndex][0],windows[windowIndex][1]-1,windows[windowIndex][2],windows[windowIndex][3],windows[windowIndex][4],windows[windowIndex][5]))
			#next window
			windowIndex += 1
			current_c = windows[windowIndex][0]
			current_b = windows[windowIndex][1]
			current_e = windows[windowIndex][2]
			current_total= 0
			neg_support = 0
			pos_support = 0 
			# deal with empty windows
			lineProcessed = False
			while lineProcessed is False:
				if c == current_c and pos >= current_b and pos <= current_e:
					current_total += depth
					lineProcessed = True
				else:
					avgDepth = float(current_total) / windowSize
					windows[windowIndex][3] = avgDepth
					for i in range(current_b,current_e + 1):
						neg_support += counts[current_c][i].count('-')
					for i in range(current_b,current_e + 1):
						pos_support += counts[current_c][i].count('+')
					windows[windowIndex][4] = neg_support
					windows[windowIndex][5] = pos_support
					if windows[windowIndex][0]!='EBV' and windows[windowIndex][0] != 'eColiK12' and windows[windowIndex][0] != 'pCC1FOS':
						outFile.write('%s\t%i\t%i\t%f\t%i\t%i\n' % (windows[windowIndex][0],windows[windowIndex][1]-1,windows[windowIndex][2],windows[windowIndex][3],windows[windowIndex][4],windows[windowIndex][5]))
					#next window
					windowIndex += 1
					if windowIndex >= len(windows):
						print 'out of range for window index'
						print line
						sys.exit()
					current_c = windows[windowIndex][0]
					current_b = windows[windowIndex][1]
					current_e = windows[windowIndex][2]
					current_total= 0
					neg_support = 0
					pos_support = 0
#outside of loop, need to print out the last window
	avgDepth = float(current_total) / windowSize
	windows[windowIndex][3] = avgDepth
	for i in range(current_b,current_e + 1):
		neg_support += counts[current_c][i].count('-')
	for i in range(current_b,current_e + 1):
		pos_support += counts[current_c][i].count('+')
	windows[windowIndex][4] = neg_support
	windows[windowIndex][5] = pos_support
	if windows[windowIndex][0]!='EBV' and windows[windowIndex][0] != 'eColiK12' and windows[windowIndex][0] != 'pCC1FOS':
		outFile.write('%s\t%i\t%i\t%f\t%i\t%i\n' % (windows[windowIndex][0],windows[windowIndex][1]-1,windows[windowIndex][2],windows[windowIndex][3],windows[windowIndex][4],windows[windowIndex][5]))	
	inFile.close()
	outFile.close()
	return windows

def read_95_coverage(index):
	cov_file = open("each_fraction_95%_cov",'r')
	threshold = []
	for i in cov_file:
		i = i.strip()
		threshold.append(float(i))
	return threshold[int(index)-1]

def calc_95_coverage(windows):
	each_cov = []
	for i in range(len(windows)):
		if windows[i][3]>0:
			each_cov.append(float(windows[i][3]))
	each_cov.sort()
	k = (len(each_cov)-1) * 0.05
	up = math.floor(k)
	low = math.ceil(k)
	if up == low:
		threshold = each_cov[int(k)]
	else:
		threshold = each_cov[int(up)]
	return threshold

def find_clones(windows,outfile,overlap_outfile,index):
	f_out = open(outfile,'w')
	f_overlap = open(overlap_outfile,'w')
	print >>f_out,"track type=bedGraph name='%s_clone' description='clone' visibility=2 priority=priority " %(index)
	print >>f_overlap,"track type=bedGraph name='%s_overlapclone' description='overlap_clone' visibility=2 priority=priority " %(index)
	windowIndex = 0
	coverage = 0
	m = 0
	current_c = windows[windowIndex][0]
	current_b = windows[windowIndex][1]
	current_e = windows[windowIndex][2]
	depth = windows[windowIndex][3]
	support_start = []  # record the index of neg support
	support_end = []   #record the index of pos_support
	neg_support = 0
	pos_support = 0
	first = False
	neg_intermediate = []
	pos_intermediate = []
	while first is False:
		if depth > 0.00:
			first = True
			start = windowIndex
			end = windowIndex
			coverage += depth
		else:
			windowIndex += 1
			current_c = windows[windowIndex][0]
			current_b = windows[windowIndex][1]
			current_e = windows[windowIndex][2]
			depth = windows[windowIndex][3]
	while windowIndex < len(windows)-1:
		windowIndex += 1
		if windows[windowIndex][3] == 0.00:
			m += 1
		else:
			if windows[windowIndex][0]==current_c and m<=3:
				end = windowIndex
				coverage += windows[windowIndex][3]
				m = 0
			else:
#				if current_c != "chr1":
#					break
				if current_c == "EBV" or current_c =='eColiK12' or current_c == 'pCC1FOS':
					break
				if current_c != 'EBV' and current_c != 'eColiK12' and current_c != 'pCC1FOS':
					for i in range(start-1,end+1):
						if windows[i][4] > 0 and windows[i][0]==current_c:
							support_start.append(i)
					for i in range(start,end+2):
						if windows[i][5] > 0 and windows[i][0]==current_c:
							support_end.append(i)
					if len(support_start) != 0:
						for i in range(0,len(support_start)):
							if support_start[i] >= start-2 and support_start[i]<=start+2:
								if neg_support==0:
									start = support_start[i]
								neg_support += windows[support_start[i]][4]
							else:
								if support_start[i]<=end:
									n_below_avg = 0
									clone_window = []
									for j in range(support_start[i],end+1):
										clone_window.append(windows[j][3])
									threshold = np.mean(clone_window)
									for j in range(start,support_start[i]):
										if windows[j][3] <= threshold:
											n_below_avg += 1
									if n_below_avg >= 2./3*(support_start[i] - start) and (support_start[i]-start)<=(end-support_start[i]+1) and neg_support==0:
										start = support_start[i]
										neg_support = windows[support_start[i]][4]
									else:
										neg_intermediate.append(support_start[i])
								else:
									neg_intermediate.append(support_start[i])
#									print "overlapping clones by neg_support%d:%s\t%i\t%i\t%f\t%d" %(windows[support_start[i]][1],current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1)
					if len(support_end) != 0:
						r = range(0,len(support_end))
						r.reverse()
						for i in r:
							if support_end[i] >= end-2 and support_end[i]<=end+2:
								if pos_support==0 and support_end[i]>= start:
									end = support_end[i]
								pos_support += windows[support_end[i]][5]
							else:
								if support_end[i]>=start:
									n_below_avg = 0
									clone_window = []
									for j in range(start,support_end[i]+1):
										clone_window.append(windows[j][3])
									threshold = np.mean(clone_window)
									for j in range(support_end[i]+1,end+1):
										if windows[j][3] <= threshold:
											n_below_avg += 1
									if n_below_avg >= 2./3*(end - support_end[i]) and (end-support_end[i])<=(support_end[i]-start+1) and pos_support==0:
										pos_support = windows[support_end[i]][5]
									else:
										pos_intermediate.append(support_end[i])
								else:
									pos_intermediate.append(support_end[i])
#									print "overlapping clones by neg_support%d:%s\t%i\t%i\t%f\t%d" %(windows[support_end[i]][1],current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1)
					if len(pos_intermediate)>=2 or len(neg_intermediate) >=2 :
						print "overlapping clones:%s\t%i\t%i\t%f\t%d\tthreshold:%f" %(current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1],threshold)
						print >>f_overlap,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1],neg_support,pos_support)
					elif len(neg_intermediate)==1 and len(pos_intermediate)==1 and neg_intermediate[0]> pos_intermediate[0]:
						if pos_intermediate[0]>=start and neg_intermediate[0]<=end:
							print "two distinct clones:%s\t%i\t%i\t%f\t%d\tthreshold:%f" %(current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1],threshold)
							print >>f_out,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1],windows[pos_intermediate[0]][2],coverage/(end-start+1),windows[pos_intermediate[0]][2]-windows[start][1],neg_support,windows[pos_intermediate[0]][5])
							print >>f_out,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[neg_intermediate[0]][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[neg_intermediate[0]][1],windows[neg_intermediate[0]][4],pos_support)
						else:
							print >>f_overlap,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1],neg_support,pos_support)
					elif (neg_support >0 and len(neg_intermediate)==1 and len(pos_intermediate)==0 and pos_support==0) or (pos_support>0 and len(pos_intermediate)==1 and len(neg_intermediate)==0 and neg_support==0):
						print "possible inversions:%s\t%i\t%i\t%f\t%d\tthreshold:%f" %(current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1],threshold)
						print >>f_overlap,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1],neg_support,pos_support)
					elif len(neg_intermediate)>=1 or len(pos_intermediate)>=1:
						print "mixed clones:%s\t%i\t%i\t%f\t%d\tthreshold:%f" %(current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1],threshold)
						print >>f_overlap,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1],neg_support,pos_support)
					elif len(neg_intermediate)==0 and len(pos_intermediate)==0:
						print >>f_out,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1],neg_support,pos_support)
				current_c = windows[windowIndex][0]
				start = windowIndex
				end = windowIndex
				coverage = 0
				m = 0
				coverage += windows[windowIndex][3]
				support_start = []
				support_end = []
				neg_intermediate = []
				pos_intermediate = []
				neg_support = 0
				pos_support = 0

def find_clones_v2(windows,outfile,overlap_outfile,index):   # This is for older versions where windows are in 1-based
	f_out = open(outfile,'w')
	f_overlap = open(overlap_outfile,'w')
	print >>f_out,"track type=bedGraph name='%s_clone' description='clone' visibility=2 priority=priority " %(index)
	print >>f_overlap,"track type=bedGraph name='%s_overlapclone' description='overlap_clone' visibility=2 priority=priority " %(index)
	windowIndex = 0
	coverage = 0
	m = 0
	current_c = windows[windowIndex][0]
	current_b = windows[windowIndex][1]
	current_e = windows[windowIndex][2]
	depth = windows[windowIndex][3]
	support_start = []  # record the index of neg support
	support_end = []   #record the index of pos_support
	neg_support = 0
	pos_support = 0
	first = False
	neg_intermediate = []
	pos_intermediate = []
	while first is False:
		if depth > 0.00:
			first = True
			start = windowIndex
			end = windowIndex
			coverage += depth
		else:
			windowIndex += 1
			current_c = windows[windowIndex][0]
			current_b = windows[windowIndex][1]
			current_e = windows[windowIndex][2]
			depth = windows[windowIndex][3]
	while windowIndex < len(windows)-1:
		windowIndex += 1
		if windows[windowIndex][3] == 0.00:
			m += 1
		else:
			if windows[windowIndex][0]==current_c and m<=3:
				end = windowIndex
				coverage += windows[windowIndex][3]
				m = 0
			else:
#				if current_c != "chr1":
#					break
				if current_c == "EBV" or current_c =='eColiK12' or current_c == 'pCC1FOS':
					break
				if current_c != 'EBV' and current_c != 'eColiK12' and current_c != 'pCC1FOS':
					for i in range(start-1,end+1):
						if windows[i][4] > 0 and windows[i][0]==current_c:
							support_start.append(i)
					for i in range(start,end+2):
						if windows[i][5] > 0 and windows[i][0]==current_c:
							support_end.append(i)
					if len(support_start) != 0:
						for i in range(0,len(support_start)):
							if support_start[i] >= start-2 and support_start[i]<=start+2:
								if neg_support==0:
									start = support_start[i]
								neg_support += windows[support_start[i]][4]
							else:
								if support_start[i]<=end:
									n_below_avg = 0
									clone_window = []
									for j in range(support_start[i],end+1):
										clone_window.append(windows[j][3])
									threshold = np.mean(clone_window)
									for j in range(start,support_start[i]):
										if windows[j][3] <= threshold:
											n_below_avg += 1
									if n_below_avg >= 2./3*(support_start[i] - start) and (support_start[i]-start)<=(end-support_start[i]+1) and neg_support==0:
										start = support_start[i]
										neg_support = windows[support_start[i]][4]
									else:
										neg_intermediate.append(support_start[i])
								else:
									neg_intermediate.append(support_start[i])
#									print "overlapping clones by neg_support%d:%s\t%i\t%i\t%f\t%d" %(windows[support_start[i]][1],current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1)
					if len(support_end) != 0:
						r = range(0,len(support_end))
						r.reverse()
						for i in r:
							if support_end[i] >= end-2 and support_end[i]<=end+2:
								if pos_support==0 and support_end[i]>= start:
									end = support_end[i]
								pos_support += windows[support_end[i]][5]
							else:
								if support_end[i]>=start:
									n_below_avg = 0
									clone_window = []
									for j in range(start,support_end[i]+1):
										clone_window.append(windows[j][3])
									threshold = np.mean(clone_window)
									for j in range(support_end[i]+1,end+1):
										if windows[j][3] <= threshold:
											n_below_avg += 1
									if n_below_avg >= 2./3*(end - support_end[i]) and (end-support_end[i])<=(support_end[i]-start+1) and pos_support==0:
										pos_support = windows[support_end[i]][5]
									else:
										pos_intermediate.append(support_end[i])
								else:
									pos_intermediate.append(support_end[i])
#									print "overlapping clones by neg_support%d:%s\t%i\t%i\t%f\t%d" %(windows[support_end[i]][1],current_c,windows[start][1],windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1)
					if len(pos_intermediate)>=2 or len(neg_intermediate) >=2 :
						print "overlapping clones:%s\t%i\t%i\t%f\t%d\tthreshold:%f" %(current_c,windows[start][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1,threshold)
						print >>f_overlap,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1,neg_support,pos_support)
					elif len(neg_intermediate)==1 and len(pos_intermediate)==1 and neg_intermediate[0]> pos_intermediate[0]:
						if pos_intermediate[0]>=start and neg_intermediate[0]<=end:
							print "two distinct clones:%s\t%i\t%i\t%f\t%d\tthreshold:%f" %(current_c,windows[start][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1,threshold)
							print >>f_out,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1]-1,windows[pos_intermediate[0]][2],coverage/(end-start+1),windows[pos_intermediate[0]][2]-windows[start][1]+1,neg_support,windows[pos_intermediate[0]][5])
							print >>f_out,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[neg_intermediate[0]][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[neg_intermediate[0]][1]+1,windows[neg_intermediate[0]][4],pos_support)
						else:
							print >>f_overlap,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1,neg_support,pos_support)
					elif (neg_support >0 and len(neg_intermediate)==1 and len(pos_intermediate)==0 and pos_support==0) or (pos_support>0 and len(pos_intermediate)==1 and len(neg_intermediate)==0 and neg_support==0):
						print "possible inversions:%s\t%i\t%i\t%f\t%d\tthreshold:%f" %(current_c,windows[start][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1,threshold)
						print >>f_overlap,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1,neg_support,pos_support)
					elif len(neg_intermediate)>=1 or len(pos_intermediate)>=1:
						print "mixed clones:%s\t%i\t%i\t%f\t%d\tthreshold:%f" %(current_c,windows[start][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1,threshold)
						print >>f_overlap,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1,neg_support,pos_support)
					elif len(neg_intermediate)==0 and len(pos_intermediate)==0:
						print >>f_out,'%s\t%d\t%d\t%f\t%d\t%i\t%i' %(current_c,windows[start][1]-1,windows[end][2],coverage/(end-start+1),windows[end][2]-windows[start][1]+1,neg_support,pos_support)
				current_c = windows[windowIndex][0]
				start = windowIndex
				end = windowIndex
				coverage = 0
				m = 0
				coverage += windows[windowIndex][3]
				support_start = []
				support_end = []
				neg_intermediate = []
				pos_intermediate = []
				neg_support = 0
				pos_support = 0

def read_in_1k_windows(windowOut):
	f = open(windowOut,"r")
	windows = []
	for line in f:
		if line[0]== "t":
			continue
		line = line.rstrip()
		line = line.split("\t")
		windows.append([line[0],int(line[1]),int(line[2]),float(line[3]),int(line[4]),int(line[5])])	  #chr,start,end,cov,size,number of neg strand, number of pos strang
	f.close()
	return windows

# read in input
mpileupFile = sys.argv[1]
bamfile = sys.argv[2]
index = sys.argv[3]
print index

# Set up 1kb windows
windowFile = '/home/jmkidd/kidd-lab-scratch/shiya-projects/indian_Jacob/script/new_genome_window.txt'	  # This one is already sorted!!!
#windowFile = '/home/jmkidd/kidd-lab/jmkidd-projects/dogs/fosmids/process/canFam3_1.1kbp.nogap.windows.bed'

"""
windowSize = 1000.0 # assume each window is 1000 unmasked bp
chromOrder = get_chromOrder()
windows = set_1k_windows(windowFile)
print "1k windows set up done"

# Set up trans_read
counts=set_1bp_list()
counts = find_trans_read(bamfile,'pCC1FOS',index,counts)
read_file = mpileupFile.replace("mpileup","trans_read.txt")
counts=read_trans_read(read_file,counts)
print "trans read set up done"

# record mpileup info to windows
windowOut = mpileupFile.replace("mpileup","1k_coverage")
windows = mpileup_to_window(mpileupFile,windows,windowSize,counts,windowOut,index)
print "mpileup to window done"
"""

windowOut = mpileupFile.replace("mpileup","1k_coverage")
windows = read_in_1k_windows(windowOut)      # You can start from this step if you already get 1k_coverage file

# find clones
outfile = mpileupFile.replace("mpileup","clone")
overlap_outfile = mpileupFile.replace("mpileup","clone_overlap")
find_clones_v2(windows,outfile,overlap_outfile,index)
print "find clones done"