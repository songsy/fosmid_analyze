#!/usr/bin/env python
# python fosmid_clone_stats.py sample_name
# Shiya Song
# 26th June 2014
# This script summarize clone statistics


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pickle
import sys
import genutils
from subprocess import *
import os
import commands

def get_pool_name(file,pool_name,index):
	inFile = open(file,'r')
	for line in inFile:
		line = line.rstrip()
		line = line.split()
		sample = line[0]
		if sample == 'runName':
			continue
		pool_name[sample]=index
	return pool_name

def find_clone_stats(pool_name,sample_name):
	writeOut = False
	sizes1 = []
	sizes2 = []
	sizes3 = []
	Size = []
	Coverage = []
	minLen = 10000
	maxLen = 50000
	minDP = 0.25
	summary_file1 = open('%s_p1.clones.%s.select.10k-50k.0.25.txt' %(sample_name,commands.getoutput('date +"%m%d%Y"')),'w')
	summary_file2 = open('%s_p2.clones.%s.select.10k-50k.0.25.txt' %(sample_name,commands.getoutput('date +"%m%d%Y"')),'w')
	summary_file3 = open('%s_p3.clones.%s.select.10k-50k.0.25.txt' %(sample_name,commands.getoutput('date +"%m%d%Y"')),'w')
	summary_file1.write('fraction\tnumClones\tmeanLen\tmedianLen\tminLen\tmaxLen\tmeanCoverage\tmedianCoverage\tminCoverage\tmaxCoverage\n')
	outDirBase= '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools/%s/' %(sample_name)
	poolNames = pool_name.keys()
	poolNames.sort()
	for pool in poolNames:
		sizes = []
		coverages = []
		baseDir = outDirBase+pool + '/'
		cloneFileName = baseDir + pool + '.markdup.clone'
		outFileName = cloneFileName + '.sel.%i.%.2f' % (minLen,minDP)
		if writeOut is True:
			outFile = open(outFileName,'w')
		f = open(cloneFileName,'r')
		for line in f:
			line = line.rstrip()
			line = line.split()
			if line[0] == 'track':
				continue
			aveDP = float(line[3])
			cLen = int(line[4])
			if aveDP >= minDP and cLen >= minLen and cLen<= maxLen:
				if pool_name[pool]==1:
					sizes1.append(cLen)
				if pool_name[pool]==2:
					sizes2.append(cLen)
				if pool_name[pool]==3:
					sizes3.append(cLen)
				sizes.append(cLen)
				coverages.append(aveDP)
				if writeOut is True:
					outFile.write('\t'.join(line)+'\n')
		if len(sizes) == 0:
			print 'no clones',cloneFileName
			res = [0,0,0,0,0,0,0,0,0]
		else:
			res = [len(sizes),np.mean(sizes),np.median(sizes),min(sizes),max(sizes),np.mean(coverages),np.median(coverages),min(coverages),max(coverages)]
		if pool_name[pool]==1:
			summary_file1.write(pool+'\t'+'\t'.join(map(str,res))+'\n')
		if pool_name[pool]==2:
			summary_file2.write(pool+'\t'+'\t'.join(map(str,res))+'\n')
		if pool_name[pool]==3:
			summary_file3.write(pool+'\t'+'\t'.join(map(str,res))+'\n')
		Coverage = Coverage+coverages
	Size = sizes1+sizes2+sizes3
	stats = [len(Size),np.mean(Size),np.median(Size),min(Size),max(Size),np.mean(Coverage),np.median(Coverage),min(Coverage),max(Coverage)]
	summary = open('clones.%s.select.10k-50k.0.25.txt' %(commands.getoutput('date +"%m%d%Y"')),'a') 
#	summary.write('sample\tnumClones\tmeanLen\tmedianLen\tminLen\tmaxLen\tmeanCoverage\tmedianCoverage\tminCoverage\tmaxCoverage\n')
	print >>summary,sample_name+'\t'+"\t".join(map(str,stats))
	fig = plt.figure()
	ax1 = fig.add_subplot(211,title=sample_name)
	ax2 = fig.add_subplot(212)
	ax1.set_ylim([0,40000])
	bins = 40
	ax1.hist(sizes1,bins,facecolor='green',alpha=0.5,histtype='stepfilled',label='library1')
	ax1.hist(sizes2,bins,facecolor='blue',alpha=0.5,histtype='stepfilled',label='library2')
	if len(sizes3)!=0:
		ax1.hist(sizes3,bins,facecolor='orange',alpha=0.5,histtype='stepfilled',label='library3')
	ax1.legend(prop={'size':10})
	ax1.set_xlabel('inferred clone insert size')
	ax1.set_ylabel('count')
#	fig.suptitle(sample_name, fontsize=10)
	ax2.hist(Coverage,bins=50,range=(0,10),facecolor='green',alpha=0.75,histtype='stepfilled')
	ax2.set_xlabel('1kb coverage of inferred clone')
	ax2.set_ylabel('count')
	fig.savefig(sample_name+"_clone_stats.pdf",format='pdf')

def find_clone_stats_NA20847(pool_name,sample_name):
	writeOut = False
	Size = []
	Coverage = []
	minLen = 10000
	maxLen = 50000
	minDP = 0.25
	outDirBase= '/home/jmkidd/kidd-lab-scratch/shiya-projects/indian_Jacob/fastq_barcode/pools/'
	poolNames = pool_name.keys()
	poolNames.sort()
	for pool in poolNames:
		sizes = []
		coverages = []
		cloneFileName = outDirBase + '_clone_v4_filter'
		f = open(cloneFileName,'r')
		for line in f:
			line = line.rstrip()
			line = line.split()
			if line[0] == 'track':
				continue
			aveDP = float(line[3])
			cLen = int(line[4])
			if aveDP >= minDP and cLen >= minLen and cLen<= maxLen:
				sizes.append(cLen)
				coverages.append(aveDP)
		if len(sizes) == 0:
			print 'no clones',cloneFileName
			res = [0,0,0,0,0,0,0,0,0]
		else:
			res = [len(sizes),np.mean(sizes),np.median(sizes),min(sizes),max(sizes),np.mean(coverages),np.median(coverages),min(coverages),max(coverages)]
		Coverage = Coverage+coverages
		Size = Size+sizes
	stats = [len(Size),np.mean(Size),np.median(Size),min(Size),max(Size),np.mean(Coverage),np.median(Coverage),min(Coverage),max(Coverage)]
	summary = open('clones.%s.select.10k-50k.0.25.txt' %(commands.getoutput('date +"%m%d%Y"')),'a') 
#	summary.write('sample\tnumClones\tmeanLen\tmedianLen\tminLen\tmaxLen\tmeanCoverage\tmedianCoverage\tminCoverage\tmaxCoverage\n')
	print >>summary,sample_name+'\t'+"\t".join(map(str,stats))
	fig = plt.figure()
	ax1 = fig.add_subplot(211,title=sample_name)
	ax2 = fig.add_subplot(212)
	ax1.set_ylim([0,40000])
	bins = 40
	ax1.hist(sizes1,bins,facecolor='green',alpha=0.5,histtype='stepfilled',label='library1')
	ax1.legend(prop={'size':10})
	ax1.set_xlabel('inferred clone insert size')
	ax1.set_ylabel('count')
#	fig.suptitle(sample_name, fontsize=10)
	ax2.hist(Coverage,bins=50,range=(0,10),facecolor='green',alpha=0.75,histtype='stepfilled')
	ax2.set_xlabel('1kb coverage of inferred clone')
	ax2.set_ylabel('count')
	fig.savefig(sample_name+"_clone_stats.pdf",format='pdf')

def find_pos(pos,list):
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

def find_clone_cov(pool_name,sample_name):
	windowFile = '/home/jmkidd/kidd-lab-scratch/shiya-projects/indian_Jacob/script/new_genome_window.txt'
	windowSize = 1000.0 # assume each window is 1000 unmasked bp
	chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
	chromLens = genutils.read_chrom_len(chromLenFile)
	coverage = {}
	counts = {}
	pos1 = {}
	pos2 = {}
	print 'Initializing chrom lists'
	for chromName in chromLens:
		counts[chromName] = []
		coverage[chromName] = []
		pos1[chromName]=[]
		pos2[chromName]=[]
	f=open(windowFile,'r')
	for line in f:
		line = line.rstrip()
		line = line.split("\t")
		c = line[0]
		if c=="EBV" or c=="eColiK12" or c=="pCC1FOS":
			continue
		b = int(line[1])
		e = int(line[2])
		counts[c].append(0)
		coverage[c].append(0)
		pos1[c].append(b)
		pos2[c].append(e)
	outDirBase= '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/pools/%s/' %(sample_name)
	poolNames = pool_name.keys()
	poolNames.sort()
	minLen = 10000
	minDP = 0.25
	for pool in poolNames:
		baseDir = outDirBase+pool + '/'
		cloneFileName = baseDir + pool + '.markdup.clone'
		inFileName = cloneFileName + '.sel.%i.%.2f' % (minLen,minDP)
		f = open(inFileName,'r')
		for line in f:
			line = line.rstrip()
			line = line.split()
			if line[0] == 'track':
				continue
			chr = line[0]
			start = int(line[1])
			end = int(line[2])
			cov = float(line[3])
			size = int(line[4])/1000
			find1 = find_pos(start,pos1[chr])
			find2 = find_pos(end,pos2[chr])
			if find1>=0:
				for j in range(find1,find2+1):
					if j < len(counts[chr]):
						counts[chr][j]+=1
						coverage[chr][j]+=cov
					else:
						print line,find1,find2,len(counts[chr])
			else:
				print chr,pos1
	
	'''
	f=open(windowFile,'r')
	i = 0
	outFile = open("wgs_clone_coverage.txt","w")
	old_c="chr1"
	for line in f:
		line = line.rstrip()
		line = line.split("\t")
		c = line[0]
		if c!=old_c:
			i = 0
			old_c = c
		if c=="EBV" or c=="eColiK12" or c=="pCC1FOS":
			continue
		b = int(line[1]) + 1 #make them all 1 based
		e = int(line[2])
		outFile.write('%s\t%i\t%i\t%f\t%i\n' % (c,b,e,coverage[c][i],counts[c][i]))
		i +=1
	dbfile=open('NA19240_wgs_clone_coverage_pickle','wb')
	pickle.dump(counts,dbfile)
	pickle.dump(coverage,dbfile)
	'''
	return counts,coverage

def find_clone_cov_plot(counts,coverage,sample):
#	f=open('clone_stats.txt','a')
	Count = []
	Coverage = []
	for chr in counts.keys():
		if chr[:4]=="chrU":
			continue
		for i in range(len(counts[chr])):
			Count.append(int(counts[chr][i]))
			Coverage.append(float(coverage[chr][i]))
	m=[np.mean(Count),np.median(Count),np.std(Count),np.mean(Coverage),np.median(Coverage),np.std(Coverage),min(Coverage),max(Coverage)]
#	print >>f,sample
#	print >>f,"\t".join(map(str,m))
	y=np.bincount(Count)
#	print >>f,"# of 1kb regions not covered",y[0], "# of 1kb regions covered by 1 clone", y[1], "total regions", len(Count)
	fig = plt.figure()
	ax1 = fig.add_subplot(211,title=sample)
	ax2 = fig.add_subplot(212)
	ax1.hist(Count,bins=30,range=(0,30),facecolor='blue',alpha=0.75,histtype='stepfilled')
	ax1.set_xlabel('number of clone per 1kb window')
	ax1.set_ylabel('count')
	ax2.hist(Coverage,bins=50,range=(0,100),facecolor='green',alpha=0.75,histtype='stepfilled')
	ax2.set_xlabel('coverage per 1kb window')
	ax2.set_ylabel('count')
	fig.savefig("%s_clone_wgs_coverage.pdf" %(sample),format='pdf')
	
if __name__=="__main__":
	'''
	sample_name = sys.argv[1]
	wk = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/'
	if sample_name =='HG02799':
		pool_name = {}
		file = wk+'HG02799-fosmids.plate4.txt.info'
		pool_name=get_pool_name(file,pool_name,2)
		file = wk+'HG02799-fosmids.plates1_2.txt.info'
		pool_name=get_pool_name(file,pool_name,1)
	if sample_name =='HG03108':
		pool_name = {}
		file = wk+'HG03108-fosmids.plate2.txt.info'
		pool_name=get_pool_name(file,pool_name,1)
		file = wk+'HG03108-fosmids.plates3_4.txt.info'
		pool_name=get_pool_name(file,pool_name,2)
	if sample_name =='NA21302':
		pool_name = {}
		file = wk+'NA21302-fosmids.plate1.txt.info'
		pool_name=get_pool_name(file,pool_name,1)
		file = wk+'NA21302-fosmids.plate2.txt.info'
		pool_name=get_pool_name(file,pool_name,2)
		file = wk+'HG21302_p3.txt.info'
		pool_name=get_pool_name(file,pool_name,3)

#	find_clone_stats(pool_name,sample_name)
	counts,coverage=find_clone_cov(pool_name,sample_name)
	find_clone_cov_plot(counts,coverage,sample_name)
	'''
	
	find_clone_stats_NA20847(pool_name,'NA20847')
	counts,coverage=find_clone_cov(pool_name,'NA20847')
	find_clone_cov_plot(counts,coverage,'NA20847')