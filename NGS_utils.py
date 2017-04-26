#!/usr/bin/env python
# python NGS_utils.py
# Shiya Song
# 19th June 2013
# grab all the link in the website and download them 
import signal
import os
import genutils
import gzip

def get_chromOrder(species):
	chromOrder = {}
	if species == "human":
		chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
	elif species == "gorilla":
		chromLenFile = "/home/jmkidd/kidd-lab-scratch/shiya-projects/primate_demography/gorGor/gorGor3.1/gorGor3.fa.fai"
	inFile = open(chromLenFile,'r')
	for i,l in enumerate(inFile):
		l = l.rstrip()
		c = l.split()[0]
		chromOrder[c] = i
	return chromOrder

def get_chromOrder_v2(chromLenFile):
	chromOrder = {}
	inFile = open(chromLenFile,'r')
	for i,l in enumerate(inFile):
		l = l.rstrip()
		c = l.split()[0]
		chromOrder[c] = i
	return chromOrder

def get_chromList():
#	chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
	chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg18/genome-index/hg18.fa.fai'
	chr_list=[]
	inFile = open(chromLenFile,'r')
	for i,l in enumerate(inFile):
		l = l.rstrip()
		c = l.split()[0]
		chr_list.append(c)
	return chr_list
chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
chromLens = genutils.read_chrom_len(chromLenFile)

class BedIterator:
	def __init__(self, filename):
		if filename[-2:]=='gz':
			self.file = gzip.open(filename, "r")
		else:
			self.file = open(filename, "r")
	def __iter__(self):
		return self
	def __next__(self):
		return self.next()
	def next(self):
		line = next(self.file)
		fields = line.rstrip().split()
		return tuple(fields)

class VcfIterator:
	def __init__(self, filename):
		if filename[-2:]=='gz':
			self.file = gzip.open(filename, "r")
		else:
			self.file = open(filename, "r")
	def __iter__(self):
		return self  
	def __next__(self):
		return self.next()  	
	def next(self):
		line = next(self.file)
		while line[0] == "#":
			line = next(self.file)
		fields = line.strip().split()
		chrom = fields[0]
		pos = int(fields[1])
		alleles = [fields[3]]
		for alt_a in fields[4].split(","):
			alleles.append(alt_a)
		geno = fields[9][:3]
		phased = geno[1]=='|'
		return (chrom, pos, tuple(alleles), (int(geno[0]), int(geno[2])), phased)

def read_tab_file(file):
	f = open(file,"r")
	locus = []
	for line in f:
		line = line.rstrip().split("\t")
		locus.append(line)
	return locus

def sort_locus(locus,chromOrder):
	locus.sort( key=lambda w: int(w[1]) )
	locus.sort( key=lambda w: chromOrder[w[0]])
	return locus

def print_locus(outfile,locus):
	f = open(outfile,"w")
	for i in locus:
		print >>f,"\t".join(map(str,i))

def print_window(outfile,locus):
	f = open(outfile,"w")
	for i in locus:
		print >>f,"\t".join(map(str,i))
		
def extract_read_from_bam(bamfile,readfile):
	signal.signal(signal.SIGPIPE, signal.SIG_DFL)
	samtools = "samtools view -h " + bamfile
	f_sam = os.popen(samtools, 'r')
	f_read = open(readfile,"r")
	samout = readfile.replace("txt","sam")
	f_out = open(samout,"w")
	read = []
	for i in f_read:
		read.append(i.strip())
	for line in f_sam:
		if line[0]=="@":
			f_out.write(line)
		elif line.strip().split("\t")[0] in read:
			f_out.write(line)
	cmd = "samtools view -bS " + samout + "> " + samout.replace("sam","bam")
	os.popen(cmd)

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

def read_vcf(file):
	f = open(file,'r')
	for line in f:
		if line[0]=='#':
			continue
		line = line.strip().split('\t')
		chr = line[0]
		pos = int(line[1])
		ref = line[3]
		alt = line[4].split(',')
		
def vcf_to_bed(file,index):
	if file[-2:]=='gz':
		f=gzip.open(file, "r")
	else:
		f=open(file,"r")
	locus = []
	for line in f:
		if line[0]=='#':
			continue
		line = line.strip().split('\t')
		if line[0][0]!='c':
			chr = 'chr'+line[0]
		else:
			chr = line[0]
		pos = line[1]
		ref=line[3]
		alt=line[4].split(',')
		GT = line[9+index].split(':')[0]
		if GT[0]!=GT[2]:
			allele = ''
			if GT[0]=='0':
				allele += ref
			else:
				allele += alt[int(GT[0])-1]
			if GT[2]=='0':
				allele += ref
			else:
				allele += alt[int(GT[2])-1]
			locus.append([chr,int(pos)-1,pos,''.join(GT),allele])
	return locus

def vcf_to_bed_by_chr(file,chrom,index):
	chromOrder = get_chromOrder("human")
	if file[-2:]=='gz':
		f=gzip.open(file, "r")
	else:
		f=open(file,"r")
	locus = []
	for line in f:
		if line[0]=='#':
			continue
		line = line.strip().split('\t')
		if line[0][0]!='c':
			chr = 'chr'+line[0]
		else:
			chr = line[0]
		if chr!=chrom and chromOrder[chr]<chromOrder[chrom]:
			continue
		elif chr==chrom:
			pos = line[1]
			ref=line[3]
			alt=line[4].split(',')
			GT = line[9+index].split(':')[0].split('|')
			if GT[0]!=GT[1]:
				allele = ''
				if GT[0]=='0':
					allele += ref
				else:
					allele += alt[int(GT[0])-1]
				if GT[1]=='0':
					allele += ref
				else:
					allele += alt[int(GT[1])-1]
				locus.append([int(pos),GT,allele])
		else:
			break
	return locus

def compare_allele(nuc1,nuc2):
	if nuc1[0]==nuc2[0] and nuc1[1]==nuc2[1]:
		return True

def compare_allele_v2(nuc1,nuc2):
	if nuc1[0]==nuc2[1] and nuc1[0]==nuc2[1]:
		return True

