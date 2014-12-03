#!/usr/bin/env python3

import sys
import gzip
import string
import copy
import argparse
import io
import re

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

def read_chrom_len(faiFileName):
    chromLens = {}
    inFile = open(faiFileName,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        chromLens[line[0]] = int(line[1])
    inFile.close()
    return chromLens    

class MaskGenerator:
    def __init__(self, filename, chr):
        self.lastCalledPos = -1
        self.lastStartPos = -1
        self.file = io.TextIOWrapper(gzip.open(filename, "w"))
        self.chr = chr
    
    # assume 1-based coordinate, output in bed format
    def addCalledPosition(self, pos):
        if self.lastCalledPos == -1:
            self.lastCalledPos = pos
            self.lastStartPos = pos
        elif pos == self.lastCalledPos + 1:
            self.lastCalledPos = pos
        else:
            self.file.write("{}\t{}\t{}\n".format(self.chr, self.lastStartPos - 1, self.lastCalledPos))
            self.lastStartPos = pos
            self.lastCalledPos = pos

class MaskIterator:
	def __init__(self, filename, negative=False):
		if filename[-3:] == ".gz":
			self.file = io.TextIOWrapper(gzip.open(filename, "r"))
		else:
			self.file = io.TextIOWrapper(open(filename, "r"))
		self.eof = False
		self.lastPos = 1
		self.negative = negative
		self.readLine()

	def readLine(self):
		try:
			line = next(self.file)
			fields = line.strip().split()
			if len(fields) == 2:
				self.start = int(fields[0])
				self.end = int(fields[1])
			else:
				self.start = int(fields[1]) + 1
				self.end = int(fields[2])
		except StopIteration:
			self.eof = True
	def getVal(self, pos):
		assert pos >= self.lastPos
		self.lastPos = pos
		while pos > self.end and not self.eof:
			self.readLine()
		if pos >= self.start and pos <= self.end:
			return True if not self.negative else False
		else:
			return False if not self.negative else True

class MergedMask:
	def __init__(self, mask_iterators):
		self.maskIterators = mask_iterators

	def getVal(self, pos):
		return all((m.getVal(pos) for m in self.maskIterators))

class VcfIterator:
	def __init__(self, filename):
		self.file = io.TextIOWrapper(gzip.open(filename, "r"))
  
	def __iter__(self):
		return self
  
	def __next__(self):
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
		m1 = re.search('DP=(\d+);',fields[7])
		m2 = re.search('MQ=(\d+).',fields[7])
		return (chrom, pos, tuple(alleles),geno,m1,m2)

class VcfPASSIterator:
	def __init__(self, filename):
		if filename[-3:] == ".gz":
			self.file = io.TextIOWrapper(gzip.open(filename, "r"))
		else:
			self.file = io.TextIOWrapper(open(filename, "r"))
  
	def __iter__(self):
		return self
  
	def __next__(self):
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
		pass_not = fields[6]=='PASS'
		return (chrom, pos, pass_not)