#!/usr/bin/env python3

import sys
import re
import gzip
import utils
import argparse
import io
from vcf_utils import *
import pickle
def chr_list():
	list =[]
	chromOrder={}
	j = 1
	for i in range(1,23):
		list.append('chr'+str(i))
		chromOrder['chr'+str(i)]=j
		j=j+1
	list.append('chrX')
	chromOrder['chrX']=j
	chromOrder['chrY']=j+1
	chromOrder['chr11_gl000202_random']=j+2
	return list,chromOrder

chr_list,chromOrder=chr_list()
def find_pos(pos,list):
	left = 0
	right = len(list)
	find = -1
	while right - left >0:
		midpoint = int((right + left)/2)
		if pos < list[midpoint]:
			right = midpoint
		elif pos > list[midpoint]:
			left = midpoint+1
		elif pos == list[midpoint]:
			left = midpoint
			find = midpoint
			break
	return find
	
def read_pass(vcf):
	snp = []
	for i in VcfPASSIterator(vcf):
		if i[0]==CHROM:
			if i[2]:
				snp.append(i[1])
	return snp

if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--depth", type=float)
	parser.add_argument("--gvcf")
	parser.add_argument("--vcf")
	parser.add_argument("--chr")
	parser.add_argument("--mask")
	parser.add_argument("--log")
	parser.add_argument("--legend_file", help="Impute2 reference panel legend file, can be gzipped or not")
	parser.add_argument("--minGQ", type=float, default=20)
	parser.add_argument("--minMapQ", type=float, default=20)
	
	args = parser.parse_args()
	sites_parser = None
	if args.legend_file is not None:
		sites_parser = utils.LegendParser(args.legend_file)
	minDepth = args.depth / 2.0
	maxDepth = args.depth * 2.0
	lastPos = 0
	line_cnt = 0
	chr_ = ""
	CHROM = args.chr
	snp = read_pass(args.vcf)
	if args.gvcf[-2:]=='gz':
		f = io.TextIOWrapper(gzip.open(args.gvcf, "r"))
	else:
		f = io.TextIOWrapper(open(args.gvcf, "r"))
	additional = 0
	callable = 0
	pass_snp = 0
	het_snp = 0
	het_snp_additional = 0
	filtered = 0
	het_snp_filtered = 0
	filter1 = 0
	filter2 = 0
	DP_list = []
	for line in f:
		if line[0] == '#':
			print(line, end="")
			continue
		fields = line.strip().split('\t')
		if chr_ == "":
			chr_ = fields[0]
			mask = utils.MaskGenerator(args.mask, chr_)
		pos = int(fields[1])
		refAllele = fields[3]
		altAllele = fields[4]
		info = fields[7]
		genotypes = fields[9].split(":")
#		if line_cnt % 10000 == 0:
#			print("parsing position {}".format(pos), file=sys.stderr)
		line_cnt += 1
		
		if sites_parser is not None:
			while not sites_parser.end and sites_parser.pos < pos:
				sites_parser.tick()

		if altAllele=="<NON_REF>":	 # reference block
			dp = int(genotypes[1])
			gq = int(genotypes[2])
			end_pos_match = re.search("END=(\d+)",info)
			if end_pos_match:
				end_pos = int(end_pos_match.group(1))
			else:
				print('End not found:',line,file=sys.stderr)
				end_pos = pos
			if dp >= minDepth and dp <= maxDepth and gq >= args.minGQ:
				for each_pos in range(pos,end_pos+1):
					if sites_parser is not None:
						while not sites_parser.end and sites_parser.pos < each_pos:
							sites_parser.tick()
					mask.addCalledPosition(each_pos)
					callable +=1
					if sites_parser is not None and sites_parser.pos == each_pos:
						fields[1] = str(each_pos)
						fields[3] = sites_parser.ref_a
						fields[4] = sites_parser.alt_a
						fields[7]="."
						fields[8]="GT:DP:GQ:PL"
						fields[9]=genotypes[0]+":"+genotypes[1]+":"+genotypes[2]+":"+genotypes[4]
						print("\t".join(fields))
		elif re.match("^[ACTGactg]$", refAllele) and re.match("^[ACTGactg],<NON_REF>$", altAllele):		# snps, exclude indels
			dp_match = re.search("DP=(\d+)", info)
			mq_match = re.search("MQ=(\d+)", info)
			if not (dp_match and mq_match):
				continue
			dp = int(dp_match.group(1))
			mq = int(mq_match.group(1))
			fields[4]=fields[4].replace(',<NON_REF>','')
			if dp >= minDepth and dp <= maxDepth and mq >= args.minMapQ:
#			if mq >= args.minMapQ:
				find = find_pos(pos,snp)
				if altAllele != "." and re.match("^[01][/|][01]", genotypes[0]):
					if find >=0:							# PASS VQSR filtering
						print("\t".join(fields))
						mask.addCalledPosition(pos)
						callable +=1
						pass_snp +=1
						if genotypes[0]=='0/1':
							het_snp +=1
					elif sites_parser is not None and sites_parser.pos == pos:			# In reference panel
						assert refAllele == sites_parser.ref_a
						print("\t".join(fields))
						mask.addCalledPosition(pos)
						callable +=1
						additional +=1
						if genotypes[0]=='0/1':
							het_snp_additional +=1
					else:
						filtered +=1
						if genotypes[0]=='0/1':
							het_snp_filtered +=1
				else:
					print('Error1',line,file=sys.stderr)
			else:
				find = find_pos(pos,snp)
				if find >=0:
					filter1 +=1
				elif sites_parser is not None and sites_parser.pos == pos:
					filter2 +=1
	f_out = open(args.log,'a')
	print(args.chr,callable,pass_snp,het_snp,additional,het_snp_additional,filtered,het_snp_filtered,filter1,filter2,file=f_out)
