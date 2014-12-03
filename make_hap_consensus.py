#!/usr/bin/env python3

import gzip
from vcf_utils import *
import argparse
import io
import re
import subprocess

chromLenFile = '/home/jmkidd/kidd-lab/genomes/hg19/hg19-EBV-fos-ecoli/hg19.fa.fai'
chromLens = read_chrom_len(chromLenFile)
chromOrder = get_chromOrder("human")
def build_chr_array():
	hap1 = []
	hap2 = []
	for i in range(0,chromLens[CHROM]+2):
		hap1.append(0)
		hap2.append(0)
	return hap1,hap2

def read_whole_genome_vcf(whole_genome_vcf,chrom):
	f = open(whole_genome_vcf,"r")
	whole_indel = {}
	whole_indel[chrom]=[]
	for line in f:
		if line[0]=="#":
			continue
		line = line.strip().split("\t")
		chr = line[0]
		pos = int(line[1])
		ref = line[3]
		alt = line[4].split(",")
		geno = line[9][:3]
		if chr!=chrom and chromOrder[chr]<chromOrder[chrom]:
			continue
		elif chr==chrom:
			whole_indel[chr].append([pos,ref,alt,geno])
		else:
			break
	return whole_indel

def find_pos(pos,list):
	left = 0
	right = len(list)
	find = -1
	while right - left >0:
		midpoint = (right + left)/2
		if pos < list[midpoint][0]:
			right = midpoint
		elif pos > list[midpoint][0]:
			left = midpoint+1
		elif pos == list[midpoint][0]:
			left = midpoint
			find = midpoint
			break
	return find

def read_phased_snps(file):
	f = open(file,"r")
	SNP = {}
	for line in f:
		line = line.strip().split("\t")
		chr = line[0]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			SNP[line[0]+" "+line[1]]=line[2]+line[3]
		else:
			break
	return SNP

def read_hap_snp_vcf(file,hap):	#This is for GATK, emit all sites call, snp only
	f = gzip.open(file, 'r')
	het =0
	concord = 0
	force_concord =0
	discord =0
	unphased_decide = 0
	extra = 0
	MIN_QUAL=20
	MIN_DP = 5
	MAX_DP = 80
	QUAL = []
	DP = []
	for line in io.TextIOWrapper(f):
		if line[0]=="#":
			continue
		line = line.strip().split('\t')
		pos = int(line[1])
		chr = line[0]
		info = line[7]
		ref = line[3]
		alt = line[4].split(',')
		GT = line[9].split(":")[0]
		m1 = re.search('DP=(\d+);',info)
		m2 = re.search('MQ=(\d+).',info)
		if GT=="0/0":
			if int(m1.group(1))>=MIN_DP and int(m1.group(1))<=MAX_DP and int(m2.group(1))>=MIN_QUAL:
				hap[pos]=ref
		elif GT=="0/1" or GT=="1/2":
			continue
		elif GT=="1/1":
			if int(m1.group(1))>=MIN_DP and int(m1.group(1))<=MAX_DP and int(m2.group(1))>=MIN_QUAL:
				hap[pos]=alt[0]
		elif GT=="2/2":
			if int(m1.group(1))>=MIN_DP and int(m1.group(1))<=MAX_DP and int(m2.group(1))>=MIN_QUAL:
				hap[pos]=alt[1]
	return hap

def read_whole_genome_all_sites_vcf(wgs_vcf,hap1,hap2,SNP,mask):
	f_out=io.TextIOWrapper(gzip.open("/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/BLOCK/{}/{}_{}_wgs_first_phased.vcf.gz".format(SAMPLE,SAMPLE,CHROM),'w'))
	f_out.write('##fileformat=VCFv4.1\n')
#	f_out.write('##fileDate={}\n'.format(subprocess.getoutput('date +"%m%d%y"')))
	f_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n')
	f_out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(SAMPLE))
	concord = 0
	opposite = 0
	discord = 0
	additional = 0
	unphased = 0
	MIN_QUAL=20
	MIN_DP = 5
	MAX_DP = 50
	a1 = None
	a2 = None
	for i in VcfIterator(wgs_vcf):
		wgs_snp = True
		call = True
		pos = i[1]
		allele = i[2]
		geno = i[3]
		m1 = i[4]
		m2 = i[5]
		if m1==None or m2==None:
			if hap1[pos]!=0 and hap2[pos]!=0:
				try:
					a1 = allele.index(hap1[pos])
					a2 = allele.index(hap2[pos])
				except ValueError:
					continue
				mask.addCalledPosition(pos)
				if a1!=0 or a2!=0:
					geno = str(a1)+'|'+str(a2)
					f_out.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t{}\n'.format(CHROM,pos,allele[0],",".join(allele[1:]),geno))

		else:
			if geno[0]==geno[2]:
				if int(m1.group(1))<MIN_DP or int(m1.group(1))>MAX_DP or int(m2.group(1))<MIN_QUAL or geno[0]=='.':
					try:
						a1 = allele.index(hap1[pos])
						a2 = allele.index(hap2[pos])
						call = True
					except ValueError:
						call = False
					if call:
						mask.addCalledPosition(pos)
						if a1!=0 or a2!=0:
							geno = str(a1)+'|'+str(a2)
							f_out.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t{}\n'.format(CHROM,pos,allele[0],",".join(allele[1:]),geno))
				else:
					mask.addCalledPosition(pos)
					if geno[0]!='0':
						geno=geno.replace('/','|')
						f_out.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t{}\n'.format(CHROM,pos,allele[0],",".join(allele[1:]),geno))
			else:
				try:
					g1 = SNP[CHROM+" "+str(pos)][0]
					g2 = SNP[CHROM+" "+str(pos)][1]
				except KeyError:
					wgs_snp = False
				if wgs_snp:
					if g1==hap1[pos] and g2==hap2[pos]:
						concord +=1
					else:
						if g1==hap2[pos] and g2==hap1[pos]:
							opposite +=1
						else:
							discord +=1
#							print(i,g1,g2,hap1[pos],hap2[pos])
					try:
						a1 = allele.index(g1)
						a2 = allele.index(g2)
					except ValueError:
						print("error1","\t","\t".join(i))
						continue
					mask.addCalledPosition(pos)
					if a1!=0 or a2!=0:
						geno = str(a1)+'|'+str(a2)
						f_out.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t{}\n'.format(CHROM,pos,allele[0],",".join(allele[1:]),geno))
				else:
					try:
						a1 = allele.index(hap1[pos])
						a2 = allele.index(hap2[pos])
					except ValueError:
						if int(m1.group(1))>=MIN_DP and int(m1.group(1))<=MAX_DP and int(m2.group(1))>=MIN_QUAL:
							unphased +=1
							f_out.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t{}\n'.format(CHROM,pos,allele[0],",".join(allele[1:]),geno))
							continue
						continue
					additional +=1
					mask.addCalledPosition(pos)
					if a1!=0 or a2!=0:
						geno = str(a1)+'|'+str(a2)
						f_out.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t{}\n'.format(CHROM,pos,allele[0],",".join(allele[1:]),geno))
	print(CHROM,concord,opposite,discord,additional,unphased)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='assign parental allele in block')
	parser.add_argument("--chr", dest='chr',help="sample name")
	parser.add_argument("--mask", dest='mask',help="sample name")
	args = parser.parse_args()
	
	SAMPLE = 'NA19240'
	CHROM = args.chr
#	whole_genome_vcf = "/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/merge/NA19240_combined_pass_vcf"
#	whole_snp= read_whole_genome_vcf(whole_genome_vcf,chr)
	phased_snp_file="/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/phased_snps_all"
	SNP=read_phased_snps(phased_snp_file)
#	not_shown_file="/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/SNPs_not_shown_April_July_v3_bed"
#	SNP_not_shown = read_phased_snps(not_shown_file,chr)
	hap1,hap2=build_chr_array()
	
	hap1_file = "NA19240_hap1."+CHROM+"_snp.vcf.gz"
	hap2_file = "NA19240_hap2."+CHROM+"_snp.vcf.gz"
	
	hap1=read_hap_snp_vcf(hap1_file,hap1)
	hap2=read_hap_snp_vcf(hap2_file,hap2)
	
	mask = MaskGenerator(args.mask, CHROM)
	all_sites_vcf="/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/merge/vcf_by_chr/NA19240_" + CHROM+"_all_sites.vcf.gz"
	read_whole_genome_all_sites_vcf(all_sites_vcf,hap1,hap2,SNP,mask)