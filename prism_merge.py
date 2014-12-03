import sys
import commands
import os
import argparse
from NGS_utils import *
import glob
import re
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab
import pickle

chromOrder = get_chromOrder("human")
def merge_phase(info_file,chr,id_to_block,block_to_snp,block_to_mec):
	f_list=[]
	f=open('%s.chr%s.prism.cmds' %(args.sample,chr),'r')
	for line in f:
		m=re.search('%s.chr%s.all.\d*.\d*.phase' %(args.sample,chr),line)
		if m is not None:
			f_list.append(m.group())
		else:
			print line
	f_out=open('%s.chr%s.all.phase' %(args.sample,chr),'w')
#	f_out2=open('%s.chr%s.all.path' %(args.sample,chr),'w')
	flip=False
	start_id = None
	end_id = None
	end = None
	first = True
	for file in f_list:
		start = None
		if not os.path.isfile('prism_all_block/'+str(file)):
			cmds = "grep '%s' %s.chr%s.prism.cmds " %(str(file),args.sample,chr)
			print cmds
#			os.popen(cmds)
			continue
		for i in BedIterator('prism_all_block/'+str(file)):
			if i[0]=='#ID':
				continue
			if start is None:
				start = i[1]
				start_id = i[0]
			end = i[1]
			end_id = i[0]
		if first is True:
			for i in BedIterator('prism_all_block/'+str(file)):
				if i[0]=='#ID':
					continue
				string = '\t'.join(block_to_mec[id_to_block[i[0]]])
				snp = block_to_snp[id_to_block[i[0]]]
				print >>f_out,'%s\t%s\t%s\t%s\t%s' %(i[0],i[1],i[2],snp,string)
			first = False
			last_end = end
			last_end_id = end_id
#			for i in BedIterator('prism_all_block/'+str(file).replace('phase','path')):
#				if i[0]=='#ID':
#					continue
#				print >>f_out2,'%s\t%s\t%s\t%s\t%s' %(i[0],i[1],i[2],i[3],i[4])
		else:
			try:
				assert last_end_id ==start_id
			except AssertionError:
				print 'redo',str(file),chr,last_end_id,start_id
			if last_end ==start:
				for i in BedIterator('prism_all_block/'+str(file)):
					if i[0]=='#ID':
						continue
#					if i[0]==start_id:
#						continue
					try:
						m=id_to_block[i[0]]
					except KeyError:
						print 'error',file,i[0]
					try:
						snp = block_to_snp[id_to_block[i[0]]]
						string = '\t'.join(block_to_mec[id_to_block[i[0]]])
					except KeyError:
						print i[0],m
						continue
					print >>f_out,'%s\t%s\t%s\t%s\t%s' %(i[0],i[1],i[2],snp,string)
				last_end = end
				last_end_id = end_id
#				for i in BedIterator('prism_all_block/'+str(file).replace('phase','path')):
#					if i[0]=='#ID':
#						continue
#					print >>f_out2,'%s\t%s\t%s\t%s\t%s' %(i[0],i[1],i[2],i[3],i[4])
			elif int(last_end)==1-int(start):
				print 'flip',start_id,end_id
				for i in BedIterator('prism_all_block/'+str(file)):
					if i[0]=='#ID':
						continue
#					if i[0]==start_id:
#						continue
					try:
						snp = block_to_snp[id_to_block[i[0]]]
						string = '\t'.join(block_to_mec[id_to_block[i[0]]])
					except KeyError:
						print i[0],m
						continue
					print >>f_out,'%s\t%s\t%s\t%s\t%s' %(i[0],str(1-int(i[1])),i[2],snp,string)
				last_end = str(1-int(end))
				last_end_id = end_id
#				for i in BedIterator('prism_all_block/'+str(file).replace('phase','path')):
#					if i[0]=='#ID':
#						continue
#					print >>f_out2,'%s\t%s\t%s\t%s\t%s' %(i[0],i[2],i[1],i[4],i[3])
	
def get_id_to_block(file):
	id_to_block={}
	for i in BedIterator(file):
		if i[4]=='None,drop' or i[4]=='interleaving':
			continue
		id_to_block[i[3]]=i[1]
	return id_to_block

def get_block_snp(file):
	block_to_snp={}
	for i in BedIterator(file):
		block_to_snp[i[1]]=len(i)-4
	return block_to_snp

def get_block_mec(file):
	block_to_mec={}
	for i in BedIterator(file):
		block_to_mec[i[1]]=i
	return block_to_mec

def phase_info(info_file,chr,id_to_block,block_to_snp,block_to_mec):
	f_list=[]
	f=open('%s.chr%s.prism.cmds' %(args.sample,chr),'r')
	for line in f:
		m=re.search('%s.chr%s.all.\d*.\d*.phase' %(args.sample,chr),line)
		if m is not None:
			f_list.append(m.group())
		else:
			print line
	f_out=open('%s.chr%s.all.phase' %(args.sample,chr),'w')
	f_redo='%s.prism.redo.v2.cmds' %(args.sample)
	flip=False
	start_id = None
	end_id = None
	end = None
	first = True
	for file in f_list:
		start = None
		if not os.path.isfile(str(file)):
			cmds = "grep '%s' %s.chr%s.prism.cmds >> %s" %(str(file),args.sample,chr,f_redo)
			print cmds
			os.popen(cmds)
			continue
		for i in BedIterator(str(file)):
			if i[0]=='#ID':
				continue
			if start is None:
				start = i[1]
				start_id = i[0]
			end = i[1]
			end_id = i[0]
			for i in BedIterator(str(file)):
				if i[0]=='#ID':
					continue
				try:
					m=id_to_block[i[0]]
				except KeyError:
					print 'error',file,i[0]
					continue
				string = '\t'.join(block_to_mec[id_to_block[i[0]]])
				snp = block_to_snp[id_to_block[i[0]]]
				print >>f_out,'%s\t%s\t%s\t%s\t%s' %(i[0],i[1],i[2],snp,string)

def draw_hist(chr,affy_posterior,None_posterior):
	posterior=[]
	f='%s.chr%s.all.phase' %(args.sample,chr)
	snp_num=[]
	for i in BedIterator(f):
		posterior.append(float(i[2]))
		snp_num.append(int(i[3]))
		if i[7]=="None":
			None_posterior.append(float(i[2]))
		else:
			affy_posterior.append(float(i[2]))
	heatmap, xedges, yedges = np.histogram2d(posterior, snp_num, bins=[[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],[0,10,100,500,1000,max(snp_num)]])
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	plt.clf()
	plt.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
#	X,Y=np.meshgrid(xedges,yedges)
#	plt.pcolormesh(X,Y,heatmap)
#	fig,ax1=plt.subplots(figsize=(7,7))
#	fig.canvs.set_window_title(args.sample+' chr'+str(chr)+' posterior & number of snp correlation')
#	xticks=pylab.setp(ax1,)
	pos = range(0,max(snp_num),max(snp_num)/5)
	if len(pos)==5:
		pos+=[max(snp_num)]
	pylab.yticks(pos,["0","10","100","500","1000",str(max(snp_num))])	
	plt.colorbar()
	plt.xlabel('posterior')
	plt.ylabel('number of snp')
	plt.title(args.sample+' chr'+str(chr)+' posterior & number of snp heatmap')
	plt.savefig(args.sample+'_'+str(chr)+'_heatmap.png',format='png',dpi=200)
	return affy_posterior,None_posterior

def calc_switch_error(chr,tot_switch_error,tot_block):
	f='%s.chr%s.all.phase' %(args.sample,chr)
	count = [0,0]
	count_snp = [0,0]
	first = True
	switch_error = 0
	phase={}
	for i in BedIterator(f):
		if i[7]!='None':
			count[int(i[1])]+=1
			count_snp[int(i[1])]+=int(i[8])
			if first:
				old = int(i[1])
				first = False
			else:
				if old!=int(i[1]):
					switch_error +=1
					old = int(i[1])
		phase[i[4]+" "+i[5]]=int(i[1])
	print count,count_snp,switch_error-1
	tot_switch_error+=switch_error-1
	tot_block +=sum(count)
	return tot_switch_error,tot_block,phase

def calc_N50(num):
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

def read_hap(BLOCK,block_list,phase):
	snp_decode={}
	for i in range(len(block_list)):
		chr = block_list[i].split(" ")[0]
		pos = block_list[i].split(" ")[1]
		try:
			phase_info = phase[block_list[i]]
		except KeyError:
			phase_info = None
		for j in BLOCK[block_list[i]].index:
			if BLOCK[block_list[i]]["hap"][j] !=0.5:
				allele1 = int(BLOCK[block_list[i]]["hap"][j])
				allele2 = 1-int(BLOCK[block_list[i]]["hap"][j])
				if phase_info==0:
					snp_decode[j]=str(allele1)+'|'+str(allele2)
				elif phase_info ==1:
					snp_decode[j]=str(allele2)+'|'+str(allele1)
			else:
				allele1 = int(BLOCK[block_list[i]]["hap"][j])
				allele2 = 1-int(BLOCK[block_list[i]]["hap"][j])
				snp_decode[j]='0/1'
	return snp_decode

def SNPs_not_shown(file):
	f = open(file,"r")
	SNP_not_shown = {}
	SNP_not_shown[CHROM]=[]
	for line in f:
		line = line.strip().split("\t")
		chr = line[1]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			SNP_not_shown[line[1]].append(int(line[2]))
		else:
			break
	return SNP_not_shown

def read_haplotype(hap_file,SNP_not_shown):
	f = open(hap_file,"r")
	BLOCK = {}
	block_list=[]
	decode = {}
	mec = {}
	start = True
	for i in range(1,2800000):
		decode[i] = 0
	for line in f:
		line = line.strip()
		if line[0]=="B":
			start = True
		elif line[0]=="*":
			if start == False:
				block_info = {}
				block_info["hap"] = pd.Series(hap,index=block)
				BLOCK[chr+" "+pos] = pd.DataFrame(block_info)
		else:
			col = line.split("\t")
			if start == True:
				pos = col[4]
				chr = col[3]
				if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
					continue
				elif chr==CHROM:
					find = find_pos(int(pos),SNP_not_shown[chr])
					if find <0:			# SNPs are shown in the matrix
						start = False			
						pos = col[4]
						block = []				# store the position in each block
						block_list.append(chr+" "+pos)		# store the first position of each block
						hap = []				# store the 1/0 in each block
						mec[chr+" "+pos] = [0,0,0,0,0,0,0,0,0,0]		# number of snp, number of snp from Affy, number of support from Affy, MEC value before, MEC value after correction, switch_error, switch_error_corrected, source type
					else:		# SNPs not shown in the matrix
						continue
				else:
					break
			find = find_pos(int(col[4]),SNP_not_shown[chr])
			if find >=0:
				continue
			if col[1] == "-":
				block.append(col[4])
				hap.append(0.5)
				if decode[int(col[0])]==0:
					decode[int(col[0])] = chr +" "+ pos
			else:
				block.append(col[4])	  # just the position of the snp
				hap.append(int(col[1]))			# "1" or "0" for the first haplotype
				decode[int(col[0])] = chr +" "+ pos		# link the refhap result (# of the snp) to the first snp of that block	
	f.close()
	return BLOCK,block_list,decode,mec
def write_vcf_v2(wgs_vcf,snp_decode):
	f = gzip.open(wgs_vcf,'r')
	f_out=gzip.open("%s%s/gVCF_calls/%s.%s.prism.v2.phased.vcf.gz" %(args.wgs_dir,SAMPLE,SAMPLE,CHROM),'w')
	for line in f:
		line = line.strip()
		if line[0] == "#":
			f_out.write(line+'\n')
			continue
		line1 = line.split('\t')
		chr = line1[0]
		if chr!=CHROM and chromOrder[chr]<chromOrder[CHROM]:
			continue
		elif chr==CHROM:
			ref = line1[3]
			alt = line1[4].split(',')
			pos = line1[1]
			info = line1[9]
			info = info.split(':')
			if info[0] == "1/2":
				if len(ref)!=1 or len(alt[0])!=1 or len(alt[1])!=1:
					continue
				else:
					try:
						info[0]=str(int(snp_decode[pos][0])+1)+snp_decode[pos][1]+str(int(snp_decode[pos][2])+1)
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
					except KeyError:
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
			elif info[0] == "0/1":
				if len(ref)!=1 or len(alt[0])!=1:
					continue
				else:
					try: 
						info[0]=snp_decode[pos]
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
					except KeyError:
						print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
			elif info[0]=="1/1":
				if len(ref)!=1 or len(alt[0])!=1:
					continue
				info[0] = "1|1"
				print >>f_out,"%s\t.\t%s\t%s" %("\t".join(line1[:7]),line1[8],":".join(info))
		else:
			break

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='run Prism')
	parser.add_argument("--sample", dest='sample',help="sample name")
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	parser.add_argument("--analysis_dir", dest='analysis_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/',help="directory for whole genome sequencing file")

	args = parser.parse_args()
	SAMPLE = args.sample

	affy_posterior=[]
	None_posterior=[]
	for i in range(1,23):
		print 'chr',i
		if(i==23):
			i='X'
		info_file = '%s.chr%s.info' %(args.sample,i)
		id_to_block=get_id_to_block(info_file)
		block_file = '%s.chr%s.local.blocks.all' %(args.sample,i)
		block_to_snp=get_block_snp(block_file)
		mec_file = '../../BLOCK/%s_gvcf/%s_refhap_track_info_chr%s.txt' %(args.sample,args.sample,i)
		block_to_mec=get_block_mec(mec_file)
		merge_phase(info_file,i,id_to_block,block_to_snp,block_to_mec)
#		phase_info(info_file,i,id_to_block,block_to_snp,block_to_mec)
		affy_posterior,None_posterior=draw_hist(i,affy_posterior,None_posterior)

	fig,ax1=plt.subplots(figsize=(7,7))
	bins = 20
#	ax1.hist(affy_posterior,bins,facecolor='green',alpha=0.5,histtype='bar',label='Affy')
	ax1.hist(None_posterior,bins,facecolor='blue',alpha=0.5,histtype='bar',label='None')
	ax1.set_xlabel('posterior')
	ax1.set_ylabel('count')
	plt.savefig(args.sample+'_heatmap.png',format='png',dpi=200)

	tot_switch_error = 0
	tot_block = 0
	for i in range(1,23):
		print 'chr',i
		CHROM='chr'+str(i)
		tot_switch_error,tot_block,phase=calc_switch_error(i,tot_switch_error,tot_block)
		snp_not_shown_file = args.analysis_dir +args.sample + "_gvcf_SNPs_not_shown"
		SNP_not_shown=SNPs_not_shown(snp_not_shown_file)
		haplotype_file = args.analysis_dir +args.sample+'_gvcf_snp_haplotype_detail'
		BLOCK,block_list,decode,mec=read_haplotype(haplotype_file,SNP_not_shown)
		snp_decode = read_hap(BLOCK,block_list,phase)
		wgs_vcf = args.wgs_dir + SAMPLE+ '/gVCF_calls/%s.%s.vcf.gz' %(SAMPLE,CHROM)
		write_vcf_v2(wgs_vcf,snp_decode)
#	print tot_switch_error,tot_block,float(tot_switch_error)/tot_block

	
		