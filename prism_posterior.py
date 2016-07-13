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

def calc_switch_error_v2(chr,tot_switch_error,tot_block_wrong,tot_snp_wrong,tot_block,tot_snp):
	f='%s.chr%s.all.phase' %(args.sample,chr)
	count = [0,0]
	count_snp = [0,0]
	posterior = [[],[]]
	first = True
	switch_error = 0
	phase={}
	last_block = '1'
	is_switch = []
	switch_position_left = []
	switch_position_right = []
	for i in BedIterator(f):
		if i[0]==last_block:
			continue
		last_block = i[0]
		if i[7]!='None':
			count[int(i[1])]+=1
			count_snp[int(i[1])]+=int(i[8])
			is_switch.append(int(i[1]))
			switch_position_left.append(int(i[5]))
			switch_position_right.append(int(i[6]))
			if first:
				old = int(i[1])
				first = False
			else:
				if old!=int(i[1]):
					switch_error +=1
					old = int(i[1])
		phase[i[4]+" "+i[5]]=int(i[1])
	print count,count_snp,min(count),min(count_snp),switch_error-1
	tot_switch_error+=switch_error-1
	tot_block_wrong += min(count)
	tot_snp_wrong += min(count_snp)
	tot_block +=sum(count)
	tot_snp += sum(count_snp)
	return tot_switch_error,tot_block_wrong,tot_snp_wrong,tot_block,tot_snp,phase,is_switch,switch_position_left,switch_position_right

def calc_switch_error_v3(chr,tot_switch_error,tot_block_wrong,tot_snp_wrong,tot_block,tot_snp):
	f='%s.chr%s.all.phase' %(args.sample,chr)
	count = [0,0]
	count_snp = [0,0]
	posterior_range = [0,0.5,0.6,0.7,0.8,0.9,1]
	posterior_block = [[0,0] for x in range(len(posterior_range))]
	posterior_snp = [[0,0] for x in range(len(posterior_range))]
	first = True
	switch_error = 0
	phase={}
	last_block = '1'
	is_switch = []
	switch_position_left = []
	switch_position_right = []
	posterior = []
	snp_num = []
	switch_posterior = []
	switch = False
	for i in BedIterator(f):
		if i[0]==last_block:
			if len(posterior)!=0:
				posterior.pop()
				snp_num.pop() 
			continue
		last_block = i[0]
		if i[7]!='None':
			count[int(i[1])]+=1
			count_snp[int(i[1])]+=int(i[8])
			posterior.append(float(i[2]))
			snp_num.append(int(i[8]))
			for j in range(len(posterior_range)):
				if float(i[2])>=posterior_range[j]:
					posterior_block[j][int(i[1])]+=1
					posterior_snp[j][int(i[1])]+=int(i[8])
			is_switch.append(int(i[1]))
			switch_position_left.append(int(i[5]))
			switch_position_right.append(int(i[6]))
			if first:
				old = int(i[1])
				old_posterior = float(i[2])
				first = False
			else:
				if old!=int(i[1]):
					switch_error +=1
					old = int(i[1])
					if chr==22:
						print old_posterior
					switch_posterior.append(old_posterior)
					switch = True
				else:
					switch = False
				old_posterior = float(i[2])
		phase[i[4]+" "+i[5]]=int(i[1])
	flip = 1 if count_snp[1]>count_snp[0] else 0
	print count,count_snp,min(count),min(count_snp),switch_error-1,flip
	tot_switch_error+=switch_error-1
	tot_block_wrong += min(count)
	tot_snp_wrong += min(count_snp)
	tot_block +=sum(count)
	tot_snp += sum(count_snp)
	return switch_posterior,posterior,snp_num,posterior_block,posterior_snp,flip,tot_switch_error,tot_block_wrong,tot_snp_wrong,tot_block,tot_snp,phase,is_switch,switch_position_left,switch_position_right

def calc_switch_error_length(is_switch,switch_position_left,switch_position_right):
	start = -1
	switch_error = 0
	switch_length = []
	phase_length = []
	last_error_pos = None
	last_end_pos = None
	extra = 0
	for k in range(len(switch_position_left)-1):
		try:
			assert switch_position_left[k+1]>=switch_position_right[k]
		except AssertionError:
			print 'Error',switch_position_left[k+1],switch_position_right[k]
	if is_switch.count(0)>=is_switch.count(1):
		for k in range(len(is_switch)):
			if is_switch[k]==1:
				if start ==-1:
					start = k
			else:
				if start !=-1:
					end = k
					switch_error +=1
					if int(switch_position_right[end-1])!=int(switch_position_left[start]):
						switch_length.append(int(switch_position_right[end-1])-int(switch_position_left[start]))
						phase_length.append(int(switch_position_right[end-1])-int(switch_position_left[start]))
						if last_end_pos is None:
							extra +=1
							phase_length.append(int(switch_position_left[start])-int(switch_position_left[0]))
						else:
							extra +=1
							phase_length.append(int(switch_position_left[start])-last_end_pos)
						last_end_pos=int(switch_position_right[end-1])
					start = -1
		if is_switch[k]==1:
			if start !=-1:
				end = k
				switch_error +=1
				if int(switch_position_right[end-1])!=int(switch_position_left[start]):
					switch_length.append(int(switch_position_right[end])-int(switch_position_left[start]))
					phase_length.append(int(switch_position_right[end])-int(switch_position_left[start]))
					if last_end_pos is None:
						extra +=1
						phase_length.append(int(switch_position_left[start])-int(switch_position_left[0]))
					else:
						extra +=1
						phase_length.append(int(switch_position_left[start])-last_end_pos)
	else:
		for k in range(len(is_switch)):
			if is_switch[k]==0:
				if start ==-1:
					start = k
			else:
				if start !=-1:
					end = k
					switch_error +=1
					if int(switch_position_right[end-1])!=int(switch_position_left[start]):
						switch_length.append(int(switch_position_right[end-1])-int(switch_position_left[start]))
						phase_length.append(int(switch_position_right[end-1])-int(switch_position_left[start]))
						if last_end_pos is None:
							extra +=1
							phase_length.append(int(switch_position_left[start])-int(switch_position_left[0]))
						else:
							extra +=1
							phase_length.append(int(switch_position_left[start])-last_end_pos)
						last_end_pos=int(switch_position_right[end-1])
					start = -1
		if is_switch[k]==0:
			if start !=-1:
				end = k
				switch_error +=1
				if int(switch_position_right[end])!=int(switch_position_left[start]):
					switch_length.append(int(switch_position_right[end])-int(switch_position_left[start]))
					phase_length.append(int(switch_position_right[end])-int(switch_position_left[start]))
					if last_end_pos is None:
						extra +=1
						phase_length.append(int(switch_position_left[start])-int(switch_position_left[0]))
					else:
						extra +=1
						phase_length.append(int(switch_position_left[start])-last_end_pos)
	return switch_error,switch_length,phase_length
	
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
	'''
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
	'''
	tot_switch_error = 0
	tot_block = 0
	tot_block_wrong = 0
	tot_snp_wrong = 0
	tot_snp = 0
	phase_length = []
	Switch_error = 0
	switch_length = []
	posterior_range = [0,0.5,0.6,0.7,0.8,0.9,1]
	Posterior_block = [[0,0] for x in range(len(posterior_range))]
	Posterior_snp = [[0,0] for x in range(len(posterior_range))]
	Posterior=[]
	Snp_num = []
	Switch_posterior = []
	for i in range(1,23):
		print 'chr',i
		CHROM='chr'+str(i)
		switch_posterior,posterior,snp_num,posterior_block,posterior_snp,flip,tot_switch_error,tot_block_wrong,tot_snp_wrong,tot_block,tot_snp,phase,is_switch,switch_position_left,switch_position_right=calc_switch_error_v3(i,tot_switch_error,tot_block_wrong,tot_snp_wrong,tot_block,tot_snp)
		Posterior += posterior
		Snp_num += snp_num
		Switch_posterior += switch_posterior
		for j in range(len(Posterior_block)):
			Posterior_block[j][0]+=posterior_block[j][0] if flip==0 else posterior_block[j][1]
			Posterior_block[j][1]+=posterior_block[j][1] if flip==0 else posterior_block[j][0]
			Posterior_snp[j][0]+=posterior_snp[j][0] if flip==0 else posterior_snp[j][1]
			Posterior_snp[j][1]+=posterior_snp[j][1] if flip==0 else posterior_snp[j][0]
		'''
		switch_error,sswitch_length,sphase_length=calc_switch_error_length(is_switch,switch_position_left,switch_position_right)
		print switch_error,np.mean(sswitch_length),np.mean(sphase_length)
		Switch_error+=switch_error
		switch_length=switch_length+sswitch_length
		phase_length=phase_length+sphase_length
		'''
#	print tot_switch_error,tot_block,float(tot_switch_error)/tot_block,tot_block_wrong,tot_block,float(tot_block_wrong)/tot_block,tot_snp_wrong,tot_snp,float(tot_snp_wrong)/tot_snp
#	print Switch_error,np.mean(switch_length),np.mean(phase_length)
	for j in range(len(Posterior_block)):
		print j,posterior_range[j],float(Posterior_block[j][0])/sum(Posterior_block[j]),float(Posterior_snp[j][0])/sum(Posterior_snp[j])
	
	heatmap, xedges, yedges = np.histogram2d(Posterior, Snp_num, bins=[[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],[0,10,100,500,1000,max(Snp_num)]])
	extents = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	plt.clf()
	plt.imshow(heatmap.T,extent=extents,aspect="auto",interpolation='nearest',origin='lower',cmap='jet')
	pos = range(0,max(Snp_num),max(Snp_num)/5)
	if len(pos)==5:
		pos+=[max(Snp_num)]
	pylab.yticks(pos,["0","10","100","500","1000",str(max(Snp_num))])	
	plt.colorbar()
	plt.xlabel('posterior')
	plt.ylabel('number of snp')
	plt.title(args.sample+' posterior & number of snp heatmap')
	plt.savefig(args.sample+'_heatmap.pdf',format='pdf',dpi=200)

	print 'num of switch', len(Switch_posterior)
	fig,ax1=plt.subplots(figsize=(7,7))
	bins = 20
#	ax1.hist(affy_posterior,bins,facecolor='green',alpha=0.5,histtype='bar',label='Affy')
	ax1.hist(Switch_posterior,bins,facecolor='blue',alpha=0.5,histtype='bar',label='None')
	ax1.set_xlabel('posterior')
	ax1.set_ylabel('count')
	plt.savefig(args.sample+'_switch_posterior_heatmap.pdf',format='pdf',dpi=200)