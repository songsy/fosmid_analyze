#!/usr/bin/env python
# python make-psmc-from-macs.py
# Shiya Song
# 11 Jan 2015

import sys
import argparse
import random

def code_to_string(code):
	string = ''
	for i in range(len(code)):
		if code[i]=='0':
			string+='A'
		else:
			string+='T'
	return string

def have_het(string):
	a=list(string)
	zero = a.count('0')
	one = a.count('1')
	if zero==0 or one ==0:
		return False
	else:
		return True

def macs_to_record_single(f):
	pos = []
	prev_pos = 0
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			substring = col[4][0]+col[4][2]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos:
					pos.append(prev_pos+1)
					prev_pos +=1
				else:
					pos.append(int(int(args.L)*float(col[2])))
					prev_pos = int(int(args.L)*float(col[2]))
	return pos

def macs_to_record_unphased(f):
	pos = [0,0,0,0]  # normal, random,all_het, skip
	prev_pos = [0,0,0,0]
	mask_pos = []
	unknown_status = False
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			index = 0
			unphase = FALSE
			col = line.split()
			substring = col[4][0]+col[4][2]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos[index]:
					pos[index].append(prev_pos[index]+1)
					prev_pos[index] +=1
				else:
					pos[index].append(int(int(args.L)*float(col[2])))
					prev_pos[index] = int(int(args.L)*float(col[2]))

			s=random.uniform(0, 1)
    		if s<=args.freq:
    			unphase = TRUE
    			index +=1
    			# random mode 
    			a = random.randint(0,1)
    			b = random.randint(0,1)
    			substring = col[4][0+a]+col[4][2+b]
    			if have_het(substring) is True:
					if int(int(args.L)*float(col[2]))<=prev_pos[index]:
						pos[index].append(prev_pos[index]+1)
						prev_pos[index] +=1
					else:
						pos[index].append(int(int(args.L)*float(col[2])))
						prev_pos[index] = int(int(args.L)*float(col[2]))
				index +=1
    			# all het mode
				het = False
				for i in range(2):
					for j in range(2):
						substring = col[4][0+i]+col[4][2+j]
						if have_het(substring) is True:
							het = True
				if het:
					if int(int(args.L)*float(col[2]))<=prev_pos[index]:
						pos[index].append(prev_pos[index]+1)
						prev_pos[index] +=1
					else:
						pos[index].append(int(int(args.L)*float(col[2])))
						prev_pos[index] = int(int(args.L)*float(col[2]))
				index +=1
				# skip mode
				if int(int(args.L)*float(col[2]))<=prev_pos[index]:
					pos[index].append(prev_pos[index]+1)
					mask_pos.append(1)
					prev_pos[index] +=1
				else:
					pos[index].append(int(int(args.L)*float(col[2])))
					prev_pos[index] = int(int(args.L)*float(col[2]))
					mask_pos.append(1)
    		else:
    			substring = col[4][0]+col[4][2]
				if have_het(substring) is True:
					for index in range(1,4):
						if int(int(args.L)*float(col[2]))<=prev_pos[index]:
							pos[index].append(prev_pos[index]+1)
							prev_pos[index] +=1
							if index==3:
								mask_pos.append(0)
						else:
							pos[index].append(int(int(args.L)*float(col[2])))
							prev_pos[index] = int(int(args.L)*float(col[2]))
							if index==3:
								mask_pos.append(0)
	assert len(mask_pos)==len(pos[3])
	return pos,mask_pos

def macs_to_record_all(f):
	pos = [[],[],[],[],[],[]]
	prev_pos = [0,0,0,0,0,0]
	for line in f:
		line = line.strip()
		if line[:4]=='SITE':
			col = line.split()
			index = 0
			for j in range(int(args.nsample)):
				for hap in range(int(args.nhap)):
#					print j,hap+int(args.nhap)
					substring = col[4][j]+col[4][hap+int(args.nhap)]
					if have_het(substring) is True:
						if int(int(args.L)*float(col[2]))<=prev_pos[index]:
							pos[index].append(prev_pos[index]+1)
							prev_pos[index] +=1
						else:
							pos[index].append(int(int(args.L)*float(col[2])))
							prev_pos[index] = int(int(args.L)*float(col[2]))
					index +=1
			substring = col[4][0]+col[4][1]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos[4]:
					pos[4].append(prev_pos[4]+1)
					prev_pos[4] +=1
				else:
					pos[4].append(int(int(args.L)*float(col[2])))
					prev_pos[4] = int(int(args.L)*float(col[2]))
			substring = col[4][2]+col[4][3]
			if have_het(substring) is True:
				if int(int(args.L)*float(col[2]))<=prev_pos[5]:
					pos[5].append(prev_pos[5]+1)
					prev_pos[5] +=1
				else:
					pos[5].append(int(int(args.L)*float(col[2])))
					prev_pos[5] = int(int(args.L)*float(col[2]))
	return pos

def write_psmcfa_unphased(pos,f_out,index):
	hetChar = 'K'
	homChar = 'T'
	failChar = 'N'
	winLen = 100
	for comparison in range(len(pos)-1):
		outFile = f_out[comparison]
		outFile.write('>%s\n' % (index))
		windowStart = 0
		windowEnd = windowStart + winLen - 1
		windows = []
		while windowStart < int(args.L):
			windows.append(homChar)
			windowStart = windowEnd + 1
			windowEnd = windowStart + winLen - 1			  
		#now, drop in the het sites into the windows
		for p in pos[comparison]:
			wID = int(p/winLen)
			try:
				windows[wID] = hetChar
			except IndexError:
				print "Error:",wID,p,winLen,len(windows)
				windows[wID-1] = hetChar
		numHet = 0
		numHom = 0
		for i in range(len(windows)):
			outFile.write('%s' % (windows[i]))
			if windows[i] == hetChar:
				numHet += 1
			else:
				numHom += 1
			if (i+1) % 50 == 0:
				outFile.write('\n')
		outFile.write('\n')
		print 'Record %i Comparison %i position %i windows %i homWind %i hetWind %i' % (index,comparison,len(pos[comparison]),len(windows),numHom,numHet)
	# last comparison: with skip mode
	outFile = f_out[-1]
	outFile.write('>%s\n' % (index))
	windowStart = 0
	windowEnd = windowStart + winLen - 1
	windows = []
	while windowStart < int(args.L):
		windows.append(homChar)
		windowStart = windowEnd + 1
		windowEnd = windowStart + winLen - 1			  
		#now, drop in the het sites into the windows
	pID = 0
	start = TRUE
	while pID < len(pos[-1]):
		if mask_pos[pID]==1:
			if start:
				start_ID = pID
				end_ID = pID
				start = False
			end_ID = pID
		else:
			if start = False:
				start = True
				if start_ID!=0:
					wID=int(pos[-1][start_ID-1]/winLen)
					leftover = pos[-1][start_ID-1]%winLen
					if leftover>10:
						try:
							windows[wID] = hetChar
						except IndexError:
							print "Error:",wID,pID,winLen,len(windows)
							windows[wID-1] = hetChar
				if end_ID != len(pos[-1])-1:
					wID=int(pos[-1][end_ID+1]/winLen)
					leftover = winLen-pos[-1][start_ID-1]%winLen
					if leftover>10:
						try:
							windows[wID] = hetChar
						except IndexError:
							print "Error:",wID,pID,winLen,len(windows)
							windows[wID-1] = hetChar
				for i in range(start_ID,end_ID+1):
					wID=int(pos[-1][i]/winLen)
					if windows[wID]==homChar:
						windows[wID]=failChar
	numHet = 0
	numHom = 0
	numN = 0
	for i in range(len(windows)):
		outFile.write('%s' % (windows[i]))
		if windows[i] == hetChar:
			numHet += 1
		elif windows[i] == failChar:
			numN += 1
		else:
			numHom +=1
		if (i+1) % 50 == 0:
			outFile.write('\n')
	outFile.write('\n')
	print 'Record %i Comparison %i position %i windows %i homWind %i hetWind %i NWind %i' % (index,comparison,len(pos[comparison]),len(windows),numHom,numHet,numN)


def write_psmcfa_multiple(pos,f_out,index):
	hetChar = 'K'
	homChar = 'T'
	failChar = 'N'
	winLen = 100
	for comparison in range(len(pos)):
		outFile = f_out[comparison]
		outFile.write('>%s\n' % (index))
		windowStart = 0
		windowEnd = windowStart + winLen - 1
		windows = []
		while windowStart < int(args.L):
			windows.append(homChar)
			windowStart = windowEnd + 1
			windowEnd = windowStart + winLen - 1			  
		#now, drop in the het sites into the windows
		for p in pos[comparison]:
			wID = int(p/winLen)
			try:
				windows[wID] = hetChar
			except IndexError:
				print "Error:",wID,p,winLen,len(windows)
				windows[wID-1] = hetChar
		numHet = 0
		numHom = 0
		for i in range(len(windows)):
			outFile.write('%s' % (windows[i]))
			if windows[i] == hetChar:
				numHet += 1
			else:
				numHom += 1
			if (i+1) % 50 == 0:
				outFile.write('\n')
		outFile.write('\n')
		print 'Record %i Comparison %i position %i windows %i homWind %i hetWind %i' % (index,comparison,len(pos[comparison]),len(windows),numHom,numHet)

def write_psmcfa_single(pos,f_out,index):
	hetChar = 'K'
	homChar = 'T'
	failChar = 'N'
	winLen = 100
	if True:
		outFile = f_out
		outFile.write('>%s\n' % (index))
		windowStart = 0
		windowEnd = windowStart + winLen - 1
		windows = []
		while windowStart < int(args.L):
			windows.append(homChar)
			windowStart = windowEnd + 1
			windowEnd = windowStart + winLen - 1			  
		#now, drop in the het sites into the windows
		for p in pos:
			wID = int(p/winLen)
			try:
				windows[wID] = hetChar
			except IndexError:
				print "Error:",wID,p,winLen,len(windows)
				windows[wID-1] = hetChar
		numHet = 0
		numHom = 0
		for i in range(len(windows)):
			outFile.write('%s' % (windows[i]))
			if windows[i] == hetChar:
				numHet += 1
			else:
				numHom += 1
			if (i+1) % 50 == 0:
				outFile.write('\n')
		outFile.write('\n')
		print 'Record %i Comparison %i position %i windows %i homWind %i hetWind %i' % (index,1,len(pos),len(windows),numHom,numHet)

def read_macs_unphase(macs_file):
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			f_out_prefix = macs_file+"_"+str(index)+"_hetFQ." 
			f_out = []
			index +=1
			for j in range(4):  # org, random, all_het,skip
				f_out.append(open(f_out_prefix+str(j),"w"))		
			f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_unphased(f)
			print 'len:',len(pos)
			write_psmcfa_unphased(pos,f_out,i % int(args.nchr))
		else:
			f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_unphased(f)
			print 'len:',len(pos)
			write_psmcfa_unphased(pos,f_out,i % int(args.nchr))

def read_macs_single(macs_file):
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			f_out_prefix = macs_file+"_"+str(index)+"_hetFQ.0" 
			f_out = open(f_out_prefix,"w")
			index +=1
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_v2(f)
			write_psmcfa_v2(pos,f_out,i % int(args.nchr))
		else:
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_v2(f)
			write_psmcfa_v2(pos,f_out,i % int(args.nchr))

def read_macs_multiple(macs_file):  #diploid+4*pseudo-diploid
	index = 0
	start = True
	for i in range(int(args.nrep)):
		if i % int(args.nchr) == 0:
			f_out_prefix = macs_file+"_"+str(index)+"_hetFQ." 
			f_out = []
			index +=1
			for j in range(6):
				f_out.append(open(f_out_prefix+str(j),"w"))
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_v3(f)
			write_psmcfa(pos,f_out,i % int(args.nchr))
		else:
			if args.version:
				f = open(macs_file+"_"+str(i)+"_sites_%s.txt" %(args.version),'r')
			else:
				f = open(macs_file+"_"+str(i)+"_sites.txt",'r')
			pos = macs_to_record_v3(f)
			write_psmcfa(pos,f_out,i % int(args.nchr))

if __name__=="__main__":

	parser = argparse.ArgumentParser(description='macs to psmcfa file')
	parser.add_argument("--file", dest='file',help="macs file prefix")
	parser.add_argument("--nrep", default=1000,dest='nrep',help="number of repetition")
	parser.add_argument("--nchr", default=100,dest='nchr',help="number of chromosome")
	parser.add_argument("--nhap", default=2,dest='nhap',help="number of haplotype per individual")
	parser.add_argument("--nsample", default=2,dest='nsample',help="number of sample")
	parser.add_argument("--L", default=3e07,dest='L',help="length of fragments")
	parser.add_argument("--version",default=None,dest='version')
	parser.add_argument("--mode",default=None,dest='mode')
	parser.add_argument("--freq",default=0.1,dest='freq')

	args = parser.parse_args()
	if args.mode=='unphase':
		read_macs_unphase(args.file)
	else:
		read_macs_v2(args.file)
	