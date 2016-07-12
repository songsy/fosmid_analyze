#!/usr/bin/env python
# python phasing_compare.py
# Shiya Song
# 23 Sept 2014
# Comparing phasing statistics in terms of switch error and mean distance between switch errors

import sys
from NGS_utils import *
import argparse
import pandas as pd
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import pickle
import os
def create_snp_list(file1):
	snp_pos = []
	snp_geno = []
	snp_phase = []
	snp_phase2 = []
	for i in VcfIterator(file1):
		if i[4] is True and i[3][0]!=i[3][1]: 
			snp_pos.append(i[1])				# store the position of the snp, used for fast search 
			snp_geno.append(i[2])		# store the ref and alt allele of the snp
			snp_phase.append(i[3])			# store the phasing information, like ["0011","TTGG","TG"]
			snp_phase2.append(())
	snp = {}
	snp["geno"] = pd.Series(snp_geno,index=snp_pos)
	snp["phase1"] = pd.Series(snp_phase,index=snp_pos)
	snp["phase2"] = pd.Series(snp_phase2,index=snp_pos)
	snp = pd.DataFrame(snp)
	return snp

def read_snp(SNP,file2):		# read in hapmap phasing after liftOver
	for i in VcfIterator(file2):
		if i[4] is True and i[3][0]!=i[3][1]: 
			try:
				info = SNP["geno"][i[1]]
			except KeyError:
				continue
			if SNP["geno"][i[1]]!=i[2]:
				phase = []
				allele = i[2][int(i[3][0])]+i[2][int(i[3][1])]
				if info[0]+info[1]==allele[0]+allele[1] or info[0]+info[1]==allele[1]+allele[0]:
					if i[3][0]<2:
						phase.append(i[3][0])
					elif i[3][0]==2:
						phase.append(1)
					if i[3][1]<2:
						phase.append(i[3][1])
					elif i[3][0]==2:
						phase.append(1)
					SNP["phase2"][i[1]]=tuple(phase)
				else:
					print i,SNP["geno"][i[1]]
			else:
				SNP["phase2"][i[1]]=i[3]
	return SNP

def compare_phasing(SNP,switch_distance,switch_length):
	list = SNP.index
	tot = 0
	opposite = 0
	is_switch = []
	switch = []
	switch_position = []
	for j in range(len(list)):
		if SNP["phase1"][list[j]]!=() and SNP["phase2"][list[j]]!=():
			tot +=1
			if SNP["phase1"][list[j]]!=SNP["phase2"][list[j]]:
#				print list[j],SNP["phase1"][list[j]],SNP["phase2"][list[j]]
				opposite +=1
				is_switch.append(1)
				switch.append(0)
				switch_position.append(list[j])
			else:
				is_switch.append(0)
				switch.append(0)
				switch_position.append(list[j])
	k = 0
	s = 0
	e = 0
	switch_error = 0
	last_error_pos = None
	if is_switch.count(0)>=is_switch.count(1):
		for k in range(len(is_switch)):
			if k==0:
				last = is_switch[k]
				s = k
			elif is_switch[k]==last:
				continue
			else:
				e = k
				if is_switch[s]==1:
					switch_error +=1
					switch_length.append(int(switch_position[e-1])-int(switch_position[s]))
					if last_error_pos is None:
						last_error_pos = switch_position[s]
					else:
						switch_distance.append(switch_position[s]-last_error_pos)
						last_error_pos = switch_position[s]
				s = k
				last = is_switch[k]
		e = k
		if is_switch[s]==1:
			switch_error +=1
			switch_length.append(int(switch_position[e-1])-int(switch_position[s]))
			if last_error_pos is None:
				last_error_pos = switch_position[s]
			else:
				switch_distance.append(switch_position[s]-last_error_pos)
				last_error_pos = switch_position[s]
	else:
		for k in range(len(is_switch)):
			if k==0:
				last = is_switch[k]
				s = k
			elif is_switch[k]==last:
				continue
			else:
				e = k
				if is_switch[s]==0:
					switch_error +=1
					switch_length.append(int(switch_position[e-1])-int(switch_position[s]))
					if last_error_pos is None:
						last_error_pos = switch_position[s]
					else:
						switch_distance.append(switch_position[s]-last_error_pos)
						last_error_pos = switch_position[s]
				s = k
				last = is_switch[k]
		e = k
		if is_switch[s]==0:
			switch_error +=1
			switch_length.append(int(switch_position[e-1])-int(switch_position[s]))
			if last_error_pos is None:
				last_error_pos = switch_position[s]
			else:
				switch_distance.append(switch_position[s]-last_error_pos)
				last_error_pos = switch_position[s]
	print chr,tot,opposite,switch_error,len(switch_distance),len(switch_length),np.mean(switch_distance),np.mean(switch_length)
	return switch_distance,switch_length,tot,opposite,switch_error

def compare_phasing_v2(SNP,switch_point_length,switch_distance,switch_length):
	list = SNP.index
	tot = 0
	opposite = 0
	is_switch = []
	switch = []
	switch_position = []
	for j in range(len(list)):
		if SNP["phase1"][list[j]]!=() and SNP["phase2"][list[j]]!=():
			tot +=1
			if SNP["phase1"][list[j]]!=SNP["phase2"][list[j]]:
#				print list[j],SNP["phase1"][list[j]],SNP["phase2"][list[j]]
				opposite +=1
				is_switch.append(1)
				switch.append(0)
				switch_position.append(list[j])
			else:
				is_switch.append(0)
				switch.append(0)
				switch_position.append(list[j])
	k = 0
	s = 0
	e = 0
	switch_error = 0
	last_error_pos = None
	start = -1
	if is_switch.count(0)>=is_switch.count(1):
		for k in range(len(is_switch)):
			if is_switch[k]==1:
				if start ==-1:
					start = k	
			else:
				if start !=-1:
					end = k
					switch_error +=1
#					print start,end-1,switch_position[end-1],switch_position[start],last_error_pos
					switch_length.append(int(switch_position[end-1])-int(switch_position[start]))
					if end==start:
						switch_point_length.append(1)
					else:
						switch_point_length.append(end-start)
					if last_error_pos is None:
						last_error_pos = int(switch_position[start])
					else:
						switch_distance.append(int(switch_position[start])-last_error_pos)
						last_error_pos = int(switch_position[start])
					start = -1
		if is_switch[k]==1:
			if start !=-1:
				end = k
				switch_error +=1
				switch_length.append(int(switch_position[end])-int(switch_position[start]))
				if end==start:
					switch_point_length.append(1)
				else:
					switch_point_length.append(end-start)
				if last_error_pos is None:
					last_error_pos = int(switch_position[start])
				else:
					switch_distance.append(int(switch_position[start])-last_error_pos)
					last_error_pos = int(switch_position[start])
	else:
		for k in range(len(is_switch)):
			if is_switch[k]==0:
				if start ==-1:
					start = k	
			else:
				if start !=-1:
					end = k
#					print start,end-1,switch_position[end-1],switch_position[start],last_error_pos
					switch_error +=1
					switch_length.append(int(switch_position[end-1])-int(switch_position[start]))
					if end==start:
						switch_point_length.append(1)
					else:
						switch_point_length.append(end-start)
					if last_error_pos is None:
						last_error_pos = int(switch_position[start])
					else:
						switch_distance.append(int(switch_position[start])-last_error_pos)
						last_error_pos = int(switch_position[start])
					start = -1
		if is_switch[k]==0:
			if start !=-1:
				end = k
				switch_error +=1
				switch_length.append(int(switch_position[end])-int(switch_position[start]))
				if end==start:
					switch_point_length.append(1)
				else:
					switch_point_length.append(end-start)
				if last_error_pos is None:
					last_error_pos = int(switch_position[start])
				else:
					switch_distance.append(int(switch_position[start])-last_error_pos)
					last_error_pos = int(switch_position[start])
	print chr,tot,opposite,switch_error,len(switch_distance),len(switch_length),len(switch_point_length),np.mean(switch_distance),np.mean(switch_length),sum(switch_length)/(len(switch_length)-switch_point_length.count(1)),np.mean(switch_point_length),switch_point_length.count(1)
	return switch_point_length,switch_distance,switch_length,tot,opposite,switch_error

def compare_phasing_v3(block_file,SNP,f_out):   # ignore flip error when calculating distance, and mean length
	list = SNP.index
	if os.path.isfile(block_file) is False:
		block_index=[0 for i in list]
		block_list_left = [0]
		block_list_right = [0]
	else:
		block_list_left = []
		block_list_right = []
		for i in BedIterator(block_file):
			if len(block_list_left)==0:
				block_list_left.append(int(i[1]))
				block_list_right.append(int(i[2]))
			elif int(i[1])>block_list_right[-1]:
				block_list_left.append(int(i[1]))
				block_list_right.append(int(i[2]))
		for i in range(1,len(block_list_left)):
			assert block_list_left[i]>block_list_left[i-1]
		list = SNP.index
		block_index = []
		for j in range(len(list)):
			find = find_interval(list[j],block_list_left)
			try:
				assert list[j]>=block_list_left[find] and list[j]<=block_list_right[find]
				block_index.append(find)
			except AssertionError:
#				print 'Assertion Error',list[j],block_list_left[find],block_list_right[find]
				block_index.append(-1)
		print '-1:',block_index.count(-1)
	list = SNP.index
	tot = 0
	opposite = 0
	is_switch = []
	label = []
	switch = []
	switch_position = []
	for j in range(len(list)):
		if block_index[j]==-1:
			continue
		if SNP["phase1"][list[j]]!=() and SNP["phase2"][list[j]]!=():
			tot +=1
			if SNP["phase1"][list[j]]!=SNP["phase2"][list[j]]:
#				print list[j],SNP["phase1"][list[j]],SNP["phase2"][list[j]]
				opposite +=1
				is_switch.append(1)
				switch.append(0)
				switch_position.append(list[j])
				label.append(block_index[j])
			else:
				is_switch.append(0)
				switch.append(0)
				switch_position.append(list[j])
				label.append(block_index[j])
	switch_error,switch_point_length,switch_distance,switch_length,phase_length=calc_phase_error(is_switch,switch_position,label,block_list_left,block_list_right,f_out)
	print chr,tot,opposite,switch_error,len(switch_distance),len(switch_length),len(switch_point_length),len(phase_length),np.mean(switch_distance),np.mean(switch_length),np.mean(switch_point_length),switch_point_length.count(1),np.mean(phase_length)
	return switch_point_length,switch_distance,switch_length,phase_length,tot,opposite,switch_error

def find_interval(pos,list):
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
		if pos<list[midpoint] and midpoint!=0:
			midpoint=midpoint-1
	return midpoint

def calc_phase_error(is_switch,switch_position,label,block_list_left,block_list_right,f_out):   # also ignore flips
	switch_distance = []   # record distance between switch errors
	switch_length = []  	# record the length of switch errors
	switch_point_length =[]   # record the number of discordance in each switch
	phase_length = []
	k = 0
	s = 0
	e = 0
	last_error_pos = None
	last_end_pos = None
	start = -1
	switch_error = 0
	switch_pos=[]
	extra = 0
	if len(is_switch)==0:
		return switch_error,switch_point_length,switch_distance,switch_length,phase_length
	if is_switch.count(0)>=is_switch.count(1):
		for k in range(len(is_switch)):
			if is_switch[k]==1:
				if start ==-1:
					cur_label=label[k]
					switch_pos=[switch_position[k]]
					start = k
					'''
					if last_end_pos is None:
						if k!=0:
							extra +=1
							print 'extra:',start,k,is_switch[k]
							phase_length.append(int(switch_position[start])-int(switch_position[0]))
					else:
						extra +=1
						print 'extra:',start,k,is_switch[k],switch_position[start],last_end_pos
						phase_length.append(int(switch_position[start])-last_end_pos)
					last_end_pos=int(switch_position[start])
					'''
				else:
					if label[k]==cur_label:
						switch_pos.append(switch_position[k])
					else:
						if len(switch_pos)!=0:
							print >>f_out,'%s\t%s\t%s\t%s' %(chr,block_list_left[cur_label],block_list_right[cur_label],",".join(map(str,switch_pos)))
						cur_label=label[k]
						switch_pos=[switch_position[k]]
			else:
				if start !=-1:
					end = k
					switch_error +=1
#					print start,end-1,switch_position[end-1],switch_position[start],last_error_pos
					if int(switch_position[end-1])-int(switch_position[start])<0:
						print start,end,is_switch,switch_position
					if end!=start+1:
						switch_point_length.append(end-start)
						switch_length.append(int(switch_position[end-1])-int(switch_position[start]))
						phase_length.append(int(switch_position[end-1])-int(switch_position[start]))
						if last_end_pos is None:
							extra +=1
							phase_length.append(int(switch_position[start])-int(switch_position[0]))
						else:
							extra +=1
							phase_length.append(int(switch_position[start])-last_end_pos)
						last_end_pos=int(switch_position[end-1])
						if last_error_pos is None:
							last_error_pos = int(switch_position[start])
						else:
							if int(switch_position[start])-last_error_pos<0:
								print start,end,is_switch,switch_position,last_error_pos
							switch_distance.append(int(switch_position[start])-last_error_pos)
							last_error_pos = int(switch_position[start])
					else:
						switch_point_length.append(end-start)
					start = -1
		if is_switch[k]==1:
			if start !=-1:
				end = k
				switch_error +=1
				if end!=start:
					switch_point_length.append(end-start+1)
					switch_length.append(int(switch_position[end])-int(switch_position[start]))
					phase_length.append(int(switch_position[end-1])-int(switch_position[start]))
					if last_end_pos is None:
						extra +=1
						phase_length.append(int(switch_position[start])-int(switch_position[0]))
					else:
						extra +=1
						phase_length.append(int(switch_position[start])-last_end_pos)
					if last_error_pos is None:
						last_error_pos = int(switch_position[start])
					else:
						if int(switch_position[start])-last_error_pos<0:
							print 'last',start,end,is_switch,switch_position,last_error_pos
						switch_distance.append(int(switch_position[start])-last_error_pos)
						last_error_pos = int(switch_position[start])
				else:
					switch_point_length.append(1)
		if len(switch_pos)!=0:
			print >>f_out,'%s\t%s\t%s\t%s' %(chr,block_list_left[cur_label],block_list_right[cur_label],",".join(map(str,switch_pos)))
	else:
		for k in range(len(is_switch)):
			if is_switch[k]==0:
				if start ==-1:
					cur_label=label[k]
					switch_pos=[switch_position[k]]
					start = k
					'''
					if last_end_pos is None:
						if k!=0:
							extra +=1
							phase_length.append(int(switch_position[start])-int(switch_position[0]))
					else:
						extra +=1
						phase_length.append(int(switch_position[start])-last_end_pos)
					last_end_pos=int(switch_position[start])
					'''
				else:
					if label[k]==cur_label:
						switch_pos.append(switch_position[k])
					else:
						if len(switch_pos)!=0:
							print >>f_out,'%s\t%s\t%s\t%s' %(chr,block_list_left[cur_label],block_list_right[cur_label],",".join(map(str,switch_pos)))
						cur_label=label[k]
						switch_pos=[switch_position[k]]
			else:
				if start !=-1:
					end = k
#					print start,end-1,switch_position[end-1],switch_position[start],last_error_pos
					switch_error +=1
					if int(switch_position[end-1])-int(switch_position[start])<0:
						print start,end,is_switch,switch_position
					if end!=start+1:
						switch_point_length.append(end-start)
						switch_length.append(int(switch_position[end-1])-int(switch_position[start]))
						phase_length.append(int(switch_position[end-1])-int(switch_position[start]))
						if last_end_pos is None:
							extra +=1
							phase_length.append(int(switch_position[start])-int(switch_position[0]))
						else:
							extra +=1
							phase_length.append(int(switch_position[start])-last_end_pos)
						last_end_pos=int(switch_position[end-1])
						if last_error_pos is None:
							last_error_pos = int(switch_position[start])
						else:
							if int(switch_position[start])-last_error_pos<0:
								print start,end,is_switch,switch_position,last_error_pos
							switch_distance.append(int(switch_position[start])-last_error_pos)
							last_error_pos = int(switch_position[start])
					else:
						switch_point_length.append(end-start)
					start = -1
		if is_switch[k]==0:
			if start !=-1:
				end = k
				switch_error +=1
				if end!=start:
					switch_point_length.append(end-start+1)
					switch_length.append(int(switch_position[end])-int(switch_position[start]))
					phase_length.append(int(switch_position[end-1])-int(switch_position[start]))
					if last_end_pos is None:
						extra +=1
						phase_length.append(int(switch_position[start])-int(switch_position[0]))
					else:
						extra +=1
						phase_length.append(int(switch_position[start])-last_end_pos)
					if last_error_pos is None:
						last_error_pos = int(switch_position[start])
					else:
						if int(switch_position[start])-last_error_pos<0:
							print 'last',start,end,is_switch,switch_position,last_error_pos
						switch_distance.append(int(switch_position[start])-last_error_pos)
						last_error_pos = int(switch_position[start])
				else:
					switch_point_length.append(1)
		if len(switch_pos)!=0:
			print >>f_out,'%s\t%s\t%s\t%s' %(chr,block_list_left[cur_label],block_list_right[cur_label],",".join(map(str,switch_pos)))
	print 'length',extra,len(switch_length),len(phase_length)
	return switch_error,switch_point_length,switch_distance,switch_length,phase_length

def compare_phasing_prism(prism_file,SNP,switch_point_length,switch_distance,switch_length,phase_length,f_out):
	block_list_left = []
	block_list_right = []
	for i in BedIterator(prism_file):
		if i[4]!='interleaving' and i[4]!='None,interleave' and i[1]!=i[2]:
			block_list_left.append(int(i[1]))
			block_list_right.append(int(i[2]))
	for i in range(1,len(block_list_left)):
		assert block_list_left[i]>block_list_left[i-1]
	list = SNP.index
	block_index = []
	for j in range(len(list)):
		find = find_interval(list[j],block_list_left)
		try:
			assert list[j]>=block_list_left[find] and list[j]<=block_list_right[find]
			block_index.append(find)
		except AssertionError:
#			print list[j],find,block_list_left[find],block_list_right[find]
			block_index.append(-1)
#	print len(block_index),len(list)
#	print block_index
	Switch_error = 0
	Tot = 0
	Opposite=0
	start = True
	for j in range(len(list)):
		if block_index[j]==-1:
			continue
		if start is True:
			tot = 0
			opposite = 0
			is_switch = []
			label = []
			switch = []
			switch_position = []
			start = False
			prev = block_index[j]
		if block_index[j]!=prev:
			if tot>1:
				switch_error,sswitch_point_length,sswitch_distance,sswitch_length,sphase_length=calc_phase_error(is_switch,switch_position,label,block_list_left,block_list_right,f_out)
				if switch_error!=0:
					print block_index[j],tot,opposite,switch_error,len(sswitch_distance),len(sswitch_length),len(sswitch_point_length),len(sphase_length),np.mean(sswitch_distance),np.mean(sswitch_length),np.mean(sswitch_point_length),sswitch_point_length.count(1),np.mean(sphase_length)
				else:
					print block_index[j],tot,opposite,switch_error
				Tot+=tot
				Opposite+=min(opposite,tot-opposite)
				Switch_error+=switch_error
				switch_distance=switch_distance+sswitch_distance
				switch_length=switch_length+sswitch_length
				switch_point_length=switch_point_length+sswitch_point_length
				phase_length=phase_length+sphase_length
			tot = 0
			opposite = 0
			is_switch = []
			switch = []
			label = []
			switch_position = []
			start = False
			prev = block_index[j]
			continue
		if SNP["phase1"][list[j]]!=() and SNP["phase2"][list[j]]!=():
			tot +=1
			if SNP["phase1"][list[j]]!=SNP["phase2"][list[j]]:
#				print list[j],SNP["phase1"][list[j]],SNP["phase2"][list[j]]
				opposite +=1
				is_switch.append(1)
				switch.append(0)
				switch_position.append(list[j])
				label.append(block_index[j])
			else:
				is_switch.append(0)
				switch.append(0)
				switch_position.append(list[j])
				label.append(block_index[j])
	if tot>1:
		switch_error,sswitch_point_length,sswitch_distance,sswitch_length,sphase_length=calc_phase_error(is_switch,switch_position,label,block_list_left,block_list_right,f_out)
		if switch_error!=0:
			print block_index[j],tot,opposite,switch_error,len(sswitch_distance),len(sswitch_length),len(sswitch_point_length),len(sphase_length),np.mean(sswitch_distance),np.mean(sswitch_length),np.mean(sswitch_point_length),sswitch_point_length.count(1),np.mean(sphase_length)
		else:
			print block_index[j],tot,opposite,switch_error
		Tot+=tot
		Opposite+=min(opposite,tot-opposite)
		Switch_error+=switch_error
		switch_distance=switch_distance+sswitch_distance
		switch_length=switch_length+sswitch_length
		switch_point_length=switch_point_length+sswitch_point_length
		phase_length=phase_length+sphase_length
	print chr,Tot,Opposite,Switch_error,len(switch_distance),len(switch_length),len(switch_point_length),len(phase_length),np.mean(switch_distance),np.mean(switch_length),np.mean(switch_point_length),switch_point_length.count(1),(sum(switch_point_length)-switch_point_length.count(1))/(len(switch_point_length)-switch_point_length.count(1)),np.mean(phase_length)
	return switch_point_length,switch_distance,switch_length,phase_length,Tot,Opposite,Switch_error


if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--name",help="Input VCF files")
	parser.add_argument("--wgs_dir", dest='wgs_dir',default='/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/',help="directory for whole genome sequencing file")
	parser.add_argument("--phase_panel", dest='phase_panel',default='1KGphase1',help="reference panel")
	parser.add_argument("--prism", dest='prism',default=0,help="prism or not")
	parser.add_argument("--phase", dest='phase',default=0,help="prism or not")
	args = parser.parse_args()
	
	Switch_distance = []   # record distance between switch errors
	Switch_length = []  	# record the length of switch errors
	Switch_point_length =[]   # record the number of discordance in each switch
	Phase_length = []
	tot = 0
	different = 0
	switch_error = 0
	if args.phase=='shapeit':
		f_out=open('%s_shapeit_switch_error.txt' %(args.name),'w')
	else:
		f_out=open('%s_1KGphase3_switch_error.txt' %(args.name),'w')
	for chr in range(1,23):
		switch_distance = []   # record distance between switch errors
		switch_length = []  	# record the length of switch errors
		switch_point_length =[]   # record the number of discordance in each switch
		phase_length = []
		chr = 'chr'+str(chr)
		if args.prism=='1':
			file1 = '%s%s/gVCF_calls/%s.%s.prism.v2.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
		elif args.name=='NA12878':
			file1 = '%s%s/all_sites/%s.%s.fosmid.phase.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
		else:
			file1 = '%s%s/gVCF_calls/%s.%s.fosmid.v2.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
		if args.phase=='shapeit':
			if args.phase_panel=='1KGphase1':
				file2 = '%s%s/gVCF_calls/%s.%s.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
				if args.name=='NA12878':
					file2 = '%s%s/all_sites/%s.%s.phased.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
			elif args.phase_panel=='1KGphase3':
				file2 = '%s%s/gVCF_calls/%s.%s.1KGphase3.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
				if args.name=='NA12878':
					file2 = '%s%s/all_sites/%s.%s.1KGphase3.vcf.gz' %(args.wgs_dir,args.name,args.name,chr)
		else:
			file2 = '/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase3/%s/%s.%s.1KGphase3.snp.vcf.gz' %(args.name,args.name,chr)
		print file1,file2
		SNP = create_snp_list(file1)
		SNP = read_snp(SNP,file2)
		if args.prism=='1':
			prism_file = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/Prism/%s/%s.%s.info' %(args.name,args.name,chr)
			switch_point_length,switch_distance,switch_length,phase_length,a,b,c=compare_phasing_prism(prism_file,SNP,switch_point_length,switch_distance,switch_length,phase_length,f_out)		
		else:
			block_file = '/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/analysis/BLOCK/%s_gvcf/%s_refhap_track_info_%s.txt' %(args.name,args.name,chr)
			switch_point_length,switch_distance,switch_length,phase_length,a,b,c=compare_phasing_v3(block_file,SNP,f_out)
		if b<a-b:
			different +=b
		else:
			different += a-b
		'''	
		fig = plt.figure()	
		ax1 = fig.add_subplot(311,title=args.name+'_'+chr)
		ax2 = fig.add_subplot(312)
		ax3 = fig.add_subplot(313)
		bins = 40
		ax1.hist(switch_distance,bins,facecolor='green',alpha=0.5,histtype='stepfilled',label='library1')
		ax1.legend(prop={'size':10})
		ax1.set_xlabel('switch distance')
		ax1.set_ylabel('count')
#		fig.suptitle(sample_name, fontsize=10)
		ax2.hist(switch_length,bins=40,facecolor='green',alpha=0.75,histtype='stepfilled')
		ax2.set_xlabel('switch length')
		ax2.set_ylabel('count')
		ax3.hist(phase_length,bins=40,facecolor='green',alpha=0.75,histtype='stepfilled')
		ax3.set_xlabel('phase length')
		ax3.set_ylabel('count')
		fig.savefig(args.name+"_hist_1KGphase_%s.pdf" %(chr),format='pdf')
		'''
		tot +=a
		switch_error += c
		Switch_distance=Switch_distance+switch_distance
		Switch_length=Switch_length+switch_length
		Switch_point_length=Switch_point_length+switch_point_length
		Phase_length=Phase_length+phase_length
		singleton=Switch_point_length.count(1)
	dbfile = open('%s_%s_%s_SNP_pickle' %(args.name,args.phase,args.phase_panel),'wb')
	pickle.dump(Switch_distance,dbfile)
	pickle.dump(Switch_length,dbfile)
	pickle.dump(Switch_point_length,dbfile)
	pickle.dump(Phase_length,dbfile)
	print switch_error,len(Switch_distance),len(Switch_length),len(Switch_point_length),len(Phase_length)
	print tot,different,1-float(different)/tot,switch_error,float(switch_error)/tot,singleton,float(singleton)/tot,np.mean(Switch_distance),np.mean(Switch_length),np.mean(Switch_point_length),Switch_point_length.count(1),(sum(Switch_point_length)-Switch_point_length.count(1))/(len(Switch_point_length)-Switch_point_length.count(1)),np.mean(Phase_length)

			
