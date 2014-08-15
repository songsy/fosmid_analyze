#!/usr/bin/env python
# python 1000G_trio_indel.py
# Shiya Song
# 20 Nov 2013
# check 1000G trio indel calls for NA19240,first convert vcf to bed format, and liftover to hg19 coordinate

import signal
import os
import pickle

def chrom_list():
	chr_list=[]
	for i in range(1,23):
		chr_list.append("chr"+str(i))
	chr_list.append("chrX")
	return chr_list

chr_list = chrom_list()

def vcf_to_bed(file):
	INDEL = []
	INDEL.append([])
	signal.signal(signal.SIGPIPE, signal.SIG_DFL)
	gc = 'gunzip -c ' + file
	f = os.popen(gc, 'r')
	line_number = 0
	for line in f:
		if line[0]=="#":
			continue
		line_number +=1
		line = line.strip().split("\t")
		chr = "chr"+line[0]
		pos = int(line[1])
		child=line[9].split(":")[0]
		father=line[10].split(":")[0]
		mother=line[11].split(":")[0]
		alt = line[4].split(",")
		if chr=="chrX":
			if child=="0/0" or child=="1/1" or child=="2/2":
				phase = child.replace("/","|")
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				if father=="0":
					if child=="0/1":
						phase="0|1"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					elif child=="0/2":
						phase=="0|2"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					elif child=="0/3":
						phase=="0|3"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					else:
						print line
				elif father =="1":
					if child=="0/1":
						phase="1|0"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					elif child=="1/2":
						phase="1|2"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					elif child=="1/3":
						phase="1|3"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					else:
						print line
				elif father == "2":
					if child=="0/2":
						phase="2|0"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					elif child=="1/2":
						phase="2|1"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					elif child=="2/3":
						phase="2|3"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					else:
						print line
				elif father == "3":
					if child=="0/3":
						phase="3|0"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					elif child=="1/3":
						phase="3|1"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					elif child=="2/3":
						phase="3|2"
						print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
						INDEL.append([line[3],line[4],phase])
					else:
						print line
				else:
					print line
			continue
		if child=="0/0" or child=="1/1" or child=="2/2" or child=="3/3" or child=="4/4":
			phase = child.replace("/","|")
			print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
			INDEL.append([line[3],line[4],phase])
		elif child=="0/1":
			if father =="0/1":
				if mother == "0/1":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="0/0" or mother=="0/2" or mother=="0/3" or mother=="0/4":
					phase="1|0"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="1/1" or mother =="1/2" or mother =="1/3" or mother=="1/4":
					phase="0|1"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="0/1":
				if father == "0/1":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="0/0" or father=="0/2" or father=="0/3" or father=="0/4":
					phase="0|1"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="1/1" or father =="1/2" or father =="1/3" or father=="1/4":
					phase="1|0"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="0/0" and (mother=="1/1" or mother=="1/2" or mother=="1/3" or mother=="1/4"):
				phase="0|1"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/1" or father=="1/2" or father=="1/3" or father=="1/4") and mother=="0/0":
				phase="1|0"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="0/2" or father=="0/3" or father=="0/4") and (mother=="1/2" or mother=="1/1"  or mother=="1/3" or mother=="1/4"):
				pahse="0|1"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/2" or father=="1/1" or father=="1/3" or father=="1/4") and (mother=="0/2" or mother=="0/3" or mother=="0/4"):
				pahse="1|0"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		elif child=="0/2":
			if father =="0/2":
				if mother == "0/2":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="0/0" or mother=="0/1" or mother=="0/3" or mother=="0/4":
					phase="2|0"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="2/2" or mother =="1/2" or mother =="2/3" or mother=="2/4":
					phase="0|2"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="0/2":
				if father == "0/2":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="0/0" or father=="0/1" or father=="0/3" or father=="0/4":
					phase="0|2"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="2/2" or father =="1/2" or father =="2/3" or father=="2/4":
					phase="2|0"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="0/0" and (mother=="1/2" or mother=="2/2" or mother=="2/3" or mother=="2/4"):
				phase="0|2"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/2" or father=="2/2" or father=="2/3" or father=="2/4") and mother=="0/0":
				phase="2|0"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="0/1" or father=="0/3" or father=="0/4") and (mother=="1/2" or mother=="2/2"  or mother=="2/3" or mother=="2/4"):
				pahse="0|2"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/2" or father=="2/2" or father=="2/3" or father=="2/4") and (mother=="0/1" or mother=="0/3" or mother=="0/4"):
				pahse="2|0"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		elif child=="0/3":
			if father =="0/3":
				if mother == "0/3":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="0/0" or mother=="0/1" or mother=="0/3" or mother=="0/4":
					phase="3|0"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="3/3" or mother =="1/3" or mother =="2/3" or mother=="3/4":
					phase="0|3"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="0/3":
				if father == "0/3":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="0/0" or father=="0/1" or father=="0/3" or father=="0/4":
					phase="0|3"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="3/3" or father =="1/3" or father =="2/3" or father=="3/4":
					phase="3|0"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="0/0" and (mother=="1/3" or mother=="3/3" or mother=="2/3" or mother=="3/4"):
				phase="0|3"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/3" or father=="3/3" or father=="2/3" or father=="3/4") and mother=="0/0":
				phase="3|0"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="0/1" or father=="0/2" or father=="0/4") and (mother=="1/3" or mother=="3/3"  or mother=="2/3" or mother=="3/4"):
				pahse="0|3"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/3" or father=="3/3" or father=="2/3" or father=="3/4") and (mother=="0/1" or mother=="0/2" or mother=="0/4"):
				pahse="3|0"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		elif child=="0/4":
			if father =="0/4":
				if mother == "0/4":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="0/0" or mother=="0/1" or mother=="0/3" or mother=="0/2":
					phase="4|0"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="4/4" or mother =="1/4" or mother =="2/4" or mother=="3/4":
					phase="0|4"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="0/4":
				if father == "0/4":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="0/0" or father=="0/1" or father=="0/3" or father=="0/2":
					phase="0|4"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="4/4" or father =="1/4" or father =="2/4" or father=="3/4":
					phase="4|0"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="0/0" and (mother=="1/4" or mother=="4/4" or mother=="2/4" or mother=="3/4"):
				phase="0|4"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/4" or father=="4/4" or father=="2/4" or father=="3/4") and mother=="0/0":
				phase="4|0"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="0/1" or father=="0/2" or father=="0/3") and (mother=="1/4" or mother=="4/4"  or mother=="2/4" or mother=="3/4"):
				pahse="0|4"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/4" or father=="4/4" or father=="2/4" or father=="3/4") and (mother=="0/1" or mother=="0/2" or mother=="0/3"):
				pahse="4|0"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		elif child=="1/2":
			if father =="1/2":
				if mother == "1/2":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="1/1" or mother=="0/1" or mother=="1/3" or mother=="1/4":
					phase="2|1"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="0/2" or mother =="2/2" or mother =="2/3" or mother=="2/4":
					phase="1|2"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="1/2":
				if father == "1/2":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="1/1" or father=="0/1" or father=="1/3" or father=="1/4":
					phase="1|2"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="0/2" or father =="2/2" or father =="2/3" or father=="2/4":
					phase="2|1"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="1/1" and (mother=="2/2" or mother=="0/2" or mother=="2/3" or mother=="2/4"):
				phase="1|2"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="2/2" or father=="0/2" or father=="2/3" or father=="2/4") and mother=="1/1":
				phase="2|1"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/3" or father=="1/4" or father=="0/1") and (mother=="2/2" or mother=="2/3"  or mother=="2/4" or mother=="0/2"):
				pahse="1|2"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="2/2" or father=="0/2" or father=="2/3" or father=="2/4") and (mother=="0/1" or mother=="1/3" or mother=="1/4"):
				pahse="2|1"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		elif child=="1/3":
			if father =="1/3":
				if mother == "1/3":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="1/1" or mother=="0/1" or mother=="1/2" or mother=="1/4":
					phase="3|1"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="0/3" or mother =="3/3" or mother =="2/3" or mother=="3/4":
					phase="1|3"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="1/3":
				if father == "1/3":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="1/1" or father=="0/1" or father=="1/2" or father=="1/4":
					phase="1|3"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="0/3" or father =="3/3" or father =="2/3" or father=="3/4":
					phase="3|1"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="1/1" and (mother=="3/3" or mother=="0/3" or mother=="2/3" or mother=="3/4"):
				phase="1|3"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="3/3" or father=="0/3" or father=="2/3" or father=="3/4") and mother=="1/1":
				phase="3|1"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="0/1" or father=="1/2" or father=="1/4") and (mother=="3/3" or mother=="2/3"  or mother=="3/4" or mother=="0/3"):
				pahse="1|3"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="3/3" or father=="2/3" or father=="3/4" or father=="0/3") and (mother=="0/1" or mother=="1/2" or mother=="1/4"):
				pahse="3|1"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		elif child=="1/4":
			if father =="1/4":
				if mother == "1/4":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="1/1" or mother=="0/1" or mother=="1/2" or mother=="1/3":
					phase="4|1"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="0/4" or mother =="4/4" or mother =="2/4" or mother=="3/4":
					phase="1|4"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="1/4":
				if father == "1/4":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="1/1" or father=="0/1" or father=="1/2" or father=="1/3":
					phase="1|4"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="0/4" or father =="4/4" or father =="2/4" or father=="3/4":
					phase="4|1"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="1/1" and (mother=="4/4" or mother=="0/4" or mother=="2/4" or mother=="3/4"):
				phase="1|4"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="4/4" or father=="0/4" or father=="2/4" or father=="3/4") and mother=="1/1":
				phase="4|1"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="0/1" or father=="1/2" or father=="1/3") and (mother=="4/4" or mother=="2/4"  or mother=="3/4" or mother=="0/4"):
				pahse="1|4"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="4/4" or father=="2/4" or father=="3/4" or father=="0/4") and (mother=="0/1" or mother=="1/2" or mother=="1/3"):
				pahse="4|1"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		elif child=="2/3":
			if father =="2/3":
				if mother == "2/3":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="2/2" or mother=="0/2" or mother=="1/2" or mother=="2/4":
					phase="3|2"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="0/3" or mother =="3/3" or mother =="1/3" or mother=="3/4":
					phase="2|3"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="2/3":
				if father == "2/3":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="2/2" or father=="0/1" or father=="1/2" or father=="1/4":
					phase="2|3"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="0/3" or father =="3/3" or father =="1/3" or father=="3/4":
					phase="3|2"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="2/2" and (mother=="3/3" or mother=="1/3" or mother=="0/3" or mother=="3/4"):
				phase="2|3"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="3/3" or father=="1/3" or father=="0/3" or father=="3/4") and mother=="2/2":
				phase="3|2"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="1/2" or father=="2/4" or father=="0/2") and (mother=="3/3" or mother=="1/3"  or mother=="3/4" or mother=="0/3"):
				pahse="2|3"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="3/3" or father=="1/3" or father=="3/4" or father=="0/3") and (mother=="1/2" or mother=="0/2" or mother=="2/4"):
				pahse="3|2"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		elif child=="2/4":
			if father =="2/4":
				if mother == "2/4":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="2/2" or mother=="0/2" or mother=="1/2" or mother=="2/3":
					phase="4|2"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="0/4" or mother =="4/4" or mother =="1/4" or mother=="3/4":
					phase="2|4"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="2/4":
				if father == "2/4":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="2/2" or father=="0/2" or father=="1/2" or father=="2/3":
					phase="2|4"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="0/4" or father =="4/4" or father =="1/4" or father=="3/4":
					phase="4|2"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="2/2" and (mother=="4/4" or mother=="1/4" or mother=="0/4" or mother=="3/4"):
				phase="2|4"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="4/4" or father=="1/4" or father=="0/4" or father=="3/4") and mother=="2/2":
				phase="4|2"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="0/2" or father=="1/2" or father=="2/3") and (mother=="4/4" or mother=="1/4"  or mother=="3/4" or mother=="0/4"):
				pahse="2|4"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="4/4" or father=="1/4" or father=="3/4" or father=="0/4") and (mother=="0/2" or mother=="1/2" or mother=="2/3"):
				pahse="4|2"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		elif child=="3/4":
			if father =="3/4":
				if mother == "3/4":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif mother =="3/3" or mother=="0/3" or mother=="1/3" or mother=="2/3":
					phase="4|3"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif mother =="0/4" or mother =="4/4" or mother =="1/4" or mother=="2/4":
					phase="3|4"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif mother =="3/4":
				if father == "3/4":
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],child)
					INDEL.append([line[3],line[4],child])
				elif father =="3/3" or father=="0/3" or father=="1/3" or father=="2/3":
					phase="3|4"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				elif father =="0/4" or father =="4/4" or father =="1/4" or father=="2/4":
					phase="4|3"
					print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
					INDEL.append([line[3],line[4],phase])
				else:
					print "bad",line
			elif father=="3/3" and (mother=="4/4" or mother=="1/4" or mother=="0/4" or mother=="2/4"):
				phase="3|4"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="4/4" or father=="1/4" or father=="0/4" or father=="2/4") and mother=="3/3":
				phase="4|3"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="0/3" or father=="1/3" or father=="2/3") and (mother=="4/4" or mother=="1/4"  or mother=="2/4" or mother=="0/4"):
				pahse="3|4"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			elif (father=="4/4" or father=="1/4" or father=="2/4" or father=="0/4") and (mother=="0/3" or mother=="1/3" or mother=="2/3"):
				pahse="4|3"
				print "%s\t%i\t%i\t%i\t%s\t%s\t%s" %(chr,pos-1,pos,line_number,line[3],line[4],phase)
				INDEL.append([line[3],line[4],phase])
			else:
				print "bad",line
		else:
			print "bad",line
	return INDEL

def read_whole_genome_vcf(whole_genome_vcf):
	f = open(whole_genome_vcf,"r")
	whole_indel = {}
	for i in chr_list:
		whole_indel[i]=[]
	for line in f:
		if line[0]=="#":
			continue
		line = line.strip().split("\t")
		chr = line[0]
		pos = int(line[1])
		ref = line[3]
		alt = line[4].split(",")
		geno = line[9][:3]
		if chr not in chr_list:
			break
		whole_indel[chr].append([pos,ref,alt,geno])
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

def compare_vcf(KG_file,INDEL,whole_indel):			# compare 1000G indel calls with our whole genome call
	not_shown = 0
	shown = 0
	same = 0
	INDEL_same = {}
	for i in range(1,23):
		INDEL_same["chr"+str(i)]=[]
	INDEL_same["chrX"]=[]
	f=open(KG_file,"r")
	for line in f:
		line = line.strip().split("\t")
		index=int(line[3])
		chr = line[0]
		pos = int(line[2])
		ref = line[4]
		alt = INDEL[index][1]
		genotype = INDEL[index][2]
		geno = [genotype[0],genotype[1],genotype[2]]
		try:
			assert ref==INDEL[index][0]
		except AssertionError:
			print "assertation",line,index,INDEL[index]
			ref==INDEL[index][0]
		if chr not in chr_list:
			continue
		find = find_pos(pos,whole_indel[chr])
		if find<0:
			not_shown += 1
		else:
			shown +=1
#			print chr,pos,ref,alt,whole_indel[chr][find]
			alt = alt.split(",")
			if whole_indel[chr][find][1]==ref:
				if len(whole_indel[chr][find][2])==1 and len(alt)==1:
					if whole_indel[chr][find][2][0]==alt[0]:
						same +=1
						print "same",chr,pos,ref,alt,geno,whole_indel[chr][find]
						INDEL_same[chr].append([pos,ref,alt,"".join(geno)])
					else:
						print chr,pos,ref,alt,geno,whole_indel[chr][find]
				elif len(alt)==1:
					if alt[0] in whole_indel[chr][find][2]: 
						same +=1
						index = whole_indel[chr][find][2].index(alt[0])
						print "before",chr,pos,ref,alt,geno,whole_indel[chr][find]
						if index!=0:
							geno[0]=geno[0].replace("1",str(index+1))
							geno[2]=geno[2].replace("1",str(index+1))
							print "after",chr,pos,ref,alt,geno,whole_indel[chr][find]
						if int(geno[0])> int(geno[2]) and geno[1]=="/":
							print "exchange",chr,pos,ref,alt,geno,whole_indel[chr][find]
							m=geno[0]
							geno[0]=geno[2]
							geno[2]=m
						INDEL_same[chr].append([pos,ref,whole_indel[chr][find][2],"".join(geno)])
					else:
						print chr,pos,ref,alt,geno,whole_indel[chr][find]
				else:
					a = int(geno[0])
					if a!=0:
						if alt[a-1] in whole_indel[chr][find][2]:
							index = whole_indel[chr][find][2].index(alt[a-1])
							a = str(index+1)
						else:
							a = "NA"
					b = int(geno[2])
					if b!=0:
						if alt[b-1] in whole_indel[chr][find][2]:
							index = whole_indel[chr][find][2].index(alt[b-1])
							b = str(index+1)
						else:
							b = "NA"
					if a !="NA" and b!="NA":
						same +=1
						geno[0]=a
						geno[2]=b
						if int(geno[0])> int(geno[2]) and geno[1]=="/":
							print "exchange",chr,pos,ref,alt,geno,whole_indel[chr][find]
							m=geno[0]
							geno[0]=geno[2]
							geno[2]=m
						print "same multi het",chr,pos,ref,alt,whole_indel[chr][find],"".join(map(str,geno))
						INDEL_same[chr].append([pos,ref,whole_indel[chr][find][2],"".join(map(str,geno))])
					else:
						print chr,pos,ref,alt,geno,whole_indel[chr][find]
			else:
				print chr,pos,ref,alt,geno,whole_indel[chr][find]
	print not_shown,shown,same
	return INDEL_same

def compare_vcf_v2(KG_file,INDEL,whole_indel):			# compare 1000G indel calls with our whole genome call
	not_shown = 0
	shown = 0
	same = 0
	f=open(KG_file,"r")
	KG_phased = 0
	KG_whole = 0
	for line in f:
		line = line.strip().split("\t")
		index=int(line[3])
		chr = line[0]
		pos = int(line[2])
		ref = line[4]
		alt = INDEL[index][1]
		genotype = INDEL[index][2]
		geno = [genotype[0],genotype[1],genotype[2]]
		try:
			assert ref==INDEL[index][0]
		except AssertionError:
			print "assertation",line,index,INDEL[index]
			ref==INDEL[index][0]
		if chr not in chr_list:
			continue
		if geno[1]=="|":
			KG_phased +=1
	for chr in chr_list:
		for i in range(len(whole_indel[chr])):
			if whole_indel[chr][i][3]==whole_indel[chr][i][5]:
				KG_whole +=1
	print KG_phased,KG_whole
	
def compare_vcf_v3(fosmid_INDEL,fosmid_INDEL_wgs):
	tot = 0
	same = 0
	fosmid_phased = {}
	for chr in chr_list:
		fosmid_phased[chr] = []
		for i in range(len(fosmid_INDEL[chr])):
			if fosmid_INDEL[chr][i][5]!="":
				pos = fosmid_INDEL[chr][i][0]
				find = find_pos(pos,fosmid_INDEL_wgs[chr])
				tot +=1
				if find>=0:
					if fosmid_INDEL[chr][i][3]==fosmid_INDEL_wgs[chr][find][5]:
						same +=1
						fosmid_phased[chr].append(fosmid_INDEL[chr][i])
#						print "\t".join(fosmid_INDEL[chr][i])					
	print "comparison between fosmid and fosmid_wgs"
	print tot,same

def compare_vcf_v4(fosmid_INDEL,fosmid_INDEL_wgs):
	fosmid_phased = {}
	a = 0
	for chr in chr_list:
		fosmid_phased[chr] = []
		for i in range(len(fosmid_INDEL[chr])):
			fosmid_phased[chr].append(fosmid_INDEL[chr][i][0:4])
		for j in range(len(fosmid_INDEL_wgs[chr])):
				pos = fosmid_INDEL_wgs[chr][j][0]
				find = find_pos(pos,fosmid_INDEL[chr])
				if find<0:
					try:
						fosmid_phased[chr].append(fosmid_INDEL_wgs[chr][j][0:3]+[fosmid_INDEL_wgs[chr][j][5]])
					except IndexError:
						print fosmid_INDEL_wgs[chr][j]
		print chr,len(fosmid_phased[chr])
		a += len(fosmid_phased[chr])
	print a
	return fosmid_phased

def intersect(KG_file,fosmid_INDEL,INDEL):
	fosmid_KG = 0
	fosmid_KG_whole = 0
	INDEL_same = {}
	for i in range(1,23):
		INDEL_same["chr"+str(i)]=[]
	INDEL_same["chrX"]=[]
	f=open(KG_file,"r")
	for line in f:
		line = line.strip().split("\t")
		index=int(line[3])
		chr = line[0]
		pos = int(line[2])
		ref = line[4]
		alt = INDEL[index][1]
		genotype = INDEL[index][2]
		geno = [genotype[0],genotype[1],genotype[2]]
		try:
			assert ref==INDEL[index][0]
		except AssertionError:
			print "assertation",line,index,INDEL[index]
			ref==INDEL[index][0]
		if chr not in chr_list:
			continue
		find = find_pos(pos,fosmid_INDEL[chr])
		if find<0:
			continue
		else:
			alt = alt.split(",")
			if fosmid_INDEL[chr][find][1]==ref:
				if len(fosmid_INDEL[chr][find][2])==1 and len(alt)==1:
					if fosmid_INDEL[chr][find][2][0]==alt[0]:
						fosmid_KG +=1
						print "same",chr,pos,ref,alt,geno,fosmid_INDEL[chr][find]
						if fosmid_INDEL[chr][find][5]!="":
							fosmid_KG_whole +=1
							INDEL_same[chr].append([pos,ref,alt,"".join(geno),fosmid_INDEL[chr][find][5],fosmid_INDEL[chr][find][3]])  #KG,whole_genome,fosmid
					else:
						print chr,pos,ref,alt,geno,fosmid_INDEL[chr][find]
				elif len(alt)==1:
					if alt[0] in fosmid_INDEL[chr][find][2]: 
						fosmid_KG +=1
						index = fosmid_INDEL[chr][find][2].index(alt[0])
						print "before",chr,pos,ref,alt,geno,fosmid_INDEL[chr][find]
						if index!=0:
							geno[0]=geno[0].replace("1",str(index+1))
							geno[2]=geno[2].replace("1",str(index+1))
							print "after",chr,pos,ref,alt,geno,fosmid_INDEL[chr][find]
						if int(geno[0])> int(geno[2]) and geno[1]=="/":
							print "exchange",chr,pos,ref,alt,geno,fosmid_INDEL[chr][find]
							m=geno[0]
							geno[0]=geno[2]
							geno[2]=m
						if fosmid_INDEL[chr][find][5]!="":
							fosmid_KG_whole +=1
							INDEL_same[chr].append([pos,ref,fosmid_INDEL[chr][find][2],"".join(geno),fosmid_INDEL[chr][find][5],fosmid_INDEL[chr][find][3]])  #KG,whole_genome,fosmid
					else:
						print chr,pos,ref,alt,geno,fosmid_INDEL[chr][find]
				else:
					a = int(geno[0])
					if a!=0:
						if alt[a-1] in fosmid_INDEL[chr][find][2]:
							index = fosmid_INDEL[chr][find][2].index(alt[a-1])
							a = str(index+1)
						else:
							a = "NA"
					b = int(geno[2])
					if b!=0:
						if alt[b-1] in fosmid_INDEL[chr][find][2]:
							index = fosmid_INDEL[chr][find][2].index(alt[b-1])
							b = str(index+1)
						else:
							b = "NA"
					if a !="NA" and b!="NA":
						fosmid_KG +=1
						geno[0]=a
						geno[2]=b
						print "same multi het",chr,pos,ref,alt,fosmid_INDEL[chr][find],"".join(map(str,geno))
						if int(geno[0])> int(geno[2]) and geno[1]=="/":
							print "exchange",chr,pos,ref,alt,geno,fosmid_INDEL[chr][find]
							m=geno[0]
							geno[0]=geno[2]
							geno[2]=m
						if fosmid_INDEL[chr][find][5]!="":
							fosmid_KG_whole +=1
							INDEL_same[chr].append([pos,ref,fosmid_INDEL[chr][find][2],"".join(map(str,geno)),fosmid_INDEL[chr][find][5],fosmid_INDEL[chr][find][3]])  #KG,whole_genome,fosmid
					else:
						print chr,pos,ref,alt,geno,fosmid_INDEL[chr][find]
			else:
				print chr,pos,ref,alt,geno,fosmid_INDEL[chr][find]
	print fosmid_KG,fosmid_KG_whole
	return INDEL_same

def intersect_phased(KG_file,fosmid_INDEL,fosmid_INDEL_wgs,INDEL):
	n1,n2,n3,n12,n13,n23,n123=[0,0,0,0,0,0,0]
	INDEL_all = {}
	for i in range(1,23):
		INDEL_all["chr"+str(i)]=[]
	INDEL_all["chrX"]=[]
	f=open(KG_file,"r")
	KG_indel_phased = 0
	for line in f:
		line = line.strip().split("\t")
		index=int(line[3])
		chr = line[0]
		pos = int(line[2])
		ref = line[4]
		alt = INDEL[index][1]
		genotype = INDEL[index][2]
		geno = [genotype[0],genotype[1],genotype[2]]
		try:
			assert ref==INDEL[index][0]
		except AssertionError:
			print "assertation",line,index,INDEL[index]
			ref==INDEL[index][0]
		if chr not in chr_list:
			continue
		if geno[1]=="|":
			KG_indel_phased +=1
		else:
			continue
		if geno[0]=='0' and geno[2]=='0':
			continue
		denovo_geno = ''
		wgs_geno = ''
		find1 = find_pos(pos,fosmid_INDEL[chr])
		find2 = find_pos(pos,fosmid_INDEL_wgs[chr])
		if find1<0 and find2<0:
			INDEL_all[chr].append([pos,ref,alt,genotype,wgs_geno,denovo_geno])			#KG,fosmid_wgs,fosmid_denovo
		elif find1>=0:
			alt = alt.split(",")
			if fosmid_INDEL[chr][find1][1]==ref:
				if len(fosmid_INDEL[chr][find1][2])==1 and len(alt)==1:
					if fosmid_INDEL[chr][find1][2][0]==alt[0]:
						denovo_geno = fosmid_INDEL[chr][find1][3]
						if fosmid_INDEL[chr][find1][5]!="":
							if find2>=0:
								wgs_geno = fosmid_INDEL_wgs[chr][find2][5]
						if denovo_geno == genotype:
							n13+=1
						INDEL_all[chr].append([pos,ref,'\t'.join(alt),genotype,wgs_geno,denovo_geno])  #KG,whole_genome,fosmid
				elif len(alt)==1:
					if alt[0] in fosmid_INDEL[chr][find1][2]: 
						index = fosmid_INDEL[chr][find1][2].index(alt[0])
						print "before",chr,pos,ref,alt,geno,fosmid_INDEL[chr][find1]
						if index!=0:
							geno[0]=geno[0].replace("1",str(index+1))
							geno[2]=geno[2].replace("1",str(index+1))
							print "after",chr,pos,ref,alt,geno,fosmid_INDEL[chr][find1]
						genotype = "".join(geno)
						denovo_geno = fosmid_INDEL[chr][find1][3]
						if fosmid_INDEL[chr][find1][5]!="":
							if find2>=0:
								wgs_geno = fosmid_INDEL_wgs[chr][find2][5]
						if denovo_geno == genotype:
							n13+=1
						INDEL_all[chr].append([pos,ref,'\t'.join(alt),genotype,wgs_geno,denovo_geno])  #KG,whole_genome,fosmid
				else:
					a = int(geno[0])
					if a!=0:
						if alt[a-1] in fosmid_INDEL[chr][find1][2]:
							index = fosmid_INDEL[chr][find1][2].index(alt[a-1])
							a = str(index+1)
						else:
							a = "NA"
					b = int(geno[2])
					if b!=0:
						if alt[b-1] in fosmid_INDEL[chr][find1][2]:
							index = fosmid_INDEL[chr][find1][2].index(alt[b-1])
							b = str(index+1)
						else:
							b = "NA"
					if a !="NA" and b!="NA":
						geno[0]=a
						geno[2]=b
						print "same multi het",chr,pos,ref,'\t'.join(alt),fosmid_INDEL[chr][find1],"".join(map(str,geno))
						genotype = "".join(map(str,geno))
						denovo_geno = fosmid_INDEL[chr][find1][3]
						if fosmid_INDEL[chr][find1][5]!="":
							if find2>=0:
								wgs_geno = fosmid_INDEL_wgs[chr][find2][5]
						if denovo_geno == genotype:
							n13+=1
						INDEL_all[chr].append([pos,ref,'\t'.join(alt),genotype,wgs_geno,denovo_geno])  #KG,whole_genome,fosmid
		elif find2>=0:
			alt = alt.split(",")
			if fosmid_INDEL_wgs[chr][find2][1]==ref:
				if len(fosmid_INDEL_wgs[chr][find2][2])==1 and len(alt)==1:
					if fosmid_INDEL_wgs[chr][find2][2][0]==alt[0]:
						wgs_geno = fosmid_INDEL_wgs[chr][find2][5]
						INDEL_all[chr].append([pos,ref,'\t'.join(alt),genotype,wgs_geno,denovo_geno])  #KG,whole_genome,fosmid
				elif len(alt)==1:
					if alt[0] in fosmid_INDEL_wgs[chr][find2][2]: 
						index = fosmid_INDEL_wgs[chr][find2][2].index(alt[0])
						print "before",chr,pos,ref,'\t'.join(alt),geno,fosmid_INDEL_wgs[chr][find2]
						if index!=0:
							geno[0]=geno[0].replace("1",str(index+1))
							geno[2]=geno[2].replace("1",str(index+1))
							print "after",chr,pos,ref,alt,geno,fosmid_INDEL_wgs[chr][find2]
						genotype = "".join(geno)
						wgs_geno = fosmid_INDEL_wgs[chr][find2][5]
						INDEL_all[chr].append([pos,ref,'\t'.join(alt),genotype,wgs_geno,denovo_geno])  #KG,whole_genome,fosmid
				else:
					a = int(geno[0])
					if a!=0:
						if alt[a-1] in fosmid_INDEL_wgs[chr][find2][2]:
							index = fosmid_INDEL_wgs[chr][find2][2].index(alt[a-1])
							a = str(index+1)
						else:
							a = "NA"
					b = int(geno[2])
					if b!=0:
						if alt[b-1] in fosmid_INDEL_wgs[chr][find2][2]:
							index = fosmid_INDEL_wgs[chr][find2][2].index(alt[b-1])
							b = str(index+1)
						else:
							b = "NA"
					if a !="NA" and b!="NA":
						geno[0]=a
						geno[2]=b
						print "same multi het",chr,pos,ref,'\t'.join(alt),fosmid_INDEL_wgs[chr][find2],"".join(map(str,geno))
						genotype = "".join(map(str,geno))
						wgs_geno = fosmid_INDEL_wgs[chr][find2][5]
						INDEL_all[chr].append([pos,ref,'\t'.join(alt),genotype,wgs_geno,denovo_geno])  #KG,whole_genome,fosmid
	for chr in chr_list:
		for i in range(len(fosmid_INDEL[chr])):
			find = find_pos(pos,INDEL_all[chr])
			wgs_geno = ''
			denovo_geno = fosmid_INDEL[chr][i][3]
			if find <0:					# not observed in 1000G
				if fosmid_INDEL[chr][i][5]!="":
					pos = fosmid_INDEL[chr][i][0]
					ref = fosmid_INDEL[chr][i][1]
					alt = fosmid_INDEL[chr][i][2]
					find2 = find_pos(pos,fosmid_INDEL_wgs[chr])
					if find2>=0:
						wgs_geno = fosmid_INDEL_wgs[chr][find2][5]
						INDEL_all[chr].append([pos,ref,'\t'.join(alt),'',wgs_geno,denovo_geno])  #KG,whole_genome,fosmid
					else:
						INDEL_all[chr].append([pos,ref,'\t'.join(alt),'',wgs_geno,denovo_geno])  #KG,whole_genome,fosmid						
				else:
#					find2 = find_pos(pos,fosmid_INDEL_wgs[chr])
#					if find2<0:
					INDEL_all[chr].append([pos,ref,'\t'.join(alt),'',wgs_geno,denovo_geno])  #KG,whole_genome,fosmid

	for chr in chr_list:
		for i in range(len(fosmid_INDEL_wgs[chr])):
			find = find_pos(pos,INDEL_all[chr])
			find1 = find_pos(pos,fosmid_INDEL[chr])
			if find<0 and find1<0:
				pos = fosmid_INDEL_wgs[chr][i][0]
				ref = fosmid_INDEL_wgs[chr][i][1]
				alt = ",".join(fosmid_INDEL_wgs[chr][i][2])
				wgs_geno = fosmid_INDEL_wgs[chr][find2][5]
				n2 +=1
				INDEL_all[chr].append([pos,ref,alt,'',wgs_geno,''])  #KG,whole_genome,fosmid
		INDEL_all[chr].sort( key=lambda w: int(w[0]) )
	print n2,n13	
	return INDEL_all
	
def KG_phased_print(KG_file,INDEL):
	KG_indel_phased = 0
	f=open(KG_file,"r")
	f_out = open("NA19240_KG_phased_indel.bed",'w')
	for line in f:
		line = line.strip().split("\t")
		index=int(line[3])
		chr = line[0]
		pos = int(line[2])
		ref = line[4]
		alt = INDEL[index][1]
		genotype = INDEL[index][2]
		geno = [genotype[0],genotype[1],genotype[2]]
		try:
			assert ref==INDEL[index][0]
		except AssertionError:
			print "assertation",line,index,INDEL[index]
			ref==INDEL[index][0]
		if chr not in chr_list:
			continue
		if geno[1]=="|":
			print >>f_out,"\t".join(line[0:3]),"\t","\t".join(INDEL[index])

def fosmid_print(fosmid_INDEL, fosmid_INDEL_wgs):
	f1 = open('fosmid_phased_indel_denovo.bed','w')
	f2 = open('fosmid_phased_indel_wgs.bed','w')
	for chr in chr_list:
		for i in range(len(fosmid_INDEL[chr])):
			if fosmid_INDEL[chr][i][3][1]=='|':
				print >>f1,chr,"\t","\t".join(map(str,fosmid_INDEL[chr][i][0:4]))
		for j in range(len(fosmid_INDEL_wgs[chr])):
			print >>f2,chr,"\t","\t".join(map(str,fosmid_INDEL_wgs[chr][j][0:3]+[fosmid_INDEL_wgs[chr][j][5]]))

def fosmid_print_combined(fosmid_INDEL, fosmid_INDEL_wgs):
	fosmid_phased = {}
	for chr in chr_list:
		fosmid_phased[chr] = []
		for i in range(len(fosmid_INDEL[chr])):
			if fosmid_INDEL[chr][i][3][1]=='|':
				if fosmid_INDEL[chr][i][5]!="":
					pos = fosmid_INDEL[chr][i][0]
					find = find_pos(pos,fosmid_INDEL_wgs[chr])
					tot +=1
					if find>=0:
						if fosmid_INDEL[chr][i][3]==fosmid_INDEL_wgs[chr][find][5]:
							fosmid_phased[chr].append(fosmid_INDEL[chr][i][0:4])
				else:
					fosmid_phased[chr].append(fosmid_INDEL[chr][i][0:4])
		for j in range(len(fosmid_INDEL_wgs[chr])):
				pos = fosmid_INDEL_wgs[chr][j][0]
				find = find_pos(pos,fosmid_INDEL[chr])
				if find<0:
					fosmid_phased[chr].append(fosmid_INDEL_wgs[chr][j][0:3]+[fosmid_INDEL_wgs[chr][j][5]])
		print chr,len(fosmid_phased[chr])
		a += len(fosmid_phased[chr])
	return fosmid_phased

def venn_compare(INDEL_same): 
	a,b,c,n12,n13,n23,n123=[0,0,0,0,0,0,0]   #KG,whole_genome,fosmid
	for chr in chr_list:
		for i in range(len(INDEL_same[chr])):
			if INDEL_same[chr][i][3]!='':
				if INDEL_same[chr][i][3]==INDEL_same[chr][i][4]:
					if INDEL_same[chr][i][4]==INDEL_same[chr][i][5]:
						n123 +=1
					elif INDEL_same[chr][i][5]!='':
						c +=1
						n12 +=1
					else:
						n12+=1
				elif INDEL_same[chr][i][3]==INDEL_same[chr][i][5]:
					n13+=1
					if INDEL_same[chr][i][4]!='':
						b+=1
				elif INDEL_same[chr][i][4]!='':
					a+=1
					b+=1
				elif INDEL_same[chr][i][5]!='':
					a+=1
					c+=1
				else:
					a+=1
			elif INDEL_same[chr][i][4]!='':
				if INDEL_same[chr][i][4]==INDEL_same[chr][i][5]:
					n23 +=1
				elif INDEL_same[chr][i][5]!='':
					b+=1
					c+=1
				else:
					b+=1
			elif INDEL_same[chr][i][5]!='':
				c +=1
	print a,b,c,n12,n13,n23,n123		

def venn_compare_v2(INDEL_same): 
	a,b,c,n12,n13,n23,n123=[0,0,0,0,0,0,0]   #KG,whole_genome,fosmid
	for chr in chr_list:
		for i in range(len(INDEL_same[chr])):
			if INDEL_same[chr][i][3]!='':
				if INDEL_same[chr][i][3]==INDEL_same[chr][i][4]:
					if INDEL_same[chr][i][4]==INDEL_same[chr][i][5]:
						n123 +=1
			if INDEL_same[chr][i][3]!='':
				if INDEL_same[chr][i][3]==INDEL_same[chr][i][4]:
					if INDEL_same[chr][i][3]!=INDEL_same[chr][i][5]:
						n12 +=1
			if INDEL_same[chr][i][3]!='':
				if INDEL_same[chr][i][3]==INDEL_same[chr][i][5]:
					if INDEL_same[chr][i][3]!=INDEL_same[chr][i][4]:
						n13 +=1
			if INDEL_same[chr][i][4]!='':
				if INDEL_same[chr][i][4]==INDEL_same[chr][i][5]:
					if INDEL_same[chr][i][3]!=INDEL_same[chr][i][4]:
						n23 +=1	
	print a,b,c,n12,n13,n23,n123		

def fosmid_print_combined(fosmid_INDEL,fosmid_INDEL_wgs):
	fosmid_phased = {}
	f1 =open('fosmid_phased_indel_combined_v2.bed','w')
	for chr in chr_list:
		fosmid_phased[chr] = []
		for i in range(len(fosmid_INDEL[chr])):
			if fosmid_INDEL[chr][i][3][1]=='|':
				if fosmid_INDEL[chr][i][5]!="":
					pos = fosmid_INDEL[chr][i][0]
					find = find_pos(pos,fosmid_INDEL_wgs[chr])
					if find>=0:
						if fosmid_INDEL[chr][i][3]==fosmid_INDEL_wgs[chr][find][5]:
							fosmid_phased[chr].append(fosmid_INDEL[chr][i][0:4])
				else:
					fosmid_phased[chr].append(fosmid_INDEL[chr][i][0:4])
		for j in range(len(fosmid_INDEL_wgs[chr])):
				pos = fosmid_INDEL_wgs[chr][j][0]
				find = find_pos(pos,fosmid_INDEL[chr])
				if find<0:
					fosmid_phased[chr].append(fosmid_INDEL_wgs[chr][j][0:3]+[fosmid_INDEL_wgs[chr][j][5]])
#		print chr,len(fosmid_phased[chr])
#		a += len(fosmid_phased[chr])
		for m in range(len(fosmid_phased[chr])):
			print >>f1,chr,"\t","\t".join(map(str,fosmid_phased[chr][m]))
	return fosmid_phased

if __name__=="__main__":
	dbfile=open('NA19240_indel_pickle','rb')   # 1000G phased indels
	KG_INDEL=pickle.load(dbfile)
	dbfile=open('fosmid_INDEL_determined_pickle','rb')   # fosmid denovo determined phased indel
	fosmid_INDEL=pickle.load(dbfile)

	dbfile=open('NA19240_fosmid_of_wgs_indel_infer_pickle','rb')    # fosmid with wgs guidance phased indel
	fosmid_INDEL_wgs=pickle.load(dbfile)
	
	whole_genome_vcf = "/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/merge/NA19240_combined_indel.vcf"
	whole_indel = read_whole_genome_vcf(whole_genome_vcf)

	KG_file ="NA19240_indel_hg19.bed"

#	compare_vcf_v3(fosmid_INDEL,fosmid_INDEL_wgs)

#	INDEL_all = intersect_phased(KG_file,fosmid_INDEL,fosmid_INDEL_wgs,KG_INDEL)
#	dbfile=open('NA19240_indel_KG_wgs_fosmid_determined_all_pickle','wb')
#	pickle.dump(INDEL_all,dbfile)
#	venn_compare_v2(INDEL_all)
	
	fosmid_phased=fosmid_print_combined(fosmid_INDEL,fosmid_INDEL_wgs)