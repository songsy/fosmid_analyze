#!/usr/bin/python

import os
import numpy as np
import argparse
import bisect
import pickle
def load_letters_by_allele(vcf_filename, chrom):
	letter_by_allele = dict()
	one_two_positions = set()
	cmd = "zcat %s | grep -P '^%s\\t' | cut -f 2,4,5,10 | cut -d ':' -f 1 |"\
		  "grep -v 1/1 | grep -v 1\|1 | grep -v 0/0 | grep -v 0\|0 | grep -v 1\|2 |"\
		  "grep -v -P '1/\.'" % (vcf_filename, chrom)
	vcf_file = os.popen(cmd)

	for j, line in enumerate(vcf_file):
		fields = line.strip().split()

		pos = int(fields[0])
		if fields[3] == "0/1" or fields[3] == "1/0"\
		   or fields[3] == "0|1" or fields[3] == "1|0":
			letter_by_allele[pos] = (fields[1], fields[2])
		elif fields[3] == "1/2":
			letters = fields[2].split(',')
			letter_by_allele[pos] = (letters[0], letters[1])
			one_two_positions.add(pos)
		else:
			print "ERROR: Could not process this line in the VCF:"
			print line
			exit(1)
	return letter_by_allele, one_two_positions
if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--vcf')
	parser.add_argument('--chr')
	parser.add_argument('--sample')
	args = parser.parse_args()
	
	letter_by_allele, one_two_positions = load_letters_by_allele(args.vcf, args.chr)
	dbfile=open(args.sample+'_'+args.chr+'_pickle','wb')
	pickle.dump(letter_by_allele,dbfile)
	pickle.dump(one_two_positions,dbfile) 