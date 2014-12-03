import os
import glob
import io

def walk_through_NA19240():
	wk= "/home/jmkidd/kidd-lab/jmkidd-projects/NA19240-pools/results/align/pools/NA19240_pool_"
	list = range(1,288)
	list.remove(120)
	list.remove(155)
	for i in list:
		if i<10:
			index = "00" + str(i)
		elif i<100:
			index = "0" + str(i)
		else:
			index = str(i)
		vcf_file = wk + index + "/NA19240_pool_" + index + ".other.bam"
		for file in glob.glob(vcf_file):
			cmd = 'rm '+file
			os.popen(cmd)
			print cmd
		vcf_file = wk + index + "/NA19240_pool_" + index + ".unknown*bam"
		for file in glob.glob(vcf_file):
			cmd = 'rm '+file
			os.popen(cmd)
			print cmd
		vcf_file = wk + index + "/NA19240_pool_" + index + "*idx"
		for file in glob.glob(vcf_file):
			cmd = 'gzip '+file
			os.popen(cmd)
			print cmd
#		vcf_file = wk + index + "/NA19240_pool_" + index + "*bai"
#		for file in glob.glob(vcf_file):
#			if file != wk + index + "/NA19240_pool_" + index + ".markdup.realign.recal.bai":
#				cmd = 'rm '+file
#				os.popen(cmd)
#				print cmd

class VcfIterator:
  def __init__(self, filename):
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

if __name__=="__main__":
	walk_through_NA19240()