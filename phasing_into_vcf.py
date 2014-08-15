import argparse
from NGS_utils import *
import gzip
import commands

def read_phased_snps(file):
	hap = {}
	for i in BedIterator(file):
		chr = i[0]
		pos = int(i[2])
		hap[pos]=i[3]
	return hap

class VcfIterator:
  def __init__(self, filename):
	self.file = gzip.open(filename, "r")
  
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
	return (chrom, pos, alleles, geno)

def read_whole_genome_vcf(whole_genome_vcf,chr,hap):
	f_out = gzip.open('NA12878.chr%s.fosmid.phase.vcf.gz' %(chr),'w')
	f_out.write('##fileformat=VCFv4.1\n')
	f_out.write('##fileDate=%s\n' %(commands.getoutput('date +"%m%d%y"')))
	f_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased Genotype">\n')
	f_out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n')
	a = 0
	for i in VcfIterator(whole_genome_vcf):
		if i[3][0]==i[3][2]:
			if i[3][0]=='1':
				f_out.write('chr%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\n' %(chr,i[1],i[2][0],",".join(i[2][1:]),i[3]))
		else:
			try:
				phase = hap[int(i[1])]
				if phase == '--':
					f_out.write('chr%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\n' %(chr,i[1],i[2][0],",".join(i[2][1:]),i[3]))
				else:
					f_out.write('chr%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\n' %(chr,i[1],i[2][0],",".join(i[2][1:]),phase[0]+'|'+phase[1]))
			except:
				a+=1
				f_out.write('chr%s\t%s\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\n' %(chr,i[1],i[2][0],",".join(i[2][1:]),i[3]))
	f_out.close()
	print chr,a


if __name__=="__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--chr", dest='chr',help="chr")
	args = parser.parse_args()
	
	phased_snp_file='/home/jmkidd/kidd-lab-scratch/shiya-projects/CEU/haplotypes/chr'+ args.chr+"_refhap_hg19.info"
	hap=read_phased_snps(phased_snp_file)
	
	all_sites_vcf="/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/NA12878/all_sites/NA12878.chr" + args.chr+".vcf.gz"
	read_whole_genome_vcf(all_sites_vcf,args.chr,hap)	