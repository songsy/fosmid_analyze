import sys
import commands
import os
import argparse
from NGS_utils import *

def make_prism_cmd(file,out_file,size,chr):
	block_file = '%s.chr%s.local.blocks.all' %(args.sample,chr)
	position_file = '%s.chr%s.positions' %(args.sample,chr)
	item = []
	interleave = []
#	os.popen("sed -i '$ d' %s" %(out_file))
#	os.popen("sed -i '$ d' %s" %(out_file.replace('.cmds','.interleave.cmds')))
	f=open('%s.prism.redo.v2.cmds' %(args.sample),'a')
	f_out1 = open(out_file,'a')
#	f_out2 = open(out_file.replace('.cmds','.redo.cmds'),'w')
	f_out3 = open(out_file.replace('.cmds','.interleave.cmds'),'a')
	REF_DIR='/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono'
	REF_PANEL="--reference-panel %s/ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz --reference-legend %s/ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz --reference-sample %s/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample" %(REF_DIR,chr,REF_DIR,chr,REF_DIR)
	vcf_file = "/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/%s/gVCF_calls/%s.chr%s.vcf.gz" %(args.sample,args.sample,chr)
	for i in BedIterator(file):
		if i[4]!='Affy' and i[4]!='None' and i[4]!='Hapmap' and i[4]!='1kG':
			if i[4]=='interleaving' or i[4]=='None,interleave' or i[4]=='Affy,interleave' or i[4]=='Hapmap,interleave' or i[4]=='1kG,interleave':
				interleave.append(i[3])
			continue
		if len(item)<3:
			item.append((i[1],i[2],i[3]))
		else:
			if int(item[-1][1])-int(item[0][0])<size:
				item.append((i[1],i[2],i[3]))
			else:
				start = item[0][0]
				end = item[-1][1]
				start_id = item[0][2]
				end_id = item[-1][2]
				path_file = '%s.chr%s.all.%s.%s.path' %(args.sample,chr,start_id,end_id)
				phase_file = '%s.chr%s.all.%s.%s.phase' %(args.sample,chr,start_id,end_id)
				cmds='prism --blocks %s --positions %s --start %s --end %s --K 100 --path %s --phase %s' %(block_file,position_file,start,end,path_file,phase_file)
				if os.path.isfile('prism_all_block/%s' %(path_file)) and os.path.isfile('prism_all_block/%s' %(phase_file)):
					cmds = cmds
#					print cmds
				else:
#					print 'no file',cmds
					print >>f_out1,cmds
#					print >>f_out2,cmds
				if len(interleave)!=0:
					while len(interleave)!=0:
						if int(interleave[-1])>int(end_id):
							interleave.pop()
						else:
							break
					if len(interleave)!=0:
						for block_id in interleave:
							cmds = 'grep "ID=%s:" %s.chr%s.local.interleaving.blocks >> %s.chr%s.all.%s.%s.interleaving.blocks' %(block_id,args.sample,chr,args.sample,chr,start_id,end_id)
							os.popen(cmds)
						unphased_block = '%s.chr%s.all.%s.%s.interleaving.blocks' %(args.sample,chr,start_id,end_id)
						phased_block = '%s.chr%s.all.%s.%s.phased.interleaving.blocks' %(args.sample,chr,start_id,end_id)
						cmds = 'prism-interleaving --unphased-blocks %s --path prism_all_block/%s --vcf %s --chr chr%s %s --phased-blocks %s' %(unphased_block,path_file,vcf_file,chr,REF_PANEL,phased_block)
						print >>f_out3,cmds
						interleave =[]
				item = [item[-1]]
	if item[-1][0]!=i[1]:
		if i[4]=='Affy' or i[4]=='None' or i[4]=='Hapmap' or i[4]=='1kG':
			item.append((i[1],i[2],i[3]))
	if len(item)!=0:
		print item
		start = item[0][0]
		end = item[-1][1]
		start_id = item[0][2]
		end_id = item[-1][2]
		path_file = '%s.chr%s.all.%s.%s.path' %(args.sample,chr,start_id,end_id)
		phase_file = '%s.chr%s.all.%s.%s.phase' %(args.sample,chr,start_id,end_id)
		cmds='prism --blocks %s --positions %s --start %s --end %s --K 100 --path %s --phase %s' %(block_file,position_file,start,end,path_file,phase_file)
		print >>f,cmds
		print >>f_out1,cmds
		print cmds
	if len(interleave)!=0:
		while len(interleave)!=0:
			if int(interleave[-1])>int(end_id):
				interleave.pop()
			else:
				break
		if len(interleave)!=0:
			for block_id in interleave:
				cmds = 'grep "ID=%s:" %s.chr%s.local.interleaving.blocks >> %s.chr%s.all.%s.%s.interleaving.blocks' %(block_id,args.sample,chr,args.sample,chr,start_id,end_id)
				os.popen(cmds)
			unphased_block = '%s.chr%s.all.%s.%s.interleaving.blocks' %(args.sample,chr,start_id,end_id)
			phased_block = '%s.chr%s.all.%s.%s.phased.interleaving.blocks' %(args.sample,chr,start_id,end_id)
			cmds = 'prism-interleaving --unphased-blocks %s --path prism_all_block/%s --vcf %s --chr chr%s %s --phased-blocks %s' %(unphased_block,path_file,vcf_file,chr,REF_PANEL,phased_block)
			print >>f_out3,cmds
			interleave =[]

def make_prism_cmd_partial(file,out_file,size,chr):
	block_file = '%s.chr%s.local.blocks' %(args.sample,chr)
	position_file = '%s.chr%s.positions' %(args.sample,chr)
	item = []
	interleave = []
	f_out1 = open(out_file.replace('.cmds','.partial.cmds'),'w')
#	f_out2 = open(out_file.replace('.cmds','.redo.cmds'),'w')
	f_out3 = open(out_file.replace('.cmds','.partial.interleave.cmds'),'w')
	REF_DIR='/home/jmkidd/kidd-lab/genomes/snp-sets/1KG/phase1/ALL.integrated_phase1_SHAPEIT_16-06-14.nomono'
	REF_PANEL="--reference-panel %s/ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz --reference-legend %s/ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz --reference-sample %s/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample" %(REF_DIR,chr,REF_DIR,chr,REF_DIR)
	vcf_file = "/home/jmkidd/kidd-lab/jmkidd-projects/additional-fosmid-pools/results/wgs-align/%s/gVCF_calls/%s.chr%s.vcf.gz" %(args.sample,args.sample,chr)
	for i in BedIterator(file):
		if i[4]!='Affy' and i[4]!='None' and i[4]!='Hapmap' and i[4]!='1kG':
			if i[4]=='interleaving' or i[4]=='None,interleave' or i[4]=='Affy,interleave' or i[4]=='Hapmap,interleave' or i[4]=='1kG,interleave':
				interleave.append(i[3])
			continue
		if len(item)<3:
			item.append((i[1],i[2],i[3]))
		else:
			if int(item[-1][1])-int(item[0][0])<size:
				item.append((i[1],i[2],i[3]))
			else:
				start = item[0][0]
				end = item[-1][1]
				start_id = item[0][2]
				end_id = item[-1][2]
				path_file = '%s.chr%s.all.%s.%s.path' %(args.sample,chr,start_id,end_id)
				phase_file = '%s.chr%s.all.%s.%s.phase' %(args.sample,chr,start_id,end_id)
				cmds='prism --blocks %s --positions %s --start %s --end %s --K 100 --path %s --phase %s' %(block_file,position_file,start,end,path_file,phase_file)
				if os.path.isfile('prism_all_block/%s' %(path_file)) and os.path.isfile('prism_all_block/%s' %(phase_file)):
					cmds = cmds
					print cmds
				else:
					print >>f_out1,cmds
#					print >>f_out2,cmds
				if len(interleave)!=0:
					while len(interleave)!=0:
						if int(interleave[-1])>int(end_id):
							interleave.pop()
						else:
							break
					if len(interleave)!=0:
						for block_id in interleave:
							cmds = 'grep "ID=%s:" %s.chr%s.local.interleaving.blocks >> %s.chr%s.all.%s.%s.interleaving.blocks' %(block_id,args.sample,chr,args.sample,chr,start_id,end_id)
							os.popen(cmds)
						unphased_block = '%s.chr%s.all.%s.%s.interleaving.blocks' %(args.sample,chr,start_id,end_id)
						phased_block = '%s.chr%s.all.%s.%s.phased.interleaving.blocks' %(args.sample,chr,start_id,end_id)
						cmds = 'prism-interleaving --unphased-blocks %s --path prism_all_block/%s --vcf %s --chr chr%s %s --phased-blocks %s' %(unphased_block,path_file,vcf_file,chr,REF_PANEL,phased_block)
						print >>f_out3,cmds
						interleave =[]
				item = [item[-1]]

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='run Prism')
	parser.add_argument("--sample", dest='sample',help="sample name")
	args = parser.parse_args()
	
	for i in range(1,23):
		if(i==23):
			i='X'
		info_file = '%s.chr%s.info' %(args.sample,i)
		out_file = '%s.chr%s.prism.cmds' %(args.sample,i)
		make_prism_cmd(info_file,out_file,1000000,i)