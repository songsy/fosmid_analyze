#!/usr/bin/env python
# Shiya Song
# 7th July 2014
# msmc_plot.py

import sys
import os
mu = 1.2e-08
G = 25

def plot_gene_flow(file,name):
	f = open(file,'r')
	outfile = file.replace("final","scale")
	f_out = open(outfile,'w')
	for line in f:
		line = line.strip()
		if line[0]=='t':
			print >>f_out,line+'\t'+'relative_coalescent_rate'
			continue
		line = line.split('\t')
		line[1]= float(line[1])/mu*G
		line.append(2*float(line[4])/(float(line[3])+float(line[5])))
		print >>f_out,'\t'.join(map(str,line))

	f_out.close()
	'''
	if os.path.isfile(name+'.gp'):
		gp_file = open(name+'.gp','a')
		gp_file.write('"%s" u 2:7 t "%s" w st ls 1,' %(outfile,name))
	else:
		gp_file = open(name+'.gp','w')
		gp_file.write('set size 1, 0.8;\n')
		gp_file.write('set xran [0:2.5e5];\n')
		gp_file.write('set mxtics 10;\n')
		gp_file.write('set mytics 10;\n')
		gp_file.write('set grid;\n')
		gp_file.write('set key right top;\n')
		gp_file.write('set xtics font "Helvetica,16";\n')
		gp_file.write('set ytics nomirror font "Helvetica,16";\n')
		gp_file.write('set xlab "Years (g=25, {/Symbol m}=1.2x10^{-8})" font "Helvetica,16";\n')
		gp_file.write('set t po eps enhance so co "Helvetica,16";\n')
		gp_file.write('set size 1, 0.8;\n')
		gp_file.write('set yran [0:1];\n')
		gp_file.write('set ylab "relative cross coalsecent rate" font "Helvetica,16";\n')
		gp_file.write('set out "%s.eps";\n' %(name))
		gp_file.write('set style line 1 lt 1 lc rgb "#FF0000" lw 4;\n')
		gp_file.write('set style line 2 lt 1 lc rgb "#00C000" lw 4;\n')
		gp_file.write('set style line 3 lt 1 lc rgb "#0080FF" lw 4;\n')
		gp_file.write('plot "%s" u 2:7 t "%s" w st ls 1;\n' %(outfile,name))
	'''
	
if __name__=="__main__":
	file = sys.argv[1]
	name = sys.argv[2]
	plot_gene_flow(file,name)
		