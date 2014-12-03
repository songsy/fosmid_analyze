import os
import sys
import commands
sample_name = sys.argv[1]
f_out = open('%s.prism.interleave.redo.cmds' %(sample_name),'w')
for i in range(1,23):
	cmd_file = '%s.%i.prism.interleave.v2.cmds' %(sample_name,i)
	f=open(cmd_file,'r')
	for line in f:
		line=line.strip().split(' ')
		file = line[-3]
		print file
		if not os.path.isfile(file):
			cmds = 'grep %s %s ' %(file,cmd_file)
			cmds = commands.getoutput(cmds)
#			cmds += ' --pickle %s_%i_pickle' %(sample_name,i)
			print cmds
			print >>f_out,cmds
		else:
			cmds = 'cat %s' %(file)
			if commands.getoutput(cmds)=='':
				print 'empty',file
				cmds = 'grep %s %s ' %(file,cmd_file)
				cmds = commands.getoutput(cmds)
#				cmds += ' --pickle %s_%i_pickle' %(sample_name,i)
#				print cmds
#				print >>f_out,cmds