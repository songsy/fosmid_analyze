#!/usr/bin/env bash

NAME=$1
for i in `seq 1 22`;
	do
		input=$NAME.chr$i.prism.cmds
		tot=$(cat $input|wc -l)
		PBS='run-steps-shiya.pbs'
		NEWPBS=run-steps-shiya-$i.pbs
		sed "s/NUM/$tot/" $PBS > $NEWPBS
		sed -i "s/NAME/$NAME/g" $NEWPBS
		sed -i "s/CHROM/chr$i/g" $NEWPBS
	done