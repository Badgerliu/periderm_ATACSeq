#!/usr/bin/env bash

mkdir bw

echo "Make sure you start virtualEnv for deeptools and go to the right folder before starting this shell, or you will have trouble!"

counter=0

for sample in *.bam
	do
		echo $sample
		describer=$(echo ${sample} | sed 's/_m_d.bam//') # define the describer with sed function
		echo $describer  
		#Normalization using deeptools
		bamCoverage --bam ${describer}_m_d.bam -o ${describer}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 3234832043 --ignoreForNormalization chrX chrM --extendReads
		echo Normalization with ${describer} is done
		let counter=counter+1
		echo ${counter} bamfiles have been finished
	done

mv *.bw ./bw
	
echo "Done with normalization! You will find the results in folder 'bw', and you can visualize the data in IGV now! --Huan"
