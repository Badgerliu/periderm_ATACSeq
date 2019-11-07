#!/usr/bin/env bash
echo "Now convert bam to bed"
echo "Make sure all bam files are ended with '_m_d.bam' "



for zyy in *_m_d.bam
	do
		echo $zyy
		describer=$(echo ${zyy} | sed 's/.bam//')
		echo $describer conversion begins
		bedtools bamtobed -i $describer.bam >$describer.bed
		echo ${describer} converstion finished
	done
	
echo "bam to bed conversion! --Huan"


mkdir ./peak

mv *_m_d.bed ./peak

cd ./peak


echo "Do peak calling"


for foo in *_m_d.bed
	do
		echo $foo
		bar=$(echo ${foo} | sed 's/_m_d.bed//') # define the describer with sed function
		echo $bar  
		#peak calling
		macs2 callpeak --nomodel -t ${bar}_m_d.bed -n ${bar} --nolambda --gsize 2.9e9 --keep-dup all --slocal 10000
		echo peak calling for ${bar} finished
	done

mv *_summits.bed ./peak
mv *.narrowPeak ./peak
mv *.xls ./peak


rm *_m.bed  # (optional) remove the chrM deleted bedfiles just to save disk space.
	
echo "Done with peak calling. You will find results in 'peak' folder.  --Huan"


