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


mkdir ./shift
for sample in *_m_d.bed
	do
		echo $sample
		describer=$(echo ${sample} | sed 's/.bed//')
		echo $describer  
		#shift the read position
		awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $3 + 4, $4, $5, $6; else print $1, $2 - 5, $3 - 5, $4, $5, $6}' ${describer}.bed >${describer}_shifted.bed
		echo ${describer} shift finished
	done
	
echo "Done with position shift, and now moving shifted files to shift folder! --Huan"

mv *_shifted.bed ./shift

cd ./shift

mkdir ./peak
echo "Do peak calling"


for foo in *_m_d_shifted.bed
	do
		echo $foo
		bar=$(echo ${foo} | sed 's/_m_d_shifted.bed//') # define the describer with sed function
		echo $bar  
		#peak calling
		macs2 callpeak --nomodel -t ${bar}_m_d_shifted.bed -n ${bar} --nolambda --gsize 2.9e9 --keep-dup all --slocal 10000
		echo peak calling for ${bar} finished
	done

mv *_summits.bed ./peak
mv *.narrowPeak ./peak
mv *.xls ./peak


rm *_m.bed  # (optional) remove the chrM deleted bedfiles just to save disk space.
	
echo "Done with peak calling. You will find results in 'peak' folder.  --Huan"


