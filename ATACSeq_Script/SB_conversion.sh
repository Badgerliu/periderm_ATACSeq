#!/usr/bin/env bash

echo 'now convert sam files from size_selection'

for sample in *.sam  
 do  
   echo $sample  
   describer=$(echo ${sample} | sed 's/.sam//')  
   echo $describer  
   
   # Convert file from SAM to BAM format  
   samtools view -bT mm10.fa $sample > ${describer}.uns.bam  
   
   # Sort BAM file  
   samtools sort -o ${describer} ${describer}.uns.bam    
   
   # index the bam file  
   samtools index ${describer} 
   
   # Remove intermediate files  
   rm ${describer}.uns.bam  
 done 
