## Periderm_ATACSeq

---------------

### Scripts for periderm_ATACSeq paper published in eLife (https://elifesciences.org/articles/51325) 



#### General workflow for ATAC-seq data process
- Environment setup (Ubuntu 18.0 or MacOS 10.12)
- QC using FastQC
- Trim the raw sequencing data using Trimmomatic (v0.36): (example paired-end raw read: C1_1.fq.gz, C1_2.fq.gz; dump unpaired output)
```java
java -jar trimmomatic-0.38.jar PE -phred33 C1_1.fq.gz C1_2.fq.gz C1_1_paired.fq.gz C1_1_unpaired.fq.gz C1_2_paired.fq.gz C1_2_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:5
```

- Mapping with Bowtie2 (example using trimmed paired reads, mapped to mm10 genome build);removal of PCR duplicates using picard tools; convert from sam output to bam output)

  ```shell
  bowtie2 -x mm10 -1 C1_1_paired.fq.gz -2 C1_2_paired.fq.gz -p 8 -S C1.sam 
  sed '/chrM/d;/random/d;/chrUn/d' C1.sam > C1_m.sam 
  samtools view -bT mm10.fa C1_m.sam > C1_m.bam 
  rm *.sam 
  samtools sort -o C1.bam C1_m.bam 
  java -jar /home/whuss/Huan/seqtool/picard/picard.jar MarkDuplicates I=C1.bam O=C1_m_d.bam M=dups.txt REMOVE_DUPLICATES=true 
  samtools index C1_m_d.bam 
  rm C1.bam 
  rm C1_m.bam
  ```

  


- Size selection using custom scripts (located in /ATACSeq_Script/)

* * 1) size_selection.py (using deeptools as environment) to generate samfiles of NucleosomeFreeRegions and monoNucleosomeRegions

* * 2) SB_conversion.sh to generate bamfiles of nucleosomeFreeRegions

* QC and plot for correlation were performed using example scrpt below:

  ```shell
  multiBamSummary bins --bamfiles *.bam --minMappingQuality 30 --region 19 -out readCounts.npz --outRawCounts readCounts.tab
  
  plotCorrelation -in readCounts.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_readCounts.png --outFileCorMatrix SpearmanCorr_readCounts.tab
  
  plotCorrelation -in readCounts.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_PearsonCorr_readCounts.png --outFileCorMatrix PearsonCorrCorr_readCounts.tab
  ```

  


- bigWig outputs using deeptools: generate_bw_mm10.sh (for mouse data), generate_bw_hg19.sh (for human data) generate_bw_zv9.sh (for zebrafish data)

- peakCalling using MACS2: shift_and_peak_calling_zv9.sh (for zebrafish data); shift_and_peak_calling_mm10.sh (for mouse data); shift_and_peak_calling_hg19.sh (for human data)

- DiffBind using EdgeR was used for detection of different accessible regions in the consensus peaks between two group, example script was used: DiffBind_mouse_palatal_epithelium_C_vs_P.R

- Motif enrichment was performed using Homer (v3) (example using C_enriched peaks)

  ```shell
  findMotifsGenome.pl C_enriched_peaks.txt mm10 C_enriched_peaks/ -size 500 -len 8,10,12 -mask -dumpFasta 
  ```

  





#### General workflow for H3K27Ac ChIP-seq data process

- Environment setup (Ubuntu 18.0)
- QC using FastQC
- Trim the raw sequencing data using Trimmomatic (v0.36): (example paired-end raw read: HIOEC_D3_H3K27Ac_1.fq.gz, HIOEC_D3_H3K27Ac_2.fq.gz; dump unpaired output)

```
java -jar trimmomatic-0.38.jar PE -phred33 HIOEC_D3_H3K27Ac_1.fq.gz HIOEC_D3_H3K27Ac_2.fq.gz HIOEC_D3_H3K27Ac_1_paired.fq.gz HIOEC_D3_H3K27Ac_1_unpaired.fq.gz HIOEC_D3_H3K27Ac_2_paired.fq.gz HIOEC_D3_H3K27Ac_2_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
```

- Mapping with Bowtie2 (example using trimmed paired reads, mapped to hg19 genome build);removal of PCR duplicates using picard tools; convert from sam output to bam output)

  ```shell
  bowtie2 -x mm10 -1 HIOEC_D3_H3K27Ac_1_paired.fq.gz -2 HIOEC_D3_H3K27Ac_2_paired.fq.gz -p 8 -S HIOEC_D3_H3K27Ac.sam 
  sed '/chrM/d;/random/d;/chrUn/d' HIOEC_D3_H3K27Ac.sam > HIOEC_D3_H3K27Ac_m.sam 
  samtools view -bT mm10.fa HIOEC_D3_H3K27Ac_m.sam > HIOEC_D3_H3K27Ac_m.bam 
  rm *.sam 
  samtools sort -o HIOEC_D3_H3K27Ac.bam HIOEC_D3_H3K27Ac_m.bam 
  java -jar /home/whuss/Huan/seqtool/picard/picard.jar MarkDuplicates I=HIOEC_D3_H3K27Ac.bam O=HIOEC_D3_H3K27Ac_m_d.bam M=dups.txt REMOVE_DUPLICATES=true 
  samtools index HIOEC_D3_H3K27Ac_m_d.bam 
  rm HIOEC_D3_H3K27Ac.bam 
  rm HIOEC_D3_H3K27Ac_m.bam
  ```

  


- bigWig outputs using deeptools: generate_bw_hg19.sh (for human data) 

- peakCalling using MACS2: shift_and_peak_calling_hg19.sh (for human data)

  

#### General workflow for gkmSVM-based machine learning process

- Scripts:  all included in **gkmSVM_script**/ folder; example script for generating vector based on GPAEs without promoters gkmSVM_fish.R.

- To generate null sequence (negative sets): fish sequences using genNullSeqs_drerio.R, mouse sequenes using genNullSeqs_mice_real.R

- Nr10mers.fa was used to probe the enriched motif(s) in the training sets.




#### General workflow for examine the pattern of motifs combination in GPAEs

- Scripts: all included in **motif_combination**/ folder.
- To perform hierarchy clustering for the motif combination pattern in GPAEs, we used hierarchy_cluster_motif_combination_pattern_in_GPAE.R
- To count the GPAEs with different sets of 3-motif combinations and 2-motif combinations, we used motif_combination_count.R



#### Example script for plot

example scripts for plot

#### Scripts for Figure 1

**Name** ./Figure_1/

Figure 1B: matrix for heatmap was generated by seqMiner and plot using R (example script: heatmap.R)

Figure 1C and D: using matrix for Figure 1B and plot by R density plot (example script density_plot.R)

Figure 1H: scatter plot using R (example script scatter_plot.R)

Figure 1I: boxplot using the normalized accessibility for each indicated group, normalization was performed using DiffBind (EdgeR algorithm)



#### Scritps for Figure 2

**Name** ./Figure_2/

Figure 2C: centipede script for footprint (Centipede_script_cldne_GRHL3.R)



#### Scripts for Figure 3

**Name** ./Figure_3/

Figure 3 B: ROC and P-R curve (ROC_and_PR.R)



#### Scripts for the density plots and heatmaps in Figure 5 are similar as those in Figure 1.


#### Other information about the deposited datasets in GEO are listed in Other data

