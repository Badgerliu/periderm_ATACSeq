### Periderm_ATACSeq

---------------　

**Scripts for periderm_ATACSeq paper published in..**



####General workflow for ATAC-seq data process

- Trim the raw sequencing data using Trimmomatic:
```
java -jar trimmomatic-0.38.jar PE -phred33 C1_1.fq.gz C1_2.fq.gz C1_1_paired.fq.gz C1_1_unpaired.fq.gz C1_2_paired.fq.gz C1_2_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:5
```


#### Scripts for Figure 1

**Name** ./Figure_1/

**Comments** Used for generating heatmap for Figure 1B


