### Periderm_ATACSeq

---------------　

**Scripts for periderm_ATACSeq paper published in..**



####General workflow for ATAC-seq data process
- Environment setup

- Trim the raw sequencing data using Trimmomatic: (example paired-end raw read: C1_1.fq.gz, C1_2.fq.gz; dump unpaired output)
```
java -jar trimmomatic-0.38.jar PE -phred33 C1_1.fq.gz C1_2.fq.gz C1_1_paired.fq.gz C1_1_unpaired.fq.gz C1_2_paired.fq.gz C1_2_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:5
```

- Mapping with Bowtie2 (example using trimmed paired reads, mapped to mm10 genome build);removal of PCR duplicates using picard tools; convert from sam output to bam output


- Size selection using custom scripts
-- size_selection.py (using deeptools as environment) to generate samfiles of NucleosomeFreeRegions and monoNucleosomeRegions
-- SB_conversion.sh to generate bamfiles of nucleosomeFreeRegions

- To get bigWig outputs using deeptools

- To 
#### Scripts for Figure 1

**Name** ./Figure_1/

**Comments** Used for generating heatmap for Figure 1B


