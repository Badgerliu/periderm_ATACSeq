flow: true
sequence:true
### Periderm_ATACSeq

---------------　

**Scripts for periderm_ATACSeq paper published in..**
[] ok
[x]

```flow
st=>start: Start
op=>operation: Your Operation
cond=>condition: Yes or No?
e=>end
st->op->cond
cond(yes)->e
cond(no)->op
```


####General workflow for ATAC-seq data process

```flow
st=>start: raw paired-end fastq files
op1=>operation: Trimmomatic
op2=>operation:bowtie2 mapping
op3=>operation:Size_selection.py
op4=>operation:SB_conversion.sh
cond=>condition:output?
op5=>operation:generate_bw.sh
op6=>operation:shift_and_peak_calling.sh
e1=>end:bw files
e2=>end:peak files
st->op1->op2->op3->op4->cond
cond(bw)->op5->e1
cond(peak)->op6->e2
```

#### Scripts for Figure 1

**Name** ./Figure_1/

**Comments** Used for generating heatmap for Figure 1B


