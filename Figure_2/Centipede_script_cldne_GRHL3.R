setwd("/Users/huanliu/Documents/ATAC_seq/re_analysis/centipede_GFP_fo/")


### Install CENTIPEDE####
install.packages("CENTIPEDE", repos="http://R-Forge.R-project.org") 
library(CENTIPEDE)  #test whether centipede is installed or not
install.packages("devtools")
devtools::install_github("slowkow/CENTIPEDE.tutorial")
#### Install Rsamtools####
source("https://bioconductor.org/biocLite.R")
biocLite("Rsamtools")


#### Run centipede####
library(CENTIPEDE)
library("Rsamtools")
library(CENTIPEDE.tutorial)
## Before running the cnetipede function, run centipede_data_2 function
cen <- centipede_data(
  bam_file = "BT2_GFP_d_m_shifted.sorted.bam",
  fimo_file = "grhl3_sites.txt.gz",
  pvalue = 1e-4,
  flank_size = 100
)

head(cen$regions)


### footprint for cldne3 locus with GRHL motif
pdf('cldne_3_grhl_motif.pdf')
plot(cen$mat[1042,], xlab = "Position", ylab = "Read Start Sites", type = "h",
     col = rep(c("blue", "red"), each = 213),ylim=c(0,35))
abline(v = c(100, 113, 313, 326) + 0.5, lty = 2)
abline(v = 213 + 0.5)
dev.off() ## edit the output using AI

