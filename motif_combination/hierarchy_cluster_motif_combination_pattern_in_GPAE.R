#### install packages####

install.packages('dendextend')
install.packages('gplots')
install.packages('NbClust')
setwd("/Users/huanliu/elife_revision/motif_combination")
motif1 <- read.csv("/Users/huanliu/elife_revision/motif_combination/GPAEs_motif_annotate.csv")


#### Determine the number of clusters ####

library(NbClust)
motif=as.matrix(motif1[,5:11])
scaled_motif <-scale(motif)
NbClust(data = motif, diss = NULL, distance = "euclidean",
        min.nc = 2, max.nc = 25, method = "ward.D2")

NbClust(data = scaled_motif, diss = NULL, distance = "euclidean",
        min.nc = 2, max.nc = 25, method = "ward.D2")

######HC of major motif########
motif=as.matrix(motif1[,5:11])
motif_major <- subset(motif, select = c("GRHL", "TEAD", "KLF",
                                           "FOS", "TFAP2", "GATA", "CEBP"))

set.seed(12)


d_motif<-dist(motif_major, method="euclidean")

hc_d_motif<-hclust(d_motif, method="complete")
min(d_motif) #0
max(d_motif) #19.62142
hist(d_motif)
####set up the color of heatmap
my_palette<- colorRampPalette(c("#ffffff", "#fcc5c0", "#fa9fb5","#f768a1","#dd3497","#ae017e","#7a0177","#49006a"))(n = 32)
color_break<-c(seq(0,0.9,length=9), seq(1,2.9, length=8), seq(3,5.9,length=8), seq(6.9,19, length=8))

dend_motif<-as.dendrogram(hc_d_motif)
library(dendextend)
####Set up the color of branches (10 branches)
cols_branches<-c("#d53e4f", "#f46d43", "#fdae61",
                 "#fee08b", "#ffffbf", "#e6f598",
                 "#abdda4", "#66c2a5","#3288bd","#542788") #10branches

dend_motif<-color_branches(dend_motif, k=10,col=cols_branches) # cuttree in to 10 clusters
source("https://raw.githubusercontent.com/talgalili/dendextend/master/R/attr_access.R")
col_labels<-get_leaves_branches_col(dend_motif)
col_labels <- col_labels[order(order.dendrogram(dend_motif))]
library(gplots)
pdf(file="h_clust_GPAE_motif_combination.pdf", height = 100, width = 100)
heatmap.2(motif_major,  
          main = "GPAE_motif_combination_hc",  
          trace="none",          
          margins =c(5,7),      
          col=my_palette,        
          breaks=color_break,     
          dendrogram="row",      
          Rowv = dend_motif,  
          Colv = "NA", 
          keysize = 0.5,
          key.xlab = "Concentration (index)",
          cexRow =0.2,
          cexCol = 0.4,
          na.rm = TRUE,
          RowSideColors = col_labels, # to add nice colored strips        
          colRow = col_labels # to add nice colored labels - only for qplots 2.17.0 and higher
) 
dev.off()
#export hc result using cutree
cutcutcut<-cutree(dend_motif, 2:10)
write.csv(cutcutcut, file="GPAE_motif_combination_hc_10_repeat.csv")

