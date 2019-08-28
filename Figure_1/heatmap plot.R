setwd("/Users/huanliu/R_script_example/")

library(gplots)
library(pheatmap)

###heatmap of mouse periderm NFR against mouse E145 H3K27Ac ChIP-seq#####
density <- read.delim("P_enriched_NFR_E145_H3K27Ac.txt")
m<- as.matrix(density[2:nrow(density), 2: ncol(density)])
range(m)   ##0, 216
range(log2(m+1))  ### 7.761551
bk =unique(c(seq(-0.1,6, length=100),seq(6,8, length=100)))
hmcols<-colorRampPalette(c("white",'#1ad853'))(length(bk)-1)
ncol(density)  ## 201

nrow(density)  ##14256

## mouse periderm NFR against mouse E145 H3K27Ac ChIP-seq
m1<-as.matrix(density[2:nrow(density),2:201])
range(m1) # (0,216)
m1<- log2(m1+1)
range(m1) #(0,10.37395)  ## 0.000000 7.761551
hist(m1)   ## check the distribution of m1



m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("mouse_palatal_periderm_NFR_E145_facial_H3K27Ac_1.png", width=301, height=1452)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off()


###heatmap of mouse periderm NFR_fo against mouse E145 H3K27Ac ChIP-seq#####
density <- read.delim("P_enriched_NFR_fo_E145_H3K27Ac.txt")
m<- as.matrix(density[2:nrow(density), 2: ncol(density)])
range(m)   ##0, 216
range(log2(m+1))  ### 7.761551
bk =unique(c(seq(-0.1,6, length=100),seq(6,8, length=100)))
hmcols<-colorRampPalette(c("white",'#c51b8a'))(length(bk)-1)
ncol(density)  ## 201

nrow(density)  ##6069

## mouse periderm NFR_fo against mouse E145 H3K27Ac ChIP-seq
m1<-as.matrix(density[2:nrow(density),2:201])
range(m1) # (0,216)
m1<- log2(m1+1)
range(m1) #(0,10.37395)  ## 0.000000 7.761551




m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("mouse_palatal_periderm_NFR_fo_E145_facial_H3K27Ac.png", width=301, height=606)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off()


###heatmap of mouse periderm NFR_no_fo against mouse E145 H3K27Ac ChIP-seq#####
density <- read.delim("P_enriched_NFR_no_fo_E145_H3K27Ac.txt")
m<- as.matrix(density[2:nrow(density), 2: ncol(density)])
range(m)   ##0, 190
range(log2(m+1))  ### 7.577429
bk =unique(c(seq(-0.1,6, length=100),seq(6,8, length=100)))
hmcols<-colorRampPalette(c("white",'#c51b8a'))(length(bk)-1)
ncol(density)  ## 201

nrow(density)  ##9640

## mouse periderm NFR_no_fo against mouse E145 H3K27Ac ChIP-seq
m1<-as.matrix(density[2:nrow(density),2:201])
range(m1) # (0,216)
m1<- log2(m1+1)
range(m1) #(0,7.577429)  ## 0.000000 7.577429




m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("mouse_palatal_periderm_NFR_no_fo_E145_facial_H3K27Ac.png", width=301, height=964)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off()





####heatmap of mouse periderm NFR_fo against mouse palatal ATAC-seq P_enriched####
density <- read.delim("P_enriched_fo_P2density.txt")
m<- as.matrix(density[2:nrow(density), 2: ncol(density)])
range(m)   ##0, 3701
range(log2(m+1))  ### 11.85409
bk =unique(c(seq(-0.1,8, length=100),seq(8,12, length=100)))
hmcols<-colorRampPalette(c("white",'#1ad853'))(length(bk)-1)
ncol(density)  ## 201

nrow(density)  ##6079

## mouse periderm NFR against mouse palatal periderm NFRS
m1<-as.matrix(density[2:nrow(density),2:ncol(density)])
range(m1) # (0,3701)
m1<- log2(m1+1)
range(m1) #(0,7.577429)  ## 0.000000 11.85409



m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("P_fo_P2_NFR.png", width=301, height=607)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off()

####heatmap of mouse periderm NFR_no_fo against mouse palatal ATAC-seq P_enriched####
density <- read.delim("P_enriched_no_fo_P2density.txt")
m<- as.matrix(density[2:nrow(density), 2: ncol(density)])
range(m)   ##0, 2708
range(log2(m+1))  ### 11.40354
bk =unique(c(seq(-0.1,8, length=100),seq(8,12, length=100)))
hmcols<-colorRampPalette(c("white",'#1ad853'))(length(bk)-1)
ncol(density)  ## 201

nrow(density)  ##9640

## mouse periderm NFR against mouse palatal periderm NFRS
m1<-as.matrix(density[2:nrow(density),2:ncol(density)])
range(m1) # (0,3701)
m1<- log2(m1+1)
range(m1) #(0,7.577429)  ## 0.000000 11.40354
hist(m1)   ## check the distribution of m1



m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("P_no_fo_P2_NFR.png", width=301, height=964)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off()














###heatmap of P enriched NFR against P2 or C2 NFR#####
density <- read.delim("P_enriched_P2_C2_density.txt")
m<- as.matrix(density[2:nrow(density), 2: ncol(density)])
range(m)
range(log2(m+1))  ### 11.85409
bk =unique(c(seq(-0.1,8, length=100),seq(8,12, length=100)))
hmcols<-colorRampPalette(c("white",'#1ad853'))(length(bk)-1)
ncol(m1)

## P_enriched_NFR against P2_NFR.bam
m1<-as.matrix(density[2:nrow(density),2:201])
range(m1) # (0,3701)
m1<- log2(m1+1)
range(m1) #(0,11.85409)
hist(m1)
nrow(m1)


m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("P_enriched_NFR_P2_NFR.png", width=301, height=1428)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off() ## notice quartz_off_screen 
dev.off()



## P_enriched_NFR against C2_NFR.bam
m1<-as.matrix(density[2:nrow(density),207:406])
range(m1) # (0,1326)
m1<- log2(m1+1)
range(m1) #(0,10.37395)
hist(m1)
nrow(m1)


m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("P_enriched_NFR_C2_NFR.png", width=301, height=1427)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off() ## notice quartz_off_screen 
dev.off()




###heatmap of P enriched NFR against P2 or C2 NFR#####
density <- read.delim("C_enriched_P2_C2_density.txt")
m<- as.matrix(density[2:nrow(density), 2: ncol(density)])
range(m)
range(log2(m+1))  ### 11.85409
bk =unique(c(seq(-0.1,8, length=100),seq(8,11, length=100)))
hmcols<-colorRampPalette(c("white",'#1ad853'))(length(bk)-1)
ncol(m1)

## C_enriched_NFR against P2_NFR.bam
m1<-as.matrix(density[2:nrow(density),2:201])
range(m1) # (0,1282)
m1<- log2(m1+1)
range(m1) #(0,10.32531)
hist(m1)
nrow(m1) #3888


m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("C_enriched_NFR_P2_NFR.png", width=301, height=389)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off() ## notice quartz_off_screen 
dev.off()



## C_enriched_NFR against C2_NFR.bam
m1<-as.matrix(density[2:nrow(density),207:406])
range(m1) # (0,1326)
m1<- log2(m1+1)
range(m1) #(0,10.37395)
hist(m1)
nrow(m1) #3888


m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("C_enriched_NFR_C2_NFR.png", width=301, height=388)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off() ## notice quartz_off_screen 
dev.off()



####Plot density for different cluster of NFR against H3K27Ac ChIP-seq####
mouse_periderm_NFR_H3K27Ac_T <- read.csv("~/R_script_example/mouse_periderm_NFR_H3K27Ac_T.csv")
density<-mouse_periderm_NFR_H3K27Ac_T
range(density$NFR_fo_H3K27Ac, density$R_no_fo_H3K27Ac)  ## 1.948, 10.479
pdf("NFR_H3K27Ac_readDensity.pdf")
plot(density$NFR_fo_H3K27Ac, type="l", ylim=c(1,12), lwd=2, col="#c51b8a")
lines(density$R_no_fo_H3K27Ac, lty=2, lwd=2, col="#c51b8a")

legend("topleft", c("Cluster 1","Cluster 2"), cex=0.8, col="#c51b8a", lty=c(1,2), lwd=2, bty="n")
dev.off()



####Plot density for different cluster of NFR against ATAC-seq####
density<- read.csv("~/R_script_example/mouse_periderm_NFR_cluster_density_ATAC.csv")

range(density$P_fo_P2, density$P_no_fo_P2)  ## 5.052282 54.721336
pdf("NFR_ATAC_readDensity.pdf")
plot(density$P_fo_P2, type="l", ylim=c(1,60), lwd=2, col="#1ad853")
lines(density$P_no_fo_P2, lty=2, lwd=2, col="#1ad853")

legend("topleft", c("Cluster 1","Cluster 2"), cex=0.8, col="#1ad853", lty=c(1,2), lwd=2, bty="n")
dev.off()


density <- read.delim("C1_CandP.txt")
m<- as.matrix(density[2:nrow(density), 2: ncol(density)])
range(m)   ##0, 3701
range(log2(m+1))  ### 11.5228
bk =unique(c(seq(-0.1,0.5, length=100),seq(0.5,11, length=100)))
hmcols<-colorRampPalette(c("white",'#1ad853'))(length(bk)-1)
ncol(density)  ## 201

nrow(density)  ##6079

## top20K_C1_vs_C_enriched_NFR
m1<-as.matrix(density[2:nrow(density),2:201])
range(m1) # (0,3701)
m1<- log2(m1+1)
range(m1) #(0,6.6)
hist(m1)
nrow(m1)


m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("top20K_C1_vs_C_enriched_NFR.png", width=301, height=718)
pheatmap(m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off() ## notice quartz_off_screen 
dev.off()



## top20K_C1_vs_P_enriched_NFR
m1<-as.matrix(density[2:nrow(density),207:406])
range(m1) # (0,1326)
m1<- log2(m1+1)
range(m1) #(0,1137395)
hist(m1)
nrow(m1)


m.row.sum<-cbind(m1, rowSums(m1))

o1<-rev(order(m.row.sum[,ncol(m.row.sum)]))
m.row.sum<- m.row.sum[o1,]
png("top20K_C1_vs_P_enriched_NFR.png", width=301, height=717)
pheatmap( m.row.sum[,1:(ncol(m.row.sum)-1)], cluster_rows = F, cluster_cols = F, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE)
dev.off() ## notice quartz_off_screen 
dev.off()


