setwd("/Users/huanliu/elife_revision/motif_combination")
mydata <- read.csv("/Users/huanliu/elife_revision/motif_combination/GPAEs_motif_annotate.csv")



try_1<- data.frame(lapply(mydata[,-(1:4)], function(x) ifelse(x!=0, mydata$peakID, x)))  ## replace the none 0 value with the name of peak_ID in each row
try_1
try_1 <- lapply(try_1,
                 function(x) {
                   x[x==0] <- NA
                   as.character(x)
                 }) 


###Calculate the combination number of 2####
cc_2 <- combn(names(try_1),2,
            FUN=function(x) {
              data.frame(matrix(x,nrow=1),
                         val=sum(try_1[[x[1]]]==try_1[[x[2]]],na.rm=TRUE))
            },
            simplify=FALSE)
motif_comb_2<-do.call(rbind,cc_2)
write.csv(motif_comb_2, file="motif_comb_2.csv", row.names = FALSE)

###Calculate the combination number of 3####
cc_3 <- combn(names(try_1),3,
              FUN=function(x) {
                data.frame(matrix(x,nrow=1),
                           val=sum((try_1[[x[1]]]==try_1[[x[2]]]) & (try_1[[x[3]]]==try_1[[x[2]]]),na.rm=TRUE))
              },
              simplify=FALSE)
motif_comb_3<-do.call(rbind,cc_3)
write.csv(motif_comb_3, file="motif_comb_3.csv", row.names = FALSE)


