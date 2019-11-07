

library(gkmSVM)
library(ROCR)
install.packages(readr)
library(readr)
setwd("/Users/huanliu/machine_learning_plot")
### Prepare the function for calculating auPRC (subfunction of gkmsvm_trainCV())
auPRC <- function(perf) {
  rec <- perf@x.values
  prec <- perf@y.values
  result <- list()
  for (i in 1:length(perf@x.values)) {
    result[i] <- list(sum((rec[[i]][2:length(rec[[i]])] - 
                             rec[[i]][2:length(rec[[i]]) - 1]) * prec[[i]][-1]))
  }
  return(result)
}
## Evaluate fish_training kernel ####
fish_training <- read_delim("/Users/huanliu/machine_learning_plot/fish_training_1x_cvpred.out", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
colnames(fish_training) <-c("sequenceid", "preds","labs","nCV")

head(fish_training)
## computing a simple ROC curve (x-axis: fpr, y-axis: tpr)
library(ROCR) # optional 

pred <- prediction( fish_training$preds, fish_training$labs)
perf <- performance(pred,"tpr","fpr")  # tpr stands for true positive rate, fpr stands for false positive rate
pdf(file="fish_training_ROC.pdf", width=5, height=5)
plot(perf, colorize=TRUE)  # right y-axis stand for the SVM score distribution
dev.off()
perf <- performance(pred, "auc") # calculate auROC
perf@y.values[[1]]   # return the auROC=0.8785235

## precision/recall curve (x-axis: recall, y-axis: precision)
perf1 <- performance(pred, "prec", "rec")
pdf(file="fish_training_PRC.pdf", width=5, height=5)
plot(perf1,colorize=TRUE)  # right y-axis stand for the SVM score distribution
dev.off()

auPRC(perf1)   # return the auPRC=0.8745557

