setwd("/Users/huanliu/Documents/periderm_profile/integration_periderm_NFR/")


accessibility <- read.csv("z9_target_G10_differential_nfr_fo_accessibility.csv")
accessibility$labels<-factor(accessibility$labels, levels = c("PSAE", "nonPSAE"))
pdf(file = "accessibility_in_periderm.pdf", height=5, width=3)
boxplot(GFP_score~labels,data = accessibility,
        varwidth=TRUE,
        col=ifelse(accessibility$labels == c("PSAE"), "green","blue"),
        outline=FALSE,
        names=c("PSAE","non-PSAE"),
        horizontal=FALSE,
        ylim=c(0,12))
dev.off()
wilcox.test(GFP_score ~ labels, data=accessibility)

pdf(file = "accessibility_in_non-periderm.pdf", height=5, width=3)
boxplot(noGFP_score~labels,data = accessibility,
        varwidth=TRUE,
        col=ifelse(accessibility$labels == c("PSAE"), "green","blue"),
        outline=FALSE,
        names=c("PSAE","non-PSAE"),
        horizontal=FALSE,
        ylim=c(0,12))
dev.off()

wilcox.test(noGFP_score ~ labels, data=accessibility)
