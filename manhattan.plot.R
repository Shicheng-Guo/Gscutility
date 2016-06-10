
library(qqman)
res <- read.table("MVP_ALL_CvsT_BHadjust.txt", header=TRUE,sep="\t")
SNP=res$probeID
CHR=res$CHR
BP=res$MAPINFO
P=res$P.Value
manhattaninput=data.frame(SNP,CHR,BP,P)
pdf("manhattan.pdf")
manhattan(manhattaninput,col = c("blue4", "orange3"),ylim = c(0,25),genomewideline=-log10(1e-17),lwd=2, suggestiveline=F)
dev.off()
