library(qqman)
res <- read.table("MVP_ALL_CvsT_BHadjust.txt", header=TRUE,sep="\t")
SNP=res$probeID
CHR=res$CHR
BP=res$MAPINFO
P=res$P.Value
manhattaninput=data.frame(SNP,CHR,BP,P)
pdf("manhattan.pdf")
manhattan(manhattaninput)
dev.off()
