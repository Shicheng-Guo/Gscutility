# use minfi library to deal with mh450k array data

library("minfi")
RG.raw <- read.450k.exp(base = slide.folder, targets = files.table)
methyl.norm <- preprocessIllumina(RG.raw, bg.correct = TRUE, normalize = "controls")
beta.table <- getBeta(methyl.norm)

