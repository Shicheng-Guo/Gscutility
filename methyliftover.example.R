liftover450k <- function(fileDir, fileKeys, fileSkipLine = 1, outputFileDir, outputFileName){
  
  # Load the Illumina 450k Methylation Array for merging
  illumData <- IlluminaHumanMethylation450kanno.ilmn12.hg19@data
  illumDataFrame <- illumData$Locations
  illumDataFrame <- cbind(illumDataFrame, illumData$Islands.UCSC )
  illumDataFrame$Type <- illumData$Manifest$Type
  illumDataFrame$probe <- rownames(illumDataFrame)
  illumDataTable <- data.table(as.data.frame(illumDataFrame))
  
  inputFile <- fread(fileDir, skip = fileSkipLine)
  WGBSliftover <- merge(illumDataTable, inputFile, by.x=c('chr', 'pos'), by.y=fileKeys)
  rownames(WGBSliftover) <- WGBSliftover$probe
  setwd(outputFileDir)
  fileName <- paste(outputFileName, '.RData', sep='')
  save(WGBSliftover, file = fileName)
}

file=list.files(pattern="*bedGraph$")
dir="/home/shg047/oasis/twin/MethylFreq/"
for(i in 1:length(file)){
f=paste(dir,file[i],sep="")
liftover450k(f,fileKeys=c("V1","V2"),fileSkipLine = 1,outputFileDir=getwd(),outputFileName=file[i])
}
