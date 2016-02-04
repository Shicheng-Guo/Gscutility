# This is updated version by Shicheng Guo with some bug fix to the popular version in the internet.

#-------------------------------------------------------------------------------------------
# Script for converting from plink (.ped) to fastPHASE (.inp) input format
#-------------------------------------------------------------------------------------------
setwd("/home/shg047/work/bak/data/plink")
#-------------------------------------------------------------------------------------------
# Read in datafile and map
#-------------------------------------------------------------------------------------------
test.ped="mydata.haploview.ped";
data <- read.table(test.ped, colClasses = "character", header=F)
dim(data)
data[1:10, 1:10]
test.map="mydata.map"
map <- read.table(test.map, header=F)
dim(map)
map[1:10, 1:4]
#-------------------------------------------------------------------------------------------
# Write IDs, genotypes and bp positions to separate new files
#-------------------------------------------------------------------------------------------
# IDs
data_IDs <- data[,2]
# Genotypes  
data_2 <- data[,-(1:6)]
dim(data_2)
data_2[1:10, 1:10]
# bp positions 
pos <- data.frame(map[,4])
pos2 <- rbind("P", pos)
dim(pos2)
pos2[1:10, 1]
# type 
type<-apply(data_2,2,function(x) length(table(x))>4)[seq(1,ncol(data_2),by=2)]
type<-gsub("TRUE","M",type)
type<-gsub("FALSE","S",type)
type<-paste(type,sep="",collapse="")
# Remove data and map (memory issues)
rm(data, map)
ls()
#-------------------------------------------------------------------------------------------
# Change missing genotype data code from "0" to "?"
#-------------------------------------------------------------------------------------------
# If memory is insufficient, write data_2 to file, replace 0's in ConText, and re-import
length(which(data_2==0))    			#Are there any zero values (i.e. missing data)?
data_2[data_2==0] <- "?"					#Replace 0's with ?'s
#-------------------------------------------------------------------------------------------
# Make header lines - note use of append=T in write.table to generate outfile (faster than rbind...etc)
#-------------------------------------------------------------------------------------------
numIDs <- nrow(data_2)
write.table(numIDs, "test.inp", quote=F, row.names=F, col.names=F)
numSNP <- ncol(data_2)/2
write.table(numSNP, "test.inp", quote=F, row.names=F, col.names=F, append=T)
P <- t(pos2)
write.table(P, "test.inp", quote=F, row.names=F, col.names=F, append=T)
write.table(type, "test.inp", quote=F, row.names=F, col.names=F, append=T)
#-------------------------------------------------------------------------------------------
# Add IDs and genotypes
#-------------------------------------------------------------------------------------------
for (i in 1:nrow(data_2)){
  ID <- paste("#", "ID", data_IDs[i],sep="")
  write.table(ID, "test.inp", quote=F, row.names=F, col.names=F, append=T)
  geno <- matrix(data_2[i,], nrow=2)
  write.table(geno, "test.inp", quote=F, row.names=F, col.names=F, append=T)
}

# test.inp can be fed into fastPHASE
