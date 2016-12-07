#### iEVORA.R
#### Description: an R-script to perform feature selection on a large omic data set using the iEVORA algorithm. The algorithm is aimed at large DNA methylation data sets, where one may wish to find features (mostly CpGs) which differ between two normal cellular phenotypes, but with one of these phenotypes representing cells which are at risk of neoplastic transformation. The algorithm uses a differential variability (DV) step to increase power/sensitivity, but uses a standard t-test to rank significant DV CpGs (DVCs). This last step is done to regularize the DV-test, which is overly sensitive to single outliers.
#### Author: Andrew E Teschendorff (a.teschendorff@ucl.ac.uk)
#### Date: 26th Nov.2015
#### Copyright 2015 Andrew Teschendorff
#### Copyright permission: iEVORA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details (see <http://www.gnu.org/licenses/> )

#### LIBRARIES NEEDED
library(qvalue);

#### AUXILIARY FUNCTIONS

doDV <- function(tmp.v,pheno.v){
    co.idx <- which(pheno.v==0);
    ca.idx <- which(pheno.v==1);
    bt.o <- bartlett.test(x=tmp.v,g=pheno.v);
    pv <- bt.o$p.value;
    logR <- log2(var(tmp.v[ca.idx])/var(tmp.v[co.idx]));
    avCA <- mean(tmp.v[ca.idx]);
    avCO <- mean(tmp.v[co.idx]);
    out.v <- c(logR,pv,avCA,avCO);
    names(out.v) <- c("log(V1/V0)","P(BT)","Av1","Av0");
    return(out.v);
}


doTT <- function(tmp.v,pheno.v){
    tt.o <- t.test(tmp.v ~ pheno.v);
    out.v <- c(-tt.o$stat,tt.o$p.val);
    names(out.v) <- c("t","P");
    return(out.v);
}

#### MAIN USER FUNCTION

#### INPUT ARGS:
#### data.m: data matrix with rows labeling features, and columns labeling samples. Rownames should be feature/probe IDs and should be provided.
#### pheno.v: a binary phenotype sample vector with entries either 0 or 1 and in the same order as the columns of data.m
#### thDV: a significance q-value (FDR) threshold for the differential variability test. By default this is 0.001.
#### thDM: a significance p-value threshold for differential means. By defaul this is 0.05.

#### OUTPUT ARGS:
#### topDVMC.m: a matrix of ranked differentially variable (DV) and differentially methylated CpGs (DVMCs), ranked according to the t-statistic P-value, but selected using the Bartlett's DV test. Columns label the t-statistic, its P-value, the mean of phenotype-1, the mean of phenotype-0, the log-ratio of the variances of phenotype-1 to phenotype-0, the Bartlett's test P-value and q-value.

iEVORA <- function(data.m,pheno.v,thDV=0.001,thDM=0.05){

    statDVC.m <- t(apply(data.m,1,doDV,pheno.v));
    print("Estimated DV statistics");
    qvDVC.v <- qvalue(statDVC.m[,2])$qval;
    dvc.idx <- which(qvDVC.v < thDV);
    nDVC <- length(dvc.idx);
    if( nDVC > 0 ){
       statDMC.m <- t(apply(data.m[dvc.idx,],1,doTT,pheno.v))
       print("Preparing output");
       tmp.s <- sort(statDMC.m[,2],decreasing=FALSE,index.return=TRUE);
       pvDMC.v <- tmp.s$x;
       ntop <- length(which(pvDMC.v < 0.05));
       if(ntop > 0){
           topDVMC.m <- cbind(statDMC.m[tmp.s$ix[1:ntop],],statDVC.m[dvc.idx[tmp.s$ix[1:ntop]],c(3:4,1:2)],qvDVC.v[dvc.idx[tmp.s$ix[1:ntop]]]);
           colnames(topDVMC.m) <- c("t","P(TT)","Av1","Av0","log[V1/V0]","P(BT)","q(BT)");
           rownames(topDVMC.m) <- rownames(statDMC.m)[tmp.s$ix[1:ntop]];
       }
       else {
         print("NO DVMCs IN THIS DATA!");
       }
     }
     else {
         print("NO DVCs! TRY RELAXING THE thDV THRESHOLD!");
     }

     return(topDVMC.m);
}

#### Example of use
#### source("iEVORA.R");
#### topDVMC.m <- iEVORA(data.m,pheno.v)
