#!/usr/bin/env Rscript
# calculate position allele frequency (A,T,G,C/sum)

suppressPackageStartupMessages(
  library(tidyverse))

#WDdir="./AF_counts"
WDdir="./AF_counts"
dir.create(WDdir)
setwd(WDdir)

#############
# functions 

# read bam_count data from file
Read_File_Fx <- function(sfilepath) {
  inrc  <- read_delim(sfilepath,
                      delim = "\t", 
                      col_names = F)
  colnames(inrc) <- c("Genome","Position","Ref","Counts","Q","A","C","G","T","N")
  inrc <- as.data.frame(inrc)
  return (inrc)
}

# calculate A,C,G,T,N counts in each position
Position_Count_ACGT_Fx <- function (inrc) {
  data <- inrc[1:4]
  data$A_count <- as.numeric(lapply(strsplit(inrc$A, ":"), "[", 2))
  data$C_count <- as.numeric(lapply(strsplit(inrc$C, ":"), "[", 2))
  data$G_count <- as.numeric(lapply(strsplit(inrc$G, ":"), "[", 2))
  data$T_count <- as.numeric(lapply(strsplit(inrc$T, ":"), "[", 2))
  # data$N_count <- as.numeric(lapply(strsplit(inrc$N, ":"), "[", 2))
  data$ACGT_count <- data$A_count+data$C_count+data$G_count+data$T_count
  
  data$A_AF <- data$A_count/(data$ACGT_count)
  data$C_AF <- data$C_count/(data$ACGT_count)
  data$G_AF <- data$G_count/(data$ACGT_count)
  data$T_AF <- data$T_count/(data$ACGT_count)
  # data$N_AF <- data$N_count/as.numeric(data$Counts)
  # data$INDEL_AF <- (1-data$A_AF-data$C_AF-data$G_AF-data$T_AF-data$N_AF)
  return (data)
}


##########
## main ##
##########

#load sample list
slist <- scan("../example/Sample_name_list.txt", what="",sep = "\n")

colnameslist <- list()
for (i in 1:3215){
  j=i*4-3
  k=i*4-2
  l=i*4-1
  m=i*4
  colnameslist[j] <- paste0(i,"_A")
  colnameslist[k] <- paste0(i,"_C")
  colnameslist[l] <- paste0(i,"_G")
  colnameslist[m] <- paste0(i,"_T")
}
colnameslist[3215*4+1] <- "max.count"
colnameslist[3215*4+2] <- "ID_new"
AF.summary <-  data.frame(matrix(ncol = 3215*4+2, 
                                 nrow = length(slist))
)
colnames(AF.summary) <- colnameslist
AF.all <-  data.frame(matrix(ncol = 3215*4+2, 
                             nrow = 0)
)
colnames(AF.all) <- colnameslist
datalist = list()

dir.create("./ind.AF")
n=1
for (sload in 1:length(slist)){
  sname <- slist[sload] #sname <- slist[1]
  print(paste0("processing: ",sname,"  (", n,"/",length(slist),")")) 
  
  sfilepath <- paste0("../minimap2_mapping/SH1212-C5_ref_splice-hq-realign/",sname,"/bam-readcount/",sname,".bam-readcount")
  
  suppressMessages(
    inrc <- Read_File_Fx(sfilepath))
  suppressMessages(
    position.out <- Position_Count_ACGT_Fx(inrc))
  
  for (pointer in 1:length(position.out$Position)){
    pos.label <- position.out$Position[pointer]
    AF.all[1,pos.label*4-3] <- position.out$A_AF[pointer]
    AF.all[1,pos.label*4-2] <- position.out$C_AF[pointer]
    AF.all[1,pos.label*4-1] <- position.out$G_AF[pointer]
    AF.all[1,pos.label*4] <- position.out$T_AF[pointer]
    AF.all[3215*4+1] <- max(position.out$Counts)
    AF.all[3215*4+2] <- sname
  }
  write.csv(AF.all, file = paste0("./ind.AF/",sname,"_ind.AF.csv"))
  datalist[[sload]] <- AF.all
  n = n + 1
}

AF.summary = do.call(rbind,datalist)
write.csv(AF.summary, file = paste0("./AF.summary.csv"), row.names = F)
save(AF.summary, file = "./AF.summary.RData")

