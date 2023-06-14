#!/usr/bin/env Rscript
# calculate Shannon entropy at each genomic position
# calculate mean Shannon entropy of HBV genomic regions

# Dr. Chih-Jen Huang
# Genomic Research Center, Academia Sinica, Taipei, Taiwan
# May.2023

suppressPackageStartupMessages(library(readr,
                                       ggplot2)
)

WDdir="./un-normalized_shannon"

dir.create(WDdir)
setwd(WDdir)
dir.create("ind")

########################################
# Position AF from bam_readcounts file #
########################################


#############
# functions

# read bam_count data from file
Read_File_Fx <- function(sfilepath) {
  inrc  <- suppressMessages(read_delim(sfilepath,
                                       delim = "\t",
                                       col_names = F))
  colnames(inrc) <- c("Genome","Position","Ref","Counts","Q","A","C","G","T","N")
  inrc <- as.data.frame(inrc)
  return (inrc)
}

# calculate A,C,G,T,N counts in each position
Position_Count_ACGT_Fx <- suppressMessages(function (inrc) {
  data <- inrc[1:4]
  data$A_count <- as.numeric(lapply(strsplit(inrc$A, ":"), "[", 2))
  data$C_count <- as.numeric(lapply(strsplit(inrc$C, ":"), "[", 2))
  data$G_count <- as.numeric(lapply(strsplit(inrc$G, ":"), "[", 2))
  data$T_count <- as.numeric(lapply(strsplit(inrc$T, ":"), "[", 2))
  data$N_count <- as.numeric(lapply(strsplit(inrc$N, ":"), "[", 2))
  
  data$A_AF <- data$A_count/as.numeric(data$Counts)
  data$C_AF <- data$C_count/as.numeric(data$Counts)
  data$G_AF <- data$G_count/as.numeric(data$Counts)
  data$T_AF <- data$T_count/as.numeric(data$Counts)
  data$N_AF <- data$N_count/as.numeric(data$Counts)
  data$INDEL_AF <- (1-data$A_AF-data$C_AF-data$G_AF-data$T_AF-data$N_AF)
  return (data)
})

# calculate A,C,G,T,N counts in each position
pos.entropy <- function(cA, cC, cG, cT) {        #input nucleotide base counts
  cA <- as.numeric(cA)
  cC <- as.numeric(cC)
  cG <- as.numeric(cG)
  cT <- as.numeric(cT)
  sumACGT <- cA +cC +cG + cT
  pA <- cA/sumACGT
  pC <- cC/sumACGT
  pG <- cG/sumACGT
  pT <- cT/sumACGT
  eA=eC=eG=eT=0
  if (log(pA) != -Inf){
    eA <- -pA*log(pA)
  }
  if (log(pC) != -Inf){
    eC <- -pC*log(pC)
  }
  if (log(pG) != -Inf){
    eG <- -pG*log(pG)
  }
  if (log(pT) != -Inf){
    eT <- -pT*log(pT)
  }
  
  posEntropy <- (eA +eC +eG +eT)
  
  
  nucl.div.f <- (1-(pA^2+pC^2+pG^2+pT^2)) # for nucleotide diversity calculation
  return (c(posEntropy, nucl.div.f))
}


##########
## main ##
##########

#load sample list

slist <- scan("../example/Sample_name_list.txt", what="",sep = "\n") # patient ID list

shannon.summary.all <- matrix(, nrow = 0,ncol = 23) #create empty matrix for final entropy report
colnames(shannon.summary.all) <- c("Sample_Name", "WG", "P", "X", "E", "preC", "C", "Large S", "Middle S", "S", "preS1", "preS2",
                                   "WG.div", "P.div", "X.div", "E.div", "preC.div","C.div", "Large S.div", "Middle S.div", "S.div", "preS1.div", "preS2.div")
shannon.table <- matrix(,nrow=0,ncol=3215) #create empty matrix
colnames(shannon.table) <- c(1:3215)


# traverse all samples
n=1
for (i in 1:length(slist)){
  slist[i]
  sname <- slist[i]
  print(paste0("processing: ",sname,"  (", n,"/",length(slist),")")) 
  #sname <- slist[1] #example for loading single sample
  (sfilepath <- paste0("../minimap2_mapping/SH1212-C5_ref_splice-hq-realign/",sname,"/bam-readcount/",sname,".bam-readcount"))
  inrc <- Read_File_Fx(sfilepath)
  position.out <- Position_Count_ACGT_Fx(inrc) # calculate AF
  ###############################
  # Shannon entropy calculation #
  ###############################
  # calculate entropy raw by raw. Use "for" loop, or it will return an error: "Error in rep("A", cA) : 'times' 引數不正???"
  for (r in 1:nrow(position.out)){
    count4base <- position.out$A_count[r]+position.out$C_count[r]+position.out$G_count[r]+position.out$T_count[r]
    position.out$count4base[r] <- sum(as.numeric(position.out$A_count[r]),as.numeric(position.out$C_count[r]),as.numeric(position.out$G_count[r]),as.numeric(position.out$T_count[r]))
    #position.out$shannon[r] <- "NA"
    if (count4base == 0){
      position.out$shannon[r] <- "NA" #if deletion (4 bases sum = 0)
      # print("NA")
      position.out$nucl.div.f[r] <- "NA"
    }else{
      # print(position.out$count4base[r],position.out$A_count[r],position.out$C_count[r],position.out$G_count[r],position.out$T_count[r])
      pos.entropy.function.out <- pos.entropy(position.out$A_count[r],position.out$C_count[r],position.out$G_count[r],position.out$T_count[r])
      position.out$shannon[r] <- pos.entropy.function.out[1]
      position.out$nucl.div.f[r] <- pos.entropy.function.out[2]
    }
  }
  
  dir.create(paste0("./ind/",sname))
  
  # SE.plot <- ggplot(position.out,aes(Position,shannon)) + geom_step() + labs(title= paste0("Shannon Entropy of ",sname),x="Position (bp)",y="Shannon Entropy")
  # ggsave(file=paste0("./",sname,"/Shannon_Entropy_",sname,".png"), plot = SE.plot, device = "png", width = 180, height = 144, unit = "mm")
  write.csv(position.out, file = paste0("./ind/",sname,"/Position_AF_",sname,".csv"))
  
  # save entropy to table
  idx <- match( colnames(shannon.table ), position.out$Position)
  shannon.table <- rbind(shannon.table,position.out$shannon[ idx ])
  rownames(shannon.table)[i] <- sname
  #shannon.mean <- mean(as.numeric(position.out$shannon), na.rm=TRUE)
  # remove outlier at position 2334, 2335 (due to mapping issue)                                                                                                         # CDS region from SH1212-C5 gff3 file
  WG.shannon <- position.out%>% filter(Position %in% (26:3200)) %>% filter(!Position %in% c(2334,2335))     # 26-3210 (remove primer region)
  P.shannon <- position.out%>% filter(Position %in% (487:3081)) %>% filter(!Position %in% c(2334,2335))     # P: 487-3018
  # X.shannon <- rbind(position.out%>% filter(Position %in% (2769:3215)),
  #                   position.out%>% filter(Position %in% (1:18)))  #actually primer region (head)         # X: 2769-3233(18)
  X.shannon <- position.out%>% filter(Position %in% (2769:3200))
  # E.shannon <- rbind(position.out%>% filter(Position %in% (3209:3215)), #actually primer region (tail)    # E: 3209-3847(632)
  #                   position.out%>% filter(Position %in% (1:632)))
  E.shannon <- position.out%>% filter(Position %in% (26:632)) # E = preC + Core
  preC.shannon <- position.out%>% filter(Position %in% (26:80)) # head 32bp in primer region, 55bp avalible
  C.shannon <- position.out%>% filter(Position %in% (81:632))
  LS.shannon <- position.out%>% filter(Position %in% (1061:2230))                                           # pre-S1/pre-S2/S: 1061-2230
  MS.shannon <- position.out%>% filter(Position %in% (1385:2230))                                           # pre-S2/S: 1385-2230
  S.shannon <- position.out%>% filter(Position %in% (1550:2230))                                            # S: 1550-2230
  preS1.shannon <- position.out%>% filter(Position %in% (1061:1384))
  preS2.shannon <- position.out%>% filter(Position %in% (1385:1549))
  
  shannon.mean.WG <- mean(as.numeric(WG.shannon$shannon), na.rm=TRUE)
  shannon.mean.P <- mean(as.numeric(P.shannon$shannon), na.rm=TRUE)
  shannon.mean.X <- mean(as.numeric(X.shannon$shannon), na.rm=TRUE)
  shannon.mean.E <- mean(as.numeric(E.shannon$shannon), na.rm=TRUE)
  shannon.mean.preC <- mean(as.numeric(preC.shannon$shannon), na.rm=TRUE)
  shannon.mean.C <- mean(as.numeric(C.shannon$shannon), na.rm=TRUE)
  shannon.mean.LS <- mean(as.numeric(LS.shannon$shannon), na.rm=TRUE)
  shannon.mean.MS <- mean(as.numeric(MS.shannon$shannon), na.rm=TRUE)
  shannon.mean.S <- mean(as.numeric(S.shannon$shannon), na.rm=TRUE)
  shannon.mean.preS1 <- mean(as.numeric(preS1.shannon$shannon), na.rm=TRUE)
  shannon.mean.preS2 <- mean(as.numeric(preS2.shannon$shannon), na.rm=TRUE)
  
  diversity.mean.WG <- mean(as.numeric(WG.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.P <- mean(as.numeric(P.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.X <- mean(as.numeric(X.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.E <- mean(as.numeric(E.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.preC <- mean(as.numeric(preC.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.C <- mean(as.numeric(C.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.LS <- mean(as.numeric(LS.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.MS <- mean(as.numeric(MS.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.S <- mean(as.numeric(S.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.preS1 <- mean(as.numeric(preS1.shannon$nucl.div.f), na.rm=TRUE)
  diversity.mean.preS2 <- mean(as.numeric(preS2.shannon$nucl.div.f), na.rm=TRUE)
  
  shannon.summary <- matrix(c(sname,shannon.mean.WG,shannon.mean.P,shannon.mean.X,shannon.mean.E,shannon.mean.preC,shannon.mean.C,
                              shannon.mean.LS,shannon.mean.MS,shannon.mean.S,shannon.mean.preS1,shannon.mean.preS2,
                              diversity.mean.WG,diversity.mean.P,diversity.mean.X,diversity.mean.E,diversity.mean.preC,diversity.mean.C,
                              diversity.mean.LS,diversity.mean.MS,diversity.mean.S,diversity.mean.preS1,diversity.mean.preS2
  ),ncol=23,byrow=TRUE)
  colnames(shannon.summary) <- c("Sample_Name", "WG", "P", "X", "E", "preC","C", "Large S", "Middle S", "S", "preS1", "preS2",
                                 "WG.div", "P.div", "X.div", "E.div", "preC.div","C.div", "Large S.div", "Middle S.div", "S.div", "preS1.div", "preS2.div")
  rownames(shannon.summary) <- sname
  write.csv(shannon.summary, file = paste0("./ind/",sname,"/ind.shannon.diversity.summary_",sname,".csv"))
  shannon.summary.all <- rbind(shannon.summary.all,shannon.summary)
  #SE.plot <- ggplot(WG.shannon,aes(Position,shannon)) + geom_step()  + labs(title= paste0("Shannon Entropy of ",sname),x="Position (bp)",y="Shannon Entropy") # + ylim(0,1)
  #ggsave(file=paste0("./",sname,"/Shannon_Entropy_",sname,".png"), plot = SE.plot, device = "png", width = 180, height = 144, unit = "mm")
  
  n = n + 1
}

write.csv(shannon.summary.all, file = "./UnNorm_shannon.diversity.summary.csv")
write.csv(shannon.table, file = paste0("./shannon.table.csv"))
as.data.frame(shannon.table) %>% select(!c(1:25,2334:2335,3201:3215)) -> shannon.table.rm.outlier
write.csv(shannon.table.rm.outlier, file = paste0("./shannon.table.rm.outlier.csv"))

save(shannon.table,shannon.table.rm.outlier,slist, file = "./shannon.table.RData")


##########################################################################
#sliding windows test run
shannon.table.v <- as.data.frame(shannon.table[,26:3200])
colnames(shannon.table.v) <- paste0("Pos_", colnames(shannon.table.v))
shannon.table.v$Pos_2334 <- "NA"
shannon.table.v$Pos_2335 <- "NA"
shannon.table.v <- as.matrix(shannon.table.v)

library(runner)

sliding_w <- function(windowsize){ 
  shannon.table.win <- matrix(,nrow=0,ncol=3175) # create empty matrix
  colnames(shannon.table.win) <- paste0("Pos_",c(26:3200))
  
  for (i in 1:length(rownames(shannon.table.v))){ 
    shannon.table.win <- rbind(shannon.table.win, 
                               runner(                    # run sliding windows
                                 as.numeric(shannon.table.v[i,]),         # row as input
                                 k = windowsize,                    # window size = 51, slide always = 1 
                                 lag = -((windowsize-1)/2),
                                 f = function(x) mean(x, na.rm = TRUE)                   # function
                               )                          # lag = 2; data shift 2 column
                               
    )
    row.names(shannon.table.win)[i] <- row.names(shannon.table.v)[i] # fill up rownames
  }
  # runner: https://cran.r-project.org/web/packages/runner/vignettes/apply_any_r_function.html
  
  write.csv(shannon.table.win, file = paste0("shannon.table.sliding_w",windowsize,".csv"))
}

sliding_w(51)
sliding_w(31)
sliding_w(15)

##########################
# Individual Shannon plot
# scale y from 0-1

#slist <- scan("../461_e-pos_reads_1000.list.txt", what="",sep = "\n")
for (i in 1:length(slist)){
  (sname <- slist[i])
  # sname <- slist[1]
  position.out <- read.csv(paste0("./ind/",sname,"/Position_AF_",sname,".csv"))
  WG.shannon <- position.out%>% filter(Position %in% (26:3200))
  SE.plot <- ggplot(WG.shannon,aes(Position,shannon)) + geom_step() + ylim(0,1.5) + labs(title= paste0("Shannon Entropy of ",sname),x="Position (bp)",y="Shannon Entropy") # + ylim(0,1)
  ggsave(file=paste0("./ind/",sname,"/Shannon_Entropy_",sname,"_1.pdf"), plot = SE.plot, width = 180, height = 144, unit = "mm")                                        
}

for (i in 1:length(slist)){
  slist[i]
  
  sname <- slist[i]
  #sname <- "6815019"  # 2810739 2823696 4911686 5820941 6815019
  
  position.out <- read.csv(paste0("./ind/",sname,"/Position_AF_",sname,".csv"))
  WG.shannon <- position.out%>% filter(Position %in% (26:3200)) %>% filter(!Position %in% (2334:2335))
  SE.plot <- ggplot(WG.shannon,aes(Position,shannon)) + geom_step() + ylim(0,1.5) + labs(title= paste0("Shannon Entropy of ",sname),x="Position (bp)",y="Shannon Entropy") # + ylim(0,1)
  ggsave(file=paste0("./ind/",sname,"/Shannon_Entropy_",sname,"_1_rm2334,2335.pdf"), plot = SE.plot, width = 180, height = 144, unit = "mm")                                        
}

