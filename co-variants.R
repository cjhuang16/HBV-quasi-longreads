#!/usr/bin/env Rscript
# extract nt sequence at multiple positions (co-variants) from bam file

library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(openxlsx)

WDdir="./co-variants"
dir.create(WDdir)
setwd(WDdir)
dir.create("./ind/")
#############
# functions #
#############

# genomic position liftover
.PosToIS <- function(pos){
  if(pos <= 536){
    IS.pos = pos+1820
  }else if(pos <= 1395){
    IS.pos = pos+1826
  }else if(pos <= 3215){
    IS.pos = pos-1395
  }else{
    print("Length over 3215")
    next
  }
  return(IS.pos)
}
.PosToSH <- function(pos){
  if (pos <= 1820){
    SH.pos <- pos + 1395
  }else if(pos <= 2356){
    SH.pos <- pos - 1820
  }else if(pos <= 2362){
    SH.pos <- "del"
  }else{
    SH.pos <- pos -1826
  }
  return(SH.pos)
}

# extract nt bases at specific site from bam file
.NtFromBam <- function(sample.name){ 
  dir.create(paste0(WDdir,"/ind/",sample.name))
  bamFile <- BamFile(paste0("../minimap2_mapping/SH1212-C5_ref_splice-hq-realign/",sample.name,"/",sample.name,"_SH1212-C5_minimap2_filter_final.clip10.sorted.bam"))
  aln = GenomicAlignments::stackStringsFromBam(bamFile, 
                                               index=bamFile,  
                                               param=param, use.names=T, D.letter="-", N.letter=".", Lpadding.letter="+", Rpadding.letter="+")
  
  as.character(aln[1:1]) #show first read
  length(aln) #total read number
  tb.bases <- setNames(data.frame(matrix(ncol = 2, nrow = length(aln))),c("read_names","bases")) #create empty df
  
  for (i in 1:length(aln)){
    tb.bases$read_names[i] <- (names(aln[i]))
    sep.bases <- substring(as.character(aln[i]),q.sites,q.sites)    # substring(x, start, end)
    sep.bases <- gsub("\\.", "-", sep.bases)                        # replace "." to "-" ("." for N in CIGAR tag and "-" for D)
    sep.bases <- gsub("\\+", "-", sep.bases)                        # replace "+" to "-" ("." for N in CIGAR tag and "-" for D)
    tb.bases$bases[i] <- paste(sep.bases,collapse="")               # paste multiple site into one string
  }
  write.table(tb.bases, file=paste0(WDdir,"/ind/",sample.name,"/",sample.name,".",paste(q.sites,collapse="-"),".nt.table.csv"), row.names = F)
  (tb.bases.count <- tb.bases %>% dplyr::count(bases))
  write.table(tb.bases.count, file=paste0(WDdir,"/ind/",sample.name,"/",sample.name,".",paste(q.sites,collapse="-"),".nt.count.table.csv"), row.names = F)
}


########### main
ori.sites <- c(1762, 1764) # 1762, 1764 co-variants for example
q.sites <- list()

# position liftover
for (i in ori.sites){ x <- .PosToSH(i)
q.sites <- append(q.sites,x)
(q.sites <- unlist(q.sites))
}

#parameters for GenomicAlignments
which <- IRangesList(JX661488.1=IRanges(1, 3215)) #,
what <- c("qname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)

sample.list <- scan(file = "../example/Sample_name_list.txt", what = "")

nu.count = 1
for (i in 1:length(sample.list)){
  print(paste0("processing: ",sample.list[i],"    ",nu.count,"/",length(sample.list)))
  .NtFromBam(sample.list[i])
  nu.count = nu.count + 1
}

# merge results

nt <- c("A","C","G","T","-")
nt2 <- nt
#nt3 <- nt
#nt4 <- nt
x <- setNames(data.frame(matrix(ncol = 1, nrow = 0)),"bases")

for (i in 1:5){
  for (j in 1:5){
    #    for (k in 1:5){
    #      for (l in 1:5){
    #y <- data.frame(bases=(paste0(nt[i],nt2[j],nt3[k],nt4[l])))
    y <- data.frame(bases=(paste0(nt[i],nt2[j])))
    x<- rbind(x,y)
  }
}
#  }
#}


for (i in 1:length(sample.list)){
  sample.name <- sample.list[i]
  tb.bases.count <- read.table(file=paste0(WDdir,"/ind/",sample.name,"/",sample.name,".",paste(q.sites,collapse="-"),".nt.count.table.csv"), sep = "",header = T)
  idx <- match(x$bases,tb.bases.count$bases)
  x$n <- tb.bases.count$n[idx]
  colnames(x)[i+1] <- sample.name
}
rownames(x) <- x$bases
x <- x[,2:(length(sample.list)+1)]
x <- as.data.frame(t(x))

y <- x[, !grepl( "-" , names( x ) ) ]
z <- x[, grepl( "-" , names( x ) ) ]
x <- cbind(y,z)
write.csv(x, file = paste0("./",paste(q.sites,collapse="-"),".1762-1764_summary.csv"),na = "")

write.csv(y, file = paste0("./",paste(q.sites,collapse="-"),".1762-1764_ACGT_summary.csv"),na = "")

write.csv(z, file = paste0("./",paste(q.sites,collapse="-"),".1762-1764_with-_summary.csv"),na = "")
write.xlsx(z, file = paste0("./",paste(q.sites,collapse="-"),".1762-1764_with-_summary.xlsx"), rowNames = TRUE)
AG <- y[, grepl( "..AG" , names( y ) ) ]

sum.x <- cbind(x, total_reads = rowSums(x,na.rm=T))
sum.x[is.na(sum.x)] <- 0
write.xlsx(sum.x, file = paste0("./",paste(q.sites,collapse="-"),".1762-1764_summary.xlsx"), rowNames = TRUE) 
sum.y <- cbind(y, total_reads = rowSums(y,na.rm=T))
sum.y[is.na(sum.y)] <- 0
write.xlsx(sum.y, file = paste0("./",paste(q.sites,collapse="-"),".1762-1764_ACGT_summary.xlsx"), rowNames = TRUE)

ratio.x <- data.frame(matrix(nrow = length(sample.list), ncol = 26))
for (i in 1:length(colnames(sum.x))){
  ratio.x[,i] <- sum.x[,i]/sum.x[,26]
}
colnames(ratio.x) <- colnames(sum.x)
rownames(ratio.x) <- rownames(sum.x)
#ratio.x[is.na(ratio.x)] <- 0
write.xlsx(ratio.x, file = paste0("./",paste(q.sites,collapse="-"),".1762-1764_ratio.xlsx"), rowNames = TRUE) 

ratio.y <- data.frame(matrix(nrow = length(sample.list), ncol = 17))
for (i in 1:length(colnames(sum.y))){
  ratio.y[,i] <- sum.y[,i]/sum.y[,17]
}
colnames(ratio.y) <- colnames(sum.y)
#rownames(ratio.y) <- rownames(sum.y)
ratio.y[is.na(ratio.y)] <- 0
write.xlsx(ratio.y, file = paste0("./",paste(q.sites,collapse="-"),".1762-1764_ACGT_ratio.xlsx"), rowNames = TRUE)

