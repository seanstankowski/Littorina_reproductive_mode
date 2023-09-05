##### Convert vcf file to sites file for ARGWEAVER analysis
##### written by Arka Pal
Sys.time() #last update -"2022-05-09 13:53:17 CEST"

#### USAGE: Rscript vcf2sites.R <sample.txt> <input.gt> <chr> <begin> <end>

##
cmd_args <- commandArgs(trailingOnly = T)
cat('Sample FilePath -', cmd_args[1])
cat('\nGenotype FilePath -', cmd_args[2])
out <- gsub(pattern = '.gt', replacement = '_RefAlt.sites', cmd_args[2])
cat('\nOutput FilePath -', out)
cat('\n')

#cmd_args <- c("sample/sample.txt","optix/VCFs/optix_variants.gt")

#setwd("~/Proj-Def-Haplo/")

## Read samples and change column names
cat('\nReading sample file')
samp <- read.table(cmd_args[1])
colnames(samp) <- "ID"
samp$hap1 <- paste0(samp$ID,"_1")
samp$hap2 <- paste0(samp$ID,"_2")
samp

## Read vcf file and change column names
cat('\nReading genotype file')
gt <- read.table(file = cmd_args[2])
#gt <- gt[,c(1,2,4,5,10:ncol(gt))]
colnames(gt) <- c("chr","pos","REF","ALT",samp$ID)
head(gt)

## Make variant sites table 
#initialize
sites <- data.frame(pos = gt$pos, seq = rep("NA", nrow(gt)))
#write genotypes for each variant site
for (posID in c(1:nrow(gt))){
  sites[posID,3] <- gt[posID,'REF']
  sites[posID,4] <- gt[posID,'ALT']
  cat(sprintf("\nWriting %i out of %i sites", posID, nrow(gt)))
  for (sampID in samp$ID){
    bp1 <- if (substr(gt[posID,sampID],1,1)==0) {gt[posID,'REF']} 
    else if (substr(gt[posID,sampID],1,1)==1){gt[posID,'ALT']}
    else {"N"}
    sites[posID,2] <- ifelse(sites[posID,2]=="NA",bp1,paste0(sites[posID,2],bp1))
    
    bp2 <- if (substr(gt[posID,sampID],3,3)==0) {gt[posID,'REF']}
    else if (substr(gt[posID,sampID],3,3)==1){gt[posID,'ALT']}
    else {"N"}
    sites[posID,2] <- ifelse(sites[posID,2]=="NA",bp2,paste0(sites[posID,2],bp2))
  }
}

## Write .sites file
cat("Writing",out,"file")
cat('NAMES', sprintf("%s\t%s", samp$hap1, samp$hap2), file = out, sep ="\t")
cat('\nREGION', cmd_args[3], cmd_args[4], cmd_args[5], file = out, sep="\t", append = T)
cat('\n', file = out, append = T)
write.table(sites, file = out, col.names = F, row.names = F, 
            sep = "\t", append = T, quote = F)
cat("\nDONE\n")
