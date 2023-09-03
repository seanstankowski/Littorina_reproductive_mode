###### PCA analysis of LGC2.1 ##########
### This script was modified version designed by Katie Hearn from the University of Shefield.
### It runs a PCA on a genotype file formates as: CHROM, POS, SAMPLES; looking for 
### three clusters that match the pattern on an inversion. 
### James Reeve - University of Gothenburg
### 15/11/2022

#### Preparation ####
### Clear the environment
rm(list = ls())
dev.off()
setwd("/Users/james")
options(stringsAsFactors = FALSE)

### Packages
library(tidyverse)
library(data.table)
library(adegenet)
library(cluster)

### Get linkage map
# Original crab map from Westram et al. 2018
Lmap <- read.table("Documents/Inversion_detection/map_V11.txt", header = TRUE)
Lmap <- Lmap[Lmap$LG == 2, ]

### Define breakpoints of LGC2.1
# Use the consensus positions in Westram et al. 2021
bp <- c(0.3, 14.2)
# Filter Lmap to positions on LGC2.1
Lmap.2 <- Lmap[Lmap$av >= bp[1] & Lmap$av <= bp[2], ]

### Get genotypes
GT <- read.table("/Volumes/PhDseqData/for_Sean/LGC2.1_genotypes.txt")[, -c(3:4)] # exclude alt and ref nucleotide-identities

# Add sample names
# Exteact from VCF file with: `gzcat LGC2.1.vcf.gz | awk '/\#/' | tail -n 1 | awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=""; print $0}' > sample_IDs.txt`
SnailIDs <- read.table("/Volumes/PhDseqData/for_Sean/sample_IDs.txt")
SnailIDs[grep("W_arc_04_La", SnailIDs)] <- "W_arc_04_La" # rename strangly named sample
colnames(GT) <- c("CHROM", "POS", SnailIDs)


#### Data wrangeling ####

### Filters
maf <- ceiling(108*0.1) # minor allele count
miss.filt <- ceiling(108*0.05) # missing individuals

### Filter genotypes
GT.2 <- apply(GT, 1, function(X){
  tmp <- X[-c(1:2)] # remove "CHROM" & "POS"
  tmp2 <- tmp[tmp != "./."] # Exclude missing sites
  
  # Allele counts
  N_ref <- sum(tmp2 == "0/0")
  N_het <- sum(tmp2 == "0/1")
  N_alt <- sum(tmp2 == "1/1")
  
  # Is allele count > min.allele.count
  if((2*N_ref + N_het) > maf & 
     (2*N_alt + N_het) > maf &
     sum(tmp == "./.") > miss.filt) {return(X)}
})
GT.2 <- do.call(rbind.data.frame, GT.2)
colnames(GT.2) <- colnames(GT)

# Convert 'POS' back to a numeric vector
GT.2$POS <- as.numeric(GT.2$POS)

# Remove any contigs with < 10 SNPs
GT.2 <- GT.2 %>% group_by(CHROM) %>% filter(n() >= 10)

### Transpose the genotype table
# Combine 'CHROM' and 'POS' into 'SNP'
tmp <- GT.2 %>% unite(SNP, c("CHROM", "POS"), sep = ":")
# Transpose
GT_trans <- transpose(tmp, keep.names = "SnailID", make.names = "SNP")
rm(tmp)


#### Run PCA ####

### Create genind object
GT_genind <- df2genind(GT_trans, ploidy = 2, sep = "/", NA.char = "_")

### Impute missing genotypes using the mean score
GT_genind <- scaleGen(GT_genind, NA.method = "mean", scale = FALSE)

### Add column names back to genotype table
rownames(GT_genind) <- SnailIDs

### Run PCA
GT_PCA <- dudi.pca(GT_genind, center = FALSE, scale = FALSE, nf = 3, scannf = FALSE)

### Extract PC scores and % variance
PCs <- GT_PCA$li
var_percent <- GT_PCA$eig / sum(GT_PCA$eig) * 100


#### Inspect PCA ####
p1 <- ggplot(PCs)+
  geom_point(aes(Axis1, Axis2))+
  labs(x = paste0("PC1 (", round(var_percent[1], 2), "%)"),
       y = paste0("PC2 (", round(var_percent[2], 2), "%)"))+
  theme_bw()

p2 <- ggplot(PCs)+
  geom_point(aes(Axis1, Axis3))+
  labs(x = paste0("PC1 (", round(var_percent[1], 2), "%)"),
       y = paste0("PC3 (", round(var_percent[3], 2), "%)"))+
  theme_bw()

p3 <- ggplot(PCs)+
  geom_point(aes(Axis2, Axis3))+
  labs(x = paste0("PC2 (", round(var_percent[2], 2), "%)"),
       y = paste0("PC3 (", round(var_percent[3], 2), "%)"))+
  theme_bw()

ggarrange(p1,p2, p3, ncol = 1)

### Add sampling information
sample_info <- read.csv("Documents/Seans_WGS_sample_information_March2020.csv", sep=";", check.names = FALSE)

PCs.2 <- merge(PCs, sample_info, by.x = "row.names", by.y = "Sample_ID")

p1 <- ggplot(PCs.2)+
  geom_point(aes(Axis1, Axis2, colour = Country))+
  labs(x = paste0("PC1 (", round(var_percent[1], 2), "%)"),
       y = paste0("PC2 (", round(var_percent[2], 2), "%)"))+
  theme_bw()

p2 <- ggplot(PCs.2)+
  geom_point(aes(Axis1, Axis3, colour = Country))+
  labs(x = paste0("PC1 (", round(var_percent[1], 2), "%)"),
       y = paste0("PC3 (", round(var_percent[3], 2), "%)"))+
  theme_bw()

p3 <- ggplot(PCs.2)+
  geom_point(aes(Axis2, Axis3, colour = Country))+
  labs(x = paste0("PC2 (", round(var_percent[2], 2), "%)"),
       y = paste0("PC3 (", round(var_percent[3], 2), "%)"))+
  theme_bw()

ggarrange(p1,p2, p3, ncol = 1, common.legend = TRUE)


#### Genotype inversion ####

### Run K means clustering
tmp <- lapply(2:6, function(K){
  # Kmeans clustering
  clust <- kmeans(PCs.2$Axis1, centers = K, iter.max = 1e5, nstart = 100)
  # Calcaulte average pairwise-difference between clusters (Silouette Score)
  sil <- mean(silhouette(clust$cluster, dist(PCs.2$Axis1))[,"sil_width"])
  # Append silhouette score to outpur
  res <- append(clust, sil)
  names(res)[10] <- "Avg_Silhouette_Weight"
  return(res)
})

# Find best K
K <- which.max(sapply(tmp, `[[`, "Avg_Silhouette_Weight")) + 1

### Name clusters
# Best cluster information
bestK <- tmp[[K - 1]]

# Identify middle cluster
min_clust <- which.min(bestK$centers) # Cluster with lowest centre
max_clust <- which.max(bestK$centers) # Cluster with highest centre
bestK$cluster[bestK$cluster != min_clust & bestK$cluster != max_clust] <- "RA"

# Add cluster scores to 'PCs.2'
PCs.2$genotype <- bestK$cluster

### Relabel cluster with most crab snails 'RR'
crab_clust <- PCs.2[PCs.2$Species == "saxatilis" & PCs.2$Ecotype %in% c("Crab", "crab", "crabish"), "genotype"]
uClust <- unique(crab_clust[crab_clust != "RA"])
# Find frequency of crab individuals in each homokaryotype cluster
pClust <- sapply(as.numeric(uClust), function(i){sum(crab_clust == i) / bestK$size[i]})
# Assign cluster with highest crab frequency as 'RR'
RR <- uClust[which.max(pClust)]

### Assign homokaryotype clusters
PCs.2$genotype[PCs.2$genotype == RR] <- "RR"
PCs.2$genotype[!(PCs.2$genotype %in% c("RR", "RA"))] <- "AA"

### Check clusters make sense
ggplot(PCs.2)+
  geom_point(aes(x = Axis1, y = Axis2, colour = genotype))+
  labs(x = paste0("PC1 (", round(var_percent[1], 2), "%)"),
       y = paste0("PC2 (", round(var_percent[2], 2), "%)"))+
  theme_bw()

### Save inversion frequencies
write.csv(PCs.2, file = "/Volumes/PhDseqData/for_Sean/LGC2.1_inversion_frequency.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)



p1 <- ggplot(PCs.2)+
  geom_point(aes(Axis1, Axis2, colour = genotype, pch = Species))+
  scale_shape_manual(values = c(15, 17, 19, 19, 19))+
  labs(x = paste0("PC1 (", round(var_percent[1], 2), "%)"),
       y = paste0("PC2 (", round(var_percent[2], 2), "%)"))+
  theme_bw()

p2 <- ggplot(PCs.2)+
  geom_point(aes(Axis1, Axis3, colour = genotype, pch = Species))+
  scale_shape_manual(values = c(15, 17, 19, 19, 19))+
  labs(x = paste0("PC1 (", round(var_percent[1], 2), "%)"),
       y = paste0("PC3 (", round(var_percent[3], 2), "%)"))+
  theme_bw()

p3 <- ggplot(PCs.2)+
  geom_point(aes(Axis2, Axis3, colour = genotype, pch = Species))+
  scale_shape_manual(values = c(15, 17, 19, 19, 19))+
  labs(x = paste0("PC2 (", round(var_percent[2], 2), "%)"),
       y = paste0("PC3 (", round(var_percent[3], 2), "%)"))+
  theme_bw()

ggarrange(p1,p2, p3, ncol = 1, common.legend = TRUE)
