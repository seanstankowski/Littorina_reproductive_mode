
#1 ----
subset_genomic_interval <- function(x, chromStart, chromEnd){x[which(x$chromEnd >= chromStart & x$chromStart <= chromEnd),]}

#2
plot_snps <- function(snps, clade, colour_name){
  nsites <- nrow(snps)
  for(i in c(1:nsites)){
    majClade <- which(snps[i,] == clade)
    count <- length(majClade)
    points(rep(snps[i,'pos'],count), majClade-2, col=adjustcolor(colour_name, alpha = 0.7), cex = 0.7)
  }
}

#3
select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-Ntip(tree)])
}

#4 ---- #
extract_edges <- function(tree_list, Ne){
  out <- as.list(tree_list)
  for (treeID in c(1:nrow(tree_list))){
    cat(treeID, tree_list[treeID,"chromStart"], tree_list[treeID,"chromEnd"],"\n")
    
    tr <- tree_list[treeID,"tree"]
    tr.egg <- pruneTree(tree = tr, seqs = as.character(eggID)) # 0-indexed
    tr.brood <- pruneTree(tree = tr, seqs = as.character(broodID)) # 0-indexed
    
    apeTree <- read.tree(text=tr)
    
    ## Tree
    nleaf <- length(apeTree$tip.label)
    timeScale <- 1/Ne
    rootage <- getTreeDepth(apeTree)*timeScale
    rootNode <- getRootNode(apeTree)
    edge <- as.data.frame(apeTree$edge)
    names(edge) <- c("anc", "dec")
    edge$anc.label <- sapply(apeTree$edge[,1], select.tip.or.node, tree = apeTree)
    edge$dec.label <- sapply(apeTree$edge[,2], select.tip.or.node, tree = apeTree)
    edge <- rbind(edge, data.frame(anc=-1, dec=rootNode, anc.label=-1, dec.label=rootNode))
    nnode <- nrow(edge)
    maxtime <- 20*timeScale
    edge$length <- c(apeTree$edge.length*timeScale, max(maxtime - rootage, 1))
    edge$depth <- rep(rootage, nnode)
    
    for (i in 1:nnode){
      j <- i
      while (1) {if (edge$dec[j] == rootNode) break
        edge$depth[i] <- edge$depth[i] - edge$length[j]
        j <- which(edge[,"dec"] == edge[j,"anc"])
      }
    }
    
    edge[which(edge$anc == -1),"dec.label"] <- edge[1,"anc.label"]
    
    for (dec.nodeID in c(1:nrow(edge))){
      decs <- Descendants(apeTree, edge$dec[dec.nodeID],"tips")[[1]]
      edge$tips.from.dec[dec.nodeID] <- list(apeTree$tip.label[decs])
      # edge$majClade.node[dec.nodeID] <- all(as.character(majClade_lat-1) %in% apeTree$tip.label[decs])
      # edge$minClade.node[dec.nodeID] <- all(as.character(minClade_lat-1) %in% apeTree$tip.label[decs])
    }
    
    out$edges[[treeID]] <- edge
  }
  return(out)
}

#5
extract_trees_at_SNPs <- function(snp, tree, spr){
  trees_at_SNPs <- data.frame()
  #SPRs_at_SNPs <- data.frame()
  for (rowID in c(1:nrow(snp))){
    cat (rowID, snp$snpID[rowID], snp$pos[rowID],"\n")
    pos <- snp$pos[rowID]
    trees_at_SNPs <- rbind(trees_at_SNPs, tree[which(tree$chromStart <= pos & tree$chromEnd >= pos),])
    #SPRs_at_SNPs <- rbind(SPRs_at_SNPs, spr[which(spr$position == trees_at_SNPs[snpID,"chromEnd"]),])
  }
  #snp <- cbind(snp[1:nrow(snp),], trees_at_SNPs[1:nrow(snp),], SPRs_at_SNPs[1:nrow(snp),])
  snp <- cbind(snp[1:nrow(snp),], trees_at_SNPs[1:nrow(snp),])
  return(snp)
}

#6 edit
subset_clades <- function(snp, alleleCount){snp[which(snp$majAlleleCount_brood==alleleCount),-c(3:102,106)]}

#7
find_clade_depths <- function(x, snp, count, inds){
  out <- data.frame()
  for (i in c(1:nrow(x))){
    snpID <- x[i,"snpID"]
    edge <- snp$tree.attributes[[which(snp$snpID == snpID)]]
    rowID <- sapply(edge$tips.from.dec, function(decs) all(inds %in% decs))
    node <- edge[which(rowID == T),]
    node <- node[which(sapply(node$tips.from.dec, length) == 20-count),]
    
    if(nrow(node) !=0){
      if(node$minClade.node != F){
        node$pos <- tr1_snps$pos[which(tr1_snps$snpID == snpID)]
        node$snpID <- snpID
        out <- rbind(out,node)
      }
    }
  }
  return(out)
}

#8
find_hapblocks <- function(inds, node_depth, singleton, trees){
  block <- data.frame()
  inds <- as.character(inds)
  # if(pop == 'all'){
  #   latTrees = trees
  # } else{
  #   latTrees = optixTreeAttributes_lat_list
  # }
  # for (treeID in c(1:length(latTrees$edges))){
  for (treeID in c(1:length(trees$edges))){
    cat (treeID,"\n")
    tr <- trees$edges[[treeID]]
    rowID <- sapply(tr$tips.from.dec, function(x) all(inds %in% x))
    node <- tr[which(rowID == T),]
    if (singleton == TRUE){
      node <- node[which(sapply(node$tips.from.dec, length) == length(inds)),]
    }else{
      node <- node[which(round(node$depth,3) == node_depth & node$length != 0),]
    }
    if(nrow(node) == 0){node[1,] = rep(NA, ncol(node))}
    node$diff <- any(setdiff(broodID, inds) %in% node$tips.from.dec[[1]])
    node$treeID <- treeID
    node$chromStart <- trees$chromStart[treeID]
    node$chromEnd <- trees$chromEnd[treeID]
    node$tree.len <- node$chromEnd-node$chromStart+1
    node$area <- node$tree.len*node$length
    block <- rbind(block,node)
  }
  return(block)
}

#9
find_disjunct_hapblocks <- function(block){
  block <- block[which(block$anc != "NA"),]
  #block <- block[which(block$diff != TRUE),]
  breaks <- which(diff(block$treeID)!=1)
  disjunct_blocks <- split(block, cumsum(1:nrow(block) %in% (breaks+1)))
  # for (i in seq(length(disjunct_blocks))){disjunct_blocks[[i]] <- subset(disjunct_blocks[[i]], diff = "FALSE")}
  return(disjunct_blocks)
}

#10
plot_hapblock <- function(hapblock, colour_name, node_depth){
  for(i in c(1:length(hapblock))){
    if(node_depth == T){
      depth = 1e-20
    }else{depth <- hapblock[[i]]$depth[1]}
    treespan <- vector()
    
    for (j in 1:nrow(hapblock[[i]])){
      treespan <- c(treespan, hapblock[[i]]$chromStart[j], hapblock[[i]]$chromEnd[j])
    }
    polygon(x = c(treespan, rev(treespan)),
            y = c(rep((hapblock[[i]]$length + depth), each = 2), rep(depth, 2*nrow(hapblock[[i]]))),
            col = adjustcolor(colour_name, alpha=0.3), 
            border = colour_name)
  }
}