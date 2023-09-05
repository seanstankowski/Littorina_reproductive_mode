##' Average color
##' @param cols A list of colors
##' @return A color created by averaging all the R, G, B values of the input colors
##' @export
##' @import grDevices
##' @import graphics
meanColor <- function(cols) {
  x <- rowSums(col2rgb(cols))/length(cols)
  rgb(x[1], x[2], x[3], maxColorValue=255)
}

getname <- function(idx, edge) {
  if (edge$name[idx] != "") return(edge$name[idx])
  node <- edge[idx,"dec"]
  decs <- which(edge[,"anc"]==node)
  if (length(decs)==0) return("")
  min(getname(decs[1], edge), getname(decs[2], edge))
}

getLeafIdx <- function(idx, edge, leafOrder) {
  node <- edge[idx,"dec"]
  decs <- which(edge[,"anc"]==node)
  if (length(decs)==0) return(which(leafOrder == edge[idx,"name"]))
  (getLeafIdx(decs[1], edge, leafOrder) +
      getLeafIdx(decs[2], edge, leafOrder))/2
}

getColorForEdge <- function(i, edge, colors, default="black") {
  leafDecs <- getAllDescendants(i, edge)
  col <- character(length(leafDecs))
  for (i in 1:length(leafDecs)) {
    if (is.null(colors[[leafDecs[i]]])) {
      col[i] <- default
    } else col[i] <- colors[[leafDecs[i]]]
  }
  meanColor(col)
}

getAllDescendants <- function(i, edge) {
  node <- edge[i,"dec"]
  decs <- which(edge[,"anc"]==node)
  if (length(decs)==0) return(edge[i,"name"])
  if (length(decs) != 2) stop("Error")
  c(getAllDescendants(decs[1], edge),
    getAllDescendants(decs[2], edge))
}

getDiscreteTime <- function(time, mod, tol=1.0e-5) {
  if (is.null(mod)) return(time)
  tidx <- which(abs(mod$times-time)/max(mod$times) < tol)
  if (length(tidx) != 1) {
    cat("time=", time, "\n")
    cat("times=", mod$times, "\n")
    stop("Error discretizing time")
  }
  tidx
}

getPathFromStr <- function(str, matchIdx, isSpr=FALSE) {
  chars <- strsplit(str, "")[[1]]
  if (!isSpr) {
    if (matchIdx >= 5 && substr(str, matchIdx-4, matchIdx-1)=="spr_")
      return(0)
    pos <- matchIdx+9
    startPos <- pos
  } else {
    pos <- matchIdx+13
    startPos <- pos
  }
  while(chars[pos] >= "0" &&
        chars[pos] <= "9") {
    pos <- pos+1
  }
  pos <- pos-1
  as.numeric(substr(str, startPos, pos))
}

getTimeFromStr <- function(str, matchIdx) {
  chars <- strsplit(str, "")[[1]]
  pos <- matchIdx
  while (chars[pos] != "=") {pos=pos+1}
  pos=pos+1
  startPos=pos
  while ((chars[pos] >= "0" && chars[pos] <= "9") || chars[pos]==".")
    pos=pos+1
  as.numeric(substr(str, startPos, pos-1))
}


## only returns branches with non-zero pop paths
## return value is data frame with columns: path, leaf, numAnc
## numAnc=0 implies leaf branch, otherwise count up trea from leaf by numAnc
getPopPaths <- function(tree) {
  matches <- as.numeric(gregexpr("pop_path=", tree)[[1]])
  if (matches[1] == -1) return(NULL)
  rv <- NULL
  chars <- strsplit(tree, "")[[1]]
  for (i in 1:length(matches)) {
    path <- getPathFromStr(tree, matches[i])
    if (path==0) next
    nodeCount <- 0
    pos <- matches[i]-1
    while (chars[pos] != '[') pos <- pos-1  # get out of NHX section
    pos <- pos-1
    while(chars[pos] != ":" || chars[pos-1] == ")") {
      if (chars[pos] == ']') {
        while (chars[pos] != '[') {
          pos <- pos-1
        }
      }                    
      if (chars[pos] == ")") {
        nodeCount <- nodeCount+1
      }
      else if (chars[pos] == ")") {
        if (nodeCount==0) stop("ERROR")
        nodeCount <- nodeCount-1
      }
      pos <- pos-1
    }
    endpos <- pos-1
    while (chars[pos] != "(" && chars[pos] != ",") pos <- pos-1
    leaf <- substr(tree, pos+1, endpos)
    rv <- rbind(rv, data.frame(path=path, leaf=leaf, numAnc=nodeCount))
  }
  rv
}

removeNHX <- function(tree) {
  gsub("\\[&&NHX[^\\]]+\\]", "", tree, perl=TRUE)
}

## call rphast prune.tree but first get rid of NHX tag
pruneTree <- function(tree, seqs, all.but=TRUE) {
  if (!requireNamespace("rphast", quietly=TRUE)) {
    stop("Need rphast package installed")
  }
  tree <- gsub("\\[&&NHX[^\\]]+\\]", "", tree, perl=TRUE)
  rphast::prune.tree(removeNHX(tree), seqs, all.but=all.but)
}


##' Node postorder
##' @inheritParams postorder.arg
##' @param idx The current idx
##' @note This is a recursive helper function for postorder.arg
##' @seealso \code{\link{postorder.arg}}
postorderRec <- function(edge, mod, idx, leafOrder=NULL) {
  node <- edge[idx,"dec"]
  w <- which(edge[,"anc"]==node)
  if (length(w) == 0) return(idx)
  if (length(w) != 2) stop("ERRROR")
  if (is.null(mod)) {
    pop <- c(0,0)
  } else {
    pop <- c(meanPop(w[1], edge, mod),
             meanPop(w[2], edge, mod))
  }
  if (pop[2] < pop[1]) {
    w <- rev(w)
  } else if (!is.null(leafOrder)) {
    if (getLeafIdx(w[2], edge, leafOrder) <
        getLeafIdx(w[1], edge, leafOrder))
      w <- rev(w)
  } else if (getname(w[2], edge) <
             getname(w[1], edge)) {
    w <- rev(w)
  }
  c(postorderRec(edge, mod, w[1], leafOrder),
    postorderRec(edge, mod, w[2], leafOrder),
    idx)
}

plotSpr <- function(tree, edgeData, mod=mod, col="black", timeScale=1, lwd=1) {
  if (is.null(mod)) stop("plotSpr currently only works if mod is provided")
  cat(tree, "\n")
  arrowLen <- par("pin")[1]*0.01
  if (timeScale != 1) mod <- scalePopModel(mod, timeScale)
  rtimePos <- as.numeric(gregexpr("recomb_time=", tree)[[1]])
  if (length(rtimePos) > 1) {
    stop("too many recomb times in a tree")
  }
  if (rtimePos == -1) return(NULL)
  ctimePos <- as.numeric(gregexpr("coal_time=", tree)[[1]])
  if (length(ctimePos) > 1 || ctimePos == -1) stop("error")
  sprPath <- 0
  sprPos <- as.numeric(gregexpr("spr_pop_path=", tree)[[1]])
  if (length(sprPos) > 1) stop("Error")
  if (sprPos != -1)
    sprPath <- getPathFromStr(tree, sprPos, isSpr=TRUE)
  matches <- c(rtimePos, ctimePos)
  chars <- strsplit(tree, "")[[1]]
  times <- numeric(2)
  
  rv <- NULL
  for (i in 1:length(matches)) {
    nodeCount <- 0
    pos <- matches[i]-1
    while (chars[pos] != '[') pos <- pos-1  # get out of NHX section
    pos <- pos-1
    while(chars[pos] != ":" || chars[pos-1] == ")") {
      if (chars[pos] == ']') {
        while (chars[pos] != '[') {
          pos <- pos-1
        }
      }                    
      if (chars[pos] == ")") {
        nodeCount <- nodeCount+1
      }
      else if (chars[pos] == ")") {
        if (nodeCount==0) stop("ERROR")
        nodeCount <- nodeCount-1
      }
      pos <- pos-1
    }
    endpos <- pos-1
    while (chars[pos] != "(" && chars[pos] != ",") pos <- pos-1
    leaf <- substr(tree, pos+1, endpos)
    rv <- rbind(rv, data.frame(leaf=leaf, numAnc=nodeCount, time=getTimeFromStr(tree, matches[i])))
  }
  
  ## rv[1,] is recomb
  rv[,"time"] <- rv[,"time"]*timeScale
  recombTime <- getDiscreteTime(rv[1,"time"], mod)
  coalTime <- getDiscreteTime(rv[2,"time"], mod)
  if (recombTime > coalTime) stop("Error")
  edgeIdx <- numeric(2)
  for (j in 1:2) {
    edgeIdx[j] <- which(edgeData$name == rv[j,"leaf"])
    if (rv[j,"numAnc"] > 0) {
      for (i in 1:rv[j,"numAnc"])
        edgeIdx[j] <- which(edgeData$dec == edgeData[edgeIdx[j],"anc"])
    }
  }
  recombIdx <- edgeIdx[1]
  coalIdx <- edgeIdx[2]
  startX <- edgeData[recombIdx,sprintf("xpos%i", recombTime)]
  endX <- edgeData[coalIdx,sprintf("xpos%i", coalTime)]
  points(x=startX, y=mod$times[recombTime], pch="X", col=col)
  if (recombTime == coalTime) {
    arrows(x0=startX, x1=endX, y0=mod$times[recombTime], length=arrowLen, col=col, lwd=lwd)
    return(invisible(NULL))
  }
  
  currPop <- mod$path[sprPath+1,recombTime]+1  # one-based pop
  allXpos <- sort(unique(c(as.numeric(unlist(edgeData[,sprintf("xpos%i", 1:mod$ntime)])),
                           as.numeric(unlist(mod$popCoords[,c("x0", "x1")])))))
  ## want to put currX between recomb branch and nearest branch OR pop boundary, in direction of coal node
  w <- which(allXpos == startX)
  if (endX > startX) {
    currX <- mean(allXpos[c(w,w+1)])
  } else {
    currX <- mean(allXpos[c(w, w-1)])
  }
  segments(x0=startX, x1=currX, y0=mod$time[recombTime], col=col, lwd=lwd)
  prevPop <- currPop
  prevT <- recombTime
  for (currT in (recombTime+1):coalTime) {
    currPop <- mod$path[sprPath+1, currT]+1
    if (currPop == prevPop) {
      segments(x0=currX, y0=mod$times[prevT], y1=mod$times[currT], col=col, lwd=lwd)
    } else {
      w1 <- which(mod$div$time > mod$times[prevT] & mod$div$time < mod$times[currT] &
                    mod$div$pop1+1 == prevPop & mod$div$pop2+1 == currPop)
      if (!is.null(mod$mig)) {
        w2 <- which(mod$mig$time > mod$times[prevT] & mod$mig$time < mod$times[currT] &
                      mod$mig$pop1+1 == prevPop & mod$mig$pop2+1 == currPop)
      } else w2 <- numeric(0)
      if (length(w1)+length(w2) != 1) {
        browser()
        stop("Error2 deciding mig/div")
      }
      if (length(w1) == 1) {  #divergence
        switchTime <- mod$div$time[w1]
      } else switchTime <- mod$mig$time[w2]
      segments(x0=currX, y0=mod$times[prevT], y1=switchTime, col=col, lwd=lwd)
      if (currX < mod$popCoords[currPop,"x0"]) {
        w <- which(allXpos == mod$popCoords[currPop,"x0"])
        newX <- mean(allXpos[c(w, w+1)])
      } else {
        if (! (currX > mod$popCoords[currPop, "x1"])) stop("Error")
        w <- which(allXpos == mod$popCoords[currPop, "x1"])
        newX <- mean(allXpos[c(w-1,w)])
      }
      segments(x0=currX, x1=newX, y0=switchTime, col=col, lwd=lwd)
      segments(x0=newX, y0=switchTime, y1=mod$times[currT], col=col, lwd=lwd)
      prevPop <- currPop
      currX <- newX
    }
    prevT <- currT
  }
  arrows(x0=currX, x1=endX, y0=mod$times[coalTime], col=col, lwd=lwd, length=arrowLen)
  invisible(NULL)
}

##' Node postorder
##' @param edge The edge matrix of a tree read by ape
##' @param mod A population model for multi-population use
##' @param leafOrder If given, this character vector determines how to sort child nodes
postorder.arg <- function(edge, leafOrder=NULL, mod=NULL) {
  root <- which(edge$anc==-1)
  postorderRec(edge, mod, root, leafOrder)
}

##' Get the index of the root node for a tree
##' @param tree An object of type "phylo" (a tree loaded in ape).
##' @return The index of the root node in the tree object
getRootNode <- function(tree) {
  if (class(tree) != "phylo") {
    stop("This function is for trees read by ape")
  }
  node <- 1
  while (1) {
    w <- which(tree$edge[,2]==node)
    if (length(w)==0) break
    if (length(w) > 1) stop("error")
    node <- tree$edge[w,1]
  }
  node
}

##' Get the population path of the parent branch of tree
##' @param tree A character vector with a newick tree
##' @return The population path encoded in the newick string for outermost branch
##' @note The population path is usually zero; it is only non-zero for the multiple-population
##' ARGweaver model
getOuterPopPath <- function(tree) {
  ## search in to last closing parentheses
  parens <- gregexpr(')', tree)[[1]]
  lastParen <- parens[length(parens)]
  lastTag <- substr(tree, lastParen+1, nchar(tree))
  matches <- as.numeric(gregexpr("pop_path=", lastTag)[[1]])
  if (matches[1] == -1) return(0)
  getPathFromStr(lastTag, matches)
}

##' Get maximum distance from a node to the leaf
##' @param tree A tree, either of type "phylo" (read in by ape), or a character
##' string with a Newick tree
##' @param idx The index of the root node. Usually this will be -1 and the true root will
##' be used. The value of idx refers to the row of the tree$edge data frame,
##' for trees read by "ape"
##' @return The maximum distance from the indexed node (or root if idx=-1) to any leaf node. 
getTreeDepth <- function(tree, idx=-1) {
  if (class(tree) != "phylo") {
    if (!requireNamespace("ape", quietly=TRUE)) {
      stop("ape package is required for getTreeDepth")
    }
    tree <- ape::read.tree(text=tree)
  }
  if (idx==-1) {
    node <- unique(tree$edge[,1][which(!is.element(tree$edge[,1], tree$edge[,2]))])
  } else {
    node <- tree$edge[idx,2]
  }
  anc <- which(tree$edge[,1]==node)
  if (length(anc)==2) {
    return(max(tree$edge.length[anc[1]] + getTreeDepth(tree, anc[1]),
               tree$edge.length[anc[2]] + getTreeDepth(tree, anc[2])))
  }
  0
}

##' Draw a tree from a newick string; add to current plot ----
##' @param tree A newick string representing a tree
##' @param col Either a character string (single value) giving color of tree. Or, a
##' list giving color to plot each leaf branch. Internal branches will
##' be plotted as "average" color among child leafs.
##' @param lwd Line width
##' @param timeScale multiply all branches by this value
##' @param call.plotSpr If TRUE and if a SPR event is encoded in the tree (with NHX tags),
##' then draw the SPR event on the tree (currently only works if mod is not NULL)
##' @param cex.leafname Character size of leaf names
##' @param leafCol If given, this can be a list assigning a color to each leaf name
##' (individual names also work if haploid names are in the form indName_1 and indName_2)
##' @param leafLabels A list indexed by the leaf names in the tree, giving the label to display
##' for each leaf name. If NULL, the leaf names are displayed.
##' @param mod (Advanced; for use with multi-population version of ARGweaver) A model
##' file read in with the function readPopModel, will draw population model underneath tree
##' @note The parameters col, leafCol, and leafLabels should all be lists indexed by leaf
##' names of the tree. However, they can also be diploid individual names, if the leafs
##' are named with the convention indName_1 and indName_2
##' @seealso \code{\link{plotTrees}} for function which creates a new plot
##' @export
drawTree <- function(tree, col="black", lwd=1, timeScale=1, 
                     call.plotSpr=FALSE, cex.leafname=2, leafCol=col, leafLabels=NULL,
                     mod=NULL) {
  if (!requireNamespace("ape", quietly=TRUE)) {
    stop("drawTree requires ape package installed")
  }
  orig.mod <- mod
  haveMod <- !is.null(mod)
  if (haveMod && timeScale != 1) mod <- scalePopModel(mod, timeScale)
  if (haveMod) maxtime <- mod$times[mod$ntimes] else maxtime <- par("usr")[4]
  apeTree <- ape::read.tree(text=tree)
  nleaf <- length(apeTree$tip.label)
  rootage <- getTreeDepth(apeTree)*timeScale
  rootPath <- getOuterPopPath(tree)
  rootNode <- getRootNode(apeTree)
  logScale <- par("ylog")
  if (logScale && haveMod) mod$times[mod$times==0] <- 1
  if (!is.null(leafCol)) {
    for (l in names(leafCol))
      for (hap in 1:2)
        leafCol[[sprintf("%s_%i", l, hap)]] <- leafCol[[l]]
  }
  if (is.list(col)) {
    for (l in names(col))
      for (hap in 1:2)
        col[[sprintf("%s_%i", l, hap)]] <- col[[l]]
  }
  if (!is.null(leafLabels)) {
    for (l in names(leafLabels))
      for (hap in 1:2)
        leafLabels[[sprintf("%s_%i", l, hap)]] <- leafLabels[[l]]
  }
  edge <- as.data.frame(apeTree$edge)
  names(edge) <- c("anc", "dec")
  edge <- rbind(edge, data.frame(anc=-1, dec=rootNode))
  nnode <- nrow(edge)
  edge$length <- c(apeTree$edge.length*timeScale, max(maxtime - rootage, 1))
  edge$depth <- rep(rootage, nnode)
  for (i in 1:nnode) {
    j <- i
    while (1) {
      if (edge$dec[j] == rootNode) break
      edge$depth[i] <- edge$depth[i] - edge$length[j]
      j <- which(edge[,"dec"] == edge[j,"anc"])
    }
  }
  
  edge$name <- rep("", nrow(edge))
  for (i in 1:length(apeTree$tip.label)){
    edge[edge$dec==i,"name"] <- apeTree$tip.label[i]}
  
  edge$path <- rep(0, nnode)
  if (haveMod) {
    paths <- getPopPaths(tree)
    if (nrow(paths) > 0) {
      for (i in 1:nrow(paths)) {
        w <- which(edge$name == paths[i,"leaf"])
        if (length(w) != 1) stop("Error")
        if (paths[i,]$numAnc > 0) {
          for (j in 1:paths[i,]$numAnc) {
            w <- which(edge[,"dec"] == edge[w,"anc"])
            if (length(w) != 1) stop("ERROR")
          }
        }
        edge$path[w] <- paths[i,"path"]
      }
    }
  }
  edge[edge$depth <= 0,"depth"] <- 1e-15
  
  o <- postorder.arg(edge, mod=mod)
  if (haveMod) {
    numBranchPerPop <- rep(0, mod$npop)
    counted <- matrix(FALSE, nrow=nrow(edge), ncol=mod$npop)
    for (i in o) {
      startTime <- edge[i,"depth"]
      endTime <- startTime + edge[i, "length"]
      startT <- getDiscreteTime(startTime, mod)
      if (length(startT) != 1) stop("Error")
      endT <- getDiscreteTime(endTime, mod)
      path <- edge[i,"path"]+1  ## one-based path
      prevPop <- -1
      for (currT in startT:endT) {
        pop <- mod$paths[path,currT]+1
        if (prevPop > 0 && pop != prevPop)
          counted[i,prevPop] <- FALSE
        if (!counted[i,pop]) {
          numBranchPerPop[pop] <- numBranchPerPop[pop]+1
          counted[i,pop] <- TRUE
        }
        prevPop <- pop
      }
      anc <- which(edge[,"dec"]==edge[i,"anc"])
      if (length(anc) == 1) {
        counted[anc,pop] <- TRUE
      }
    }
  } else {
    numBranchPerPop <- (nnode+1)/2
    popStep <- (par("usr")[2]-par("usr")[1])/numBranchPerPop
    popX <- par("usr")[1]+popStep/2
  }
  
  edge$xpos <- -1
  if (haveMod) {
    tnames <- sprintf("xpos%i", 1:mod$ntime)
    for (i in 1:mod$ntime) {
      edge[,tnames[i]] <- -1
    }
    popStep <- ( mod$popCoords$x1 - mod$popCoords$x0 ) / (numBranchPerPop)
    popX <- mod$popCoords$x0 + popStep/2
    
    mod$div$total <- 0
    mod$div$count <- 0
    for (i in 1:nrow(mod$div)) {
      t1 <- which(mod$times == max(mod$times[mod$times < mod$div$time[i]]))
      t2 <- which(mod$times == min(mod$times[mod$times > mod$div$time[i]]))
      paths <- which(mod$paths[,t1]==mod$div$pop1[i] &
                       mod$paths[,t2]==mod$div$pop2[i])
      mod$div$total[i] <- sum(edge$depth < mod$div$time[i] &
                                edge$depth + edge$length > mod$div$time[i] &
                                is.element(edge$path+1, paths))
    }
    if (!is.null(mod$mig)) {
      mod$mig$total <- 0
      mod$mig$count <- 0
      for (i in 1:nrow(mod$mig)) {
        t1 <- which(mod$times == max(mod$times[mod$times < mod$mig$time[i]]))
        t2 <- which(mod$times == min(mod$times[mod$times > mod$mig$time[i]]))
        paths <- which(mod$paths[,t1]==mod$mig$pop1[i] &
                         mod$paths[,t2]==mod$mig$pop2[i])
        mod$mig$total[i] <- sum(edge$depth < mod$mig$time[i] &
                                  edge$depth + edge$length > mod$mig$time[i] &
                                  is.element(edge$path+1, paths))
      }
    }
  }
  for (i in o) {
    startTime <- edge[i,"depth"]
    endTime <- startTime + edge[i,"length"]
    if (haveMod) {
      startT <- getDiscreteTime(startTime, mod)
      if (length(startT) != 1) stop("Error")
      endT <- getDiscreteTime(endTime, mod)
      path <- edge[i,"path"]+1  ## one-based path
      startPop <- mod$paths[path,startT]+1  ## one-based population
    } else {
      startPop <- 1
      startT <- 1
      endT <- 2
    }
    node = edge$dec[i]
    anc <- which(edge$anc == node)
    if (!is.list(col)) {
      useCol <- col
    } else {
      useCol <- getColorForEdge(i, edge, col)
    }
    if (length(anc) > 0) {
      if (length(anc) != 2) stop("Error")
      edge$xpos[i] <- mean(edge$xpos[anc])
      ##            edge$xpos.orig[i] <- edge$xpos[i]
      ##            if (edge$xpos[anc[1]] < 0 || edge$xpos[anc[2]] < 0) browser()
      segments(y0=startTime, x0=edge$xpos[anc[1]], x1=edge$xpos[anc[2]], col=useCol, lwd=lwd)
    } else {
      displayName <- edge[i,"name"]
      tipColor <- ifelse(is.null(leafCol) || !is.element(displayName, names(leafCol)), "black",
                         leafCol[[displayName]])
      if (!is.null(leafLabels)) {
        if (is.element(displayName, names(leafLabels)))
          displayName <- leafLabels[[displayName]]
      }
      node <- edge[i,"dec"]
      edge$xpos[i] <- popX[startPop]
      ##            edge$xpos.orig[i] <- edge$xpos[i]
      popX[startPop] <- popX[startPop] + popStep[startPop]
      mtext(text=displayName, side=1, line=0.25, at=edge$xpos[i], las=2, cex=cex.leafname,
            col=tipColor)
    }
    
    prevPop <- startPop
    prevT <- startT
    if (startT==endT) next
    currPop <- startPop
    for (currT in (startT+1):endT) {
      if (haveMod) {
        currPop <- mod$path[path,currT]+1
      }
      if (prevPop == currPop) {
        if (haveMod) {
          segments(y0=mod$times[prevT], y1=mod$times[currT], x0=edge$xpos[i], col=useCol, lwd=lwd)
          for (j in prevT:currT) edge[i,tnames[j]] <- edge$xpos[i]
        } else {
          #                    browser()
          segments(y0=startTime, y1=endTime, x0=edge$xpos[i], col=useCol, lwd=lwd)
        }
      } else {
        ## need to get time of switch
        w1 <- which(mod$div$time > mod$times[prevT] & mod$div$time < mod$times[currT] &
                      mod$div$pop1+1==prevPop & mod$div$pop2+1==currPop)
        if (!is.null(mod$mig)) {
          w2 <- which(mod$mig$time > mod$times[prevT] & mod$mig$time < mod$times[currT] &
                        mod$mig$pop1+1==prevPop & mod$mig$pop2+1 == currPop)
        } else w2 <- numeric(0)
        if (length(w1)+length(w2) != 1) {
          browser()
          stop("Error deciding mig/div")
        }
        if (length(w1) == 1) {
          if (mod$div$total[w1] == 1) {
            switchTime <- (mod$div$y0[w1]+mod$div$y1[w1])/2
          } else {
            switchTime <- mod$div$y0[w1] + (mod$div$y1[w1]-mod$div$y0[w1])*mod$div$count[w1]/(mod$div$total[w1]-1)
          }
          mod$div$count[w1] <- mod$div$count[w1]+1
        } else {
          switchTime <- mod$mig$time[w2]
        }
        segments(y0=mod$times[prevT], y1=switchTime, x0=edge$xpos[i], col=useCol, lwd=lwd)
        if (haveMod) for (j in prevT:(currT-1)) edge[i,tnames[j]] <- edge$xpos[i]
        oldpos <- edge$xpos[i]
        edge$xpos[i] <- popX[currPop]
        popX[currPop] <- popX[currPop]+popStep[currPop]
        if (oldpos < 0 || edge$xpos[i] < 0) browser()
        segments(y0=switchTime, x0=oldpos, x1=edge$xpos[i], col=useCol, lwd=lwd)
        segments(y0=switchTime, y1=mod$times[currT], x0=edge$xpos[i], col=useCol, lwd=lwd)
        if (haveMod) edge[i,tnames[currT]] <- edge$xpos[i]
      }
      prevT <- currT
      prevPop <- currPop
    }
  }
  rv <- edge[o,]
  if (call.plotSpr) plotSpr(tree, rv, orig.mod, timeScale=timeScale)
  invisible(rv)
}

if (FALSE) { ## Code in progress
  resolveAllele <- function(tree, alleles, nodeAlleles) {
    
  }
  
  getNodeAlleles <- function(site, tree, node=NULL, nodeAlleles=list()) {
    if (is.null(node)) {
      node <- which(tree$anc==-1)
    }
    decs <- which(tree$anc==tree[node,"dec"])
    if (length(decs) == 0) {
      nodeAlleles[[node]] <- site[node, "name"]
      if (nodeAlleles[[node]] == "N")
        nodeAlleles[[node]] <- c("A", "C", "G", "T")
    } else {
      if (length(decs) != 2) stop("error")
      nodeAlleles <- getNodeAlleles(site, tree, decs[1], nodeAlleles)
      nodeAlleles <- getNodeAlleles(site, tree, decs[2], nodeAlleles)
      a1 <- nodeAlleles[[decs[1]]]
      a2 <- nodeAlleles[[decs[2]]]
      if (a1 == "N") {
        nodeAlleles[[node]] <- a2
      } else if (a2 == "N") {
        nodeAlleles[[node]] <- a1
      } else {
        possible <- intersect(a1, a2)
        if (length(possible) >= 1) {
          nodeAlleles[[node]] <- possible
          if (length(a1) > 1) {
            resolveAllele(tree, possible, decs[1], nodeAlleles)
          }
          if (length(a2) > 1) {
            resolveAllele(tree, possible, desc[2], nodeAlleles)
          }
        } else {
          nodeAlleles[[node]] <- union(a1, a2)
        }
      }
    }
  }
  
  
  getNodeAlleles <- function(site, tree, node=NULL, alleles=list()) {
    tree$node <- as.integer(rownames(tre))
    if (is.null(node)) {
      node <- tree[which(tree$anc==-1),"node"]
    }
    decs <- tree[which(tree$anc==node),"node"]
    if (length(decs) == 0) {
      alleles[[node]] <- site[node, "name"]
    } else {
      alleles[[node]] <- unique(c(alleles[[decs]]))
      for (i in 1:length(decs)) {
        alleles <- getNodeAlleles(site, tree, node=decs[i], alleles)
      }
    }}
} ## IF FALSE

drawSiteOnTree <- function(site, tree)  {
  ## first check if site is variant
  tmp <- table(as.character(site[,2:ncol(site)]))
  w <- which(names(tmp)=="N")
  if (length(w) == 1) tmp <- tmp[-w]
  if (length(tmp) == 1) return()  # no variation here
  
  
  
  
}


##' Plot a tree ----
##' @inheritParams drawTree
##' @param tree A newick string containing a tree
##' @param prune A list of leafs to prune from the tree before plotting (Null means prune none)
##' @param keepSeqs A list of leafs to keep in the tree (NULL means keep all) 
##' @param drawSpr If TRUE, draw the SPR event that turns each tree into the next one
##' @param ylab label for y axis
##' @param logScale If TRUE, plot in log scale
##' @param ylim Range for y axis (default; use range of tree)
##' @param mar Margins for plot
##' @param add If TRUE, do not create a new plot
##' @param sites If given, draw mutation events on tree
##' @param chromStart The start coordinate of the tree region (1-based; only used if SITES not null)
##' @param chromEnd The end coordinate of the tree region
##' @param ... Passed to plot function
##' @note This creates a new plot for each tree. If plotting to the screen, probably
##' want to call par(ask=TRUE) first.
##' @export
plotTree <- function(tree, prune=NULL, keepSeqs=NULL, 
                     col="black", leafCol=col, leafLabels=NULL, timeScale=1,
                     drawSpr=FALSE, sites=NULL, chromStart=NULL, chromEnd=NULL,
                     ylab="Generations", logScale=FALSE, ylim=NULL, add=FALSE,
                     mar=c(8,4,1,1), mod=NULL, cex.leafname=0.8, ...) {
  if (!is.null(mar)) {
    if (length(mar) != 4)
      stop("mar should be numeric vector of length 4")
    par(mar=mar)
  }
  if (is.null(ylim)) ylim <- c(0, getTreeDepth(tree)*timeScale)
  if (!add) {
    if (!logScale) {
      plot(c(1), c(1), ylim=ylim, xaxt="n", xlab="", yaxs="i", bty="n",type="n", ylab=ylab, ...)
    } else {
      if (ylim[1] == 0) ylim[1] <- 1
      plot(c(1), c(1), ylim=ylim, xaxt="n", xlab="", yaxs="i", bty="n",type="n", log="y", ylab=ylab, ...)
    }
  }
  if (!is.null(prune)) tree <- pruneTree(tree, seqs=prune, all.but=FALSE)
  if (!is.null(keepSeqs)) tree <- pruneTree(tree, seqs=keepSeqs, all.but=TRUE)
  rv <- drawTree(tree, call.plotSpr=drawSpr, col=col, leafCol=leafCol, cex.leafname=cex.leafname,
                 leafLabels=leafLabels, timeScale=timeScale,  mod=mod)
  
  if (!is.null(sites)) {
    sites <- sites[sites$pos >= chromStart & sites$pos <=chromEnd,]
    if (nrow(sites) > 0) {
      for (i in 1:nrow(sites)) {
        drawSiteOnTree(sites[i,], rv)
      }
    }
  }
  invisible(rv)
}


##' Plot trees from newick string ----
##' @inheritParams plotTree
##' @param trees A character vector, each element is a newick tree to plot
##' @param treeInfo If given, should be a character vector of same length as trees
##' @param regionSide Print treeInfo to the margins on this side
##' (1=bottom, 2=left, 3=top, 4=right, anything else = don't print)
##' @param regionLine Print region of eaach tree on this line of the margin
## #' @param ... Passed to plotTree function
##' @param popwidth If mod is not null, popwidth can be a numeric vector
##' of length equal to the number of populations, giving relative width
##' of each
##' @param xlab Label for x axis
##' @note This creates a new plot for each tree. If plotting to the screen, probably
##' want to call par(ask=TRUE) first.
##' @export
plotTrees <- function(trees, prune=NULL, keepSeqs=NULL, treeInfo=NULL,
                      col="black", leafCol=col, leafLabels=NULL, timeScale=1,
                      drawSpr=FALSE, ylab="Generations", xlab="", logScale=FALSE, ylim=NULL, add=FALSE,
                      mar=c(8,4,1,1),
                      #                      sites=NULL,
                      regionSide=1, regionLine=4, mod=NULL, popwidth=NULL, ...) {
  if (!is.null(treeInfo)) {
    if (length(treeInfo) != length(trees))
      warning("treeInfo should be same length as trees ", length(treeInfo), ", ", length(trees))
  }
  rv <- list()
  if (add == FALSE && !is.null(mod)) add <- TRUE
  mod2 <- NULL
  for (i in 1:length(trees)) {
    if (!is.null(mod)) mod2 <- drawPopModel(mod, popwidth=popwidth, timescale=timeScale, ylim=ylim, ylab=ylab, xlab=xlab)
    rv[[i]] <- plotTree(trees[i], prune=prune, keepSeqs=keepSeqs, col=col, leafCol=leafCol,
                        leafLabels=leafLabels, timeScale=timeScale, drawSpr=drawSpr, ylab=ylab,
                        logScale=logScale, ylim=ylim, add=add, mar=mar, mod=mod2, ...)
    if (is.element(regionSide, 1:4))
      mtext(treeInfo[i], side=regionSide, line=regionLine, ...)
    if (!is.null(treeInfo)) cat(treeInfo[i], "\n")
  }
  invisible(rv)
}


##' Plot trees from bed file produced by smc2bed ----
##' @inheritParams plotTrees
##' @param file A bed file produced by smc2bed
##' @param iter The MCMC iteration to use. "max" means use the last iteration in the file.
##' NULL or a value < 0 means to plot all iterations.
##' @param chrom The chromosome of the region where plots are desired. See note below.
##' @param start The start coordinate of region where plots are desired. See note below.
##' @param end The end coordinate of region where plots are desired. See note below.
##' @param interval Only plot trees separated by this interval in base pairs
##' @param regionSide Print region of each tree to the margins on this side
##' (1=bottom, 2=left, 3=top, 4=right, anything else = don't print)
##' @param regionLine Print region of each tree on this line of the margin
##' @param regionRep If TRUE, include MCMC rep in region string
##' @param xlab label for x axis
##' @param popwidth If mod is not NULL, relative widths for each population
##' @note If the input file is bgzipp'ed and tabixed (which is done automatically when
##' created with the script smc2bed-all), then tabix can be used to read in the file. This
##' will be much more efficient if the region chrom:start-end covers a subset of the region
##' covered by the entire file.
##' @note This creates a new plot for each tree. If plotting to the screen, probably
##' want to call par(ask=TRUE) first.
##' @export
plotTreesFromBed <- function(file=NULL, iter="max", chrom=NULL, start=-1, end=-1,
                             prune=NULL, keepSeqs=NULL, col="black", leafCol=col, 
                             leafLabels=NULL, interval=1, timeScale=1, drawSpr=FALSE,
                             ylab="Generations", xlab="", logScale=FALSE, ylim=NULL, mar=c(8,4,1,1),
                             add=FALSE,  #sitesFile=NULL,
                             regionSide=1, regionLine=4, regionRep=TRUE,
                             treeInfo=NULL, mod=NULL, popwidth=NULL,
                             ...) {
  #    if (!is.null(sitesFile)) {
  #        sites <- readSites(sitesFile)
  #    } else sites <- NULL
  if (file.exists(paste0(file, ".tbi")) && start >= 0 && end >= start) {
    if (is.null(chrom)) {
      chrom <- scan(pipe(sprintf("tabix -l %s", file)), what=character())
      if (length(chrom) > 1) {
        stop("Need to give chromosome")
      }
    }
    x <- read.table(pipe(sprintf("tabix %s %s:%.0f-%.0f", file, chrom,start, end)), header=FALSE, stringsAsFactors=FALSE)
    names(x) <- c("chrom", "chromStart", "chromEnd", "MCMC", "tree")
  } else {
    x <- scan(file, what=character(), nmax=1, quiet=TRUE)
    if (substr(x, 1, 1) == "#") {
      ## If there is a header, then this is likely output from arg-summarize.
      x <- readArgSummary(file)
      if (!is.element("tree", names(x)))
        stop("File does not have a column named tree")
    } else {
      x <- read.table(file, header=FALSE, stringsAsFactors=FALSE)
      names(x) <- c("chrom", "chromStart", "chromEnd", "MCMC", "tree")
    }
  }
  if (!is.null(iter)) {
    if (iter == "max") {
      iter <- max(x$MCMC)
      cat("using iter ", iter, "\n")
    }
    if (iter >= 0) x <- x[x$MCMC==iter,]
  }
  if (start >= 0) x <- x[x$chromEnd >= start,]
  if (end >= 0) x <- x[x$chromStart <= end,]
  if (interval > 1) {
    if (start < 0) start <- min(x$chromStart)
    if (end < 0) end <- max(x$chromEnd)
    keep <- rep(FALSE, nrow(x))
    for (i in seq(start, end, interval)) {
      w <- which(x$chromStart <= i & x$chromEnd > i)
      if (length(w) > 0) keep[w] <- TRUE
    }
    x <- x[keep,]
  }
  if (nrow(x) == 0) {
    warning("No trees")
    return(invisible(NULL))
  }
  rv <- list()
  if (is.null(treeInfo))
    treeInfo <- sprintf("%s:%.0f-%.0f", x[,1], x[,2]+1,x[,3])
  if (regionRep) treeInfo <- sprintf("%s rep=%i", treeInfo, x[,4])
  rv$plots <- plotTrees(x$tree, prune=prune, keepSeqs=keepSeqs,
                        treeInfo=treeInfo,
                        col=col, leafCol=leafCol, leafLabels=leafLabels, timeScale=timeScale,
                        drawSpr=drawSpr,
                        ylab=ylab, logScale=logScale, ylim=ylim, add=add, mar=mar, xlab=xlab,
                        #              sites=if(is.null(sites)) {NULL} else {sites[sites$chromStart > x[,2] & sites$chromEnd < x[,3],]},
                        regionSide=regionSide, regionLine=regionLine, mod=mod, popwidth=popwidth,
                        ...)
  rv$input <- x
  invisible(rv)
}