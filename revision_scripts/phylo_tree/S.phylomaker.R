S.PhyloMaker<-function (tree, spList, nodes, output.spList = T, scenarios = c("S1", "S2", "S3")) 
{
  options(scipen=999)
  tree0 <- tree
  spList[sapply(spList, is.factor)] <- lapply(spList[sapply(spList, is.factor)], as.character)
  if (any(duplicated(spList$species))) 
  {
    warning("Duplicated species detected and removed.")
    print(spList$species[duplicated(spList$species)])
  }
  spList <- spList[!duplicated(spList$species), ]
  spList.original <- spList
  spList$species <- gsub(" ", "_", spList$species)
  spList$species <- gsub("(^[[:alpha:]])", "\\U\\1", spList$species, perl = TRUE)
  spList$genus <- gsub("(^[[:alpha:]])", "\\U\\1", spList$genus, perl = TRUE)
  spList$family <- gsub("(^[[:alpha:]])", "\\U\\1", spList$family, perl = TRUE)
  rnN <- data.frame(node.label = paste("N", 1:length(tree$node.label), sep = ""), oriN = tree$node.label, stringsAsFactors = FALSE)
  nodes[,c("level","family","genus","rn","bn","taxa")]<-lapply(nodes[,c("level","family","genus","rn","bn","taxa")], as.character)
  tree$node.label <- paste("N", 1:length(tree$node.label), sep = "")
  kk<-c()
  for (i in 1:length(tree$tip.label)) {
    kk<-c(kk,substring(tree$tip.label[i],1,gregexpr("_",tree$tip.label[i])[[1]][1]-1))
  }
  m<-data.frame(num=1:length(kk),genus=kk,species=tree$tip.label)
  m<-merge(m,nodes[,c("genus","family")])
  mX <- m
  m <- m[,c("genus","family")]
  m <- m[!duplicated(m$genus),]
  dimnames(m)[[2]][2] <- "family_in_PhytoPhylo"
  m<-m[,c("genus","family_in_PhytoPhylo")]
  m0 <- spList[!duplicated(spList$genus), c("genus","family")]
  dimnames(m0)[[2]][2] <- "family_in_spList"
  mm<-merge(m0, m)
  g<-mm[which(is.na(match(paste(mm$genus,mm$family_in_spList,sep="_"),paste(mm$genus,mm$family_in_PhytoPhylo,sep="_")))),]
  if (dim(g)[1]>0)
  {
    print("Taxonomic classification not consistent between spList and PhytoPhylo.")
    print(g) 
  }
  add.tip <- spList[which(is.na(match(spList$species, tree$tip.label))), ]
  status <- rep("match(prune)", dim(spList)[1])
  status[which(is.na(match(spList$species, tree$tip.label)))] <- "match(add)"
  if (dim(add.tip)[1] == 0 & length(na.omit(match(spList$species, tree$tip.label))) == 0)
    stop("Incorrect format of species list.")
  if (length(setdiff(spList$species, tree0$tip.label)) == 0 & length(na.omit(match(spList$species, tree$tip.label))) > 0)
  {
    print("There is no species needs to be added, all the species are pruned from PhytoPhylo.")
    splis <- spList.original
    treeX <- drop.tip(tree0, setdiff(tree0$tip.label, splis$species))
    splis$status <- "match(prune)"
    phylo0 <- list(Scenario.1 = NULL, Scenario.2 = NULL, Scenario.3 = NULL, Species.list = splis)
    if ("S1" %in% scenarios) {phylo0$Scenario.1 <- treeX}
    if ("S2" %in% scenarios) {phylo0$Scenario.2 <- treeX}
    if ("S3" %in% scenarios) {phylo0$Scenario.3 <- treeX}
    phylo0[sapply(phylo0, is.null)] <- NULL
    return(phylo0)
    stop()
  }
  add.tip$sort <- ""
  add.tip$sort[which(!is.na(match(add.tip$genus, nodes[nodes$level == "G", ]$genus)))] <- "G1"
  add.tip$sort[which(is.na(match(add.tip$genus, nodes[nodes$level == "G", ]$genus)) & !is.na(match(add.tip$family, nodes[nodes$level == "F", ]$family)))] <- "F1"
  add.tip$sort[add.tip$sort == "F1"][duplicated(add.tip[add.tip$sort == "F1", ]$genus)] <- "F2"
  a <- which(add.tip$sort == "")
  if (length(a) > 0)
  {  
    print(paste("Note:", length(a), "taxa unmatch:",sep=" "))
    print(add.tip$species[a])
    status[match(add.tip$species[a], spList$species)] <- "unmatch"         
  }  
  spList.original$status <- status
  if ("S1" %in% scenarios) {
    t1 <- tree
    rnN1<-rnN
    nG <- nodes[nodes$level == "G", ]
    nF <- nodes[nodes$level == "F", ]
    data <- add.tip[add.tip$sort == "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- match(data$family[i], nF$family)
        g <- nF$gen.n[n]
        s <- nF$sp.n[n]
        if (g == 1 & s == 1) {                                                                          
          num <- grep(nF$taxa[n], t1$tip.label)
          len <- t1$edge.length[match(num, t1$edge[, 2])]
          t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num, position = len)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- grep(nF$taxa[n], t1$tip.label)
          t1$node.label[match(t1$edge[match(num, t1$edge[, 2]), 1], unique(t1$edge[, 1]))] <- paste("NN", t1$Nnode + 1, sep = "")
          rnN1$node.label[match(nF$bn[n], rnN1$node.label)]<- paste("NN", t1$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t1$Nnode + 1, sep = "")
          nF$bn.bl[n] <- len
        }             
        else {
          num <- unique(t1$edge[, 1])[match(nF$bn[n], t1$node.label)]
          len <- nF$bn.bl[n]
          t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num)
        }
      }
    }
    data <- add.tip[add.tip$sort != "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- grep(paste(data$genus[i], "_", sep = ""), t1$tip.label)
        if (length(n) == 1) {
          num <- n
          len <- t1$edge.length[match(num, t1$edge[, 2])]
          t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num, position = len)
        }
        if (length(n) > 1) {
          num <- fastMRCA(t1, t1$tip.label[min(n)], t1$tip.label[max(n)])
          len <- fastDist(t1, t1$tip.label[min(n)], t1$tip.label[max(n)])/2
          t1 <- bind.tip(t1, tip.label = data$species[i], edge.length = len, where = num)
        }   
      }
    }
    toDrop <- setdiff(1:length(t1$tip.label), which(!is.na(match(t1$tip.label, spList$species))))
    t1 <- drop.tip(t1, tip = toDrop)
    Re <- which(!is.na(match(t1$node.label, rnN1$node.label)))
    noRe <- which(is.na(match(t1$node.label, rnN1$node.label)))
    t1$node.label[Re] <- rnN1$oriN[match(t1$node.label, rnN1$node.label)[Re]]
    t1$node.label[noRe] <- ""
  }
  else {
    t1 <- NULL
  }
  if ("S2" %in% scenarios) {
    t2 <- tree
    rnN2<-rnN
    nG <- nodes[nodes$level == "G", ]
    nF <- nodes[nodes$level == "F", ]
    data <- add.tip[add.tip$sort == "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- match(data$family[i], nF$family)
        g <- nF$gen.n[n]
        s <- nF$sp.n[n]
        if (g == 1 & s == 1) {
          num <- grep(nF$taxa[n], t2$tip.label)
          len <- t2$edge.length[match(num, t2$edge[,2])] * sample((1:99)/100,1)
          t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)    
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- grep(data$species[i], t2$tip.label)
          t2$node.label[match(t2$edge[match(num, t2$edge[, 2]), 1], unique(t2$edge[,1]))] <- paste("NN", t2$Nnode + 1, sep = "")
          rnN2$node.label[match(nF$bn[n], rnN2$node.label)]<- paste("NN", t2$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t2$Nnode + 1, sep = "")
          nF$bn.bl[n] <- len
        }
        else {
          num <- unique(t2$edge[, 1])[match(nF$bn[n], t2$node.label)]
          len <- t2$edge.length[match(num,t2$edge[,2])] * sample((1:99)/100,1)
          t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- grep(data$species[i], t2$tip.label)
          t2$node.label[match(t2$edge[match(num, t2$edge[, 2]), 1], unique(t2$edge[, 1]))] <- paste("NN", t2$Nnode + 1, sep = "")
          rnN2$node.label[match(nF$bn[n], rnN2$node.label)]<- paste("NN", t2$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t2$Nnode + 1, sep = "")
          nF$bn.bl[n] <- nF$bn.bl[n]+len
        }
      }
    }
    data <- add.tip[add.tip$sort != "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- grep(paste(data$genus[i], "_", sep = ""), t2$tip.label)
        if (length(n) == 1) {
          num <- n
          len <- t2$edge.length[match(num, t2$edge[, 2])] * sample((1:99)/100,1)
          t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)
        }
        if (length(n) > 1) {
          num <- sample(n,1)
          len <- t2$edge.length[match(num, t2$edge[, 2])] * sample((1:99)/100,1)
          t2 <- bind.tip(t2, tip.label = data$species[i], edge.length = len, where = num, position = len)
        }   
      }
    }
    toDrop <- setdiff(1:length(t2$tip.label), which(!is.na(match(t2$tip.label, spList$species))))
    t2 <- drop.tip(t2, tip = toDrop)
    Re <- which(!is.na(match(t2$node.label, rnN2$node.label)))
    noRe <- which(is.na(match(t2$node.label, rnN2$node.label)))
    t2$node.label[Re] <- rnN2$oriN[match(t2$node.label, rnN2$node.label)[Re]]
    t2$node.label[noRe] <- ""
  }
  else {
    t2 <- NULL
  }
  if ("S3" %in% scenarios) {
    t3 <- tree
    rnN3<-rnN
    nG <- nodes[nodes$level == "G", ]
    nF <- nodes[nodes$level == "F", ]
    data <- add.tip[add.tip$sort == "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- match(data$family[i], nF$family)
        g <- nF$gen.n[n]
        s <- nF$sp.n[n]
        if (g == 1 & s == 1) {                                                                          
          num <- grep(nF$taxa[n], t3$tip.label)
          len <- t3$edge.length[match(num, t3$edge[, 2])] * (2/3)
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num, position = len)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- grep(nF$taxa[n], t3$tip.label)
          t3$node.label[match(t3$edge[match(num, t3$edge[, 2]), 1], unique(t3$edge[, 1]))] <- paste("NN", t3$Nnode + 1, sep = "")
          rnN3$node.label[match(nF$bn[n], rnN3$node.label)]<- paste("NN", t3$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t3$Nnode + 1, sep = "")
          nF$bn.bl[n] <- len
        }
        if (g == 1 & s > 1) {
          num <- unique(t3$edge[, 1])[match(nF$bn[n], t3$node.label)]
          if ((2/3)*nF$rn.bl[n] <= nF$bn.bl[n])  { len <-(nF$rn.bl[n]-nF$bn.bl[n])/2 }
          if ((2/3)*nF$rn.bl[n] > nF$bn.bl[n])   { len <-nF$rn.bl[n]*2/3-nF$bn.bl[n] }
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num, position = len)
          nF$gen.n[n] <- g + 1
          nF$sp.n[n] <- s + 1
          num <- unique(t3$edge[, 1])[match(nF$bn[n], t3$node.label)]
          t3$node.label[match(t3$edge[match(num, t3$edge[, 2]), 1], unique(t3$edge[, 1]))] <- paste("NN", t3$Nnode + 1, sep = "")
          rnN3$node.label[match(nF$bn[n], rnN3$node.label)]<- paste("NN", t3$Nnode + 1, sep = "")
          nF$bn[n] <- paste("NN", t3$Nnode + 1, sep = "")
          nF$bn.bl[n] <- len+nF$bn.bl[n]
        }
        if (g > 1) {
          num <- unique(t3$edge[, 1])[match(nF$bn[n], t3$node.label)]
          len <- nF$bn.bl[n]
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num)
        }
      }
    }
    data <- add.tip[add.tip$sort != "F1", ]
    if (dim(data)[1] > 0) {
      for (i in 1:dim(data)[1]) {
        n <- grep(paste(data$genus[i], "_", sep = ""), t3$tip.label)
        if (length(n) == 1) {
          len <- t3$edge.length[match(n, t3$edge[, 2])]/2
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = n, position = len)
          nG$sp.n[match(data$genus[i], nG$genus)] <- length(n) + 1
          nu <- grep(paste(data$genus[i], "_", sep = ""), t3$tip.label)
          num <- fastMRCA(t3, t3$tip.label[nu[1]], t3$tip.label[nu[2]])
          t3$node.label[match(num, unique(t3$edge[, 1]))] <- paste("NN", t3$Nnode + 1,  sep = "")
          rnN3$node.label[match(nG$bn[n], rnN3$node.label)]<- paste("NN", t3$Nnode + 1, sep = "")
          nG$bn[match(data$genus[i], nG$genus)] <- paste("NN", t3$Nnode + 1, sep = "")
          nG$bn.bl[match(data$genus[i], nG$genus)] <- len
        }
        if (length(n) > 1) {
          num <- fastMRCA(t3, t3$tip.label[min(n)], t3$tip.label[max(n)])
          len <- as.numeric(branching.times(t3))[match(num, unique(t3$edge[, 1]))]
          t3 <- bind.tip(t3, tip.label = data$species[i], edge.length = len, where = num)
        }
      }
    }
    toDrop <- setdiff(1:length(t3$tip.label), which(!is.na(match(t3$tip.label, spList$species))))
    t3 <- drop.tip(t3, tip = toDrop)
    Re <- which(!is.na(match(t3$node.label, rnN3$node.label)))
    noRe <- which(is.na(match(t3$node.label, rnN3$node.label)))
    t3$node.label[Re] <- rnN3$oriN[match(t3$node.label, rnN3$node.label)[Re]]
    t3$node.label[noRe] <- ""
  }
  else {
    t3 <- NULL
  }
  if (output.spList == FALSE) 
    spList <- NULL
  phylo <- list(Scenario.1 = t1, Scenario.2 = t2, Scenario.3 = t3, Species.list = spList.original)
  phylo[sapply(phylo, is.null)] <- NULL
  return(phylo)
}
