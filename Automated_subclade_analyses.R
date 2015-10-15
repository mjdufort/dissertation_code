## This is is one of several files containing scripts and functions used in processing and analysis of data for Matthew Dufort's Ph.D. dissertation at the University of Minnesota, titled "Coexistence, Ecomorphology, and Diversification in the Avian Family Picidae (Woodpeckers and Allies)."

## this file contains scripts and functions to calculate variables at the subclade level, and to test for relationships between subclade variables


### load packages and data

## load necessary packages
library(ape)
library(geiger)
library(phytools)
library(nlme)
library(laser)
library(DDD)

load(file="Picidae_data_for_distribution_morphology_evolution.RData") # load data needed from morphology and distribution analyses (from Morphology_data_processing.R)
load(file="Picidae_BAMM_data_for_automated_subclade_analyses.RData")  # load data objects needed from BAMM analyses (from BAMM_data_prep_and_processing.R)


### generate necessary functions for automated subclade analyses

## the function extractSubclades.all() extracts all subclades with at least min.taxa and at most max.taxa tips from the tree
# as input, it takes phy (a phylogenetic tree of class phylo), min.taxa (the minimum number of taxa for a subclade to be included), and max.taxa (the maximum number of taxa for a subclade to be include)
# it returns a list of phy objects, one for each subclade, with the list elements named with the node numbers from the original tree
extractSubclades.all <- function(phy, min.taxa=4, max.taxa=Inf) {
  require(ape)
  subclades <- list()  # initialize a list to store subclades
  ntaxa <- length(phy$tip.label)  # get the number of taxa in the tree
  j <- 0
  for (i in (ntaxa + 1):(ntaxa + phy$Nnode)) {  # loop over internal nodes
    clade.tmp <- extract.clade(phy, node=i)  # extract current subclade (the subclade descending from the current node)
    
    # if current subclade meets specifications, add it to the list
    if ((length(clade.tmp$tip.label) >= min.taxa) & (length(clade.tmp$tip.label) <= max.taxa)) {
      j <- j+1
      subclades[[j]] <- clade.tmp
      names(subclades)[j] <- as.character(i)
    }
  }
  return(subclades)
}

## the function getSubclades.withData() extracts subclades from a tree, including only subclades that have sufficient taxa with data in a vector or matrix of trait data
# as input, it takes phylist (a list of subclades, each a phylo object), taxondata (the vector of matrix of trait data, with names or rownames corresponding to taxon names), inc.data (boolean to return a treedata object for each subclade; if FALSE, returns a phylo object for each subclade), and min.taxa (the minimum number of taxa for a subclade to be included)
# it returns a list of subclades, either as phylo objects or treedata objects, with the list elements named by the node numbers in the original tree 
getSubclades.withData <- function(phylist, taxondata, inc.data=TRUE, min.taxa=4, quiet=TRUE) {
  require(geiger)
  subclades.new <- list()
  
  # get taxon names from trait data vector/matrix
  if (is.matrix(taxondata)) {
    taxon.names <- rownames(taxondata)
  } else if (is.vector(taxondata)) {
    taxon.names <- names(taxondata)
  }
  
  # loop over subclades in phylist, testing if each subclade has more than min.taxa in the data vector/matrix
  j <- 0
  for (i in 1:length(phylist)) {
    if (!quiet) print(i)
    if (!quiet) print(phylist[[i]]$tip.label %in% taxon.names)
    if (sum(phylist[[i]]$tip.label %in% taxon.names) >= min.taxa) {
      j <- j + 1
      if (inc.data) {
        subclades.new[[j]] <- treedata(phy = phylist[[i]], data=taxondata, warnings=FALSE)
      } else {
        subclades.new[[j]] <- phylist[[i]]
      }
      names(subclades.new)[j] <- names(phylist)[i]
    }
  }
  return(subclades.new)
}

## the function areOverlappingSubclades() tests a list of two or more subclades (or trees) to determine if there is any overlap in the tips included
# as input, it takes phylist (a list of subclades, each a phylo object), and getoverlaps (boolean to return overlapping taxa)
# it returns TRUE if any of the subclades in phylist share taxa, and FALSE if there are no shared taxa among them; if getoverlaps=TRUE, it returns a list containing the test value (TRUE or FALSE), and a vector of the taxa shard among subclades
areOverlappingSubclades <- function(phylist, getoverlaps=FALSE) {
  
  # generate a character vector containing the concatenated taxa from each subclade
  taxnames <- character()
  if (!is.null(phylist[[1]]$data)) {  # checks if they're treedata objects
    taxnames <- unlist(lapply(phylist, FUN = function(x) x$phy$tip.label))
  } else {
    taxnames <- unlist(lapply(phylist, FUN = function(x) x$tip.label))
  }
  
  # check for duplicates in the vector of taxon names
  duplicates <- duplicated(taxnames)
  if (!any(duplicates)) {
    return(FALSE)
  } else if (!getoverlaps) {
    return(TRUE)
  } else {
    return(list(test=TRUE, overlaps=taxnames[which(duplicates)]))}
}

## the function subcladeCombinations.all() determines all sets of reciprocally monophyletic subclades meeting a set of criteria, and returns them as a list of lists of phylo objects
# this sped-up version generates a pairwise matrix of overlapping clades, then checks if any of the subclades in the combination are TRUE in the matrix (and therefore takes advantage of speed-ups with vectorization)
# for large trees, there is a VERY large number of possible combinations, and using this function is not advisable
# as input, it takes phylist (a list of subclades, each a phylo object), min.clades (the minimum number of clades to include in a combination), and max.clades (the maximum number of clades to include in a combination)
# it returns a list of subclade combinations, each a list of phylo objects
subcladeCombinations.all <- function(phylist, min.clades=3, max.clades=Inf) {
  if (max.clades > length(phylist)) max.clades <- length(phylist)
  
  # generate matrix of pairwise subclade overlaps
  subclade.overlap.pairwise <- matrix(nrow=length(phylist), ncol=length(phylist))
  for (i in 1:nrow(subclade.overlap.pairwise)) {
    for (j in 1:ncol(subclade.overlap.pairwise)) {
      subclade.overlap.pairwise[i,j] <- areOverlappingSubclades(list(phylist[[i]], phylist[[j]]))
    }
  }
  
  subclade.names <- names(phylist)  # get the subclade names
  combinations <- list()  # initialize list to store subclade combinations
  
  complete <- FALSE  # boolean to end search
  k <- 0
  for (nclades in (min.clades:max.clades)) {  # loop over number of subclades to include in set
    if (!complete) {
      length.last <- length(combinations)
      combinations.to.test <- combn(x=(1:length(subclade.names)), m=nclades, simplify=TRUE)  # generate a matrix of combinations to test
      print(paste("Testing ", ncol(combinations.to.test), " combinations for ", nclades, " clades.", sep=""))
      
      # test each proposed combination for reciprocal monophyly; if they are reciprocally monophyletic, add to list
      for (i in 1:ncol(combinations.to.test)) {
        if ((i %% 10000) == 0) print(paste("Testing combination ",i, sep=""))
        pairwise.combinations.temp <- combn(x=combinations.to.test[,i], m=2, simplify=TRUE)
        if (!any(subclade.overlap.pairwise[cbind(pairwise.combinations.temp[1,],pairwise.combinations.temp[2,])])) {
          k <- k+1
          combinations[[k]] <- subclade.names[combinations.to.test[,i]]
        }
      }
      
      # test if any combinations were added for this number of subclades, and terminate if none were
      if (length(combinations)==length.last) {
        complete <- TRUE
        print(paste("No successful combinations for ", nclades, " clades; stopping search.", sep=""))
      }
    }
  }
  return(combinations)
}

## the function subcladeCombinations.random() generates a random sample of combinations of reciprocally monophyletic subclades meeting a set of criteria
# this samples by selecting a subclade at random, then selecting another from all the possible subclades that don't overlap the first, and continuing doing that until there aren't any more possibilities; this approach probably leads to the same subclades being selected repeatedly, as certain isolated subclades are almost always going to be suitable
# as input, it takes phylist (a list of subclades, each a phylo object), ncombs (the maximum number of combinations to return), min.clades (the minimum number of subclades to include in a combination), max.clades (the maximum number of subclades to include in a combination), min.taxa (the minimum number of taxa for a subclade to be considered for inclusion), max.fails (the maximum number of failures before halting the search), and report (boolean to output status updates to console)
# it returns a list of subclade combinations, each a list of phylo objects
subcladeCombinations.random <- function(phylist, ncombs=1000, min.clades=5, max.clades=Inf, min.taxa=4, max.fails=1e6, report=TRUE) {
  
  # check if the objects are phylo objects or treedata objects; also drop subclades with fewer taxa than the minimum
  for (i in names(phylist)) {
    if (class(phylist[[i]]) != "phylo") {
      if (class(phylist[[i]]$phy) == "phylo") {
        phylist[[i]] <- phylist[[i]]$phy
      } else {
        cat("\nError: item ", i, " in phylist is not a phylo or treedata object.\n", sep="")
        return()
      }
    }
    
    if (length(phylist[[i]]$tip.label) < min.taxa) phylist[[i]] <- NULL  # drop subclades with too few taxa
  }
  
  if (max.clades > length(phylist)) max.clades <- length(phylist)
  subclade.names <- names(phylist)  # extract the subclade names
  
  # generate matrix of pairwise subclade overlaps
  subclade.overlap.pairwise <- matrix(nrow=length(phylist), ncol=length(phylist), dimnames=list(subclade.names, subclade.names))
  for (i in 1:nrow(subclade.overlap.pairwise)) {
    for (j in 1:ncol(subclade.overlap.pairwise)) {
      subclade.overlap.pairwise[i,j] <- areOverlappingSubclades(list(phylist[[i]], phylist[[j]]))
    }
  }
  
  combinations <- list() # the combinations that will be returned
  all.done <- FALSE
  z <- 1
  fails <- 0
  while ((length(combinations) < ncombs) & (!all.done)) {
    combination.done <- FALSE
    combination.temp <- sample(x=subclade.names, size=1) # pick the first subclade in the possible combination
    q <- 1
    while ((length(combination.temp) < max.clades) & (!combination.done)) {
      subclades.possible.additions <- colnames(subclade.overlap.pairwise)[which(rowSums(as.matrix(subclade.overlap.pairwise[,combination.temp]))==0)] # this finds all subclades that don't overlap with any of the subclades already in the combination
      if (length(subclades.possible.additions) == 0) {
        combination.done <- TRUE
      } else {
        q <- q + 1
        combination.temp[q] <- sample(x=subclades.possible.additions, size=1)
      }
    }
    combination.temp <- sort(combination.temp)
    if ((length(combination.temp) >= min.clades) & (length(which(sapply(combinations, identical, combination.temp, simplify=TRUE)==TRUE)) < 1)) {
      combinations[[z]] <- combination.temp
      cat("Found combination ", z, "\n", sep="")
      z <- z + 1
    } else {
      fails <- fails+1
    }
    if (fails == max.fails) {
      all.done <- TRUE
      print(paste("Reached maximum failures. Returning", length(combinations), "combinations"))
    }
  }
  return(combinations)
}

## the function subcladeCombinations.sequential() determines a set of combinations of reciprocally monophyletic subclades by working its way down the tree; it slices the tree at each node and determines all valid subclades below that slice
# as input, it takes phy (a tree as a phylo object), min.taxa (the minimum number of taxa for a subclade to be included), min.clades (the minimum number of subclades to include in a combination), and max.clades (the maximum number of subclades to include in a combination)
# it returns a list of subclade combinations, each a list of phylo objects
subcladeCombinations.sequential <- function(phy, min.taxa=4, min.clades=5, max.clades=Inf) {
  require(ape)
  combinations <- list()
  
  phy.nodedepth.sorted <- sort((max(branching.times(phy)) - branching.times(phy)), decreasing=FALSE) # generate a vector of node depths
  
  l <- 0
  for (i in 1:length(phy.nodedepth.sorted)) {
    candidate.nodes <- phy$edge[,2][(node.depth.edgelength(phy)[phy$edge[,1]] <= phy.nodedepth.sorted[i]) & (node.depth.edgelength(phy)[phy$edge[,2]] > phy.nodedepth.sorted[i]) & (phy$edge[,2] > length(phy$tip.label))]  # find all the descendant nodes from edges cut at current step in phy.nodedepth.sorted
    
    # identify nodes just below the branching point I'm examining
    candidate.nodes <- candidate.nodes[candidate.nodes > length(phy$tip.label)]
    
    # extract combination (if possible) from list of descendant subclades
    if (length(candidate.nodes) >= min.clades) {
      candidate.combination <- character()
      for (j in 1:length(candidate.nodes)) {
        
        if (length(extract.clade(phy, node=candidate.nodes[j], root.edge=0)$tip.label) >= min.taxa) {
          candidate.combination <- c(candidate.combination, candidate.nodes[j])
        }
      }
      if ((length(candidate.combination) >= min.clades) & (length(candidate.combination) <= max.clades)) {
        l <- l + 1
        combinations[[l]] <- candidate.combination
      }
    }
  }
  combinations <- combinations[!duplicated(combinations)]
  return(combinations)
}


## this function determines all members of each subclade, including those not in the tree; it uses a list of taxon proxies, and checks these proxies against the actual taxa in the subclade; it has several options for returning these taxa
# as input, it takes phylist (a list of subclades, each a phylo object or treedata object), taxon.proxies (a list containing a vector of proxies for each taxon), and to_return (a switching variable, which allows the user to select whether to return the missing taxa ("missing"), all taxa ("full"), or a list of the included and missing taxa ("split"))
# it returns a list containing vectors with the set of taxa specified by to_return
subclades.fulltaxlist <- function(phylist, taxon.proxies, to_return="full") {
  subclades.taxa_to_include <- list()
  for (i in 1:length(phylist)) {
    subclades.taxa_to_include.temp <- character()
    for (j in 1:length(taxon.proxies)) {  # loop over list of taxa
      
      # if all proxies are included in the subclade, add the current taxon to the list of included taxa
      if (all(taxon.proxies[[j]] %in% phylist[[i]]$tip.label) | all(taxon.proxies[[j]] %in% phylist[[i]]$phy$tip.label)) {
        subclades.taxa_to_include.temp <- c(subclades.taxa_to_include.temp, names(taxon.proxies)[j])
      }
    }
    subclades.taxa_to_include[[i]] <- switch(to_return, 
           missing = subclades.taxa_to_include.temp,
           full = c(phylist[[i]]$tip.label, subclades.taxa_to_include.temp),
           split = list(included=phylist[[i]]$tip.label, missing=subclades.taxa_to_include.temp))
  }
  names(subclades.taxa_to_include) <- names(phylist)
  return(subclades.taxa_to_include)
}

# the function subclade.fulltaxlist() is the same as subclades.fulltaxlist(), but it acts only on a single subclade; this allows looping or applying over a list of treedata objects and adding the full membership to the treedata object
# as input, it takes phy (a subclades, either a phylo object or treedata object), taxon.proxies (a list containing a vector of proxies for each taxon), and to_return (a switching variable, which allows the user to select whether to return the missing taxa ("missing"), all taxa ("full"), or a list of the included and missing taxa ("split"))
# it returns a vector with the set of taxa specified by to_return
subclade.fulltaxlist <- function(phy, taxon.proxies, to_return="full") {
  taxa_to_include.tmp <- character()
  for (j in 1:length(taxon.proxies)) {
    if (all(taxon.proxies[[j]] %in% phy$tip.label)) {
      taxa_to_include.tmp <- c(taxa_to_include.tmp, names(taxon.proxies)[j])
      }
    }
  taxa.to_include <- switch(to_return, 
                            missing = taxa_to_include.tmp,
                            full = c(phy$tip.label, taxa_to_include.tmp),
                            split = list(included=phy$tip.label, missing=taxa_to_include.tmp))
  return(taxa.to_include)
}
## the function getTreedata.subclades() extracts the backbone tree with subclades, and builds a treedata object including the subclade data
# as input, it takes phy (the full tree as a phylo object), subclade.combination (a vector containing the node numbers of the subclades), and subclade.data (the data for the subclades, as a matrix with node numbers as the rownames)
# it returns a treedata object, where the returned tree has only the subclades as tips, with the backbone of those nodes retained
getTreedata.subclades <- function(phy, subclade.combination, subclade.data) {
  subclade.data.selected <- subset(subclade.data, row.names(subclade.data) %in% subclade.combination)
  subclades.temp <- list()
  subclades.edge.length.temp <- numeric()
  
  # get the stem edge length for each subclade, and rename one tip in each subclade with teh subclade name
  for (i in 1:length(subclade.combination)) {
    subclades.temp[[i]] <- extract.clade(phy, node=as.numeric(subclade.combination[i]))
    subclades.edge.length.temp[i] <- phy$edge.length[which(phy$edge[,2]==as.numeric(subclade.combination[i]))] # find the stem edge length for the subclade
    phy$tip.label[phy$tip.label==subclades.temp[[i]]$tip.label[1]] <- subclade.combination[i] # rename one tip with the name of the subclade
  }
  
  # loop over subclades, dropping all tips but the one set to the subclade name above; this is done separately, as dropping tips could change the node numbers and make the step above not work properly
  for (i in 1:length(subclade.combination)) {
    phy <- drop.tip(phy, tip=subclades.temp[[i]]$tip.label[-1])  # drop the remaining tips from the subclade
    phy$edge.length[which(phy$edge[,2]==which(phy$tip.label==subclade.combination[i]))] <- subclades.edge.length.temp[i]  # finds the edge that has the subclade name as its descendant node, and changes the length
  }
  
  phy.treedata <- treedata(phy, data=subclade.data.selected, warnings=FALSE)  # generate treedata object with backbone tree and subclade data
  return(phy.treedata)
}


### generate full taxon lists (to match taxa not on tree with subclades)

## the function read.taxon.proxy.list() reads a file of taxon proxies and formats them for later use
# as input, it takes filename (the location of the file containing the taxon proxies, with each taxon name followed by all the proxy taxa that must be present for the focal taxon to be included)
# it returns a list of character vectors, where each list element is named with the focal taxon name, and the vector contains all the proxy taxa that must be present for the focal taxon to be included
read.taxon.proxy.list <- function(filename) {
  taxon.proxy.list <- strsplit(scan(file=filename, what="character", sep="\n"), split=",") # read in file as a list of character vectors
  names(taxon.proxy.list) <- sapply(taxon.proxy.list, function(x) x[1]) # set the first element in the character vector to be the name
  taxon.proxy.list <- lapply(taxon.proxy.list, function(x) x[-1]) # remove that first element
  for (i in names(taxon.proxy.list)) {
    if (length(taxon.proxy.list[[i]]) == 0) taxon.proxy.list[[i]] <- NULL
  }  # this drops empty lists, so that the only ones retained are ones that actually have proxies
  return(taxon.proxy.list)  
}

## read in files of taxon proxies
picidae.RAxML.taxon.subclade.proxies <- list()
picidae.RAxML.taxon.subclade.proxies[["full_tree"]] <- read.taxon.proxy.list(filename="picidae_taxon_proxies_for_automated_subclade_analyses_full_tree.csv")
picidae.RAxML.taxon.subclade.proxies[["morph_tree"]][["all_inds"]] <- read.taxon.proxy.list(filename="picidae_taxon_proxies_for_automated_subclade_analyses_morph_tree_all_inds.csv")
picidae.RAxML.taxon.subclade.proxies[["morph_tree"]][["complete_ind_only"]] <- read.taxon.proxy.list(filename="picidae_taxon_proxies_for_automated_subclade_analyses_morph_tree_complete_ind_only.csv")

### extract subclades from full tree and morph trees

## extract subclades from full trees
picidae.RAxML.all.BEAST_calibrated.with_proxies.subclades <- extractSubclades.all(picidae.RAxML.all.BEAST_calibrated.with_proxies)

# extract subclades from morph trees
picidae.morph.log.fully_reduced.treedata.subclades <- list()
for (i in c("all_inds", "complete_ind_only")) {
  picidae.morph.log.fully_reduced.treedata.subclades[[i]] <- extractSubclades.all(picidae.morph.log.fully_reduced.treedata[[i]]$phy)
}
rm(i)


### generate treedata-like objects for each subclade, for the data variants I'm using

## combine all the data into a treedata-like object that has a bunch of different sets of data for each subclade, for picidae
picidae.morph.log.fully_reduced.subclades.treedata <- list()
for (i in names(picidae.morph.log.fully_reduced.treedata.subclades)) {  # loop over individual inclusion
  for (j in names(picidae.morph.log.fully_reduced.treedata.subclades[[i]])) {  # loop over subclades
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]] <- list()
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy"]] <- picidae.morph.log.fully_reduced.treedata.subclades[[i]][[j]]
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["missing_taxa"]] <- subclade.fulltaxlist(phy=picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy"]], taxon.proxies=picidae.RAxML.taxon.subclade.proxies[["morph_tree"]][[i]], to_return="missing")
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["missing_count"]] <- length(picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["missing_taxa"]])
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["node.full_tree"]] <- getMRCA(phy=picidae.RAxML.all.BEAST_calibrated.with_proxies, tip=picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy"]]$tip.label)
    
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy.full_tree"]] <- extract.clade(phy=picidae.RAxML.all.BEAST_calibrated.with_proxies, node=picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["node.full_tree"]])
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["missing_taxa.full_tree"]] <- subclade.fulltaxlist(phy=picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy.full_tree"]], taxon.proxies=picidae.RAxML.taxon.subclade.proxies[["full_tree"]], to_return="missing")
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["missing_count.full_tree"]] <- length(picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["missing_taxa.full_tree"]])
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["geomean"]] <- picidae.morph.log.fully_reduced.geomean[[i]][[j]][picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy"]]$tip.label] # pull in the geomean data
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phyl_pca"]] <- picidae.morph.log.fully_reduced.phyl_pca[[i]][[j]]$pca$S[picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy"]]$tip.label,] # pull in the unscaled PCA-rotated data
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["geomean_scaled.phyl_pca"]] <- picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]][[j]]$pca$S[picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy"]]$tip.label,] # pull in the geomean-scaled PCA-rotated data
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["overlaps.scaled"]] <- picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0[["mytax"]][["migratory"]][["overlaps.scaled"]][c(picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy"]]$tip.label, picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["missing_taxa"]])]
    
    # get overlaps scaled by focal taxon range and similarity in geomean, unscaled PCA, and geomean-scaled PCA
    picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["overlaps.euclidean_scaled"]] <- list()
    for (q in c("geomean", "phyl_pca", "geomean_scaled.phyl_pca")) {
      picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["overlaps.euclidean_scaled"]][[q]] <- picidae.summed_overlaps.shp.BirdLife.UnaryUnion.buffer0.euclidean_scaled[["migratory"]][[q]][[i]][[j]][["inc_no_phylo"]][c(picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["phy"]]$tip.label, picidae.morph.log.fully_reduced.subclades.treedata[[i]][[j]][["missing_taxa"]])]
    }
  }
}
rm(i,j,q)


### fit models of diversification and morphological evolution to the entire data set

## fit diversification models to full tree
picidae.divrate.models <- list()
picidae.divrate.models[["full_tree"]] <- list()
picidae.divrate.models[["full_tree"]][["constant"]] <- try(bd_ML(branching.times(picidae.RAxML.all.BEAST_calibrated.with_proxies), missnumspec=237-length(picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label), tdmodel=0))
picidae.divrate.models[["full_tree"]][["time_dependent"]] <- try(bd_ML(branching.times(picidae.RAxML.all.BEAST_calibrated.with_proxies), missnumspec=237-length(picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label), tdmodel=1, idparsopt=1:3, initparsopt=c(0.1, 0.05, 0.1)))
picidae.divrate.models[["full_tree"]][["diversity_dependent"]] <- try(dd_ML(branching.times(picidae.RAxML.all.BEAST_calibrated.with_proxies), missnumspec=237-length(picidae.RAxML.all.BEAST_calibrated.with_proxies$tip.label), ddmodel=1))

## calculate AICc for divrate models of full tree
for (i in names(picidae.divrate.models[["full_tree"]])) {
  picidae.divrate.models[["full_tree"]][[i]][["AICc"]] <- (-2 * picidae.divrate.models[["full_tree"]][[i]]$loglik) + (2 * picidae.divrate.models[["full_tree"]][[i]]$df) + (((2 * picidae.divrate.models[["full_tree"]][[i]]$df) * (picidae.divrate.models[["full_tree"]][[i]]$df + 1)) / (picidae.RAxML.all.BEAST_calibrated.with_proxies$Nnode - picidae.divrate.models[["full_tree"]][[i]]$df - 1))
}

## fit diversification models to morph tree
picidae.divrate.models[["morph_tree"]] <- list()
picidae.divrate.models[["morph_tree"]][["constant"]] <- try(bd_ML(branching.times(picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy), missnumspec=237-length(picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy$tip.label), tdmodel=0))
picidae.divrate.models[["morph_tree"]][["time_dependent"]] <- try(bd_ML(branching.times(picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy), missnumspec=237-length(picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy$tip.label), tdmodel=1, idparsopt=1:3, initparsopt=c(0.1, 0.05, 0.1)))
picidae.divrate.models[["morph_tree"]][["diversity_dependent"]] <- try(dd_ML(branching.times(picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy), missnumspec=237-length(picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy$tip.label), ddmodel=1))

## calculate AICc for divrate models of morph tree
for (i in names(picidae.divrate.models[["morph_tree"]])) {
picidae.divrate.models[["morph_tree"]][[i]][["AICc"]] <- (-2 * picidae.divrate.models[["morph_tree"]][[i]]$loglik) + (2 * picidae.divrate.models[["morph_tree"]][[i]]$df) + (((2 * picidae.divrate.models[["morph_tree"]][[i]]$df) * (picidae.divrate.models[["morph_tree"]][[i]]$df + 1)) / (picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy$Nnode - picidae.divrate.models[["morph_tree"]][[i]]$df - 1))
}

## summarize divrate model results with aicc
for (i in names(picidae.divrate.models)) {
  for (q in names(picidae.divrate.models[[i]])) {
    cat(i, q, picidae.divrate.models[[i]][[q]]$AICc, "\n", sep=" ")
  }
}
rm(i,q)

## fit morphological evolution models to morph tree with geomean, phyl_pca, and geomean_scaled.phyl_pca (with the same models I used below)
picidae.morphrate.models <- list()
for (i in c("all_inds", "complete_ind_only")) {
  cat("\nStarting model fitting for", i, "\n", sep=" ")
    # for geomean
  cat("\nStarting geomean models.\n")
  for (q in c("BM","OU","trend","EB")) {
    cat("Starting ", q, " model\n", sep="")
    picidae.morphrate.models[[i]][["geomean"]][[q]] <- fitContinuous(phy=picidae.morph.log.fully_reduced.treedata[[i]]$phy, dat=picidae.morph.log.fully_reduced.geomean[[i]][picidae.morph.log.fully_reduced.treedata[[i]]$phy$tip.label], model=q)
  }
  
  # for phyl_pca
  cat("\nStarting phyl_pca models.\n")
  for (q in c("BM","OU","trend","EB")) {
    cat("Starting ", q, " model\n", sep="")
    picidae.morphrate.models[[i]][["phyl_pca"]][[q]] <- fitContinuous(phy=picidae.morph.log.fully_reduced.treedata[[i]]$phy, dat=picidae.morph.log.fully_reduced.phyl_pca[[i]]$pca$S, model=q)
  }
  
  # for geomean_scaled.phyl_pca
  cat("\nStarting geomean_scaled phyl_pca models.\n")
  for (q in c("BM","OU","trend","EB")) {
    cat("Starting ", q, " model\n", sep="")
    picidae.morphrate.models[[i]][["geomean_scaled.phyl_pca"]][[q]] <- fitContinuous(phy=picidae.morph.log.fully_reduced.treedata[[i]]$phy, dat=picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[[i]]$pca$S, model=q)
  }
}

# summarize morphological evolution model fits (with AICc)
for (i in names(picidae.morphrate.models)) {
  for (q in names(picidae.morphrate.models[[i]])) {
    for (r in names(picidae.morphrate.models[[i]][[q]])) {
      if ("gfit" %in% class(picidae.morphrate.models[[i]][[q]][[r]])) {
        cat(i, q, r, picidae.morphrate.models[[i]][[q]][[r]]$opt$aicc, "\n", sep=" ")
      } else if ("gfits" %in% class(picidae.morphrate.models[[i]][[q]][[r]])) {
        cat(i, q, r, picidae.morphrate.models[[i]][[q]][[r]][[1]]$opt$aicc, "\n", sep=" ")
      }
    }
  }
}
rm(i,q,r)


### calculate subclade metrics

## the function calcMetrics.subclades() calculates a huge range of metrics for a list of subclades, including fitting models of diversification and trait evolution to the subclade
# as input, it takes subclades.treedata (a list of treedata-like objects, each containing a phy and other data objects for a single subclade), BAMM_divrates (the subclade average diversification rates from BAMM), BAMM_morphrates (the subclade average trait evolution rates from BAMM), metrics (a character vector containing the metrics to calculate), return_format (format of object to be returned; can be "matrix" or "list"), and quiet (boolean to output status to console)
# it returns either a matrix or list of metrics by subclade
calcMetrics.subclades <- function(subclades.treedata, BAMM_divrates=NULL, BAMM_morphrates=NULL, metrics=c("ntaxa", "ntaxa.on_morph_tree", "total_div", "crown_age", "divrate.ms.e10", "divrate.ms.e50", "divrate.ms.e90", "divrate.ML.constant.rate", "divrate.ML.constant.AICc", "divrate.ML.constant.AIC", "divrate.ML.time_dependent.rate", "divrate.ML.time_dependent.lambda1", "divrate.ML.time_dependent.mu1", "divrate.ML.time_dependent.AICc", "divrate.ML.time_dependent.AIC", "divrate.ML.diversity_dependent.rate", "divrate.ML.diversity_dependent.K", "divrate.ML.diversity_dependent.AICc", "divrate.ML.diversity_dependent.AIC", "divrate.BAMM", "divrate.BAMM.morph_tree", "gamma", "morphrate.geomean.BM.sigsq", "morphrate.geomean.BM.AICc", "morphrate.geomean.BM.AIC", "morphrate.geomean.OU.sigsq", "morphrate.geomean.OU.alpha", "morphrate.geomean.OU.AICc", "morphrate.geomean.OU.AIC", "morphrate.geomean.trend.slope", "morphrate.geomean.trend.sigsq", "morphrate.geomean.trend.AICc", "morphrate.geomean.trend.AIC", "morphrate.geomean.EB.alpha", "morphrate.geomean.EB.sigsq", "morphrate.geomean.EB.AICc", "morphrate.geomean.EB.AIC", "morphrate.geomean.delta.delta", "morphrate.geomean.delta.sigsq", "morphrate.geomean.delta.AICc", "morphrate.geomean.delta.AIC", "morphrate.geomean.BAMM", "morphrate.phyl_pca.BM.sigsq", "morphrate.phyl_pca.PC1.BM.AICc", "morphrate.phyl_pca.PC1.BM.AIC", "morphrate.phyl_pca.PC1.OU.sigsq", "morphrate.phyl_pca.PC1.OU.alpha", "morphrate.phyl_pca.PC1.OU.AICc", "morphrate.phyl_pca.PC1.OU.AIC", "morphrate.phyl_pca.PC1.trend.slope", "morphrate.phyl_pca.PC1.trend.sigsq", "morphrate.phyl_pca.PC1.trend.AICc", "morphrate.phyl_pca.PC1.trend.AIC", "morphrate.phyl_pca.PC1.EB.alpha", "morphrate.phyl_pca.PC1.EB.sigsq", "morphrate.phyl_pca.PC1.EB.AICc", "morphrate.phyl_pca.PC1.EB.AIC", "morphrate.phyl_pca.PC1.delta.delta", "morphrate.phyl_pca.PC1.delta.sigsq", "morphrate.phyl_pca.PC1.delta.AICc", "morphrate.phyl_pca.PC1.delta.AIC", "morphrate.phyl_pca.PC13.BAMM", "morphrate.geomean_scaled.phyl_pca.BM.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.BM.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.BM.AIC", "morphrate.geomean_scaled.phyl_pca.PC1.OU.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.OU.alpha", "morphrate.geomean_scaled.phyl_pca.PC1.OU.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.OU.AIC", "morphrate.geomean_scaled.phyl_pca.PC1.trend.slope", "morphrate.geomean_scaled.phyl_pca.PC1.trend.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.trend.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.trend.AIC", "morphrate.geomean_scaled.phyl_pca.PC1.EB.alpha", "morphrate.geomean_scaled.phyl_pca.PC1.EB.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.EB.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.EB.AIC", "morphrate.geomean_scaled.phyl_pca.PC1.delta.delta", "morphrate.geomean_scaled.phyl_pca.PC1.delta.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.delta.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.delta.AIC", "morphrate.geomean_scaled.phyl_pca.PC13.BAMM", "avg_overlaps.rangesize_scaled", "avg_overlaps.euclidean_scaled.geomean", "avg_overlaps.euclidean_scaled.phyl_pca", "avg_overlaps.euclidean_scaled.geomean_scaled.phyl_pca"), return_format="matrix", quiet=TRUE) {
  # build a list of vectors to store the various metrics (list above is NOT complete)
  # loop over subclades.treedata, which is a list of treedata-like objects, each with a phy, and a bunch of different data
  # calculate metrics for each subclade, and store them in the vectors
  # reformat the list of vectors if necessary (e.g. to matrix)
  # return the reformatted subclade metrics
  require(geiger)
  require(laser)
  
  for (metric in metrics) {  # set up vectors to store subclade data in
    assign(metric, value=numeric())
    if (!quiet) print(metric)
  }
  
  for (i in names(subclades.treedata)) {  # loop over subclades, calculating metrics
    
    cat("\nStarting clade ", i, ", ", which(names(subclades.treedata)==i), " of ",  length(subclades.treedata), " total subclades.\n\n", sep="")
    
    ## diversification and tree-shape stuff
    cat("Starting diversification analyses.\n")
    
    # calculate total number of taxa
    if ("ntaxa" %in% metrics) ntaxa[[i]] <- length(subclades.treedata[[i]]$phy$tip.label) + subclades.treedata[[i]]$missing_count
    
    if ("ntaxa.on_morph_tree" %in% metrics) ntaxa.on_morph_tree[[i]] <- length(subclades.treedata[[i]]$phy$tip.label)
    
    # calculate total diversification
    if ("total_div" %in% metrics) total_div[[i]] <- log(length(subclades.treedata[[i]]$phy$tip.label) + subclades.treedata[[i]]$missing_count)
    
    # calculate clade age
    if ("crown_age" %in% metrics) crown_age[[i]] <- max(node.depth.edgelength(subclades.treedata[[i]]$phy))
    
    # calculate Magallon-Sanderson diversification rates
    if (length(intersect(c("divrate.ms.e10","divrate.ms.e50","divrate.ms.e90"), metrics)) > 0) cat("Calculating Magallon-Sanderson diversification rates.\n")
    if ("divrate.ms.e10" %in% metrics) divrate.ms.e10[[i]] <- geiger::bd.ms(phy=subclades.treedata[[i]]$phy.full_tree, missing=subclades.treedata[[i]]$missing_count.full_tree, crown=TRUE, epsilon=0.10) # calculate Magallon-Sanderson diversification rate with extinction fraction 0.10
    if ("divrate.ms.e50" %in% metrics) divrate.ms.e50[[i]] <- geiger::bd.ms(phy=subclades.treedata[[i]]$phy.full_tree, missing=subclades.treedata[[i]]$missing_count.full_tree, crown=TRUE, epsilon=0.50) # calculate Magallon-Sanderson diversification rate with extinction fraction 0.50
    if ("divrate.ms.e90" %in% metrics) divrate.ms.e90[[i]] <- geiger::bd.ms(phy=subclades.treedata[[i]]$phy.full_tree, missing=subclades.treedata[[i]]$missing_count.full_tree, crown=TRUE, epsilon=0.90) # calculate Magallon-Sanderson diversification rate with extinction fraction 0.90
    
    # calculate diversification rate using laser
    if ("divrate.laser" %in% metrics) {
      cat("Calculating laser diversification rate.\n")
      divrate.laser[[i]] <- laser::bd(subclades.treedata[[i]]$phy.full_tree)$r # get diversification rate from laser model fitting
    } 
    
    # fit constant-rate model, and return diversification rate (lambda-mu) and/or AICc
    
    if (length(intersect(c("divrate.ML.constant.rate","divrate.ML.constant.AICc", "divrate.ML.constant.AIC"), metrics)) > 1) {
      cat("Fitting constant-rate diversification model.\n")
      sink("/dev/null")
      divmodel.tmp <- try(bd_ML(branching.times(subclades.treedata[[i]]$phy.full_tree), missnumspec=subclades.treedata[[i]]$missing_count.full_tree, tdmodel=0))  # fit a constant-rate model
      sink()
      if (class(divmodel.tmp) == "try-error") {
        if ("divrate.ML.constant.rate" %in% metrics) divrate.ML.constant.rate[[i]] <- NA
        if ("divrate.ML.constant.AICc" %in% metrics) divrate.ML.constant.AICc[[i]] <- NA
        if ("divrate.ML.constant.AIC" %in% metrics) divrate.ML.constant.AIC[[i]] <- NA
      } else {
        if ("divrate.ML.constant.rate" %in% metrics) divrate.ML.constant.rate[[i]] <- with(divmodel.tmp, lambda0-mu0)  # extract diversification rate (lambda - mu) from the constant-rate model
        if ("divrate.ML.constant.AICc" %in% metrics) divrate.ML.constant.AICc[[i]] <- (-2 * divmodel.tmp$loglik) + (2 * divmodel.tmp$df) + (((2 * divmodel.tmp$df) * (divmodel.tmp$df + 1)) / (subclades.treedata[[i]]$phy$Nnode - divmodel.tmp$df - 1))  # calculate AICc for the constant-rate model
        if ("divrate.ML.constant.AIC" %in% metrics) divrate.ML.constant.AIC[[i]] <- (-2 * divmodel.tmp$loglik) + (2 * divmodel.tmp$df)  # calculate AIC for the constant-rate model
      }
      rm(divmodel.tmp)  # remove the temporary model
    }
    
    # fit time-dependent-rate model, and return diversification rate (lambda-mu), lambda1, mu1, and/or AICc
    if (length(intersect(c("divrate.ML.time_dependent.rate", "divrate.ML.time_dependent.lambda1", "divrate.ML.time_dependent.mu1", "divrate.ML.time_dependent.AICc", "divrate.ML.time_dependent.AIC"), metrics)) > 1) {
      cat("Fitting time-dependent diversification model.\n")
      sink("/dev/null")
      divmodel.tmp <- try(bd_ML(branching.times(subclades.treedata[[i]]$phy.full_tree), missnumspec=subclades.treedata[[i]]$missing_count.full_tree, tdmodel=1, idparsopt=1:3, initparsopt=c(0.1, 0.05, 0.1)))  # fit a time-dependent-rate model
      sink()
      if (class(divmodel.tmp) == "try-error") {
        if ("divrate.ML.time_dependent.rate" %in% metrics) divrate.ML.time_dependent.rate[[i]] <- NA
        if ("divrate.ML.time_dependent.lambda1" %in% metrics) divrate.ML.time_dependent.lambda1[[i]] <- NA
        if ("divrate.ML.time_dependent.mu1" %in% metrics) divrate.ML.time_dependent.mu1[[i]] <- NA
        if ("divrate.ML.time_dependent.AICc" %in% metrics) divrate.ML.time_dependent.AICc[[i]] <- NA
        if ("divrate.ML.time_dependent.AIC" %in% metrics) divrate.ML.time_dependent.AIC[[i]] <- NA
      } else {
        if ("divrate.ML.time_dependent.rate" %in% metrics) divrate.ML.time_dependent.rate[[i]] <- with(divmodel.tmp, lambda0-mu0)  # extract diversification rate (lambda - mu) from the time-dependent-rate model
        if ("divrate.ML.time_dependent.lambda1" %in% metrics) divrate.ML.time_dependent.lambda1[[i]] <- with(divmodel.tmp, lambda1)  # extract diversification rate (lambda - mu) from the time-dependent-rate model
        if ("divrate.ML.time_dependent.mu1" %in% metrics) divrate.ML.time_dependent.mu1[[i]] <- with(divmodel.tmp, mu1)  # extract diversification rate (lambda - mu) from the time-dependent-rate model
        if ("divrate.ML.time_dependent.AICc" %in% metrics) divrate.ML.time_dependent.AICc[[i]] <- (-2 * divmodel.tmp$loglik) + (2 * divmodel.tmp$df) + (((2 * divmodel.tmp$df) * (divmodel.tmp$df + 1)) / (subclades.treedata[[i]]$phy$Nnode - divmodel.tmp$df - 1))  # calculate AICc for the time-dependent-rate model
        if ("divrate.ML.time_dependent.AIC" %in% metrics) divrate.ML.time_dependent.AIC[[i]] <- (-2 * divmodel.tmp$loglik) + (2 * divmodel.tmp$df)  # calculate AIC for the time-dependent-rate model
      }
      rm(divmodel.tmp)  # remove the temporary model
    }
    
    # fit diversity-dependent-rate model, and return diversification rate (lambda-mu), K, and/or AICc
    if (length(intersect(c("divrate.ML.diversity_dependent.rate", "divrate.ML.diversity_dependent.K", "divrate.ML.diversity_dependent.AICc", "divrate.ML.diversity_dependent.AIC"), metrics)) > 1) {
      cat("Fitting diversity-dependent diversification model.\n")
      sink("/dev/null")
      divmodel.tmp <- try(dd_ML(branching.times(subclades.treedata[[i]]$phy.full_tree), missnumspec=subclades.treedata[[i]]$missing_count.full_tree, ddmodel=1))  # fit a diversity-dependent-rate model, with exponential dependence in speciation rate
      sink()
      if (class(divmodel.tmp) == "try-error") {
        if ("divrate.ML.diversity_dependent.rate" %in% metrics) divrate.ML.diversity_dependent.rate[[i]] <- NA
        if ("divrate.ML.diversity_dependent.K" %in% metrics) divrate.ML.diversity_dependent.K[[i]] <- NA
        if ("divrate.ML.diversity_dependent.AICc" %in% metrics) divrate.ML.diversity_dependent.AICc[[i]] <- NA
        if ("divrate.ML.diversity_dependent.AIC" %in% metrics) divrate.ML.diversity_dependent.AIC[[i]] <- NA
      } else {
        if ("divrate.ML.diversity_dependent.rate" %in% metrics) divrate.ML.diversity_dependent.rate[[i]] <- with(divmodel.tmp, lambda-mu)  # extract diversification rate (lambda - mu) from the diversity-dependent-rate model
        if ("divrate.ML.diversity_dependent.K" %in% metrics) divrate.ML.diversity_dependent.K[[i]] <- with(divmodel.tmp, K)  # extract diversification rate (lambda - mu) from the diversity-dependent-rate model
        if ("divrate.ML.diversity_dependent.AICc" %in% metrics) divrate.ML.diversity_dependent.AICc[[i]] <- (-2 * divmodel.tmp$loglik) + (2 * divmodel.tmp$df) + (((2 * divmodel.tmp$df) * (divmodel.tmp$df + 1)) / (subclades.treedata[[i]]$phy$Nnode - divmodel.tmp$df - 1))  # calculate AICc for the time-dependent-rate model
        if ("divrate.ML.diversity_dependent.AIC" %in% metrics) divrate.ML.diversity_dependent.AIC[[i]] <- (-2 * divmodel.tmp$loglik) + (2 * divmodel.tmp$df)  # calculate AIC for the time-dependent-rate model
      }
      rm(divmodel.tmp)  # remove the temporary model
    }
    
    # extract average diversification rate from BAMM
    if ("divrate.BAMM" %in% metrics) divrate.BAMM[[i]] <- BAMM_divrates$full_tree[as.character(subclades.treedata[[i]]$node.full_tree)] # get average subclade diversification rate from BAMM
    
    # extract average diversification rate from BAMM
    if ("divrate.BAMM.morph_tree" %in% metrics) divrate.BAMM.morph_tree[[i]] <- BAMM_divrates$morph_tree[i] # get average subclade diversification rate from BAMM
    
    # calculate gamma
    if ("gamma" %in% metrics) gamma[[i]] <- gammaStat(subclades.treedata[[i]]$phy.full_tree)
    
    
    ## morphological evolution stuff; I use fitContinuous because the functions in the mvMORPH package are really slow with more than a few variables
    cat("Starting morphological evolution analyses.\n")
    
    # fit BM model to geomean data, and extract sigma-squared and/or AICc
    if (length(intersect(c("morphrate.geomean.BM.sigsq", "morphrate.geomean.BM.AICc", "morphrate.geomean.BM.AIC"), metrics)) > 0) {
      cat("Fitting BM model to geomean data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean, model="BM"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean.BM.sigsq" %in% metrics) morphrate.geomean.BM.sigsq[[i]] <- NA
        if ("morphrate.geomean.BM.AICc" %in% metrics) morphrate.geomean.BM.AICc[[i]] <- NA
        if ("morphrate.geomean.BM.AIC" %in% metrics) morphrate.geomean.BM.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean.BM.sigsq" %in% metrics) morphrate.geomean.BM.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.geomean.BM.AICc" %in% metrics) morphrate.geomean.BM.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.geomean.BM.AIC" %in% metrics) morphrate.geomean.BM.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit Ornstein-Uhlenbeck (OU) model to geomean data, and and extract sigsq, alpha (the stable attractor parameter) and/or AICc
    if (length(intersect(c("morphrate.geomean.OU.alpha", "morphrate.geomean.OU.sigsq", "morphrate.geomean.OU.AICc", "morphrate.geomean.OU.AIC"), metrics)) > 0) {
      cat("Fitting OU model to geomean data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean, model="OU"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean.OU.alpha" %in% metrics) morphrate.geomean.OU.alpha[[i]] <- NA
        if ("morphrate.geomean.OU.sigsq" %in% metrics) morphrate.geomean.OU.sigsq[[i]] <- NA
        if ("morphrate.geomean.OU.AICc" %in% metrics) morphrate.geomean.OU.AICc[[i]] <- NA
        if ("morphrate.geomean.OU.AIC" %in% metrics) morphrate.geomean.OU.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean.OU.alpha" %in% metrics) morphrate.geomean.OU.alpha[[i]] <- morphmodel.tmp$opt$alpha
        if ("morphrate.geomean.OU.sigsq" %in% metrics) morphrate.geomean.OU.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.geomean.OU.AICc" %in% metrics) morphrate.geomean.OU.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.geomean.OU.AIC" %in% metrics) morphrate.geomean.OU.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    
    # fit trend model to geomean data, and extract sigma-squared and/or AICc
    if (length(intersect(c("morphrate.geomean.trend.slope", "morphrate.geomean.trend.sigsq", "morphrate.geomean.trend.AICc", "morphrate.geomean.trend.AIC"), metrics)) > 0) {
      cat("Fitting trend model to geomean data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean, model="trend"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean.trend.slope" %in% metrics) morphrate.geomean.trend.slope[[i]] <- NA
        if ("morphrate.geomean.trend.sigsq" %in% metrics) morphrate.geomean.trend.sigsq[[i]] <- NA
        if ("morphrate.geomean.trend.AICc" %in% metrics) morphrate.geomean.trend.AICc[[i]] <- NA
        if ("morphrate.geomean.trend.AIC" %in% metrics) morphrate.geomean.trend.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean.trend.slope" %in% metrics) morphrate.geomean.trend.slope[[i]] <- morphmodel.tmp$opt$slope
        if ("morphrate.geomean.trend.sigsq" %in% metrics) morphrate.geomean.trend.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.geomean.trend.AICc" %in% metrics) morphrate.geomean.trend.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.geomean.trend.AIC" %in% metrics) morphrate.geomean.trend.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit early burst (EB) model to geomean data, and extract alpha and/or AICc
    if (length(intersect(c("morphrate.geomean.EB.alpha", "morphrate.geomean.EB.sigsq", "morphrate.geomean.EB.AICc", "morphrate.geomean.EB.AIC"), metrics)) > 0) {
      cat("Fitting EB model to geomean data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean, model="EB"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean.EB.alpha" %in% metrics) morphrate.geomean.EB.alpha[[i]] <- NA
        if ("morphrate.geomean.EB.sigsq" %in% metrics) morphrate.geomean.EB.sigsq[[i]] <- NA
        if ("morphrate.geomean.EB.AICc" %in% metrics) morphrate.geomean.EB.AICc[[i]] <- NA
        if ("morphrate.geomean.EB.AIC" %in% metrics) morphrate.geomean.EB.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean.EB.alpha" %in% metrics) morphrate.geomean.EB.alpha[[i]] <- morphmodel.tmp$opt$a
        if ("morphrate.geomean.EB.sigsq" %in% metrics) morphrate.geomean.EB.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.geomean.EB.AICc" %in% metrics) morphrate.geomean.EB.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.geomean.EB.AIC" %in% metrics) morphrate.geomean.EB.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit delta model to geomean data, and extract delta and/or AICc
    if (length(intersect(c("morphrate.geomean.delta.delta", "morphrate.geomean.delta.sigsq", "morphrate.geomean.delta.AICc", "morphrate.geomean.delta.AIC"), metrics)) > 0) {
      cat("Fitting delta model to geomean data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean, model="delta"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean.delta.delta" %in% metrics) morphrate.geomean.delta.delta[[i]] <- NA
        if ("morphrate.geomean.delta.sigsq" %in% metrics) morphrate.geomean.delta.sigsq[[i]] <- NA
        if ("morphrate.geomean.delta.AICc" %in% metrics) morphrate.geomean.delta.AICc[[i]] <- NA
        if ("morphrate.geomean.delta.AIC" %in% metrics) morphrate.geomean.delta.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean.delta.delta" %in% metrics) morphrate.geomean.delta.delta[[i]] <- morphmodel.tmp$opt$delta
        if ("morphrate.geomean.delta.sigsq" %in% metrics) morphrate.geomean.delta.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.geomean.delta.AICc" %in% metrics) morphrate.geomean.delta.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.geomean.delta.AIC" %in% metrics) morphrate.geomean.delta.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # extract average geomean morphological evolution rate from BAMM
    if ("morphrate.geomean.BAMM" %in% metrics) morphrate.geomean.BAMM[[i]] <- BAMM_morphrates[[grep("geomean(?!_scaled)", names(BAMM_morphrates), value=TRUE, perl=TRUE)]][i]
    
    # fit BM model to phyl_pca data, and extract sigma-squared and/or AICc
    if (length(intersect(c("morphrate.phyl_pca.BM.sigsq", "morphrate.phyl_pca.PC1.BM.AICc", "morphrate.phyl_pca.PC1.BM.AIC"), metrics)) > 0) {
      cat("Fitting BM model to phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$phyl_pca, model="BM"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp[["PC1"]]$res[,"convergence"])) > 0)) {
        if ("morphrate.phyl_pca.BM.sigsq" %in% metrics) morphrate.phyl_pca.BM.sigsq[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.BM.AICc" %in% metrics) morphrate.phyl_pca.PC1.BM.AICc[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.BM.AIC" %in% metrics) morphrate.phyl_pca.PC1.BM.AIC[[i]] <- NA
      } else {
        if ("morphrate.phyl_pca.BM.sigsq" %in% metrics) morphrate.phyl_pca.BM.sigsq[[i]] <- sum(sapply(morphmodel.tmp, function(x) x$opt$sigsq))
        if ("morphrate.phyl_pca.PC1.BM.AICc" %in% metrics) morphrate.phyl_pca.PC1.BM.AICc[[i]] <- morphmodel.tmp[["PC1"]]$opt$aicc
        if ("morphrate.phyl_pca.PC1.BM.AIC" %in% metrics) morphrate.phyl_pca.PC1.BM.AIC[[i]] <- morphmodel.tmp[["PC1"]]$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit Ornstein-Uhlenbeck (OU) model to phyl_pca data, and extract sigsq, alpha (the stable attractor parameter) and/or AICc
    if (length(intersect(c("morphrate.phyl_pca.PC1.OU.alpha", "morphrate.phyl_pca.PC1.OU.sigsq", "morphrate.phyl_pca.PC1.OU.AICc"), metrics)) > 0) {
      cat("Fitting OU model to phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$phyl_pca[,"PC1"], model="OU"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.phyl_pca.PC1.OU.alpha" %in% metrics) morphrate.phyl_pca.PC1.OU.alpha[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.OU.sigsq" %in% metrics) morphrate.phyl_pca.PC1.OU.sigsq[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.OU.AICc" %in% metrics) morphrate.phyl_pca.PC1.OU.AICc[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.OU.AIC" %in% metrics) morphrate.phyl_pca.PC1.OU.AIC[[i]] <- NA
      } else {
        if ("morphrate.phyl_pca.PC1.OU.alpha" %in% metrics) morphrate.phyl_pca.PC1.OU.alpha[[i]] <- morphmodel.tmp$opt$alpha
        if ("morphrate.phyl_pca.PC1.OU.sigsq" %in% metrics) morphrate.phyl_pca.PC1.OU.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.phyl_pca.PC1.OU.AICc" %in% metrics) morphrate.phyl_pca.PC1.OU.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.phyl_pca.PC1.OU.AIC" %in% metrics) morphrate.phyl_pca.PC1.OU.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit trend model to phyl_pca data, and extract slope and/or AICc
    if (length(intersect(c("morphrate.phyl_pca.PC1.trend.slope", "morphrate.phyl_pca.PC1.trend.sigsq", "morphrate.phyl_pca.PC1.trend.AICc", "morphrate.phyl_pca.PC1.trend.AIC"), metrics)) > 0) {
      cat("Fitting trend model to phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$phyl_pca[,"PC1"], model="trend"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.phyl_pca.PC1.trend.slope" %in% metrics) morphrate.phyl_pca.PC1.trend.slope[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.trend.sigsq" %in% metrics) morphrate.phyl_pca.PC1.trend.sigsq[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.trend.AICc" %in% metrics) morphrate.phyl_pca.PC1.trend.AICc[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.trend.AIC" %in% metrics) morphrate.phyl_pca.PC1.trend.AIC[[i]] <- NA
      } else {
        if ("morphrate.phyl_pca.PC1.trend.slope" %in% metrics) morphrate.phyl_pca.PC1.trend.slope[[i]] <- morphmodel.tmp$opt$slope
        if ("morphrate.phyl_pca.PC1.trend.sigsq" %in% metrics) morphrate.phyl_pca.PC1.trend.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.phyl_pca.PC1.trend.AICc" %in% metrics) morphrate.phyl_pca.PC1.trend.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.phyl_pca.PC1.trend.AIC" %in% metrics) morphrate.phyl_pca.PC1.trend.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit early burst (EB) model to phyl_pca data, and extract alpha (the rate decline parameter) and/or AICc
    if (length(intersect(c("morphrate.phyl_pca.PC1.EB.alpha", "morphrate.phyl_pca.PC1.EB.sigsq", "morphrate.phyl_pca.PC1.EB.AICc", "morphrate.phyl_pca.PC1.EB.AIC"), metrics)) > 0) {
      cat("Fitting EB model to phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$phyl_pca[,"PC1"], model="EB"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.phyl_pca.PC1.EB.alpha" %in% metrics) morphrate.phyl_pca.PC1.EB.alpha[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.EB.sigsq" %in% metrics) morphrate.phyl_pca.PC1.EB.sigsq[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.EB.AICc" %in% metrics) morphrate.phyl_pca.PC1.EB.AICc[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.EB.AIC" %in% metrics) morphrate.phyl_pca.PC1.EB.AIC[[i]] <- NA
      } else {
        if ("morphrate.phyl_pca.PC1.EB.alpha" %in% metrics) morphrate.phyl_pca.PC1.EB.alpha[[i]] <- morphmodel.tmp$opt$a
        if ("morphrate.phyl_pca.PC1.EB.sigsq" %in% metrics) morphrate.phyl_pca.PC1.EB.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.phyl_pca.PC1.EB.AICc" %in% metrics) morphrate.phyl_pca.PC1.EB.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.phyl_pca.PC1.EB.AIC" %in% metrics) morphrate.phyl_pca.PC1.EB.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit delta model to phyl_pca data, and extract alpha (the rate decline parameter) and/or AICc
    if (length(intersect(c("morphrate.phyl_pca.PC1.delta.delta", "morphrate.phyl_pca.PC1.delta.sigsq", "morphrate.phyl_pca.PC1.delta.AICc", "morphrate.phyl_pca.PC1.delta.AIC"), metrics)) > 0) {
      cat("Fitting delta model to phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$phyl_pca[,"PC1"], model="delta"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.phyl_pca.PC1.delta.delta" %in% metrics) morphrate.phyl_pca.PC1.delta.delta[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.delta.sigsq" %in% metrics) morphrate.phyl_pca.PC1.delta.sigsq[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.delta.AICc" %in% metrics) morphrate.phyl_pca.PC1.delta.AICc[[i]] <- NA
        if ("morphrate.phyl_pca.PC1.delta.AIC" %in% metrics) morphrate.phyl_pca.PC1.delta.AIC[[i]] <- NA
      } else {
        if ("morphrate.phyl_pca.PC1.delta.delta" %in% metrics) morphrate.phyl_pca.PC1.delta.delta[[i]] <- morphmodel.tmp$opt$delta
        if ("morphrate.phyl_pca.PC1.delta.sigsq" %in% metrics) morphrate.phyl_pca.PC1.delta.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.phyl_pca.PC1.delta.AICc" %in% metrics) morphrate.phyl_pca.PC1.delta.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.phyl_pca.PC1.delta.AIC" %in% metrics) morphrate.phyl_pca.PC1.delta.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # extract average phyl_pca morphological evolution rate from BAMM
    if ("morphrate.phyl_pca.PC13.BAMM" %in% metrics) morphrate.phyl_pca.PC13.BAMM[[i]] <- BAMM_morphrates[[grep("(?<!scaled_)phyl_pca_PC1", names(BAMM_morphrates), value=TRUE, perl=TRUE)]][[i]]
    
    
    # fit BM model to geomean_scaled.phyl_pca data, and extract sigma-squared and/or AICc
    if (length(intersect(c("morphrate.geomean_scaled.phyl_pca.BM.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.BM.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.BM.AIC"), metrics)) > 0) {
      cat("Fitting BM model to geomean-scaled phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean_scaled.phyl_pca, model="BM"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp[["PC1"]]$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean_scaled.phyl_pca.BM.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.BM.sigsq[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.BM.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.BM.AICc[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.BM.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.BM.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean_scaled.phyl_pca.BM.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.BM.sigsq[[i]] <- sum(sapply(morphmodel.tmp, function(x) x$opt$sigsq))
        if ("morphrate.geomean_scaled.phyl_pca.PC1.BM.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.BM.AICc[[i]] <- morphmodel.tmp[["PC1"]]$opt$aicc
        if ("morphrate.geomean_scaled.phyl_pca.PC1.BM.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.BM.AIC[[i]] <- morphmodel.tmp[["PC1"]]$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit Ornstein-Uhlenbeck (OU) model to geomean_scaled.phyl_pca data, and extract sigsq, alpha (the stable attractor parameter) and/or AICc
    if (length(intersect(c("morphrate.geomean_scaled.phyl_pca.PC1.OU.alpha", "morphrate.geomean_scaled.phyl_pca.PC1.OU.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.OU.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.OU.AIC"), metrics)) > 0) {
      cat("Fitting OU model to geomean-scaled phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean_scaled.phyl_pca[,"PC1"], model="OU"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean_scaled.phyl_pca.PC1.OU.alpha" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.OU.alpha[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.OU.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.OU.sigsq[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.OU.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.OU.AICc[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.OU.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.OU.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean_scaled.phyl_pca.PC1.OU.alpha" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.OU.alpha[[i]] <- morphmodel.tmp$opt$alpha
        if ("morphrate.geomean_scaled.phyl_pca.PC1.OU.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.OU.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.geomean_scaled.phyl_pca.PC1.OU.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.OU.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.geomean_scaled.phyl_pca.PC1.OU.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.OU.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit trend model to geomean_scaled.phyl_pca data, and extract slope and/or AICc
    if (length(intersect(c("morphrate.geomean_scaled.phyl_pca.PC1.trend.slope", "morphrate.geomean_scaled.phyl_pca.PC1.trend.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.trend.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.trend.AIC"), metrics)) > 0) {
      cat("Fitting trend model to geomean-scaled phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean_scaled.phyl_pca[,"PC1"], model="trend"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean_scaled.phyl_pca.PC1.trend.slope" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.trend.slope[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.trend.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.trend.sigsq[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.trend.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.trend.AICc[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.trend.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.trend.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean_scaled.phyl_pca.PC1.trend.slope" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.trend.slope[[i]] <- morphmodel.tmp$opt$slope
        if ("morphrate.geomean_scaled.phyl_pca.PC1.trend.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.trend.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.geomean_scaled.phyl_pca.PC1.trend.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.trend.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.geomean_scaled.phyl_pca.PC1.trend.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.trend.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit early burst (EB) model to geomean_scaled.phyl_pca data, and extract alpha (the rate decline parameter) and/or AICc
    if (length(intersect(c("morphrate.geomean_scaled.phyl_pca.PC1.EB.alpha", "morphrate.geomean_scaled.phyl_pca.PC1.EB.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.EB.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.EB.AIC"), metrics)) > 0) {
      cat("Fitting EB model to geomean-scaled phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean_scaled.phyl_pca[,"PC1"], model="EB"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean_scaled.phyl_pca.PC1.EB.alpha" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.EB.alpha[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.EB.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.EB.sigsq[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.EB.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.EB.AICc[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.EB.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.EB.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean_scaled.phyl_pca.PC1.EB.alpha" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.EB.alpha[[i]] <- morphmodel.tmp$opt$a
        if ("morphrate.geomean_scaled.phyl_pca.PC1.EB.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.EB.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.geomean_scaled.phyl_pca.PC1.EB.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.EB.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.geomean_scaled.phyl_pca.PC1.EB.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.EB.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # fit delta model to geomean_scaled.phyl_pca data, and extract alpha (the rate decline parameter) and/or AICc
    if (length(intersect(c("morphrate.geomean_scaled.phyl_pca.PC1.delta.delta", "morphrate.geomean_scaled.phyl_pca.PC1.delta.sigsq", "morphrate.geomean_scaled.phyl_pca.PC1.delta.AICc", "morphrate.geomean_scaled.phyl_pca.PC1.delta.AIC"), metrics)) > 0) {
      cat("Fitting delta model to geomean-scaled phyl_pca data.\n")
      morphmodel.tmp <- try(fitContinuous(phy=subclades.treedata[[i]]$phy, dat=subclades.treedata[[i]]$geomean_scaled.phyl_pca[,"PC1"], model="delta"))
      if ((class(morphmodel.tmp) == "try-error") | (min(abs(morphmodel.tmp$res[,"convergence"])) > 0)) {
        if ("morphrate.geomean_scaled.phyl_pca.PC1.delta.delta" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.delta.delta[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.delta.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.delta.sigsq[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.delta.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.delta.AICc[[i]] <- NA
        if ("morphrate.geomean_scaled.phyl_pca.PC1.delta.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.delta.AIC[[i]] <- NA
      } else {
        if ("morphrate.geomean_scaled.phyl_pca.PC1.delta.delta" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.delta.delta[[i]] <- morphmodel.tmp$opt$delta
        if ("morphrate.geomean_scaled.phyl_pca.PC1.delta.sigsq" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.delta.sigsq[[i]] <- morphmodel.tmp$opt$sigsq
        if ("morphrate.geomean_scaled.phyl_pca.PC1.delta.AICc" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.delta.AICc[[i]] <- morphmodel.tmp$opt$aicc
        if ("morphrate.geomean_scaled.phyl_pca.PC1.delta.AIC" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC1.delta.AIC[[i]] <- morphmodel.tmp$opt$aic
      }
      rm(morphmodel.tmp)
    }
    
    # extract average geomean_scaled.phyl_pca morphological evolution rate from BAMM
    if ("morphrate.geomean_scaled.phyl_pca.PC13.BAMM" %in% metrics) morphrate.geomean_scaled.phyl_pca.PC13.BAMM[[i]] <- BAMM_morphrates[[grep("geomean_scaled_phyl_pca_PC1", names(BAMM_morphrates), value=TRUE)]][[i]]
    
    
    ## overlap metrics
    cat("Starting overlap metrics.\n")
    
    # calculate average of summed overlaps scaled by focal taxon range
    if ("avg_overlaps.rangesize_scaled" %in% metrics) avg_overlaps.rangesize_scaled[[i]] <- mean(subclades.treedata[[i]]$overlaps.scaled)
    
    if ("avg_overlaps.euclidean_scaled.geomean" %in% metrics) avg_overlaps.euclidean_scaled.geomean[[i]] <- mean(subclades.treedata[[i]]$overlaps.euclidean_scaled$geomean)
    
    if ("avg_overlaps.euclidean_scaled.phyl_pca" %in% metrics) avg_overlaps.euclidean_scaled.phyl_pca[[i]] <- mean(subclades.treedata[[i]]$overlaps.euclidean_scaled$phyl_pca)
    
    if ("avg_overlaps.euclidean_scaled.geomean_scaled.phyl_pca" %in% metrics) avg_overlaps.euclidean_scaled.geomean_scaled.phyl_pca[[i]] <- mean(subclades.treedata[[i]]$overlaps.euclidean_scaled$geomean_scaled.phyl_pca)

  }
  
  if (return_format == "matrix") {
    subclade_data <- matrix(nrow=length(subclades.treedata), ncol=length(metrics), dimnames=list(names(subclades.treedata), metrics))
    for (i in 1:length(metrics)) {
      if (!quiet) cat(metrics[i], ": ", get(metrics[i]), "\n", sep="")
      subclade_data[,i] <- get(metrics[i])
    }
  }
  
  if (return_format == "list") {
    subclade_data <- list()
    for (metric in metrics) subclade_data[[metric]] <- get(metric)
  }
  
  return(subclade_data)
}

## loop over the list of subclades, calculating metrics for each
picidae.subclade.data <- list()
for (i in names(picidae.morph.log.fully_reduced.subclades.treedata)) {  # loop over individual inclusion
  picidae.subclade.data[[i]] <- calc.subclade.metrics(subclades.treedata=picidae.morph.log.fully_reduced.subclades.treedata[[i]], BAMM_divrates=picidae.BAMM.divrates_by_node, BAMM_morphrates=picidae.BAMM.morphrates_by_node)
}
rm(i)

## calculate delta_aicc for all models vs. basic models (e.g time-dependent and diversity-dependent vs. constant-rate diversification, OU and trend and EB vs. BM for all morph variables)
for (i in names(picidae.subclade.data)) {  # loop over individual inclusion
  for (m in grep("divrate(?!.ML.constant)[a-zA-Z0-9._]+AICc", colnames(picidae.subclade.data[[i]]), value=TRUE, perl=TRUE)) {
    picidae.subclade.data[[i]] <- cbind(picidae.subclade.data[[i]], picidae.subclade.data[[i]][,m] - picidae.subclade.data[[i]][,"divrate.ML.constant.AICc"])
    colnames(picidae.subclade.data[[i]])[ncol(picidae.subclade.data[[i]])] <- sub("AICc", "delta_AICc", m)
  }
  for (m in grep("morphrate[a-zA-Z0-9._]+(?<!BM.)AICc", colnames(picidae.subclade.data[[i]]), value=TRUE, perl=TRUE)) {
    picidae.subclade.data[[i]] <- cbind(picidae.subclade.data[[i]], picidae.subclade.data[[i]][,m] - picidae.subclade.data[[i]][,sub("[a-zA-Z0-9]+(?=.AICc)", "BM", m, perl=TRUE)])
    colnames(picidae.subclade.data[[i]])[ncol(picidae.subclade.data[[i]])] <- sub("AICc", "delta_AICc", m)
  }
}
rm(i,m)


### generate subclade combinations for my subclade regressions

## generate subclade combinations using the random method (100 iterations), and one set using the sequential selection method; for each, set a minimum of 5 clades and 6 taxa per clade
picidae.subclade.combinations.6sp.random <- list()
picidae.subclade.combinations.6sp.sequential <- list()
for (i in names(picidae.morph.log.fully_reduced.subclades.treedata)) {
  picidae.subclade.combinations.6sp.random[[i]] <- subcladeCombinations.random(phylist=picidae.morph.log.fully_reduced.subclades.treedata[[i]], ncombs=100, min.clades=5, min.taxa=6)
  picidae.subclade.combinations.6sp.sequential[[i]] <- subcladeCombinations.sequential(phy=picidae.morph.log.fully_reduced.treedata[[i]]$phy, min.clades=5, min.taxa=6)
}
rm(i)


### generate backbone trees, subclade treedata objects, and fit models to subclade combinations

## general process
## loop over subclade combinations
# generate backbone tree for each subclade combination
# generate treedata object for each subclade combination (using the backbone tree and the subclade data)
# fit models to subclade data using pgls, with the treedata$phy and treedata$data
# save off the relevant bits from the models (slopes and intercepts, R^2)
# examine distributions of r_squared, p_values, slopes, etc. 


## the function fit.subcladeModels.bivariate() fits regression models to subclade data, iterating over various combinations of subclades
# as input, it takes phy (the phylogenetic tree of taxa, as phylo object), subclade.combinations (a list of subclade combinations, each containing a vector of node numbers in phy), subclade.data (a matrix of subclade metrics, with rownames as the subclade numbers (the node numbers)), models (an optional character vector containing models to test, formatted as "var1_vs_var2"), models_filename (the name of an optional text file with models to test, with one model on each line, formatted as "var1_vs_var2"), return_format (format to return results in; can be "matrix" or "list"), model_fitting (method for fitting models; either "pgls" or "lm"), quiet.subclade.combinations (boolean to output to console when starting the next subclade combination), quiet.models (boolean to output to console when starting fitting the next model)
# it returns a matrix or list, with the columns of the matrix or the elements of the list containig the parameter values from the specified models fit to each subclade combination 
fit.subcladeModels.bivariate <- function(phy, subclade.combinations, subclade.data, models=NULL, models_filename="Picidae_subclade_models_bivariate.txt", return_format="matrix", model_fitting="pgls", quiet.subclade.combinations=TRUE, quiet.models=TRUE) {
  # if models not provided as argument, read them from file
  if (is.null(models)) {
    models <- read.table(file=models_filename, header=F, stringsAsFactors=F)[,1]
  }
  
  # create an empty set of vectors for storing parameter values from each model
  for (model in models) {
    for (measure in c("r_squared","p_value","slope","intercept")) {
      assign(paste(model, measure, sep="."), value=numeric())
    }
  }
  
  for (i in 1:length(subclade.combinations)) {  # loop over subclade combinations
    
    # generate backbone tree and treedata object
    if (!quiet.subclade.combinations) cat("\nStarting model fitting for combination ", i, " of ", length(subclade.combinations), ".\n", sep="")
    subclade.treedata.tmp <- getTreedata.subclades(phy=phy, subclade.combination=subclade.combinations[[i]], subclade.data=subclade.data)
    
    for (model in models) {  # loop over models
      if (!quiet.models) cat("Starting model ", model, "\n", sep="")
      model_split <- strsplit(model, "_vs_")[[1]]  # split the model name into the two component variables
      y_var <- model_split[1]  # extract variable names
      x_var <- model_split[2]  # extract variable names
      
      
      if (model_fitting=="pgls") {
        # fit model using pgls
        model.tmp <- try(gls(data = data.frame(y = subclade.treedata.tmp$data[,y_var], x = subclade.treedata.tmp$data[,x_var]), model = y ~ x, na.action=na.exclude, correlation=corPagel(value=1, phy=subclade.treedata.tmp$phy), method="REML"))  # model with correlation structure based on tree, with lambda estimated
        if (class(model.tmp)=="try-error") {
          model.tmp <- try(gls(data = data.frame(y = subclade.treedata.tmp$data[,y_var], x = subclade.treedata.tmp$data[,x_var]), model = y ~ x, na.action=na.exclude, correlation=corBrownian(value=1, phy=subclade.treedata.tmp$phy), method="REML"))  # model with correlation structure based on tree, assuming Brownian Motion (if lambda estimation fails)
        }
        
        # if model still fails, set parameter values to NA
        if (class(model.tmp)=="try-error") {  
          assign(paste(model, "r_squared", sep="."), value=c(get(paste(model, "r_squared", sep=".")), NA))
          assign(paste(model, "p_value", sep="."), value=c(get(paste(model, "p_value", sep=".")), NA))
        } else {
          assign(paste(model, "r_squared", sep="."), value=c(get(paste(model, "r_squared", sep=".")), cor(subclade.treedata.tmp$data[,y_var], model.tmp$fitted)^2))
          assign(paste(model, "p_value", sep="."), value=c(get(paste(model, "p_value", sep=".")), summary(model.tmp)$tTable["x","p-value"]))
        }
        
        
      } else if (model_fitting=="lm") {
        # fit model using lm
        model.tmp <- try(lm(subclade.treedata.tmp$data[,y_var] ~ subclade.treedata.tmp$data[,x_var]))
        if (class(model.tmp)=="try-error") {
          assign(paste(model, "r_squared", sep="."), value=c(get(paste(model, "r_squared", sep=".")), NA))
          assign(paste(model, "p_value", sep="."), value=c(get(paste(model, "p_value", sep=".")), NA))
        } else {
          assign(paste(model, "r_squared", sep="."), value=c(get(paste(model, "r_squared", sep=".")), summary(model.tmp)$r.squared))
          assign(paste(model, "p_value", sep="."), value=c(get(paste(model, "p_value", sep=".")), summary(model.tmp.lm)$coefficients["x","Pr(>|t|)"]))
        }
      }
      
      # if model fitting fails, set parameter values to NA
      if (class(model.tmp)=="try-error") {
        assign(paste(model, "slope", sep="."), value=c(get(paste(model, "slope", sep=".")), NA))
        assign(paste(model, "intercept", sep="."), value=c(get(paste(model, "intercept", sep=".")), NA))
      } else {
        assign(paste(model, "slope", sep="."), value=c(get(paste(model, "slope", sep=".")), model.tmp$coefficients["x"]))
        assign(paste(model, "intercept", sep="."), value=c(get(paste(model, "intercept", sep=".")), model.tmp$coefficients["(Intercept)"]))
      }
    }
  }
  
  if (return_format == "matrix") {
    # generate a matrix with the values of r_squared, slope, and intercept for all models (in columns), and the subclade combinations as rows
    subclade_combination_model_results <- matrix(nrow=length(subclade.combinations), ncol=length(models)*4, dimnames=list(as.character(1:length(subclade.combinations)), as.vector(sapply(models, function(y) paste(y, c("r_squared","p_value","slope","intercept"), sep=".")))))
    for (i in colnames(subclade_combination_model_results)) {
      subclade_combination_model_results[,i] <- get(i)
    }
  } else if (return_format == "list") {
    # generate a list of vectors with the values of r_squared, slope, and intercept for all models (in separate list items), and the subclade combinations as elements of vectors
    subclade_combination_model_results <- list()
    for (i in as.vector(sapply(models, function(y) paste(y, c("r_squared","p_value","slope","intercept", sep="."))))) subclade_combination_model_results[[i]] <- get(i)
  } else if (return_format == "array") {
    # generate an array with all models in one dimension, the values of r_squared, slope, and intercept in another dimension, and all subclade combinations in another dimension
    subclade_combination_model_results <- array(dim=c(length(models), 4, length(subclade.combinations)), dimnames=list(models, c("r_squared","p_value","slope","intercept"), as.character(1:length(subclade.combinations))))
    for (i in models) {
      for (j in c("r_squared","p_value","slope","intercept")) {
        subclade_combination_model_results[i,j,] <- get(paste(i,j, sep="."))
      }
    }
  }
  
  return(subclade_combination_model_results)
}

picidae.subclade.combinations.6sp.random.model_params <- list()
picidae.subclade.combinations.6sp.sequential.model_params <- list()
for (i in names(picidae.subclade.combinations.6sp.random)) {  # loop over individual inclusions
  cat("\nStarting subclade model fitting for random combinations of", i, "\n", sep=" ")
  picidae.subclade.combinations.6sp.random.model_params[[i]] <- fit.subcladeModels.bivariate(phy=picidae.morph.log.fully_reduced.treedata[[i]]$phy, subclade.combinations=picidae.subclade.combinations.6sp.random[[i]], subclade.data=picidae.subclade.data[[i]], quiet.subclade.combinations=FALSE)
  cat("\nStarting subclade model fitting for sequential combinations of", i, "\n", sep=" ")
  picidae.subclade.combinations.6sp.sequential.model_params[[i]] <- fit.subcladeModels.bivariate(phy=picidae.morph.log.fully_reduced.treedata[[i]]$phy, subclade.combinations=picidae.subclade.combinations.6sp.sequential[[i]], subclade.data=picidae.subclade.data[[i]], quiet.subclade.combinations=FALSE)
}
rm(i)


### summarize results

## histograms of important models (across the different random subclade combinations)
for (i in names(picidae.subclade.combinations.6sp.random.model_params)) {
  for (j in sub(".r_squared", "", grep("[a-zA-z0-9_.]+.r_squared", colnames(picidae.subclade.combinations.6sp.random.model_params[[i]]), perl=TRUE, value=TRUE))) {
    pdf(file=paste("picidae_6sp_random", i , j, "histogram.pdf", sep="_"), height=10, width=10)
    par(mfrow=c(2,2))
    for (k in c("r_squared","p_value","slope","intercept")) {
      if (k == "p_value") {
        hist(picidae.subclade.combinations.6sp.random.model_params[[i]][,paste(j, k, sep=".")], xlab=NULL, main=k, col="gray", breaks=seq(0,1, by=0.05))
        abline(v=0.05, lwd=2, col="red")
      } else {
        hist(picidae.subclade.combinations.6sp.random.model_params[[i]][,paste(j, k, sep=".")], xlab=NULL, main=k, col="gray", breaks=10)
      }
      if (k == "slope") abline(v=0, lwd=2, col="red")
    }
    dev.off()
  }
}
rm(i,j,k)

## histograms and line plots of important models (across the different sequential subclade combinations)
for (i in names(picidae.subclade.combinations.6sp.sequential.model_params)) {
  for (j in sub(".r_squared", "", grep("[a-zA-z0-9_.]+.r_squared", colnames(picidae.subclade.combinations.6sp.sequential.model_params[[i]]), perl=TRUE, value=TRUE))) {
    pdf(file=paste("picidae_6sp_sequential.", j, ".histogram_lineplot.pdf", sep=""), height=10, width=20)
    par(mfcol=c(2,4))
    for (k in c("r_squared","p_value","slope","intercept")) {
      if (k == "p_value") {
        hist(picidae.subclade.combinations.6sp.sequential.model_params[[i]][,paste(j, k, sep=".")], xlab=NULL, main=k, col="gray", breaks=seq(0,1, by=0.05))
        abline(v=0.05, lwd=2, col="red")
      } else {
        hist(picidae.subclade.combinations.6sp.sequential.model_params[[i]][,paste(j, k, sep=".")], xlab=NULL, main=k, col="gray", breaks=10)
      }
      if (k == "slope") abline(v=0, lwd=2, col="red")
      plot(picidae.subclade.combinations.6sp.sequential.model_params[[i]][,paste(j, k, sep=".")] ~ rownames(picidae.subclade.combinations.6sp.sequential.model_params[[i]]), ylab=k, xlab="Tree slice (starting at root)", main=k, type="l")
      if (k == "p_value") abline(h=0.05, lty="dashed")
      if (k == "slope") abline(h=0, lty="dashed")
    }
    dev.off()
  }
}
rm(i,j,k)

## capture median of parameter values
picidae.subclade.combinations.6sp.random.model_params.median <- list()
for (i in names(picidae.subclade.combinations.6sp.random.model_params)) {
  medians.tmp <- apply(picidae.subclade.combinations.6sp.random.model_params[[i]], MARGIN=2, median, na.rm=TRUE)
  picidae.subclade.combinations.6sp.random.model_params.median[[i]] <- matrix(nrow=ncol(picidae.subclade.combinations.6sp.random.model_params[[i]])/4, ncol=4, dimnames=list(sub(".r_squared", "", grep("[a-zA-z0-9_.]+.r_squared", colnames(picidae.subclade.combinations.6sp.sequential.model_params[[i]]), perl=TRUE, value=TRUE)), c("r_squared", "p_value", "slope", "intercept")))
  picidae.subclade.combinations.6sp.random.model_params.median[[i]][,1] <- medians.tmp[grep("r_squared", names(medians.tmp), perl=TRUE, value=TRUE)]
  picidae.subclade.combinations.6sp.random.model_params.median[[i]][,2] <- medians.tmp[grep("p_value", names(medians.tmp), perl=TRUE, value=TRUE)]
  picidae.subclade.combinations.6sp.random.model_params.median[[i]][,3] <- medians.tmp[grep("(?<!trend.)slope", names(medians.tmp), perl=TRUE, value=TRUE)]
  picidae.subclade.combinations.6sp.random.model_params.median[[i]][,4] <- medians.tmp[grep("intercept", names(medians.tmp), perl=TRUE, value=TRUE)]
  
  ## output median values to a table
  write.csv(picidae.subclade.combinations.6sp.random.model_params.median[[i]], file=paste("picidae_subclade_combinations.6sp_random", i, "model_params.median.csv", sep="."))
}
rm(i,medians.tmp)

## capture median of parameter values without outliers (Picidae clades 234, 235; Picinae clades 208, 209)
picidae.subclade.combinations.6sp.random.model_params.no_outliers.median <- list()
for (i in names(picidae.subclade.combinations.6sp.random.model_params)) {
  medians.tmp <- apply(picidae.subclade.combinations.6sp.random.model_params[[i]][!sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) length(intersect(x, c("234","235"))) > 0),], MARGIN=2, median, na.rm=TRUE)
  picidae.subclade.combinations.6sp.random.model_params.no_outliers.median[[i]] <- matrix(nrow=ncol(picidae.subclade.combinations.6sp.random.model_params[[i]])/4, ncol=4, dimnames=list(sub(".r_squared", "", grep("[a-zA-z0-9_.]+.r_squared", colnames(picidae.subclade.combinations.6sp.sequential.model_params[[i]]), perl=TRUE, value=TRUE)), c("r_squared", "p_value", "slope", "intercept")))
  picidae.subclade.combinations.6sp.random.model_params.no_outliers.median[[i]][,1] <- medians.tmp[grep("r_squared", names(medians.tmp), perl=TRUE, value=TRUE)]
  picidae.subclade.combinations.6sp.random.model_params.no_outliers.median[[i]][,2] <- medians.tmp[grep("p_value", names(medians.tmp), perl=TRUE, value=TRUE)]
  picidae.subclade.combinations.6sp.random.model_params.no_outliers.median[[i]][,3] <- medians.tmp[grep("(?<!trend.)slope", names(medians.tmp), perl=TRUE, value=TRUE)]
  picidae.subclade.combinations.6sp.random.model_params.no_outliers.median[[i]][,4] <- medians.tmp[grep("intercept", names(medians.tmp), perl=TRUE, value=TRUE)]
}
rm(i,medians.tmp)


### output scatterplots of clade-level metrics for all subclades 

i <- "all_inds"
picidae.subclade.data.6sp.main_variant <- picidae.subclade.data[[i]][picidae.subclade.data[[i]][,"ntaxa.on_morph_tree"] >= 6,] # trim subclade data to only include clades with at least 6 taxa on the tree, as those were the ones used in model fitting

## plots of diversification rate vs. average range overlap and rate of shape evolution
pdf(file="Picidae_diversification_rates_vs_overlap_and_shape_evolution_rate.pdf", width=10, height=5, useDingbats=FALSE)
par(mfrow=c(1,2))

plot(picidae.subclade.data.6sp.main_variant[,"divrate.ML.constant.rate"] ~ picidae.subclade.data.6sp.main_variant[,"avg_overlaps.rangesize_scaled"], xlab="Average Range Overlap", ylab="Diversification Rate", pch=19)
abline(a=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["divrate.ML.constant.rate_vs_avg_overlaps.rangesize_scaled","intercept"], b=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["divrate.ML.constant.rate_vs_avg_overlaps.rangesize_scaled","slope"])
text(x=8.5,y=0.13,labels=bquote(paste("median ", R^2, " = ", .(format(round(picidae.subclade.combinations.6sp.random.model_params.median[[i]]["divrate.ML.constant.rate_vs_avg_overlaps.rangesize_scaled","r_squared"], 2), nsmall=2)), sep="")))

plot(picidae.subclade.data.6sp.main_variant[,"divrate.ML.constant.rate"] ~ picidae.subclade.data.6sp.main_variant[,"morphrate.geomean_scaled.phyl_pca.BM.sigsq"], xlab="Rate of Size-scaled Shape Evolution", ylab="", pch=19)
abline(a=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["divrate.ML.constant.rate_vs_morphrate.geomean_scaled.phyl_pca.BM.sigsq","intercept"], b=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["divrate.ML.constant.rate_vs_morphrate.geomean_scaled.phyl_pca.BM.sigsq","slope"])
text(x=0.0105,y=0.13,labels=bquote(paste("median ", R^2, " = ", .(format(round(picidae.subclade.combinations.6sp.random.model_params.median[[i]]["divrate.ML.constant.rate_vs_morphrate.geomean_scaled.phyl_pca.BM.sigsq","r_squared"], 2), nsmall=2)), sep="")))

par(mfrow=c(1,1))
dev.off()


## plots of the three morphological evolution rates vs. overlaps
pdf(file="Morpholopgical_evolution_rates_vs_overlap.pdf", width=10.5, height=3.5, useDingbats=FALSE)
par(mfrow=c(1,3), mar=c(5,5,4,2)+0.1)

plot(picidae.subclade.data.6sp.main_variant[,"morphrate.geomean.BM.sigsq"] ~ picidae.subclade.data.6sp.main_variant[,"avg_overlaps.rangesize_scaled"], xlab="Average Range Overlap", ylab=expression("Rate " ~ (sigma^2)), pch=19, main="Size Evolution")
abline(a=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["morphrate.geomean.BM.sigsq_vs_avg_overlaps.rangesize_scaled","intercept"], b=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["morphrate.geomean.BM.sigsq_vs_avg_overlaps.rangesize_scaled","slope"])
text(x=8.5,y=0.006,labels=bquote(paste("median ", R^2, " = ", .(format(round(picidae.subclade.combinations.6sp.random.model_params.median[[i]]["morphrate.geomean.BM.sigsq_vs_avg_overlaps.rangesize_scaled","r_squared"], 2), nsmall=2)), sep="")))

plot(picidae.subclade.data.6sp.main_variant[,"morphrate.phyl_pca.BM.sigsq"] ~ picidae.subclade.data.6sp.main_variant[,"avg_overlaps.rangesize_scaled"], xlab="Average Range Overlap", ylab="", main="Overall Morphological Evolution", pch=19)
abline(a=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["morphrate.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled","intercept"], b=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["morphrate.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled","slope"])
text(x=8.5,y=0.07,labels=bquote(paste("median ", R^2, " = ", .(format(round(picidae.subclade.combinations.6sp.random.model_params.median[[i]]["morphrate.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled","r_squared"], 2), nsmall=2)), sep="")))

plot(picidae.subclade.data.6sp.main_variant[,"morphrate.geomean_scaled.phyl_pca.BM.sigsq"] ~ picidae.subclade.data.6sp.main_variant[,"avg_overlaps.rangesize_scaled"], xlab="Average Range Overlap", ylab="", main="Size-scaled Shape Evolution", pch=19)
abline(a=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["morphrate.geomean_scaled.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled","intercept"], b=picidae.subclade.combinations.6sp.random.model_params.median[[i]]["morphrate.geomean_scaled.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled","slope"])
text(x=8.5,y=0.002,labels=bquote(paste("median ", R^2, " = ", .(format(round(picidae.subclade.combinations.6sp.random.model_params.median[[i]]["morphrate.geomean_scaled.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled","r_squared"], 2), nsmall=2)), sep="")))

par(mfrow=c(1,1))
dev.off()


### quantifying the inclusion of subclades in the random combinations and taxa in subclades
mean(sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) length(x)))  # calculate the average number of subclades included in the random combinations
mean(sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) sum(sapply(x, function(y) length(picidae.morph.log.fully_reduced.subclades.treedata[[i]][[y]]$phy$tip.label)))))  # calculate the average number of taxa from the morph tree included in the subclade
mean(sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) sum(sapply(x, function(y) length(picidae.morph.log.fully_reduced.subclades.treedata[[i]][[y]]$phy.full_tree$tip.label)))))  # calculate the average number of taxa from the full tree included in the subclade

## checking median values of variables from subclade combinations with and without the two outlier clades (234 and 235 in picidae analyses, 208 and 209 in picinae analyses)
sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) length(intersect(x, c("234","235"))) > 0) # identify subclade combinations including one of those 2 outlier subclades
apply(picidae.subclade.combinations.6sp.random.model_params[[i]][sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) length(intersect(x, c("234","235"))) > 0),], MARGIN=2, FUN = median) # median values for subclades combinations with the two outlier clades
apply(picidae.subclade.combinations.6sp.random.model_params[[i]][!sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) length(intersect(x, c("234","235"))) > 0),], MARGIN=2, FUN = median) # median values for subclades combinations without the two outlier clades
(apply(picidae.subclade.combinations.6sp.random.model_params[[i]][sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) length(intersect(x, c("234","235"))) > 0),], MARGIN=2, FUN = median) - apply(picidae.subclade.combinations.6sp.random.model_params[[i]][!sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) length(intersect(x, c("234","235"))) > 0),], MARGIN=2, FUN = median)) / (apply(picidae.subclade.combinations.6sp.random.model_params[[i]][sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) length(intersect(x, c("234","235"))) > 0),], MARGIN=2, FUN = max) - apply(picidae.subclade.combinations.6sp.random.model_params[[i]][sapply(picidae.subclade.combinations.6sp.random[[i]], function(x) length(intersect(x, c("234","235"))) > 0),], MARGIN=2, FUN = min))  # quantify the difference between the parameters from the model fits to subclade combinations with and without the outlier clades, as a percentage of the range of values across all subclade combinations


### output table of median AICc and delta-AICc values from fitting diversification models and morphological evolution models by subclade

i <- "all_inds"

## output diversification models: AICc for constant-rate; delta AICc for time-dependent, diversity-dependent
for (m in c("divrate.ML.constant.AICc","divrate.ML.time_dependent.delta_AICc", "divrate.ML.diversity_dependent.delta_AICc")) {
  cat(m, ": ", median(picidae.subclade.data[[i]][picidae.subclade.data[[i]][,"ntaxa.on_morph_tree"]>=6,m]), sep="")
}

## generate table of morph evolution models: AICc for BM; delta AICc for OU, EB, trend
morph_vars <- c("geomean", "phyl_pca.PC1", "geomean_scaled.phyl_pca.PC1")
models <- c("OU", "EB", "trend")
picidae.morph_models.delta_AICc <- matrix(nrow=length(morph_vars), ncol=length(models), dimnames=list(morph_vars, models))
for (m in morph_vars) {
  for (n in models) {
    picidae.morph_models.delta_AICc.table[m,n] <- median(picidae.subclade.data[[i]][picidae.subclade.data[[i]][,"ntaxa.on_morph_tree"]>=6,paste("morphrate.", m, ".", n, ".delta_AICc", sep="")])
  }
}
rm(i,m,n)


### output table of median parameter values (slope, pseudo-R^2, and p-value) from model fits to subclade combinations

i <- "all_inds"

subclade_models_to_output <- c("total_div_vs_crown_age", "divrate.ML.constant.rate_vs_morphrate.geomean.BM.sigsq", "divrate.ML.constant.rate_vs_morphrate.phyl_pca.BM.sigsq", "divrate.ML.constant.rate_vs_morphrate.geomean_scaled.phyl_pca.BM.sigsq", "divrate.ML.constant.rate_vs_avg_overlaps.rangesize_scaled", "divrate.ML.constant.rate_vs_avg_overlaps.euclidean_scaled.geomean", "divrate.ML.constant.rate_vs_avg_overlaps.euclidean_scaled.phyl_pca", "divrate.ML.constant.rate_vs_avg_overlaps.euclidean_scaled.geomean_scaled.phyl_pca", "morphrate.geomean.BM.sigsq_vs_avg_overlaps.rangesize_scaled", "morphrate.geomean.BM.sigsq_vs_avg_overlaps.euclidean_scaled.geomean", "morphrate.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled", "morphrate.phyl_pca.BM.sigsq_vs_avg_overlaps.euclidean_scaled.phyl_pca", "morphrate.geomean_scaled.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled", "morphrate.geomean_scaled.phyl_pca.BM.sigsq_vs_avg_overlaps.euclidean_scaled.geomean_scaled.phyl_pca", "gamma_vs_avg_overlaps.rangesize_scaled", "gamma_vs_avg_overlaps.euclidean_scaled.geomean", "gamma_vs_avg_overlaps.euclidean_scaled.phyl_pca", "gamma_vs_avg_overlaps.euclidean_scaled.geomean_scaled.phyl_pca", "divrate.ML.time_dependent.delta_AICc_vs_avg_overlaps.rangesize_scaled", "divrate.ML.time_dependent.lambda1_vs_avg_overlaps.rangesize_scaled", "divrate.ML.diversity_dependent.delta_AICc_vs_avg_overlaps.rangesize_scaled", "morphrate.geomean.EB.delta_AICc_vs_avg_overlaps.rangesize_scaled", "morphrate.geomean.EB.alpha_vs_avg_overlaps.rangesize_scaled", "morphrate.phyl_pca.PC1.EB.delta_AICc_vs_avg_overlaps.rangesize_scaled", "morphrate.phyl_pca.PC1.EB.alpha_vs_avg_overlaps.rangesize_scaled", "morphrate.geomean_scaled.phyl_pca.PC1.EB.delta_AICc_vs_avg_overlaps.rangesize_scaled", "morphrate.geomean_scaled.phyl_pca.PC1.EB.alpha_vs_avg_overlaps.rangesize_scaled")
params_to_output <- c("r_squared", "p_value", "slope")

## generate table to store median parameter values
picidae.subclade_models.params.table <- matrix(nrow=length(subclade_models_to_output), ncol=length(params_to_output), dimnames=list(subclade_models_to_output,params_to_output))
for (m in subclade_models_to_output) {
  for (n in params_to_output) {
    picidae.subclade_models.params.table[m,n] <- picidae.subclade.combinations.6sp.random.model_params.median[[i]][m,n]
  }
}
rm(m,n)

write.csv(picidae.subclade_models.params.table, file="picidae.subclade_models.params.median.csv")  # output table to file

## generate table to store median parameter values, without outliers
picidae.subclade_models.params.no_outliers.table <- matrix(nrow=length(subclade_models_to_output), ncol=length(params_to_output), dimnames=list(subclade_models_to_output,params_to_output))
for (m in subclade_models_to_output) {
  for (n in params_to_output) {
    picidae.subclade_models.params.no_outliers.table[m,n] <- picidae.subclade.combinations.6sp.random.model_params.no_outliers.median[[i]][m,n]
  }
}
rm(m,n)

write.csv(picidae.subclade_models.params.no_outliers.table, file="picidae.subclade_models.params.no_outliers.median.csv")  # output table to file


### output boxplots of slope, R^2, and p_value for the most important models

## generate vectors of names of models to output
subclade_models_to_output.divrates_morphrates <-c("divrate.ML.constant.rate_vs_morphrate.phyl_pca.BM.sigsq", "divrate.ML.constant.rate_vs_morphrate.geomean.BM.sigsq", "divrate.ML.constant.rate_vs_morphrate.geomean_scaled.phyl_pca.BM.sigsq")
subclade_models_to_output.divrates_overlaps <- c("divrate.ML.constant.rate_vs_avg_overlaps.rangesize_scaled", "divrate.ML.constant.rate_vs_avg_overlaps.euclidean_scaled.phyl_pca", "divrate.ML.constant.rate_vs_avg_overlaps.euclidean_scaled.geomean", "divrate.ML.constant.rate_vs_avg_overlaps.euclidean_scaled.geomean_scaled.phyl_pca")
subclade_models_to_output.morphrates_overlaps <- c("morphrate.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled", "morphrate.geomean.BM.sigsq_vs_avg_overlaps.rangesize_scaled", "morphrate.geomean_scaled.phyl_pca.BM.sigsq_vs_avg_overlaps.rangesize_scaled")
subclade_models_to_output.divmodels_overlaps <- c("divrate.ML.time_dependent.delta_AICc_vs_avg_overlaps.rangesize_scaled", "divrate.ML.diversity_dependent.delta_AICc_vs_avg_overlaps.rangesize_scaled")
subclade_models_to_output.morphmodels_overlaps <- c("morphrate.geomean.EB.delta_AICc_vs_avg_overlaps.rangesize_scaled", "morphrate.phyl_pca.PC1.EB.delta_AICc_vs_avg_overlaps.rangesize_scaled", "morphrate.geomean_scaled.phyl_pca.PC1.EB.delta_AICc_vs_avg_overlaps.rangesize_scaled")

## output boxplots for divrates vs. morphrates
for (m in params_to_output) {
  pdf(file=paste("picidae_subclade_combinations.6sp_random.model_params.boxplots.divrates_morphrates.", m, ".pdf", sep=""), height=6, width=4)
  data.tmp <- numeric()
  names.tmp <- character()
  
  # loop over models to output, storing data and the name of the data
  for (n in subclade_models_to_output.divrates_morphrates) {
    data.tmp <- c(data.tmp, picidae.subclade.combinations.6sp.random.model_params[[i]][,grep(paste(n, m, sep="."), colnames(picidae.subclade.combinations.6sp.sequential.model_params[[i]]), value=TRUE)])  # get data for current model and parameter
    names.tmp <- c(names.tmp, rep(n, times=nrow(picidae.subclade.combinations.6sp.random.model_params[[i]])))
  }
  names.tmp <- factor(x=names.tmp, levels=subclade_models_to_output.divrates_morphrates)  # set the names of the data to a factor, so that they plot in the correct order
  
  boxplot(data.tmp ~ names.tmp, names=c("Overall", "Size", "Shape"), ylim=switch(m, r_squared = c(0,1), p_value = c(0,1), slope = range(data.tmp)))  # create boxplot
  
  # add horizontal lines to boxplots as needed
  switch(m,
         r_squared = NULL,
         p_value = abline(h=0.05, lty=2),
         slope = abline(h=0, lty=2)
  )
  
  dev.off()
}
rm(m,n, data.tmp, names.tmp)


## output boxplots for divrates vs. overlaps
for (m in params_to_output) {
  pdf(file=paste("picidae_subclade_combinations.6sp_random.model_params.boxplots.divrates_overlaps.", m, ".pdf", sep=""), height=6, width=4)
  data.tmp <- numeric()
  names.tmp <- character()
  
  # loop over models to output, storing data and the name of the data
  for (n in subclade_models_to_output.divrates_overlaps) {
    data.tmp <- c(data.tmp, picidae.subclade.combinations.6sp.random.model_params[[i]][,grep(paste(n, m, sep="."), colnames(picidae.subclade.combinations.6sp.sequential.model_params[[i]]), value=TRUE)])
    names.tmp <- c(names.tmp, rep(n, times=nrow(picidae.subclade.combinations.6sp.random.model_params[[i]])))
  }
  names.tmp <- factor(x=names.tmp, levels=subclade_models_to_output.divrates_overlaps)
  
  boxplot(data.tmp ~ names.tmp, names=c("Unscaled", "Overall", "Size", "Shape"), ylim=switch(m, r_squared = c(0,1), p_value = c(0,1), slope = range(data.tmp)))  # create boxplot
  
  # add horizontal lines to boxplots as needed
  switch(m,
         r_squared = NULL,
         p_value = abline(h=0.05, lty=2),
         slope = abline(h=0, lty=2)
  )
  dev.off()
}
rm(m,n, data.tmp, names.tmp)

## output boxplots for morphrates vs. overlaps
for (m in params_to_output) {
  pdf(file=paste("picidae_subclade_combinations.6sp_random.model_params.boxplots.morphrates_overlaps.", m, ".pdf", sep=""), height=6, width=4)
  data.tmp <- numeric()
  names.tmp <- character()
  
  # loop over models to output, storing data and the name of the data
  for (n in subclade_models_to_output.morphrates_overlaps) {
    data.tmp <- c(data.tmp, picidae.subclade.combinations.6sp.random.model_params[[i]][,grep(paste(n, m, sep="."), colnames(picidae.subclade.combinations.6sp.sequential.model_params[[i]]), value=TRUE)])
    names.tmp <- c(names.tmp, rep(n, times=nrow(picidae.subclade.combinations.6sp.random.model_params[[i]])))
  }
  names.tmp <- factor(x=names.tmp, levels=subclade_models_to_output.morphrates_overlaps)
  
  boxplot(data.tmp ~ names.tmp, names=c("Overall", "Size", "Shape"), ylim=switch(m, r_squared = c(0,1), p_value = c(0,1), slope = range(data.tmp)))  # create boxplot
  
  # add horizontal lines to boxplots as needed
  switch(m,
         r_squared = NULL,
         p_value = abline(h=0.05, lty=2),
         slope = abline(h=0, lty=2)
  )
  dev.off()
}
rm(m,n, data.tmp, names.tmp)



###### run the analyses for Picinae ######

### extract subclades from full tree and morph trees

## extract subclades from full trees
picinae.RAxML.all.BEAST_calibrated.with_proxies.subclades <- extractSubclades.all(picinae.RAxML.all.BEAST_calibrated.with_proxies)

# extract subclades from morph trees
picinae.morph.log.fully_reduced.treedata.subclades <- list()
for (i in c("all_inds", "complete_ind_only")) {
  picinae.morph.log.fully_reduced.treedata.subclades[[i]] <- extractSubclades.all(picinae.morph.log.fully_reduced.treedata[[i]]$phy)
}
rm(i)

### to complete analyses for Picinae, use same code as above, replacing picidae with picinae for variable and file names
