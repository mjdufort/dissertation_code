## This is is one of several files containing scripts and functions used in processing and analysis of data for Matthew Dufort's Ph.D. dissertation at the University of Minnesota, titled "Coexistence, Ecomorphology, and Diversification in the Avian Family Picidae (Woodpeckers and Allies)."

## this file contains scripts for preparing data for analysis using BAMM (Bayesian Analysis of Macroevolutionary Mixture), and for processing the results of BAMM runs


### load packages and data

library(ape)
library(BAMMtools)
library(phytools)
library(geiger)
library(coda)

load(file="Picidae_data_for_distribution_morphology_evolution.RData")  # load morphology data (from Morphology_data_processing.R)
load(file="R_trees_for_combined_analyses.RData")  # load phylogenetic trees (from Tree_manipulation.R) 


### list of the main trees and data objects I'm using

## Picidae trees
picidae.RAxML.all.BEAST_calibrated  # full tree of Picidae, BEAST calibrated
picidae.RAxML.all.r8s_calibrated.ladderized  # full tree of Picidae, r8s calibrated
picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy  # morph tree of Picidae, BEAST_calibrated, for the main data set (with or without pygo)

## Picinae trees
picinae.RAxML.all.BEAST_calibrated  # full tree of Picinae, BEAST calibrated
picinae.RAxML.all.r8s_calibrated.ladderized  # full tree of Picinae, r8s calibrated
picinae.morph.log.fully_reduced.treedata[["all_inds"]]$phy  # morph tree of Picinae, BEAST_calibrated, for the main data set (with pygo)

## ladderized versions of the trees
picidae.RAxML.all.BEAST_calibrated.ladderized <- ladderize(picidae.RAxML.all.BEAST_calibrated)  # full tree of Picidae, BEAST calibrated
picidae.RAxML.all.r8s_calibrated.ladderized <- ladderize(picidae.RAxML.all.r8s_calibrated.ladderized)  # full tree of Picidae, r8s calibrated
# Picinae trees
picinae.RAxML.all.BEAST_calibrated.ladderized <- ladderize(picinae.RAxML.all.BEAST_calibrated)  # full tree of Picinae, BEAST calibrated
picinae.RAxML.all.r8s_calibrated.ladderized <- ladderize(picinae.RAxML.all.r8s_calibrated.ladderized)  # full tree of Picinae, r8s calibrated

## convert timescale of r8s trees to millions of years (as it is for the BEAST trees); r8s trees caused issues because edge lengths were in years, which messes up the scaling of the rates
picidae.RAxML.all.r8s_calibrated.ladderized.scaled <- picidae.RAxML.all.r8s_calibrated.ladderized
picidae.RAxML.all.r8s_calibrated.ladderized.scaled$edge.length <- picidae.RAxML.all.r8s_calibrated.ladderized.scaled$edge.length / 1000000  # scaled the r8s tree to Ma (instead of years)
picidae.RAxML.all.r8s_calibrated.ladderized.scaled.ladderized <- ladderize(picidae.RAxML.all.r8s_calibrated.ladderized.scaled)
picinae.RAxML.all.r8s_calibrated.ladderized.scaled <- picinae.RAxML.all.r8s_calibrated.ladderized
picinae.RAxML.all.r8s_calibrated.ladderized.scaled$edge.length <- picinae.RAxML.all.r8s_calibrated.ladderized.scaled$edge.length / 1000000  # scaled the r8s tree to Ma (instead of years)
picinae.RAxML.all.r8s_calibrated.ladderized.scaled.ladderized <- ladderize(picinae.RAxML.all.r8s_calibrated.ladderized.scaled)

# store vector of tree names
treenames <- c("picidae.RAxML.all.BEAST_calibrated", "picidae.RAxML.all.r8s_calibrated.ladderized.scaled", "picidae.morph.log.fully_reduced.treedata", "picinae.RAxML.all.BEAST_calibrated", "picinae.RAxML.all.r8s_calibrated.ladderized.scaled", "picinae.morph.log.fully_reduced.treedata")

## check that trees are valid for BAMM analyses (ultrametric, binary, with no 0-length or negative-length edges)
for (tree in treenames) {
  if (class(get(tree)) == "phylo") {
    tree.tmp <- get(tree)
    cat("\nTree: ", tree, "\n", sep="")
    cat("Ultrametric: ", is.ultrametric(tree.tmp), "\n", sep="")
    cat("Binary: ", is.binary.tree(tree.tmp), "\n", sep="")
    cat("Min edge length: ", min(tree.tmp$edge.length), "\n", sep="")
  } else if (class(get(tree)) == "list") {
    tree.tmp <- get(tree)[["all_inds"]]$phy
    
    cat("\nTree: ", tree, "\n", sep="")
    cat("Ultrametric: ", is.ultrametric(tree.tmp), "\n", sep="")
    cat("Binary: ", is.binary.tree(tree.tmp), "\n", sep="")
    cat("Min edge length: ", min(tree.tmp$edge.length), "\n", sep="")
    
  }  
}
rm(list=c("tree", grep("tmp", ls(), value=TRUE)))

## Picidae morph data
picidae.morph.log.fully_reduced.geomean[["all_inds"]]  # geomean of Picidae, with all individuals included
picidae.morph.log.fully_reduced.phyl_pca[["all_inds"]]$pca$S  # PCA of unscaled data for Picidae, with all individuals included
picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca[["all_inds"]]$pca$S  # PCA of size-scaled shape data for Picidae, with all individuals included

## Picinae morph data
picinae.morph.log.fully_reduced.geomean[["all_inds"]]  # geomean of picinae, with all individuals included
picinae.morph.log.fully_reduced.phyl_pca[["all_inds"]]$pca$S  # PCA of unscaled data for picinae, with all individuals included
picinae.morph.log.fully_reduced.geomean_scaled.phyl_pca[["all_inds"]]$pca$S  # PCA of size-scaled shape data for picinae, with all individuals included

# store vector of data object names
datanames <- c("picidae.morph.log.fully_reduced.geomean", "picidae.morph.log.fully_reduced.phyl_pca", "picidae.morph.log.fully_reduced.geomean_scaled.phyl_pca", "picinae.morph.log.fully_reduced.geomean", "picinae.morph.log.fully_reduced.phyl_pca", "picinae.morph.log.fully_reduced.geomean_scaled.phyl_pca")


### write out trees, prior files, and morphological data for BAMM analyses

## write out the trees to be used in BAMM analyses
for (tree in treenames) {
  if (class(get(tree)) == "phylo") {
    tree.tmp <- get(tree)
    write.tree(phy=tree.tmp, file=paste(tree, ".tre", sep=""))
  } else if (class(get(tree)) == "list") {
    tree.tmp <- get(tree)[["all_inds"]]$phy
    write.tree(phy=tree.tmp, file=paste(tree, ".tre", sep=""))
  }
}
rm(tree, tree.tmp)

## write out ladderized versions of the trees to be used in BAMM analyses
for (tree in treenames) {
  if (class(get(tree)) == "phylo") {
    tree.tmp <- get(tree)
    write.tree(phy=ladderize(tree.tmp), file=paste(tree, "_ladderized.tre", sep=""))
  } else if (class(get(tree)) == "list") {
    tree.tmp <- get(tree)[["all_inds"]]$phy
    write.tree(phy=ladderize(tree.tmp), file=paste(tree, "_ladderized.tre", sep=""))
  }
}
rm(tree, tree.tmp)

## output speciation-extinction priors for all trees to be used in BAMM analyses
for (tree in treenames) {
  if (length(grep("picidae", tree)) > 0) {
    taxa.tmp <- 237
  } else if (length(grep("picinae", tree)) > 0) {
    taxa.tmp <- 203
  }
  
  if (class(get(tree)) == "phylo") {
    tree.tmp <- get(tree)
    setBAMMpriors(phy=tree.tmp, total.taxa=taxa.tmp, outfile=paste(tree, "_se_priors.txt", sep=""))
  } else if (class(get(tree)) == "list") {
    tree.tmp <- get(tree)[["all_inds"]]$phy
    setBAMMpriors(phy=tree.tmp, total.taxa=taxa.tmp, outfile=paste(tree, "_se_priors.txt", sep=""))
  }
}
rm(tree, tree.tmp, taxa.tmp)

## output morphological data for all data sets to be used in BAMM analyses (need to do this before the trait evolution priors, because setBAMMpriors uses a data file, rather than a data object, for the traits)
for (dataname in datanames) {
  data.tmp <- get(dataname)[["all_inds"]]
  
  if (length(grep("picidae", dataname)) > 0) {
    tree.tmp <- picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy
  } else if (length(grep("picinae", dataname)) > 0) {
    tree.tmp <- picinae.morph.log.fully_reduced.treedata[["all_inds"]]$phy
  }
  
  if (length(grep("pca", dataname)) > 0) {
    for (i in 1:3) {
      write.table(data.tmp$pca$S[,i], file=paste(dataname, "_PC", i, "_morph_data.txt", sep=""), sep="\t")
    }
  } else if (length(grep("geomean", dataname)) > 0) {
    write.table(data.tmp[names(data.tmp) %in% tree.tmp$tip.label], file=paste(dataname, "_morph_data.txt", sep=""), sep="\t")
  }
}
rm(list=c("i", "dataname", grep("tmp", ls(), value=TRUE)))

## output morphological evolution priors for all data sets to be used in BAMM analyses
for (dataname in datanames) {
  data.tmp <- get(dataname)[["all_inds"]]
  
  if (length(grep("picidae", dataname)) > 0) {
    tree.tmp <- picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy
  } else if (length(grep("picinae", dataname)) > 0) {
    tree.tmp <- picinae.morph.log.fully_reduced.treedata[["all_inds"]]$phy
  }
    
  if (length(grep("pca", dataname)) > 0) {
    for (i in 1:3) {
      setBAMMpriors(phy=tree.tmp, traits=paste(dataname, "_PC", i, "_morph_data.txt", sep=""), outfile=paste(dataname, "_PC", i, "_morph_priors.txt", sep=""))
    }
  } else if (length(grep("geomean", dataname)) > 0) {
    setBAMMpriors(phy=tree.tmp, traits=paste(dataname, "_morph_data.txt", sep=""), outfile=paste(dataname, "_morph_priors.txt", sep=""))
  }
}
rm(list=c("i","tree","dataname", grep("tmp", ls(), value=TRUE)))


### automate the BAMM post-hoc analyses, using a "BAMM.run_name" and a "BAMM.run_tree" so that I don't have to edit the variables over and over and over

BAMM.run_names <- c("picidae_BAMM_full_tree_se_run1", "picidae_BAMM_full_tree_r8s_scaled_se_run1", "picidae_BAMM_morph_tree_se_run1", "picidae_BAMM_morph_tree_trait_geomean_run1", "picidae_BAMM_morph_tree_trait_unscaled_PC1_run1", "picidae_BAMM_morph_tree_trait_unscaled_PC2_run1", "picidae_BAMM_morph_tree_trait_unscaled_PC3_run1", "picidae_BAMM_morph_tree_trait_geomean_scaled_PC1_run1", "picidae_BAMM_morph_tree_trait_geomean_scaled_PC2_run1", "picidae_BAMM_morph_tree_trait_geomean_scaled_PC3_run1", "picinae_BAMM_full_tree_se_run1", "picinae_BAMM_full_tree_r8s_scaled_se_run1", "picinae_BAMM_morph_tree_se_run1", "picinae_BAMM_morph_tree_trait_geomean_run1", "picinae_BAMM_morph_tree_trait_unscaled_PC1_run1", "picinae_BAMM_morph_tree_trait_unscaled_PC2_run1", "picinae_BAMM_morph_tree_trait_unscaled_PC3_run1", "picinae_BAMM_morph_tree_trait_geomean_scaled_PC1_run1", "picinae_BAMM_morph_tree_trait_geomean_scaled_PC2_run1", "picinae_BAMM_morph_tree_trait_geomean_scaled_PC3_run1")

# get trees for use in post-hoc analyses
picidae.morph_tree_for_BAMM <- picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy
picinae.morph_tree_for_BAMM <- picinae.morph.log.fully_reduced.treedata[["all_inds"]]$phy

# get ladderized versions of trees for use in post-hoc analyses
picidae.morph_tree_for_BAMM.ladderized <- ladderize(picidae.morph.log.fully_reduced.treedata[["all_inds"]]$phy)
picinae.morph_tree_for_BAMM.ladderized <- ladderize(picinae.morph.log.fully_reduced.treedata[["all_inds"]]$phy)

BAMM.run_trees <- c("picidae.RAxML.all.BEAST_calibrated.ladderized", "picidae.RAxML.all.r8s_calibrated.ladderized.scaled.ladderized", rep(x="picidae.morph_tree_for_BAMM.ladderized", times=8), "picinae.RAxML.all.BEAST_calibrated.ladderized", "picinae.RAxML.all.r8s_calibrated.ladderized.scaled.ladderized", rep(x="picinae.morph_tree_for_BAMM.ladderized", times=8))

BAMM.run_types <- rep(x=c(rep("diversification", times=3), rep("trait", times=10)), times=2) # 3 diversification analyses and 10 trait analyses, for picidae and picinae

devAskNewPage(ask=TRUE)
# devAskNewPage(ask=FALSE)
opar <- par()

# create a list to store the various objects and derivatives from the BAMM analyses
BAMM.data <- list()

# iterate over the BAMM runs, importing and processing the data
for (i in 1:length(BAMM.run_names)) {
  cat("\nStarting processing of ", BAMM.run_names[i], "\n", sep="")
  BAMM.data[[i]] <- list()
  names(BAMM.data)[i] <- BAMM.run_names[i]
  
  # Read BAMM output
  BAMM.data[[i]]$eventdata <- getEventData(phy = eval(parse(text=BAMM.run_trees[i])), type=BAMM.run_types[i], eventdata=paste(BAMM.run_names[i], "_event_data.txt", sep=""), burnin=0.1, nsamples=1000)
  BAMM.data[[i]]$priordata <- read.csv(file=paste(BAMM.run_names[i], "_prior_probs.txt", sep=""), header=T)
  
  # Assessing MCMC Convergence
  BAMM.data[[i]]$mcmcout <- read.csv(paste(BAMM.run_names[i], "_mcmc_out.txt", sep=""), header=T)
  BAMM.data[[i]]$burnstart <- floor(0.1 * nrow(BAMM.data[[i]]$mcmcout))
  BAMM.data[[i]]$postburn <- BAMM.data[[i]]$mcmcout[BAMM.data[[i]]$burnstart:nrow(BAMM.data[[i]]$mcmcout),]
  
  # How Many Rate Shifts?
  BAMM.data[[i]]$post_probs <- table(BAMM.data[[i]]$postburn$N_shifts) / nrow(BAMM.data[[i]]$postburn)
  BAMM.data[[i]]$shift_probs <- summary(BAMM.data[[i]]$eventdata, print=F)
  
  # Bayes Factors for Model Comparison
  BAMM.data[[i]]$bfmat <- computeBayesFactors(postdata=BAMM.data[[i]]$mcmcout, priordata=BAMM.data[[i]]$priordata, burnin=0.1)
  
  # Bayesian credible sets of shift configurations
  BAMM.data[[i]]$priorshifts <- getBranchShiftPriors(phy=eval(parse(text=BAMM.run_trees[i])), priordata=BAMM.data[[i]]$priordata)
  BAMM.data[[i]]$css <- credibleShiftSet(BAMM.data[[i]]$eventdata, BAMM.data[[i]]$priorshifts, set.limit=0.95)
  
#   # Finding the single BEST shift configuration; continue from here (need to fix this, as it's throwing an error; could be due to changes in BAMMtools - need to check, but even the examples don't work)
#   cat("\nTrying to get best shift configuration\n")
#   BAMM.data[[i]]$best <- try(getBestShiftConfiguration(BAMM.data[[i]]$eventdata, prior=BAMM.data[[i]]$priorshifts)) # for now, I just test it, and move on if it doesn't work
#   ## can do more with this, pull out samples for different shift configurations, using BAMM.data[[i]]$css$indices
  
  # Maximum Shift Credibility configuration
  BAMM.data[[i]]$msc.set <- maximumShiftCredibility(BAMM.data[[i]]$eventdata, maximize='product')
  BAMM.data[[i]]$msc.config <- subsetEventData(BAMM.data[[i]]$eventdata, index = BAMM.data[[i]]$msc.set$sampleindex)
  
  # Plot marginal shift probabilities (which can be difficult to interpret)
  BAMM.data[[i]]$marg_probs <- marginalShiftProbsTree(BAMM.data[[i]]$eventdata)
  
  # Bayes factors for rate shifts
  BAMM.data[[i]]$bftree <- bayesFactorBranches(BAMM.data[[i]]$eventdata, BAMM.data[[i]]$priorshifts) # generates copy of tree, with branch lengths scaled by Bayes factor evidence for a rate shift on that branch
  BAMM.data[[i]]$edgemax <- which(BAMM.data[[i]]$bftree$edge.length == max(BAMM.data[[i]]$bftree$edge.length))

  # Average rates across entire tree
  cat("\nCalculating evolutionary rates across entire tree\n")
  BAMM.data[[i]]$allrates <- getCladeRates(BAMM.data[[i]]$eventdata)
  
  ## Clade-specific evolutionary rates
  ## get rates for internal nodes in tree (to use this with my tree-breaking approach)
  BAMM.data[[i]]$rates_distribution <- list()
  BAMM.data[[i]]$div_rate_avg <- numeric()
  cat("\nCalculating evolutionary rates by clade for ", length(unique(BAMM.data[[i]]$eventdata$edge[,1])), " edges:\n", sep="")
  print(unique(BAMM.data[[i]]$eventdata$edge[,1]))
  for (j in unique(BAMM.data[[i]]$eventdata$edge[,1])) {
    print(j)
    BAMM.data[[i]]$rates_distribution[[as.character(j)]] <- getCladeRates(BAMM.data[[i]]$eventdata, node=j)
    if (BAMM.run_types[i] == "diversification") {
      BAMM.data[[i]]$div_rate_avg[[as.character(j)]] <- mean(BAMM.data[[i]]$rates_distribution[[as.character(j)]]$lambda - BAMM.data[[i]]$rates_distribution[[as.character(j)]]$mu)
    } else if (BAMM.run_types[i] == "trait") {
      BAMM.data[[i]]$morph_rate_avg[[as.character(j)]] <- mean(BAMM.data[[i]]$rates_distribution[[as.character(j)]]$beta)
    }  
  }
  
  # Macroevolutionary cohort analysis
  BAMM.data[[i]]$cmat <- getCohortMatrix(BAMM.data[[i]]$eventdata)
  
  # Cumulative shift probabilities
  BAMM.data[[i]]$cst <- cumulativeShiftProbsTree(BAMM.data[[i]]$eventdata)
  
  # plot cumulative shift probabilities >= 0.95 in red
  BAMM.data[[i]]$edgecols <- rep('black', length(eval(parse(text=paste(BAMM.run_trees[i], "$edge.length", sep="")))))
  BAMM.data[[i]]$is_highprobshift <- BAMM.data[[i]]$cst$edge.length >= 0.95
  BAMM.data[[i]]$edgecols[BAMM.data[[i]]$is_highprobshift] <- 'red'
}
rm(i,j)

devAskNewPage(ask=TRUE)

## now plot the interesting stuff (including diagnostics and summaries of output)
for (i in 1:length(BAMM.data)) {
  cat("\nStarting plots for ", names(BAMM.data)[i], "\n", sep="")
  
  # Plotting event data
  par(opar)
  cat("\nPlotting event data\n")
  plot.bammdata(BAMM.data[[i]]$eventdata)
  
  # Assessing MCMC Convergence
  par(opar)
  cat("\nPlotting MCMC convergence\n")
  plot(BAMM.data[[i]]$mcmcout$logLik ~ BAMM.data[[i]]$mcmcout$generation)
  cat("\nEffective sample size of N_shifts: ", effectiveSize(BAMM.data[[i]]$postburn$N_shifts), sep="")
  cat("\nEffective sample size of logLike: ", effectiveSize(BAMM.data[[i]]$postburn$logLik), sep="")

  # How Many Rate Shifts?
  par(opar)
  cat("\nDetails on rate shifts\n")
  print(BAMM.data[[i]]$post_probs)
  summary(BAMM.data[[i]]$eventdata, print=T)
  
  # Bayes Factors for Model Comparison
  cat("\nBayes Factors for model comparisons\n")
  BAMM.data[[i]]$bfmat
  
  # Mean phylorate plot
  par(opar)
  cat("\nMean phylorate plot\n")
  plot.bammdata(BAMM.data[[i]]$eventdata, lwd=2, legend=T, spex="netdiv")
  
  # Bayesian credible sets of shift configurations
  par(opar)
  cat("\nNumber of distinct credible shift configurations\n")
  BAMM.data[[i]]$css$number.distinct
  cat("\nDetails on credible shift configurations\n")
  summary(BAMM.data[[i]]$css)
  plot.credibleshiftset(BAMM.data[[i]]$css)
  
  # Finding the single BEST shift configuration; continue from here (need to fix this, as it's throwing an error; could be due to changes in BAMMtools - need to check )
#   par(opar)
#   if (class(BAMM.data[[i]]$best) == "bammdata") {
#     cat("\nPlotting best shift configuration\n")
#     plot.bammdata(BAMM.data[[i]]$best)
#     addBAMMshifts(BAMM.data[[i]]$best)
#   }
  ## can do more with this, pull out samples for different shift configurations, using BAMM.data[[i]]$css$indices
  
  # Maximum Shift Credibility configuration
  par(opar)
  cat("\nPlotting maximum shift credibility configuration\n")
  plot.bammdata(BAMM.data[[i]]$msc.config)
  addBAMMshifts(BAMM.data[[i]]$msc.config)
  
  # Plot marginal shift probabilities (which can be difficult to interpret)
  par(opar)
  cat("\nPlotting marginal shift probabilities\n")
  plot.phylo(BAMM.data[[i]]$marg_probs, show.tip.label=FALSE)
  
  # Bayes factors for rate shifts
  par(opar)
  cat("\nPlotting Bayes Factors for rate shifts\n")
  plot.phylo(BAMM.data[[i]]$bftree, show.tip.label=FALSE) # plot this tree
  
  # Average rates across entire tree
  par(opar)
  cat("\nOutputting evolutionary rates across entire tree\n")
  if (BAMM.run_types[i] == "diversification") {
    cat("\nAverage lambda: ", mean(BAMM.data[[i]]$allrates$lambda), "\n", sep="")
    cat("\nQuantiles of lambda:", quantile(BAMM.data[[i]]$allrates$lambda, c(0.025, 0.975)), "\n", sep=" ") # 95% HPD interval   
  } else if (BAMM.run_types[i] == "trait") {
    cat("\nAverage beta: ", mean(BAMM.data[[i]]$allrates$beta), sep="")
    cat("\nQuantiles of beta:", quantile(BAMM.data[[i]]$allrates$beta, c(0.025, 0.975)), "\n", sep=" ") # 95% HPD interval   
  }
  
  # Rate-throough-time analysis
  par(opar)
  cat("\nPlotting rate through time\n")
  if (BAMM.run_types[i] == "diversification") {
    plotRateThroughTime(BAMM.data[[i]]$eventdata, ratetype='netdiv')
  } else if (BAMM.run_types[i] == "trait") {
    plotRateThroughTime(BAMM.data[[i]]$eventdata, ratetype='auto')
  } 
  
  # Macroevolutionary cohort analysis
  par(opar)
  cat("\nPlotting macroevolutionary cohorts\n")
  cohorts(BAMM.data[[i]]$cmat, BAMM.data[[i]]$eventdata)
  
  # Cumulative shift probabilities
  par(opar)
  cat("\nPlotting cumulative shift probabilities\n")
  plot.phylo(BAMM.data[[i]]$cst, show.tip.label=FALSE)
  
  # plot cumulative shift probabilities >= 0.95 in red
  par(opar)
  cat("Plotting cumulative shift probabilities >= 0.95 in red\n")
  plot.phylo(eval(parse(text=BAMM.run_trees[i])), edge.color = BAMM.data[[i]]$edgecols, show.tip.label=FALSE)
  par(opar)
}

devAskNewPage(ask=FALSE)
rm(i)


### output mean phylorate plots for net diversification to pdf files
for (i in names(BAMM.data)) {
  pdf(file=paste(i, ".mean_phylorate_plot.pdf", sep=""), width=6, height=10)
  plot.bammdata(BAMM.data[[i]]$eventdata, lwd=2, legend=T, spex="netdiv")
  if (length(BAMM.data[[i]]$css$shiftnodes[[i]]) > 0) addBAMMshifts(BAMM.data[[i]]$eventdata, shiftnodes=BAMM.data[[i]]$css$shiftnodes[[1]])
  dev.off()
}
rm(i)


### extract the average morph rates (beta) and diversification rates (lambda-mu) by node, so that I can pass them to my automated subclade analysis file

## get picidae average diversification rates
picidae.BAMM.divrates_by_node <- list()
for (i in grep("picidae", grep("se", names(BAMM.data), value=TRUE), value=TRUE)) {
  picidae.BAMM.divrates_by_node[[gsub("picidae_BAMM_", "", gsub("_se_run1", "", i))]] <- BAMM.data[[i]][["div_rate_avg"]]
}
rm(i)

## get picinae average diversification rates
picinae.BAMM.divrates_by_node <- list()
for (i in grep("picinae", grep("se", names(BAMM.data), value=TRUE), value=TRUE)) {
  picinae.BAMM.divrates_by_node[[gsub("picinae_BAMM_", "", gsub("_se_run1", "", i))]] <- BAMM.data[[i]][["div_rate_avg"]]
}
rm(i)

## get picidae average morph rates
picidae.BAMM.morphrates_by_node <- list()
picidae.BAMM.morphrates_by_node[["geomean"]] <- BAMM.data[["picidae_BAMM_morph_tree_trait_geomean_run1"]][["morph_rate_avg"]]
picidae.BAMM.morphrates_by_node[["phyl_pca_PC13"]] <- BAMM.data[["picidae_BAMM_morph_tree_trait_unscaled_PC1_run1"]][["morph_rate_avg"]] + BAMM.data[["picidae_BAMM_morph_tree_trait_unscaled_PC2_run1"]][["morph_rate_avg"]] + BAMM.data[["picidae_BAMM_morph_tree_trait_unscaled_PC3_run1"]][["morph_rate_avg"]]
picidae.BAMM.morphrates_by_node[["geomean_scaled_phyl_pca_PC13"]] <- BAMM.data[["picidae_BAMM_morph_tree_trait_geomean_scaled_PC1_run1"]][["morph_rate_avg"]] + BAMM.data[["picidae_BAMM_morph_tree_trait_geomean_scaled_PC2_run1"]][["morph_rate_avg"]] + BAMM.data[["picidae_BAMM_morph_tree_trait_geomean_scaled_PC3_run1"]][["morph_rate_avg"]]

## get picinae average morph rates
picinae.BAMM.morphrates_by_node <- list()
picinae.BAMM.morphrates_by_node[["geomean"]] <- BAMM.data[["picinae_BAMM_morph_tree_trait_geomean_run1"]][["morph_rate_avg"]]
picinae.BAMM.morphrates_by_node[["phyl_pca_PC13"]] <- BAMM.data[["picinae_BAMM_morph_tree_trait_unscaled_PC1_run1"]][["morph_rate_avg"]] + BAMM.data[["picinae_BAMM_morph_tree_trait_unscaled_PC2_run1"]][["morph_rate_avg"]] + BAMM.data[["picinae_BAMM_morph_tree_trait_unscaled_PC3_run1"]][["morph_rate_avg"]]
picinae.BAMM.morphrates_by_node[["geomean_scaled_phyl_pca_PC13"]] <- BAMM.data[["picinae_BAMM_morph_tree_trait_geomean_scaled_PC1_run1"]][["morph_rate_avg"]] + BAMM.data[["picinae_BAMM_morph_tree_trait_geomean_scaled_PC2_run1"]][["morph_rate_avg"]] + BAMM.data[["picinae_BAMM_morph_tree_trait_geomean_scaled_PC3_run1"]][["morph_rate_avg"]]


### output figures
#####
pdf(file="picidae.RAxML.bamm_div.pdf", width=4, height=10) # open pdf device
plot.bammdata(BAMM.data[[1]]$eventdata)
dev.off()

pdf(file="picidae.RAxML.bamm_geomean.pdf", width=4, height=10) # open pdf device
plot.bammdata(picidae.morph_tree_trait_geomean_run2.BAMM)
dev.off()

pdf(file="picidae.RAxML.bamm_geomean_scaled_pca1.pdf", width=4, height=10) # open pdf device
plot.bammdata(picidae.morph_tree_trait_geomean_scaled_pca1_run2.BAMM)
dev.off()


### plot maximum shift credibility trees

pdf(file="picidae.RAxML.bamm_geomean.msc.pdf", width=4, height=10) # open pdf device
plot.bammdata(picidae.morph_tree_trait_geomean_run2.BAMM.msc.config)
addBAMMshifts(picidae.morph_tree_trait_geomean_run2.BAMM.msc.config)
dev.off()

pdf(file="picidae.RAxML.bamm_geomean_scaled_pca1.msc.pdf", width=4, height=10) # open pdf device
plot.bammdata(picidae.morph_tree_trait_geomean_scaled_pca1_run2.BAMM.msc.config)
addBAMMshifts(picidae.morph_tree_trait_geomean_scaled_pca1_run2.BAMM.msc.config)
dev.off()


### save off objects I want to use in the automated subclade analyses (the diversification (lambda) and trait evolution (beta) rates by node)
save(file="Picidae_BAMM_data_for_automated_subclade_analyses.RData", list=c(grep("rates_by_node", ls(), value=TRUE)))
