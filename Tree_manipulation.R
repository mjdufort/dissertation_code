# this file contains scripts for reading, modifying, and plotting the various phylogenetic trees used in Chapter 1 and in downstream analyses
# it also contains scripts for calculating differences among trees, comparing trees both within this project, and with the results of previous supertree and supermatrix analyses of birds


## load packages
library(ape)
library(stringr)
library(phyloch)
library(phangorn)
library(phytools)

### read in, trim, and root the trees

taxa.remove.final_trees <- c("Ramphastidae_spp","Indicatoridae_spp") # taxa to remove from the trees (outgroups, taxa that aren't species)

picidae.RAxML.all <- read.tree("picidae_RAxML_all.tre")  # read RAxML tree from concatenated analysis of all loci
picidae.RAxML.all_exFGBi7 <- read.tree("picidae_RAxML_all_ex-FGBi7.tre")  # read RAxML tree from concatenated analysis of all loci minus FGB-I7
picidae.RAxML.mito <- read.tree("picidae_RAxML_mito.tre")  # read RAxML tree from concatenated analysis of all mitochondrial loci

picidae.RAxML.nuclear <- read.tree("picidae_RAxML_nuclear.tre")  # read RAxML tree from concatenated analysis of all nuclear loci
picidae.RAxML.all.PartitionFinder <- read.tree("picidae_RAxML_all_PartitionFinder.tre")  # read RAxML tree from concatenated analysis of all loci, with partitions as determined by PartitionFinder
picidae.RAxML.all.unpartitioned <- read.tree("picidae_RAxML_all_unpartitioned.tre")  # read RAxML tree from concatenated analysis of all loci, with a single partition

## read in, root, and write out the RAxML locus_by_locus trees
picidae.RAxML.locus_by_locus.locus_list <- c("ACO1i9","FGBi5","FGBi7","BRMi15","CMOS","GAPDHi11","HMGN2","IRBP","LDH","mito12S","mito16S","mitoATP6_Picus_rabieri_problem","mitoATP6","mitoATP8","mitoCOI","mitoCOIII","mitoCR","mitoCYTB_inc-Mstri","mitoCYTB_inc-Xiph-Mstri","mitoCYTB","mitoND2","mitoND3","MUSKi4","MYOi2","PEPCKi9","PERi9","RAG1","TGFBi5")
picidae.RAxML.locus_by_locus <- list()
for (i in picidae.RAxML.locus_by_locus.locus_list) {
  picidae.RAxML.locus_by_locus[[i]] <- read.tree(paste("picidae_", i, ".tre", sep=""))
  if ("Ramphastidae_spp" %in% picidae.RAxML.locus_by_locus[[i]]$tip.label) {
    picidae.RAxML.locus_by_locus[[i]] <- root(picidae.RAxML.locus_by_locus[[i]], "Ramphastidae_spp")
  } else if ("Indicatoridae_spp" %in% picidae.RAxML.locus_by_locus[[i]]$tip.label) {
    picidae.RAxML.locus_by_locus[[i]] <- root(picidae.RAxML.locus_by_locus[[i]], "Indicatoridae_spp")
  }
  write.tree(picidae.RAxML.locus_by_locus[[i]], file=paste("picidae_", i, "_rooted.tre", sep=""))
}



picidae.STAR.all <- read.tree("picidae_STAR_all.tre")  # read STAR tree from analysis of ML gene trees from RAxML analysis of individual nuclear loci and concatenated mtDNA
picidae.STAR.all.bootstrap <- read.tree("picidae_STAR_all_bootstrap.tre")  # read STAR tree from analysis of bootstrap trees from RAxML analysis of individual nuclear loci and concatenated mtDNA
picidae.STAR.all.bootstrap.rooted <- root(picidae.STAR.all.bootstrap, outgroup="Ramphastidae_spp")  # root the STAR bootstrap tree on Ramphastidae outgroup

picidae.STAR.all_exFGBi7 <- read.tree("picidae_STAR_all_ex-FGBi7.tre")
picidae.STAR.all_exFGBi7.bootstrap <- read.tree("picidae_STAR_all_ex-FGBi7_bootstrap.tre")

picidae.RAxML.all.r8s_calibrated <- read.tree("picidae_RAxML_all_r8s_calibrated_inc_calibs.tre")

picidae.RAxML.all.r8s_calibrated.no_calibs <- read.tree("picidae_RAxML_all_r8s_calibrated_no_calibs.tre")

picidae.mrbayes.all <- read.nexus("picidae_mrbayes_all.nex")  # read MrBayes tree from concatenated analysis of all loci 

# reroot trees on Ramphastidae and ladderize them for ease of viewing
picidae.RAxML.all.rooted.ladderized <- ladderize(root(picidae.RAxML.all, outgroup="Ramphastidae_spp"), right=FALSE)
picidae.RAxML.all_exFGBi7.rooted.ladderized <- ladderize(root(picidae.RAxML.all_exFGBi7, outgroup="Ramphastidae_spp"), right=FALSE)
picidae.RAxML.mito.rooted.ladderized <- ladderize(root(picidae.RAxML.mito, outgroup="Ramphastidae_spp"), right=FALSE)
picidae.RAxML.nuclear.rooted.ladderized <- ladderize(root(picidae.RAxML.nuclear, outgroup="Ramphastidae_spp"), right=FALSE)
picidae.STAR.all.rooted.ladderized <- ladderize(root(picidae.STAR.all, outgroup="Ramphastidae_spp"), right=FALSE)
picidae.STAR.all_exFGBi7.rooted.ladderized <- ladderize(root(picidae.STAR.all_exFGBi7, outgroup="Ramphastidae_spp"), right=FALSE)

picidae.STAR.all.bootstrap.rooted.ladderized <- ladderize(root(picidae.STAR.all.bootstrap, outgroup="Ramphastidae_spp"), right=FALSE)
picidae.STAR.all.bootstrap.rooted.ladderized <- rotate(picidae.STAR.all.bootstrap.rooted.ladderized, node=179, polytom=c(1,2))  # ladderizing left one node oddly positioned
picidae.STAR.all.bootstrap.rooted.ladderized$node.label <- str_replace(picidae.STAR.all.bootstrap.rooted.ladderized$node.label, pattern="\\.00", replacement="") # remove the unnecessary ".00" in the bootstrap support values for the STAR tree
picidae.STAR.all_exFGBi7.bootstrap.rooted.ladderized <- ladderize(root(picidae.STAR.all_exFGBi7.bootstrap, outgroup="Ramphastidae_spp"), right=FALSE)
picidae.RAxML.all.r8s_calibrated.ladderized <- ladderize(picidae.RAxML.all.r8s_calibrated, right=FALSE)
picidae.RAxML.all.r8s_calibrated.no_calibs.ladderized <- ladderize(drop.tip(picidae.RAxML.all.r8s_calibrated.no_calibs, tip="Ramphastidae_spp"), right=FALSE)
picidae.mrbayes.all.ladderized <- ladderize(picidae.mrbayes.all, right=FALSE)

write.tree(picidae.RAxML.all.rooted.ladderized, file="picidae_RAxML_all_rooted_ladderized.tre")  # write out the rerooted, ladderized tree


## compare trees

picidae.RAxML.all_ex_Crob_Cchr.rooted.ladderized <- ladderize(drop.tip(picidae.RAxML.all.rooted, tip=c("Campephilus_robustus", "Colaptes_chrysoides")), right=FALSE)  # drop the tips from the RAxML tree that are not present in STAR tree, to allow comparison of the two trees

picidae.STAR.all.rooted.ladderized.brlen1 <- compute.brlen(picidae.STAR.all.rooted.ladderized, 1)  # create a copy of the STAR tree with all branch lengths set to 1 (for ease of visualization)
picidae.STAR.all_exFGBi7.rooted.ladderized.brlen1 <- compute.brlen(picidae.STAR.all_exFGBi7.rooted.ladderized, 1)  # create a copy of the STAR bootstrap tree with all branch lengths set to 1 (for ease of visualization)

cophyloplot(picidae.STAR.all.rooted.ladderized.brlen1, picidae.STAR.all_exFGBi7.rooted.ladderized.brlen1, assoc=cbind(sort(picidae.STAR.all.rooted.ladderized.brlen1$tip.label),sort(picidae.STAR.all_exFGBi7.rooted.ladderized.brlen1$tip.label)), cex=0.5, show.tip.label=T, space=100, gap=4, use.edge.length=F, node.depth=2)  # plot the STAR trees with and without FGBi7 opposite each other to visualize differences

cophyloplot(picidae.STAR.all.rooted.ladderized, picidae.RAxML.all_ex_Crob_Cchr.rooted.ladderized, assoc=cbind(sort(picidae.STAR.all.rooted.ladderized$tip.label),sort(picidae.RAxML.all_ex_Crob_Cchr.rooted.ladderized$tip.label)), cex=0.5, show.tip.label=T, space=100, gap=4, use.edge.length=F, node.depth=2)  # plot the STAR tree and RAxML tree opposite each other to visualize differences

cophyloplot(picidae.RAxML.all.rooted.ladderized, picidae.RAxML.all_exFGBi7.rooted.ladderized, assoc=cbind(sort(picidae.RAxML.all.rooted.ladderized$tip.label),sort(picidae.RAxML.all_exFGBi7.rooted.ladderized$tip.label)), cex=0.5, show.tip.label=T, space=100, gap=4, use.edge.length=F, node.depth=2)  # plot the RAxML trees with and without FGBi7 opposite each other to visualize differences

cophyloplot(picidae.RAxML.all.rooted.ladderized, picidae.mrbayes.all.ladderized, assoc=cbind(sort(picidae.RAxML.all.rooted.ladderized$tip.label),sort(picidae.mrbayes.all.ladderized$tip.label)), cex=0.5, show.tip.label=T, space=100, gap=4, use.edge.length=F, node.depth=2)  # plot the RAxML tree and MrBayes tree opposite each other to visualize differences


### output the trees

## plot the RAxML tree
pdf(file="picidae_RAxML_all_tree_rooted_ladderized.pdf", width=10, height=4.5)
plot(picidae.RAxML.all.rooted.ladderized, direction="up", cex=0.35, no.margin=T, y.lim=c(0,0.35))
nodelabels(text=picidae.RAxML.all.rooted.ladderized$node.label, col="gray1", frame="none", cex=0.3, srt=90)
add.scale.bar(col="gray2", cex=0.6)
dev.off()

## remove outgroups from the RAxML tree
picidae.RAxML.all.rooted.ladderized.no_outgroups <- ladderize(drop.tip(root(picidae.RAxML.all, outgroup="Ramphastidae_spp"), tip=c("Indicatoridae_spp", "Ramphastidae_spp")), right=TRUE)  # reroot on Ramphastidae, drop outgroups, and ladderize
picidae.RAxML.all.rooted.ladderized.no_outgroups$root.edge <- max(node.depth.edgelength(picidae.RAxML.all.rooted.ladderized.no_outgroups)) * 0.02 # add a short root edge to the tree to show where the root is

## plot the RAxML tree without outgroups
pdf(file="picidae_RAxML_all_tree_rooted_ladderized_no_outgroups.pdf", width=4.5, height=10)
plot(picidae.RAxML.all.rooted.ladderized.no_outgroups, direction="rightwards", cex=0.35, no.margin=T, x.lim=c(0,0.18), root.edge=TRUE)
nodelabels(text=picidae.RAxML.all.rooted.ladderized.no_outgroups$node.label, col="gray1", frame="none", cex=0.3, srt=0)
add.scale.bar(col="gray2", cex=0.6)
dev.off()


## plot the STAR bootstrap tree (for figure 1.3)
pdf(file="picidae_STAR_all_bootstrap_tree_rooted_ladderized.pdf", width=8.1, height=5.0)
plot(picidae.STAR.all.bootstrap.rooted.ladderized, use.edge.length=F, node.depth=2, direction="up", cex=0.29, no.margin=T, label.offset=0.05)
nodelabels(text=picidae.STAR.all.bootstrap.rooted.ladderized$node.label, col="gray1", frame="none", cex=0.25, srt=90, pos=3, offset=0.5)
dev.off()

## remove outgroups from the STAR bootstrap tree
picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups <- ladderize(drop.tip(root(picidae.STAR.all.bootstrap, outgroup="Ramphastidae_spp"), tip=c("Indicatoridae_spp", "Ramphastidae_spp")), right=TRUE) # reroot the tree on Ramphastidae and ladderize it
picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups$node.label <- str_replace(picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups$node.label, pattern="\\.00", replacement="") # remove the weird ".00" in the bootstrap support values for the STAR tree
picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups$root.edge <- max(node.depth.edgelength(picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups)) * 0.02 # add a short root edge to the tree so that I can show that it's rooted (though this isn't used with the STAR tree)

## plot the STAR bootstrap tree with no outgroups
pdf(file="picidae_STAR_all_bootstrap_tree_rooted_ladderized_no_outgroups.pdf", width=3.54, height=8.1)
plot(picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups, use.edge.length=F, node.depth=2, direction="rightwards", cex=0.29, no.margin=T, label.offset=0.05)
nodelabels(text=picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups$node.label, col="gray1", frame="none", cex=0.25, srt=0, adj=1)
dev.off()



## plot the RAxML tree (uncalibrated) with RAxML and MrBayes support values (for Figure 1.2)

picidae.RAxML.all.mrbayes_support <- read.tree(file="picidae_RAxML_all_tree_mrbayes_support.tre")  # read in the RAxML tree with MrBayes support values (determined using a program outside R)
picidae.RAxML.all.mrbayes_support.rooted.ladderized.no_outgroups <- ladderize(drop.tip(root(picidae.RAxML.all.mrbayes_support, outgroup="Ramphastidae_spp"), tip=c("Indicatoridae_spp", "Ramphastidae_spp")), right=TRUE) # reroot on Ramphastidae, ladderize, and drop the outgroups

pdf(file="RAxML_all_tree_rooted_ladderized_no_outgroups_MrBayes_support.pdf", width=5.0, height=8.1)
plot(picidae.RAxML.all.rooted.ladderized.no_outgroups, direction="rightwards", cex=0.29, no.margin=T, x.lim=c(0,0.16), root.edge=TRUE, label.offset=0.0005)  # plot RAxML tree
nodelabels(text=paste(picidae.RAxML.all.rooted.ladderized.no_outgroups$node.label, picidae.RAxML.all.mrbayes_support.rooted.ladderized.no_outgroups$node.label, sep="/"), col="gray1", frame="none", cex=0.21, srt=0, adj=1)  # add RAxML support values
nodelabels(text=picidae.RAxML.all.mrbayes_support.rooted.ladderized.no_outgroups$node.label, col="gray2", frame="none", cex=0.25, srt=0)  # add MrBayes support values
add.scale.bar(col="gray3", cex=0.5)
dev.off()


picidae.mrbayes.all.ladderized.no_outgroups <- ladderize(drop.tip(picidae.mrbayes.all, tip=c("Indicatoridae_spp", "Ramphastidae_spp")), right=TRUE)  # ladderize the MrBayes tree and drop outgroups
picidae.mrbayes.all.ladderized.no_outgroups$root.edge <- max(node.depth.edgelength(picidae.mrbayes.all.ladderized.no_outgroups)) * 0.02 # add a short root edge to the tree so that I can show that it's rooted

## output the MrBayes tree
pdf(file="picidae_mrbayes_all_tree_ladderized_no_outgroups.pdf", width=4.5, height=10)
plot(picidae.mrbayes.all.ladderized.no_outgroups, direction="rightwards", cex=0.35, no.margin=T, x.lim=c(0,0.09), root.edge=TRUE)
nodelabels(text=as.character(round(as.numeric(picidae.mrbayes.all.ladderized.no_outgroups$node.label),digits=2)), col="gray1", frame="none", cex=0.3, srt=0)  # add support values
add.scale.bar(col="gray2", cex=0.6)
dev.off()


### calculate node support correlations
## this method uses only the nodes present in both trees, using matchNodes function of
## this may overestimate the correlation, as it ignores the effect of nodes that are not found in both trees (which are likely to have lower correlation than those found in both trees)

matchNodes.RAxML_mrbayes <- matchNodes(tr1=picidae.RAxML.all.rooted.ladderized.no_outgroups, tr2=picidae.mrbayes.all.ladderized.no_outgroups, method="descendants")  # match node numbers in the RAxML and MrBayes trees
support.correlation.RAxML_mrbayes <- cor.test(as.numeric(picidae.RAxML.all.rooted.ladderized.no_outgroups$node.label[matchNodes.RAxML_mrbayes[,1]-length(picidae.RAxML.all.rooted.ladderized.no_outgroups$tip.label)]), as.numeric(picidae.mrbayes.all.ladderized.no_outgroups$node.label[matchNodes.RAxML_mrbayes[,2]-length(picidae.mrbayes.all.ladderized.no_outgroups$tip.label)]), alternative="greater", method="pearson", use="pairwise.complete.obs")  # calculate the correlation between node support values in the RAxML and MrBayes trees

matchNodes.RAxML_STAR <- matchNodes(tr1=picidae.RAxML.all.rooted.ladderized.no_outgroups, tr2=picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups, method="descendants")  # match node numbers in the RAxML and STAR trees
support.correlation.RAxML_STAR <- cor.test(as.numeric(picidae.RAxML.all.rooted.ladderized.no_outgroups$node.label[matchNodes.RAxML_STAR[,1]-length(picidae.RAxML.all.rooted.ladderized.no_outgroups$tip.label)]), as.numeric(picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups$node.label[matchNodes.RAxML_STAR[,2]-length(picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups$tip.label)]), alternative="greater", method="pearson", use="pairwise.complete.obs")  # calculate the correlation between node support values in the RAxML and STAR trees

matchNodes.mrbayes_STAR <- matchNodes(tr1=picidae.mrbayes.all.ladderized.no_outgroups, tr2=picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups, method="descendants")  # match node numbers in the MrBayes and STAR trees
support.correlation.mrbayes_STAR <- cor.test(as.numeric(picidae.mrbayes.all.ladderized.no_outgroups$node.label[matchNodes.mrbayes_STAR[,1]-length(picidae.mrbayes.all.ladderized.no_outgroups$tip.label)]), as.numeric(picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups$node.label[matchNodes.mrbayes_STAR[,2]-length(picidae.STAR.all.bootstrap.rooted.ladderized.no_outgroups$tip.label)]), alternative="greater", method="pearson", use="pairwise.complete.obs")  # calculate the correlation between node support values in the MrBayes and STAR trees


### calculate correlations of branch lengths and node depths between time-calibrated  trees using different calibration methods

matchNodes.r8s.calibs.no_calibs <- matchNodes(tr1=picidae.RAxML.all.r8s_calibrated, tr2=picidae.RAxML.all.r8s_calibrated.no_calibs, method="descendants")  # match nodes between the RAxML trees with diverence times estimated in r8s with and without internal calibration points (this shouldn't be necessary, but it's safer to not assume that the nodes have the same numbers)
brlen.correlations.r8s.calibs.no_calibs <- cor.test(picidae.RAxML.all.r8s_calibrated$edge.length[match(matchNodes.r8s.calibs.no_calibs[,1],picidae.RAxML.all.r8s_calibrated$edge[,2])], picidae.RAxML.all.r8s_calibrated.no_calibs$edge.length[match(matchNodes.r8s.calibs.no_calibs[,2],picidae.RAxML.all.r8s_calibrated.no_calibs$edge[,2])], alternative="greater", method="pearson", use="pairwise.complete.obs")  # calculate correlation of branch lengths between the RAxML trees with divergence times estimated in r8s with and without internal calibration points
node.depth.correlations.r8s.calibs.no_calibs <- cor.test(node.depth.edgelength(picidae.RAxML.all.r8s_calibrated.no_calibs)[matchNodes.r8s.calibs.no_calibs[,2]], node.depth.edgelength(picidae.RAxML.all.r8s_calibrated)[matchNodes.r8s.calibs.no_calibs[,1]], alternative="greater", method="pearson", use="pairwise.complete.obs")  # calculate correlation of node depths between the RAxML trees with divergence times estimated in r8s with and without internal calibration points

picidae.RAxML.all.BEAST_calibrated <- read.beast("picidae_RAxML_all_BEAST_calibrated_MCC.tre")  # read in the RAxML tree with divergence times estimated using BEAST, with internal calibration points
write.tree(picidae.RAxML.all.BEAST_calibrated, file="picidae_RAxML_all_BEAST_calibrated_MCC.newick.tre")  # output the RAxML tree with divergence times estimated using BEAST in simple newick format

picidae.RAxML.all.BEAST_calibrated.no_calibs <- read.beast("picidae_RAxML_all_BEAST_calibrated_no_calibs_MCC.tre")  # read in the RAxML tree with divergence times estimated using BEAST, without internal calibration points
write.tree(picidae.RAxML.all.BEAST_calibrated.no_calibs, file="picidae_RAxML_all_BEAST_calibrated_no_calibs_MCC.newick.tre")  # output the RAxML tree with divergence times estimated using BEAST (without internal calibration points) in simple newick format

matchNodes.BEAST.calibs.no_calibs <- matchNodes(tr1=picidae.RAxML.all.BEAST_calibrated, tr2=picidae.RAxML.all.BEAST_calibrated.no_calibs, method="descendants")  # match nodes between the RAxML trees with diverence times estimated in BEAST with and without internal calibration points (this shouldn't be necessary, but it's safer to not assume that the nodes have the same numbers)
node.depth.correlations.BEAST.calibs.no_calibs <- cor.test(node.depth.edgelength(picidae.RAxML.all.BEAST_calibrated)[matchNodes.BEAST.calibs.no_calibs[,1]], node.depth.edgelength(picidae.RAxML.all.BEAST_calibrated.no_calibs)[matchNodes.BEAST.calibs.no_calibs[,2]], alternative="greater", method="pearson", use="pairwise.complete.obs")  # calculate correlation of node depths between the RAxML trees with divergence times estimated in BEAST with and without internal calibration points

matchNodes.r8s.BEAST.calibs <- matchNodes(tr1=picidae.RAxML.all.r8s_calibrated, tr2=picidae.RAxML.all.BEAST_calibrated, method="descendants")  # match nodes between the RAxML trees with diverence times estimated in r8s and BEAST, with internal calibration points
node.depth.correlations.r8s.BEAST.calibs <- cor.test(node.depth.edgelength(picidae.RAxML.all.r8s_calibrated)[matchNodes.r8s.BEAST.calibs[,1]], node.depth.edgelength(picidae.RAxML.all.BEAST_calibrated)[matchNodes.r8s.BEAST.calibs[,2]], alternative="greater", method="pearson", use="pairwise.complete.obs")  # calculate correlation of node depths between the RAxML trees with divergence times estimated in r8s and BEAST, with internal calibration points


### plot RAxML tree with divergence times estimated in BEAST, with HPD on node ages

picidae.RAxML.all.BEAST_calibrated.ladderized <- ladderize(picidae.RAxML.all.BEAST_calibrated, right=TRUE)  # ladderize the tree

## output the RAxML tree with divergence times estimated in BEAST
pdf(file="/Users/MattDufort/Documents/Grad_school/Dissertation_stuff/Dissertation_chapters/Chapter_1_Supermatrix_analyses/RAxML_all_tree_BEAST_calibrated.rooted.ladderized.pdf", width=3.54, height=8.1)
par(mar=c(1.5,0,0,0))
plot.phylo(picidae.RAxML.all.BEAST_calibrated.ladderized, direction="rightwards", cex=0.29, x.lim=c(-7, 72), label.offset=0.2)  # plot the tree
HPDbars(picidae.RAxML.all.BEAST_calibrated.ladderized, label="height_95%_HPD", col="gray60", lwd=2.5)  # add 95% HPD on node ages to the nodes
axisPhylo(side=1, cex.axis=0.7, line=-1, col="gray1", tck=-0.02)  # add time axis
dev.off()

write.tree(ladderize(picidae.RAxML.all.BEAST_calibrated.ladderized, right=FALSE), file = "picidae_RAxML_all_BEAST_calibrated_no_calibs_MCC_ladderized.newick.tre")  # output the RAxML tree with divergence times estimated using BEAST in simple newick format



### read in trees from previous supermatrix and supertree analyses, for comparison

## read in Davis and Page supertree
supertree.davis_page <- read.nexus("5.BirdSupertree.tre")
supertree.davis_page.picidae <- extract.clade(supertree.davis_page, node=getMRCA(phy=supertree.davis_page, tip=c("Jynx_torquilla", "Blythipicus_pyrrhotis")))  # extract the clade containing all Picidae (using taxa whose relationship spans the basal split)
write.tree(phy=supertree.davis_page.picidae, file="5.BirdSupertree.picidae.tre")  # write out the Picidae portion of the tree for use in other programs

## read in Jetz et al supermatrix tree
supertree.jetz.picidae <- read.tree("Jetz_Picidae.MCC.newick.tre")

## read in Burleigh supermatrix tree
supertree.burleigh <- read.tree("BigBird.All.newNames.7000Taxa.ML.tre")
supertree.burleigh.picidae <- extract.clade(phy=supertree.burleigh, node=getMRCA(supertree.burleigh, tip=grep('Picidae', supertree.burleigh$tip.label, value=TRUE)), root.edge=1)  # extract the clade containing all Picidae (using the naming convention that includes the family before the genus and species)
supertree.burleigh.picidae$tip.label <- gsub("Piciformes_Picidae_", "", supertree.burleigh.picidae$tip.label)  # remove the order and family from the tip names to simplify comparison


### calculate topological distance between trees

## drop outgroups from trees to work with
picidae.RAxML.all.no_outgroups <- drop.tip(picidae.RAxML.all, tip=c("Indicatoridae_spp", "Ramphastidae_spp"))  # drop outgroups from RAxML tree
picidae.STAR.all.no_outgroups <- drop.tip(picidae.STAR.all, tip=c("Indicatoridae_spp", "Ramphastidae_spp"))  # drop outgroups from STAR tree
picidae.mrbayes.all.no_outgroups <- drop.tip(picidae.mrbayes.all, tip=c("Indicatoridae_spp", "Ramphastidae_spp"))  # drop outgroups from RAxML tree

## calculate tree distance metrics for pairs of trees, after dropping tips not shared between them
treedist.RAxML_STAR <- treedist(drop.tip(picidae.RAxML.all.no_outgroups, tip=setdiff(picidae.RAxML.all.no_outgroups$tip.label, picidae.STAR.all.no_outgroups$tip.label)), drop.tip(picidae.STAR.all.no_outgroups, tip=setdiff(picidae.STAR.all.no_outgroups$tip.label, picidae.RAxML.all.no_outgroups$tip.label)))
treedist.RAxML_mrbayes <- treedist(drop.tip(picidae.RAxML.all.no_outgroups, tip=setdiff(picidae.RAxML.all.no_outgroups$tip.label, picidae.mrbayes.all.no_outgroups$tip.label)), drop.tip(picidae.mrbayes.all.no_outgroups, tip=setdiff(picidae.mrbayes.all.no_outgroups$tip.label, picidae.RAxML.all.no_outgroups$tip.label)))
treedist.mrbayes_STAR <- treedist(drop.tip(picidae.mrbayes.all.no_outgroups, tip=setdiff(picidae.mrbayes.all.no_outgroups$tip.label, picidae.STAR.all.no_outgroups$tip.label)), drop.tip(picidae.STAR.all.no_outgroups, tip=setdiff(picidae.STAR.all.no_outgroups$tip.label, picidae.mrbayes.all.no_outgroups$tip.label)))

## update tip names in Davis and Page tree to match the genera in my taxonomy (note that this checks primarily for changes in genus name, but I did NOT check if they had errors in their naming)
supertree.davis_page.picidae.taxnames_corrected <- supertree.davis_page.picidae  # create a copy of the Davis & Page supertree
supertree.davis_page.picidae.taxnames_corrected$tip.label[!(supertree.davis_page.picidae.taxnames_corrected$tip.label %in% picidae.RAxML.all.no_outgroups$tip.label)]  # identify tips in the Davis & Page tree that are not in my RAxML tree
supertree.davis_page.picidae.taxnames_corrected$tip.label[!(supertree.davis_page.picidae.taxnames_corrected$tip.label %in% picidae.RAxML.all.no_outgroups$tip.label)] <- c("Chrysocolaptes_validus", "Chrysophlegma_miniaceum", "Chrysophlegma_flavinucha", "Chrysophlegma_mentale", "Colaptes_rubiginosus", "Colaptes_rivolii", "Piculus_leucolaemus", "Dryocopus_funebris", "Micropternus_brachyurus", "Dendrocopus_leucotos", "Dendrocopos_noguchii", "Dryobates_minor", "Dryobates_pubescens", "Dryobates_nuttallii", "Dryobates_scalaris", "Leuconotopicus_borealis", "Leuconotopicus_fumigatus", "Leuconotopicus_villosus", "Leuconotopicus_albolarvatus", "Leuconotopicus_stricklandi", "Veniliornis_lignarius", "Veniliornis_mixtus", "Leiopicus_mahrattensis", "Leiopicus_medius", "Picoides_kizuki", "Picoides_canicapillus", "Melanerpes_portoricensis", "Verreauxia_africana")  # update names of tips in the Davis & Page tree that don't match names in my taxonomy

treedist.RAxML_davis_page <- treedist(drop.tip(picidae.RAxML.all.no_outgroups, tip=setdiff(picidae.RAxML.all.no_outgroups$tip.label, supertree.davis_page.picidae.taxnames_corrected$tip.label)), drop.tip(supertree.davis_page.picidae.taxnames_corrected, tip=setdiff(supertree.davis_page.picidae.taxnames_corrected$tip.label, picidae.RAxML.all.no_outgroups$tip.label)))  # calculate tree distance metrics between Davis & Page supertree and my RAxML tree

## update tip names in Jetz et al. tree to match the genera in my taxonomy (note that this checks primarily for changes in genus name, but I did NOT check if they had errors in their naming)
supertree.jetz.picidae.taxnames_corrected <- supertree.jetz.picidae  # create a copy of the Jetz et al. tree
supertree.jetz.picidae.taxnames_corrected$tip.label[!(supertree.jetz.picidae.taxnames_corrected$tip.label %in% picidae.RAxML.all.no_outgroups$tip.label)]  # identify tips in the Jetz et al. tree that are not in my RAxML tree
supertree.jetz.picidae.taxnames_corrected$tip.label[!(supertree.jetz.picidae.taxnames_corrected$tip.label %in% picidae.RAxML.all.no_outgroups$tip.label)] <- c("Pteroglossus_castanotis", "Indicator_minor", "Verreauxia_africana", "Xiphidiopicus_percussus", "Picoides_kizuki", "Picoides_canicapillus", "Picoides_moluccensis", "Picoides_maculatus", "Leiopicus_mahrattensis", "Leiopicus_medius", "Dendropicos_elliotii", "Dendropicos_pyrrhogaster", "Dendropicos_goertae", "Dendropicos_griseocephalus", "Dryobates_scalaris", "Dryobates_nuttallii", "Dryobates_pubescens", "Dryobates_minor", "Leuconotopicus_borealis", "Leuconotopicus_fumigatus", "Leuconotopicus_villosus", "Leuconotopicus_albolarvatus", "Leuconotopicus_stricklandi", "Leuconotopicus_arizonae", "Chrysocolaptes_validus", "Micropternus_brachyurus", "Chrysophlegma_miniaceum", "Chrysophlegma_flavinucha", "Chrysophlegma_mentale", "Dryocopus_funebris", "Piculus_leucolaemus")  # update names of tips in the Jetz et al. tree that don't match names in my taxonomy

treedist.RAxML_jetz <- treedist(drop.tip(picidae.RAxML.all.no_outgroups, tip=setdiff(picidae.RAxML.all.no_outgroups$tip.label, supertree.jetz.picidae.taxnames_corrected$tip.label)), drop.tip(supertree.jetz.picidae.taxnames_corrected, tip=setdiff(supertree.jetz.picidae.taxnames_corrected$tip.label, picidae.RAxML.all.no_outgroups$tip.label)))  # calculate tree distance metrics between Jetz et al. supermatrix tree and my RAxML tree

## update tip names in Burleight et al. tree to match the genera in my taxonomy (note that this checks primarily for changes in genus name, but I did NOT check if they had errors in their naming)
supertree.burleigh.picidae.taxnames_corrected <- supertree.burleigh.picidae  # create a copy of the Burleigh et al. tree
supertree.burleigh.picidae.taxnames_corrected$tip.label[!(supertree.burleigh.picidae.taxnames_corrected$tip.label %in% picidae.RAxML.all.no_outgroups$tip.label)]  # identify tips in the Burleigh et al. tree that are not in my RAxML tree
supertree.burleigh.picidae.taxnames_corrected$tip.label[!(supertree.burleigh.picidae.taxnames_corrected$tip.label %in% picidae.RAxML.all.no_outgroups$tip.label)] <- c("Piculus_leucolaemus", "Dryocopus_funebris", "Chrysophlegma_flavinucha", "Chrysophlegma_mentale", "Chrysophlegma_miniaceum", "Micropternus_brachyurus", "Dendrocopos_leucotos_SYN", "Dendrocopos_noguchii", "Dendrocopos_major_SYN", "Dryobates_pubescens", "Dryobates_minor", "Dendrocopos_minor_SYN", "Dryobates_nuttallii", "Dryobates_scalaris", "Xiphidiopicus_percussus", "Leuconotopicus_stricklandi", "Leuconotopicus_arizonae", "Leuconotopicus_albolarvatus", "Leuconotopicus_villosus", "Leuconotopicus_fumigatus", "Leuconotopicus_borealis", "Dendropicos_fuscescens_SYN", "Leiopicus_mahrattensis", "Leiopicus_medius", "Dendrocopos_canicapillus_SYN", "Picoides_canicapillus", "Picoides_moluccensis", "Picoides_maculatus", "Picoides_kizuki", "Chrysocolaptes_validus", "Verreauxia_africana")  # update names of tips in the Burleigh et al. tree that don't match names in my taxonomy

treedist.RAxML_burleigh <- treedist(drop.tip(picidae.RAxML.all.no_outgroups, tip=setdiff(picidae.RAxML.all.no_outgroups$tip.label, supertree.burleigh.picidae.taxnames_corrected$tip.label)), drop.tip(supertree.burleigh.picidae.taxnames_corrected, tip=setdiff(supertree.burleigh.picidae.taxnames_corrected$tip.label, picidae.RAxML.all.no_outgroups$tip.label)))  # calculate tree distance metrics between Burleigh et al. supermatrix tree and my RAxML tree


## save an R data object with the trees I need for later use
save(file="R_trees_for_combined_analyses.RData", list=c("picidae.RAxML.all.rooted.ladderized", "picidae.RAxML.all.r8s_calibrated.ladderized", "picidae.RAxML.all.r8s_calibrated.no_calibs.ladderized", "picidae.mrbayes.all.ladderized", "picidae.RAxML.all.BEAST_calibrated", "picidae.RAxML.all.BEAST_calibrated.no_calibs"))
