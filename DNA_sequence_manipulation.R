## This is is one of several files containing scripts and functions used in processing and analysis of data for Matthew Dufort's Ph.D. dissertation at the University of Minnesota, titled "Coexistence, Ecomorphology, and Diversification in the Avian Family Picidae (Woodpeckers and Allies)."

# this file contains functions and scripts for trimming and concatenating sequence alignments
# it also contains scripts for testing for base compositional heterogeneity

library(ape)

### generate functions for later use

## the function trimAlignmentToReference() trims an alignment so that the sequences do not extend beyond a reference sequence; this ensures that loci don't overlap due to sequences extending across loci
# as input, it takes: x (an alignment as a list of sequences), reference (either the name of a taxon (in which case ref.as.taxon should be TRUE), or the name of a sequence (ref.as.taxon should be FALSE)), and optional arguments to change the characters used for gaps and missing data
# it returns a list containing the trimmed alignment as a list of sequences, and a vector of the alignment positions removed
trimAlignmentToReference <- function(x, reference, ref.as.taxon=TRUE, gap="-", missing="?") {
  
  # find the reference taxon in the list of sequences
  if (!ref.as.taxon) {
    refnum <- which(names(x)==reference)
  } else {
    x.taxon.names <- trimTaxnamesAlign(x)
    refnum <- which(names(x.taxon.names)==reference)
  }
  if (length(refnum) != 1) {
    print(paste("Error: reference name matches",length(refnum),"sequences in the file."), sep="")
  }
  
  bounds <- range(which(!(x[[refnum]] %in% c(gap, missing)))) # find the first and last non-gap non-missing characters in the reference sequence
  
  # create a vector containing the positions of to be trimmed
  trim.head <- integer(0)
  if (bounds[1] > 1) {
    trim.head <- 1:(bounds[1]-1)
  }
  trim.tail <- integer(0)
  if (bounds[2] < length(x[[refnum]])) {
    trim.tail <- (bounds[2]+1):length(x[[refnum]])
  }
  exset <- c(trim.head, trim.tail)
  
  # remove the portions of the alignment that extend beyond the reference sequence
  x.trimmed <- removeExsetAlign(x, exset)
  
  return(list(data=x.trimmed, chars.removed=exset))
}

## the function trimTaxnamesAlign() trims sequence names in alignments (as lists of sequences) down to just the taxon name, without locus or GenBank number info
# it assumes that the taxon names are two words long, with underscores in between them; this can be modified using optional arguments sep and nwords
# as input, it takes the nexus alignment (as a list of sequences)
trimTaxnamesAlign <- function(x, sep="_", nwords=2) {
  find.underscores <- gregexpr(sep, names(x)) # find word boundaries in taxon names
  
  for (i in 1:length(names(x))) {  # loop over taxon names
    if (!is.na(find.underscores[[i]][nwords])) names(x)[i] <- strtrim(names(x)[i],find.underscores[[i]][nwords]-1)  # trim extra words from taxon names as necessary
  }
  return(x)  # return the modified alignment (as a list of sequences)
}

## the function namesEvenAlign makes the names of the taxa in a nexus alignment equal in length by adding spaces, so that they can be output using my writeNexusAlign function
# as input, it takes: x (an alignment as a list of sequences), buffer (the number of additional spaces to put between the end of the name and the beginning of the sequence)
# it returns a modified alignment
namesEvenAlign <- function(x, buffer=5) {

  max.length <- max(nchar(names(x)))  # get the length of the longest sequence name
  for (i in 1:length(names(x))) {  # loop over sequence names
    names(x)[i] <- paste(names(x)[i],paste(rep.int(' ',max.length-nchar(names(x)[i])+buffer),sep='',collapse=''))  # add spaces to end of sequence to match length of longest sequence plus a buffer
  }
  
  return(x)
}

## the function writeNexusAlign() outputs an alignment to a nexus-formatted file
# as input, it takes an alignemnt (as a list of sequences), a filename for output, a buffer length (for spacing between sequence names and sequences), and characters used for gaps and missing data (necessary as these are included in the nexus format)
# it checks for partitions in the alignment object, and includes them in the output file if present
writeNexusAlign <- function(x, filename='output.nex', buffer=5, gap="-", missing="?") {
  require(stringr)
  
  # check for partitions in the alignment; if present, split the partitions and sequence data
  if (names(x)[2]=="partitions") (partitioned <- TRUE) else (partitioned <- FALSE)
  if (partitioned) {
    x.partitions <- x$partitions
    x <- x$data
  }
  
  # trim and replace any excess whitespace after the names
  names(x) <- str_trim(names(x),side="right")
  x <- namesEvenAlign(x, buffer=buffer)  # homogenize the length of the sequence names to simplify output 
  
  # output sequences to file
  sink(filename)
  cat("#NEXUS\r\rBegin DATA;\r")
  cat(paste("\tDimensions ntax=",length(names(x))," nchar=",length(x[[1]]),";\r",sep=""))
  cat(paste("\tFormat datatype=NUCLEOTIDE gap=",gap," missing=",missing,";\r\tMatrix\r\r",sep=""))
  for (i in 1:length(x)) {
    cat(paste("\t",names(x)[i],paste(x[[i]],sep="",collapse=""),"\r",sep=""))
  }
  cat(";\rEnd;\r")
  
  # if alignment includes partitions, write partitions to a sets block
  if (partitioned) {
    # code to write out sets block
    cat("\rBegin sets;\r")
    for (i in 1:nrow(x.partitions)) {
      cat(paste("\tcharset ",rownames(x.partitions)[i]," = ",x.partitions[i,1],"-",x.partitions[i,2],";\r",sep=""))
    }
    cat("End;\r")
  }
  
  sink()
}

## the function trimLowOverlapAlign() removes positions from an alignment where fewer than a minimum number of taxa have data at that position
# as input, it takes x: (the alignment as a list of sequences), overlap.min (the minimum number of taxa required to retain a position in the alignment), and characters for gaps, missing data, and ambiguous bases
# it returns a modified alignment (as a list) and a vector of the positions removed
trimLowOverlapAlign <- function(x, overlap.min=4, gap="-", missing="?", anybase=c("n","N")) {
  require(plyr)
  require(matrixStats)
  
  x.matrix <- do.call(rbind,x)  # convert list of sequences to a matrix for efficiency
  all.remove <- c(gap,missing,anybase)  # generate a complete vector of characters not counted as data
  x.dims <- c(nrow(x.matrix), ncol(x.matrix)) # store the dimensions of the matrix (necessary because %in% converts matrix to vector)
  
  counts <- colCounts(matrix(!(x.matrix %in% all.remove), nrow=x.dims[1], ncol=x.dims[2], byrow=FALSE))  # count the number of sequences with data at each alignment position
  exset <- which(counts < overlap.min)  # identify the positions where the count of sequences with data is less than overlap.min
  x.matrix <- x.matrix[,-exset]  # remove positions (columns) where the count of sequences with data is less than overlap.min
  
  x <- alply(x.matrix,1)  # convert data matrix back to a list
  names(x) <- rownames(x.matrix)  # add sequence names to the alignment
  
  return(list(data=x, removed=exset))
}

## the function removeExsetAlign() removes positions from an alignment (as a list of sequences) and returns a new alignment without those characters
# as input, it takes: x (the alignment as a list of sequences), exset (a vector of the positions to remove)
removeExsetAlign <- function(x, exset) {
  require(plyr)
  if (length(exset) > 0) {  # check for null exset
    x.matrix <- do.call(rbind,x)  # convert list to a matrix for efficiency
    x.matrix <- x.matrix[,-exset]  # strip out columns that match the exclusion set
    x <- alply(x.matrix,1)  # convert data matrix back to a list
    names(x) <- rownames(x.matrix)  # add sequence names to the alignment
  }
  return(x)
}


## the function subsetAlign() trims an alignment to include only the positions within a specified range
# as input, it takes: x (an alignment as a list), range (a numeric vector containing the positions to be removed), and optional arguments rm.empty (boolean determining whether rows with no data remaining are removed), and characters for gaps and missing data
# it returns a modified alignment
subsetAlign <- function(x, range, rm.empty=TRUE, gap="-", missing="?") {
  require(matrixStats)
  require(plyr)
  
  x.matrix <- do.call(rbind,x) # convert list to a matrix efficiency
  x.matrix <- x.matrix[,range] # reduce matrix to only include specified range
  
  # remove sequences with only gaps and/or missing data
  x.dims <- c(nrow(x.matrix), ncol(x.matrix))  # store the dimensions of the matrix (necessary because %in% converts matrix to vector)
  counts.data <- rowCounts(matrix(!(x.matrix %in% c(gap, missing)), nrow=x.dims[1], ncol=x.dims[2], byrow=FALSE))  # calculate counts of data for each row
  x.matrix <- x.matrix[counts.data>0,]  # strip out rows with no non-gap non-missing data
  
  x <- alply(x.matrix,1) # convert data matrix back to a list
  names(x) <- rownames(x.matrix) # add sequence names to the alignment
  
  return(x)
}

## the function concatenateAlign() merges two or more alignments with some overlapping taxa (partial or full taxon overlap), replacing any empty data with the specific character for missing data
# as input, it takes x (a list of alignments, each one a list of sequences), and missing (a character to be used for missing data)
# it returns a list containing a single concatenated alignment and a vector with the start and end points of each alignment included in the concatenated data set
concatenateAlign <- function(x, missing="?") {
  require(stringr)
  
  x.matrices <- lapply(x,do.call,what=rbind)  # convert alignments to matrices
  
  # start with first alignment in the list
  x.matrices[[1]] <- apply(x.matrices[[1]],1,str_c,sep="",collapse="")  # collapse each sequence from a vector to a character string
  
  x.merged <- data.frame(taxon=names(x.matrices[[1]]), sequence=x.matrices[[1]],row.names=NULL,stringsAsFactors=FALSE)  # create a data frame to store the sequences, starting with the first alignment; it's easier to manipulate the elements in the data frame
  
  # intialize vectors to store the start and end of each partition
  starts <- rep(NA,length(x)); ends <- rep(NA,length(x))
  starts[1] <- 1
  ends[1] <- nchar(x.merged$sequence[1])
  
  for (i in 2:length(x)) {  # loop over remaining alignments
    x.matrices[[i]] <- apply(x.matrices[[i]],1,str_c,sep="",collapse="")  # collapse each sequence from a vector to a character string
    x.df.temp <- data.frame(taxon=names(x.matrices[[i]]),sequence=x.matrices[[i]],row.names=NULL,stringsAsFactors=FALSE)  # create a temporary data frame with the current alignment
    x.merged <- merge(x.merged, x.df.temp, by="taxon", all=TRUE)  # merges the two data frames, matching the elements by taxon; this adds a sequence.x (concatenated alignment to this point) and a sequence.y (current alignment to be added)
    x.merged$sequence.x[is.na(x.merged$sequence.x)] <- paste(str_c(rep("?",max(nchar(x.merged$sequence.x))),sep="",collapse=""))  # replace NAs in concatenated aligment with string of ?
    x.merged$sequence.y[is.na(x.merged$sequence.y)] <- paste(str_c(rep("?",max(nchar(x.merged$sequence.y))),sep="",collapse=""))  # replace NAs in to-be-added alignment with string of ?
    x.merged$sequence <- with(x.merged, paste(sequence.x,sequence.y,sep=""))  # concatenate sequence.x and sequence.y into sequence for each row
    x.merged <- subset(x.merged,select=c(taxon,sequence))  # remove sequence.x and sequence.y
    starts[i] <- ends[i-1] + 1  # record start of partition for current alignment to be added
    ends[i] <- nchar(x.merged$sequence[1])  # record start of partition for current alignment to be added
  }
  
  x.merged.list <- data.matrix(x.merged$sequence)  # convert the data frame to a matrix
  x.merged.list <- strsplit(x.merged.list,split="")  # split the matrix into a list
  names(x.merged.list) <- x.merged$taxon  # add the taxon names to the elements in the list
  partitions <- cbind(starts,ends)
  rownames(partitions) <- names(x)  # add partition names
  
  return(list(data=x.merged.list, partitions=partitions))
}



## the function quantifyBasesNonmissing quantifies the number of non-missing bases in a list of sequences
# as input, it takes x (an alignment as a list), and characters for gaps, missing data, and ambiguous data
# it returns a list containing a vector with the number of positions with data for each taxon (with the items in the vector named by the taxon), and the total length of the alignment
quantifyBasesNonmissing <- function(x,gap="-",missing="?",anybase=c("n","N")) {
  require(matrixStats)
  
  if (names(x)[2]=="partitions") x <- x$data  # remove partitions block if present in the alignment
  x.matrix <- do.call(rbind,x)  # convert list to a matrix for efficiency
  x.counts <- numeric(length=length(x))  # create a vector to store counts by taxon
  all.exclude <- c(gap, missing, anybase)  # create a character vector with all the characters to ignore in counting
  
  x.counts <- rowCounts(matrix(!(x.matrix %in% all.exclude),nrow=nrow(x.matrix),ncol=ncol(x.matrix),byrow=FALSE))  # count the positions in each row that don't match any of the characters to exclude
  names(x.counts) <- names(x)  # add sequence names to the elements in the vector
  
  return(list(counts=x.counts, bp.total=ncol(x.matrix)))
}

## the function printAlignmentMatrix() outputs an alignment as a figure with taxa as rows and subsets of the alignment as columns; filled pixels indicate data present, unfilled indicate missing data
# as input, it takes x (the alignment as a list of sequences), numpixels (the width of the plot in number of pixels), partitions (an optional matrix of partition start and end positions) and missing (characters to consider as data not present)
printAlignmentMatrix <- function(alignment, numpixels, partitions=NA, missing=c("-","?","n","N")) {
  if ((names(x)[2]=="partitions") & (!is.na(partitions))) {
    print("Error: partitions argument is non-null and partitions were found in alignment.")
    return()
  }
  
  if ((names(x)[2]=="partitions") & (is.na(partitions))) {
    alignment <- alignment$data
    partitions <- alignment$partitions
  }
  
  alignment.matrix <- do.call(rbind, alignment)  # convert alignment to matrix
  bin.size <- ncol(alignment.matrix) %/% numpixels  # calculate the rounded bin size (in characters of the alignment)
  bin.remainder <- ncol(alignment.matrix) %% numpixels  # calculate the remainder after the bins (in characters of the alignment)
  ends <- rep(NA, numpixels)  # set up a vector for the end-point of each bin
  bin.is.large <- rep(FALSE, numpixels)  # set up a vector to store whether a given bin includes an additional position (so the remainder gets included)
  bin.is.large[floor((1:bin.remainder) * numpixels / bin.remainder)] <- TRUE  # set an equally spaced number of bins to have an additional position included
  
  for (i in 1:numpixels) {  # loop over the bins
    larges <- sum(bin.is.large[1:i])  # calculate the number of "oversized" bins
    ends[i] <- (bin.size * (i - larges)) + ((bin.size +1 ) * larges)  # calculate the end position for ech bin
  }
  
  alignment.matrix.numeric <- matrix(as.numeric(rep(NA, nrow(alignment.matrix) * ncol(alignment.matrix))), nrow=nrow(alignment.matrix), ncol=ncol(alignment.matrix))  # initialize empty matrix to store data present/absent information by position for entire alignment
  alignment.matrix.numeric[!(alignment.matrix %in% missing)] <- 1  # set positions with data present in alignment to 1 in numeric matrix
  alignment.matrix.numeric[alignment.matrix %in% missing] <- 0  # set positions with data absent in alignment to 0 in numeric matrix
  
  alignment.matrix.reduced <- matrix(rep(NA, nrow(alignment.matrix) * numpixels), nrow=nrow(alignment.matrix), ncol=numpixels)  # initialize empty matrix to store data present/absent information by bin
  for (i in 1:nrow(alignment.matrix.reduced)) {  # loop over rows
    alignment.matrix.reduced[i,1] <- max(alignment.matrix.numeric[i,1:ends[1]])  # determine data present/absent for first bin
    for (j in 2:ncol(alignment.matrix.reduced)) {  # loop over bins
      alignment.matrix.reduced[i,j] <- max(alignment.matrix.numeric[i,((ends[j-1]+1):ends[j])])  # determine data present/absent for current bin
    }
  }
  
  image(t(alignment.matrix.reduced)[,nrow(alignment.matrix.reduced):1], col=c("white","black"), xaxt="n", yaxt="n")  # plot the reduced matrix
  
  # if partitions are present, add an axis with tick marks at the boundaries of the partitions
  if (!is.na(partitions)) {
    ticks <- c(0, ((1 / ncol(alignment.matrix)) * alignment$partitions[,2]))
    axis(1, at=ticks, labels=FALSE)
  }
}


  
### read in locus-level alignments and concatenate them into a single alignment

locus.filenames <- read.table(file="locus_filenames_for_import.txt", header=T, sep=",")  # read in the locations of files containing the locus-level alignments, as Nexus alignment files

mito.nums <- grep("mito.*", locus.filenames$locus)  # identify the mitochondrial loci (using the fact that I used a mito- prefix when naming those loci)
autosomal.names <- c("CMOS","BFIBi5","FGBi7","GAPDHi11","HMGN2","IRBP","LDH","MYOi2","PEPCKi9","PERi9","RAG1","TGFBi5")  # store the names of autosomal loci
zlinked.names <- c("ACO1i9","BRMi15","MUSKi4")  # store the names of autosomal loci

Picidae.sequences <- list()  # initialize an object to store alignments
for (i in 1:nrow(locus.filenames)) {  # loop over alignment filenames
  Picidae.sequences[[i]] <- read.nexus.data(file=locus.filenames$filename[i])  # read the current alignment file
  names(Picidae.sequences)[i] <- as.character(locus.filenames$locus[i])  # name the list element containing the current alignment
}

Picidae.sequences.trimmed <- Picidae.sequences  # create a copy of the list of alignments, to trim the names

# loop over alignments, trimming the sequence names so that they only include genus_species
for(i in 1:length(Picidae.sequences)) {
  Picidae.sequences.trimmed[[i]] <- trimTaxnamesAlign(trimLowOverlapAlign(Picidae.sequences[[i]])$data)
}

# loop over mitochondrial alignments, trimming the alignments to the Dryocopus_pileatus reference
for (i in mito.nums) {
  Picidae.sequences.trimmed[[i]] <- trimAlignmentToReference(Picidae.sequences.trimmed[[i]], reference="Dryocopus_pileatus")$data
}

## trim problematic alignment regions (identified by visual inspection of alignments)
exset.FGBi7 <- c(970:979)
Picidae.sequences.trimmed$FGBi7 <- removeExsetAlign(Picidae.sequences.trimmed$FGBi7, exset.FGBi7)
exset.mito12S <- c(79:86, 122:127, 222:226, 317:325, 374:394, 483:495, 625:629, 792:797, 892:934)
Picidae.sequences.trimmed$mito12S <- removeExsetAlign(Picidae.sequences.trimmed$mito12S, exset.mito12S)
exset.mito16S <- c(171:175)
Picidae.sequences.trimmed$mito16S <- removeExsetAlign(Picidae.sequences.trimmed$mito16S, exset.mito16S)
exset.PERi9 <- c(59:73)
Picidae.sequences.trimmed$PERi9 <- removeExsetAlign(Picidae.sequences.trimmed$PERi9, exset.PERi9)

# output trimmed alignments to Nexus alignment files
for (i in 1:length(Picidae.sequences.trimmed)) {
  writeNexusAlign(Picidae.sequences.trimmed[[i]],filename=paste(names(Picidae.sequences.trimmed)[i], ".nex", sep=""))
}


Picidae.mito.all.concatenated <- concatenateAlign(Picidae.sequences.trimmed[mito.nums])  # generate concatenated alignment of mitochondrial loci
Picidae.nuclear.all.concatenated <- concatenateAlign(Picidae.sequences.trimmed[-mito.nums])  # generate concatenated alignment of nuclear loci
Picidae.nuclear.all_ex_FGBi7.concatenated <- concatenateAlign(Picidae.sequences.trimmed[-c(mito.nums, which(names(Picidae.sequences.trimmed)=="FGBi7"))])  # generate concatenated alignment of nuclear loci excluding FGB intron 7
Picidae.autosomal.all.concatenated <- concatenateAlign(Picidae.sequences.trimmed[autosomal.names])  # generate concatenated alignment of autosomal loci
Picidae.zlinked.all.concatenated <- concatenateAlign(Picidae.sequences.trimmed[zlinked.names])  # generate concatenated alignment of z-linked loci
Picidae.all.concatenated <- concatenateAlign(c(Picidae.sequences.trimmed[mito.nums],Picidae.sequences.trimmed[-mito.nums]))  # generate concatenated alignment of all loci, with mitochondrial first, then nuclear
Picidae.all.concatenated.chromosome_ordered <- concatenateAlign(c(Picidae.sequences.trimmed[mito.nums],Picidae.sequences.trimmed[autosomal.names],Picidae.sequences.trimmed[zlinked.names]))  # generate concatenated alignment of all loci, with mitochondrial first, then autosomal, then z-linked
Picidae.all_ex_FGBi7.concatenated <- concatenateAlign(c(Picidae.sequences.trimmed[mito.nums],Picidae.sequences.trimmed[-c(mito.nums, which(names(Picidae.sequences.trimmed)=="FGBi7"))]))  # generate concatenated alignment of all loci excluding FGB intron 7

## output concatenated alignments to Nexus alignment files
writeNexusAlign(Picidae.mito.all.concatenated, filename="Picidae_align_mito_all_concatenated.nex")
writeNexusAlign(Picidae.nuclear.all.concatenated, filename="Picidae_align_nuclear_all_concatenated.nex")
writeNexusAlign(Picidae.nuclear.all_ex_FGBi7.concatenated, filename="Picidae_align_nuclear_all_ex-FGBi7_concatenated.nex")
writeNexusAlign(Picidae.all.concatenated, filename="Picidae_align_all_concatenated.nex")
writeNexusAlign(Picidae.all_ex_FGBi7.concatenated, filename="Picidae_align_all_ex_FGBi7_concatenated.nex")

## calculate and summarize data by alignment
coverage.Picidae.mito.all.concatenated <- quantifyBasesNonmissing(Picidae.mito.all.concatenated)  # quantify data present by taxon
summary(coverage.Picidae.mito.all.concatenated$counts / coverage.Picidae.mito.all.concatenated$bp.total)  # summarize coverage
coverage.Picidae.nuclear.all.concatenated <- quantifyBasesNonmissing(Picidae.nuclear.all.concatenated)  # quantify data present by taxon
summary(coverage.Picidae.nuclear.all.concatenated$counts / coverage.Picidae.nuclear.all.concatenated$bp.total)  # summarize coverage
coverage.Picidae.all.concatenated <- quantifyBasesNonmissing(Picidae.all.concatenated)  # quantify data present by taxon
summary(coverage.Picidae.all.concatenated$counts / coverage.Picidae.all.concatenated$bp.total)  # summarize coverage


### plot coverage by taxon

printAlignmentMatrix(Picidae.all.concatenated, numpixels=1000)
printAlignmentMatrix(Picidae.all.concatenated.chromosome_ordered, numpixels=1000)


# in order to get the taxa to print in the order they are on the tree
picidae.RAxML <- read.tree("picidae_RAxML.tre")
Picidae.all.concatenated.sorted_RAxML_tree <- Picidae.all.concatenated
Picidae.all.concatenated.sorted_RAxML_tree$data <- Picidae.all.concatenated$data[picidae.RAxML$tip.label[picidae.RAxML$edge[picidae.RAxML$edge[,2]<=(length(picidae.RAxML$tip.label)),2]]] # this re-orders the data matrix to match the order of the tips as they're displayed on the tree

# in order to get the taxa to print in the order they are on the tree
Picidae.all.concatenated.chromosome_ordered.sorted_RAxML_tree <- Picidae.all.concatenated.chromosome_ordered
Picidae.all.concatenated.chromosome_ordered.sorted_RAxML_tree$data <- Picidae.all.concatenated.chromosome_ordered$data[picidae.RAxML$tip.label[picidae.RAxML$edge[picidae.RAxML$edge[,2]<=(length(picidae.RAxML$tip.label)),2]]] # this re-orders the data matrix to match the order of the tips as they're displayed on the tree


pdf(file="Picidae_locus_coverage_plot_all_chromosome_ordered_RAxML_sorted.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(picidae.RAxML, show.tip.label=FALSE, x.lim=c(0,120000000), edge.width=0.75)
printAlignmentMatrix(Picidae.all.concatenated.chromosome_ordered.sorted_RAxML_tree, numpixels=1000)
dev.off()


### test for base compositional heterogeneity

## the function basefreqsChisqTest() conducts a chi-square test of each taxon's base frequencies against the average base frequencies for that locus
# as input, it takes data (the frequencies by taxon) and props (the average frequencies for that locus)
# it returns a list containing a list of the chi square test results and a vector with the p values from those test
basefreqsChisqTest <- function(data, props) {
  # takes a matrix as the data
  test.results <- list()  # initialize a list to store test results
  p.values <- numeric(nrow(data))  # create a vector to store p values
  for (i in 1:nrow(data)) {  # loop over taxa, running the test and storing the p value
    test.results[[i]] <- chisq.test(x=round(data[i,]), p=props)
    p.values[i] <- test.results[[i]]$p.value
  }
  return(list(test.results=test.results, p.values=p.values))
}

## loop through each locus and each taxon within each locus, running an exact multinomial test (or a chi-square test if the sample sizes are larger or if the multinomial test takes too long) and storing the results somewhere; could store these as a list, but I should also store just the chi-square statistic values and p-values so that I can access them more easily

basefreqs.filenames <- read.table(file="basefreqs_filenames.txt", header=T, sep=",", stringsAsFactors=FALSE)  # read list of names of files containing base frequencies by taxon for each locus

Picidae.basefreqs <- list()  # initialize a list to store base frequencies
for (i in 1:nrow(basefreqs.filenames)) {  # loop over loci
  Picidae.basefreqs[[i]] <- list()
  names(Picidae.basefreqs)[i] <- as.character(basefreqs.filenames$locus[i])
  Picidae.basefreqs[[i]][[1]] <- read.csv(file=basefreqs.filenames$filename[i], stringsAsFactors=FALSE)  # read the base frequencies by taxon for each locus
  Picidae.basefreqs[[i]][[2]] <- read.csv(file=paste("/Users/MattDufort/Documents/Grad_school/Dissertation_stuff/Picidae_supermatrix_analyses/Picidae_base_composition_tests/Picidae_locus-by-locus/", basefreqs.filenames$basefreqs_filename[i], sep=""), header=FALSE)  # read the average base frequencies for each locus
  names(Picidae.basefreqs[[i]]) <- c("observed", "averages")
}

# loop over loci, doing a chi square test of the empirical frequencies vs. the average
for (i in 1:length(Picidae.basefreqs)) {
  chisq.tmp <- basefreqsChisqTest(as.matrix(Picidae.basefreqs[[i]]$observed[,2:5]), (unlist(Picidae.basefreqs[[i]]$averages) / sum(unlist(Picidae.basefreqs[[i]]$averages))))
  Picidae.basefreqs[[i]][[3]] <- chisq.tmp$test.results
  Picidae.basefreqs[[i]][[4]] <- chisq.tmp$p.values
}

for (i in 1:length(Picidae.basefreqs)) {
  print(min(Picidae.basefreqs[[i]][[4]]))
}

# do an exact multinomial test for Jynx torquilla COI
multinomial.test(c(388,438,241,446), p=c(0.25615, 0.33515, 0.16199, 0.24671))
