# This file is part of RNAseqWrapper.
# 
# RNAseqWrapper is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# RNAseqWrapper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# See <http://www.gnu.org/licenses/> for a a copy of the GNU General Public License.

#'@title a global variable for plotting control.
#'@note using quartz on MacOSX disables the lzw compression of tiffs - this variable changes the default to "Xlib" (only on MacOS)
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}
GLOBAL_VARIABLE_TIFF_LIB <- getOption("bitmapType"); if (Sys.info()['sysname'] == "Darwin") { GLOBAL_VARIABLE_TIFF_LIB <- "Xlib" }

#'@title a global variable for plotting control.
#'@note set to true if you wish some plots to be saved as svg instead of tiff/png.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}
GLOBAL_VARIABLE_USE_SVG <- TRUE; if (Sys.info()['sysname'] == "Darwin") { GLOBAL_VARIABLE_USE_SVG <- FALSE }

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title read tables with identical columns from a directory
#'@param aDir the path to the files in "aFiles"
#'@param aExp a regular expression required within the file name
#'@param aExt a file extension required for the file to be read
#'@param aExc a regular expression which should be omitted in the file names
#'@param aHead a vector with the column names
#'@param aSep the separator used in the files
#'@param useHeader set to TRUE if tables contain column names. WARNING: aHead is still required.
#'@return a list of matrices containing values from all files - so for every header/datatype one matrix
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.read.files.from.directory <- function(aDir, aExp, aExt, aExc, aHead, aSep, useHeader = FALSE) {
  # WARNING useHeader still neads aHead (it just tells if the first row shall be removed or not)
  # read all the files from a given directory that contain a regular expression aExp and the ending aExt. All with aExc will be excluded
  # aHead is a vector with the column names, aSep is the separator
  fList <- grep(aExp, list.files(aDir, aExt), value = TRUE)
  if (aExc != "") { fList <- fList[-grep(aExc,fList)] }
  sList <- sapply(fList, function(x) substr(x, 1, nchar(x)-nchar(aExt)))
  # read first all the files in a list and get a list of all the rownames (RID)
  RID <- c()
  tList <- list()
  for (sampleName in sList) {
    if (useHeader) {tList[[sampleName]] <- read.table(file.path(aDir,paste(sampleName, aExt, sep="")), header = TRUE, sep = aSep, quote = "", row.names = 1)}
    else {tList[[sampleName]] <- read.table(file.path(aDir,paste(sampleName, aExt, sep="")), header = FALSE, sep = aSep, quote = "", row.names = 1)}
    colnames(tList[[sampleName]]) <- aHead
    RID <- union(RID, rownames(tList[[sampleName]]))
  }
  # create a list of matrices containing values from all samples - for every datatype one matrix
  out <- list()
  for (h in aHead) {
    out[[h]] <- matrix(0, length(RID), length(sList), FALSE, list(RID, sList))
    for (sampleName in sList) {
      out[[h]][rownames(tList[[sampleName]]), sampleName] <- tList[[sampleName]][,h]
    }
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title read tables with identical columns from a given list of files
#'@param aDir the path to the files in "aFiles"
#'@param aFiles a vector with file names
#'@param aHead a vector with the column names
#'@param aSep the separator used in the files
#'@param useHeader set to TRUE if tables contain column names. WARNING: aHead is still required.
#'@return a list of matrices containing values from all files - so for every header/datatype one matrix
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.read.files.from.given.list <- function(aDir, aFiles, aHead, aNames = c(), aSep = '\t', useHeader = FALSE) {
  # WARNING useHeader still neads aHead (it just tells if the first row shall be removed or not)
  # read all the given files from a given directory
  # aHead is a vector with the column names, aSep is the separator
  # read first all the files in a list and get a list of all the rownames (RID)
  RID <- c()
  temp <- list()
  for (fn in aFiles) {
    if (useHeader) { temp[[fn]] <- read.table(file.path(aDir, fn), header = TRUE, sep = aSep, quote = "", row.names = 1) }
    else { temp[[fn]] <- read.table(file.path(aDir, fn), header = FALSE, sep = aSep, quote = "", row.names = 1) }
    colnames(temp[[fn]]) <- aHead
    RID <- union(RID, rownames(temp[[fn]]))
  }
  # create a list of matrices containing values from all files - so for every header/datatype one matrix
  out <- list()
  for (h in aHead) {
    out[[h]] <- matrix(0, length(RID), length(aFiles), FALSE, list(RID, aFiles))
    for (fn in aFiles) {
      out[[h]][rownames(temp[[fn]]), fn] <- temp[[fn]][,h]
    }
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title read simple "list-files" (one column) from a directory
#'@param aDir a directory with the files
#'@param aExp a regular expression required within the file name
#'@param aExt a file extension required for the file to be read
#'@param aExc a regular expression which should be omitted in the file names
#'@return a list of vectors read from the files
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.read.lists.from.directory <- function(aDir, aExp, aExt, aExc) {
  # 
  fList <- grep(aExp, list.files(aDir, aExt), value = TRUE)
  if (aExc != "") { fList <- fList[-grep(aExc,fList)] }
  sList <- sapply(fList, function(x) substr(x, 1, nchar(x)-nchar(aExt)))
  # read all the files in a list with the fileName as key to the list
  out <- list()
  for (fn in sList) {
    out[[fn]] <- scan(file.path(aDir, paste(fn, aExt, sep = '')), what = "character")
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title read Rcount tables
#'@param aDir a directory with the files
#'@param aExp a regular expression required within the file name
#'@param aExc a regular expression which should be omitted in the file names
#'@param toReturn the column of interest (normally TH)
#'@return a data.frame with Rcount gene expression values
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.read.Rcount <- function(aDir, aExp = "", aExc = "", toReturn = "TH") {
  fullData <- f.read.files.from.directory(aDir, aExp, ".txt", aExc, c("sumUnAmb", "sumAmb", "sumAllo", "sumDistUnAmb", "sumDistAmb", "sumDistAllo", "TH", "VAL"), "\t")
  out <- round(fullData[[toReturn]])
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title read featureCounts tables
#'@param aDir a directory with the files
#'@param aExp a regular expression required within the file name
#'@param aExc a regular expression which should be omitted in the file names
#'@return a data.frame with featureCounts expression values
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.read.featureCounts <- function(aDir, aExp = "", aExc = "") {
  fullData <- f.read.files.from.directory(aDir, aExp, ".txt", aExc, c("counts"), "\t", TRUE)
  out <- fullData$counts
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Check whether a FASTQ file has a PHRED offset of 33 or 64
#'@param fastq a fastq file (may be gzipped)
#'@return the PHRED offset: either 33 or 64
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.check.phred.offset <- function(fastq) {
  # open a connection to the file
  isZipped <- (length(grep("*.gz$", fastq)) > 0)
  if (isZipped) {
    fq <- gzcon(file(fastq, "rb"))
  } else {
    fq <- file(fastq, open = "rt")
  }
  # the different formats and their ranges:
  # Sanger: 0 to 93 (rarely above 60) with ASCII 33 to 126
  # SOLiD: like Sanger ?
  # Illumina 1.8: like Sanger ?
  # Solexa/Illumina 1.0: -5 to 62 (normally only -5 to 40) with ASCII 59 to 126
  # Illumina 1.3: 0 to 62 (normally only 0 to 40) with ASCII 64 to 126
  # Illumina 1.5: Phred scores 0 to 2 have a slightly different meaning. The values 0 and 1 are no longer used and the value 2,
  # encoded by ASCII 66 "B", is used also at the end of reads as a Read Segment Quality Control Indicator.
  # check the quality range of the first 10'000 reads
  lines <- readLines(fq, 4e4)
  qualStrings <- lines[seq(4, 4e4, by = 4)]
  qualStrings <- gsub('\\', '\\\\', qualStrings, fixed = TRUE) # otherwise the backslash is an escape character
  ranges <- t(sapply(qualStrings, function(x) range(as.integer(charToRaw(x))), USE.NAMES = FALSE))
  lowest <- min(ranges[,1])
  highest <- max(ranges[,2])
  cat("FILE: ", fastq, "has quality values from", lowest, "to", highest, "\n")
  if (lowest < 58) {
    out <- 33
  } else {
    out <- 64
  }
  # close the file connection and return the PHRED offset
  close(fq)
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title summarize a table according to groups and names
#'@param x <dataframe, matrix>: numeric, colnames as in byTab
#'@param byTab <dataframe>: at least two columns: sample and group
#'@param summaryFunction: a function like mean, median, sum
#'@return the summarized table (a matrix)
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.summarize.columns <- function(x, byTab, summaryFunction) {
  byTab$sample <- as.character(byTab$sample)
  byTab$group <- as.character(byTab$group)
  groupnames <- unique(byTab$group)
  out <- matrix(0, nrow = nrow(x), ncol = length(groupnames), dimnames = list(rownames(x), groupnames))
  for (gn in groupnames) { 
    toSummarize <- byTab$sample[byTab$group == gn]
    if (length(toSummarize) > 1) {out[,gn] <- apply(x[,toSummarize], 1, summaryFunction)}
    else { out[,gn] <- x[,toSummarize] }
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title simple comparison function which makes either >= x or > 0
#'@param x values
#'@param rt threshold
#'@return logical vector
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.comp <- function(x,rt){
  if(rt==0){y <- x>0}
  else{y <- x>=rt}
  return(y)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title trim a table according to a column
#'@param x table
#'@param cn the column to trim
#'@param ext fraction of extreme values to remove (two-sided! 0.05 removes in total 10\% of all values)
#'@param viaRanks filter based on ranks - I don't remember anymore why I did this
#'@return table without the "ext" highest and the "ext" lowest values in the specified column
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.trim <- function(x, cn, ext = 0.05, viaRanks = FALSE) {
  if (viaRanks) {
    tot <- nrow(x)
    num_remove <- floor(tot*ext)
    min_rank <- num_remove
    max_rank <- tot-num_remove
    y <- rank(x[[cn]])
    out <- rownames(x)[(y > min_rank) & (y < max_rank)]
  } else {
    lowerBound <- quantile(data[,cn], ext)
    upperBound <- quantile(data[,cn], 1-ext)
    out <- rownames(data)(data[,cn] >= lowerBound & data[,cn] <= upperBound)
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title merge two matrices with different columns (rownames don't need to match completely)
#'@param x a matrix
#'@param y another matrix
#'@param rt a threshold - remove rows which are all below it (in case of 0, it tests for >, not >=)
#'@return merged matrix
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.merge.two.matrices <- function(x,y,rt){
  RID <- union(rownames(x), rownames(y))
  CID <- c(colnames(x), colnames(y))
  com <- matrix(0, nrow = length(RID), ncol = length(CID), dimnames = list(RID, CID))
  for (cname in colnames(x)) {
    com[rownames(x),cname] <- x[,cname]
  }
  for (cname in colnames(y)) {
    com[rownames(y),cname] <- y[,cname]
  }
  com.use <- rownames(com)[which(apply(f.comp(com,rt),1,sum)>0)]
  com <- com[com.use,]
  return(com)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title take the logarithm of the absolute value and append the original sign
#'@param x a vector
#'@param logFun a logarithm funtion
#'@return a vector with the log values
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.log.twosided <- function(x, logFun) {
  if (sum(is.na(x))>0) {
    cat("NA values replaced by 0\n")
    x[is.na(x)] <- 0
  }
  x[x>0] <- logFun(x[x>0])
  x[x<0] <- -logFun(-x[x<0])
  return(x)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title alternative transformation for count data
#'@param x a vector with cound data
#'@return a vector with the transformed values
#'@references 
#'TODO: ADD THE REFERENCE
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.stahel.trafo <- function(x) {
  med <- median(x, na.rm = TRUE)
  qua <- quantile(x, 0.25, na.rm = TRUE)
  if (qua == 0) {
    cat("f.stahel.trafo(): the first quartile equals to zero. Removing all zero counts to avoid -Inf.\n")
    temp <- x
    temp[temp==0] <- NA
    med <- median(temp, na.rm = TRUE)
    qua <- quantile(temp, 0.25, na.rm = TRUE)
  }
  con <- med/((med/qua)^2.9)
  out <- log(x + con)
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title calculate a GSEA-like cumulative sum
#'@param sortedGenes a named (genes) vector (e.g. logFC) - SORTED!
#'@param myGOIs genes of interest
#'@return cumulative sum (like in the GSEA plot)
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.calculate.GSEA.sum <- function(sortedGenes, myGOIs) {
  numGOIs <- length(myGOIs)
  numGenes <- length(sortedGenes)
  observed <- sortedGenes %in% myGOIs
  posStep <- (numGenes - numGOIs)/numGOIs
  forSum <- rep(-1, numGenes)
  forSum[observed] <- posStep
  if (abs(sum(forSum)) > 1) {cat(posStep*length(myGOIs), posStep, length(myGOIs))}
  cumSum <- cumsum(forSum)
  return(cumSum)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title calculate the maximal enrichment for a given gene set of interest
#'@param sortedGenes a named (genes) vector (e.g. logFC) - SORTED!
#'@param myGOIs genes of interest
#'@return maximal enrichment value (used to compare to random sampling)
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.calculate.max.GSEA.enrichment <- function(sortedGenes, myGOIs) {
  cumSum <- f.calculate.GSEA.sum(sortedGenes, myGOIs)
  cumSumMax <- cumSum[which.max(abs(cumSum))]
  return(cumSumMax)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Do a GSEA-like test
#'@param namedValues a sorted vector (e.g. logFC) with names (genes) vector (e.g. logFC) - SORTED!
#'@param GOIs a vector with the genes of interest
#'@param rDir directory where the plot is saved
#'@param filePrefix prefix for the name of the plot
#'@param mainTitle title in the plot
#'@param sortName a label vor the values (e.g., "logFC")
#'@param numReps number of random sets
#'@return a plot (saved as svg) and list with the pValue and enrichment scores
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.do.GSEA.like.test <- function(namedValues, GOIs, rDir, filePrefix = "", mainTitle = "", sortName = "LFC", numReps = 10000) {
  sortedGenes <- names(namedValues)
  ## check for existence of genes
  notInSorted <- setdiff(GOIs, sortedGenes)
  if (length(notInSorted) > 0) {f.print.message(filePrefix, "- removed", length(notInSorted), "/", length(GOIs), "GOIs - they were not in the set of sorted genes")}
  GOIs <- intersect(GOIs, sortedGenes)
  observedES <- f.calculate.max.GSEA.enrichment(sortedGenes, GOIs)
  observedSum <- f.calculate.GSEA.sum(sortedGenes, GOIs)
  sampledValues <- rep(0, numReps)
  for (i in 1:numReps) {
    sampledGenes <- sample(sortedGenes, length(GOIs))
    sampledValues[i] <- f.calculate.max.GSEA.enrichment(sortedGenes, sampledGenes)
  }
  ylimits <- c(-1000, 2000)
  ylimits <- c(min(observedES, sampledValues), max(observedES, sampledValues))
  pValue <- mean(abs(observedES) < abs(sampledValues))
  svgOutfile <- file.path(rDir, paste(filePrefix, "GSEA.svg", sep = ''))
  #pngOutfile <- file.path(rDir, paste(filePrefix, "GSEA.png", sep = ''))
  svg(svgOutfile, height = 7, width = 5)
  par(oma = c(2,2,1,1))
  layout(matrix(1:3, nrow = 3), heights = c(1, 0.2, 0.2))
  par(mar = c(0, 2, 1, 1))
  plot(observedSum, type = "l", bty = "n", xaxs = "r", xaxt = "n", yaxs = "r", xlab = "", ylab = "enrichment score", las = 1, main = mainTitle, ylim = ylimits)
  text(par("usr")[1] + par("usr")[2]/15, par("usr")[4] - par("usr")[4]/15, adj = c(0,1),labels = paste(c("P-value:", "n:"),c(round(pValue, digits = 3), length(GOIs)), collapse = "\n", sep=" "))
  plot(1:length(sortedGenes),  as.numeric(sortedGenes %in% GOIs), type = "h", bty = "n", xaxt = "n", yaxt = "n", xaxs = "r", yaxs = "r", xlab = "", ylab = "GOI-density", las = 1, main = "", ylim = c(0,1), col = "#00000064")
  par(mar = c(2, 2, 0, 1))
  plot(namedValues, type = "h", bty = "n", xaxs = "r", yaxs = "r", xlab = paste("genes sorted accoring to", sortName), ylab = sortName, las = 1, main = "")
  dev.off()
  #system(paste("rsvg-convert -a -d 300 -p 300 ", svgOutfile," > ", pngOutfile, sep = ''))
  out <- list(pValue = pValue, observedSum = observedSum, observedES = observedES, sampledES = sampledValues, observedESat = which.max(abs(observedSum)))
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title join a list of overlapping or adjacent genomic regions
#'@param data a data frame with at least three columns: chrom, start, and end
#'@param otherCols a character vector with column names which should be summarized as well using the summaryFunction
#'@param summaryFunction see "otherCols"
#'@param gap the maximal size of the gap between two fragments to be merged
#'@return a data frame with the joined fragments (cols are chrom, start, end, count, otherCols)
#'@note WARNING: THIS FUNCTION DOES NOT WORK IF FRAGMENTS ARE ENTIRELY WITHIN ANOTHER FRAGMENT - IT ASSUMES PARTIAL OVERLAP FROM ONE TO THE NEXT
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.join.overlapping.fragments <- function(data, otherCols, summaryFunction, gap = 0) {
  data <- data[with(data, order(chrom, start, end)),]
  data$regNum <- 1
  offset <- 1
  for (chr in unique(data$chrom)) {
    subData <- subset(data, chrom == chr)
    subData$regNum[1] <- offset
    if (nrow(subData) > 1) {
      temp <- (subData$start[2:nrow(subData)] - subData$end[1:(nrow(subData)-1)]) > gap
      subData$regNum[2:nrow(subData)] <- offset + cumsum(temp)
    } 
    data[data$chrom == chr, ] <- subData
    offset <- max(data$regNum) + 1
  }
  out <- data.frame(
    chrom = aggregate(data$chrom, by=list(region=data$regNum), unique, simplify = TRUE)$x,
    start = aggregate(data$start, by=list(region=data$regNum), min, simplify = TRUE)$x,
    end = aggregate(data$end, by=list(region=data$regNum), max, simplify = TRUE)$x,
    count = aggregate(data$chrom, by=list(region=data$regNum), length, simplify = TRUE)$x,
    stringsAsFactors = FALSE
  )
  for (cn in otherCols) {
    out[[cn]] = aggregate(data[[cn]], by=list(region=data$regNum), summaryFunction, simplify = TRUE)$x
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title join a list of overlapping or adjacent genomic regions - only if an additional variable has the same sign in both
#'@param data a data frame with at least three columns: chrom, start, and end
#'@param otherCols a character vector with column names which should be summarized as well using the summaryFunction
#'@param summaryFunction see "otherCols"
#'@param dirCol the name of the column which is used to decide wheter two fragments should be merged or not (based on the sign of the values)
#'@param gap the maximal size of the gap between two fragments to be merged
#'@return a data frame with the joined fragments (cols are chrom, start, end, count, otherCols)
#'@note WARNING: THIS FUNCTION DOES NOT WORK IF FRAGMENTS ARE ENTIRELY WITHIN ANOTHER FRAGMENT - IT ASSUMES PARTIAL OVERLAP FROM ONE TO THE NEXT
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.join.overlapping.fragments.directional <- function(data, otherCols, summaryFunction, dirCol, gap = 0) {
  data <- data[with(data, order(chrom, start, end)),]
  data$regNum <- 1
  offset <- 1
  for (chr in unique(data$chrom)) {
    subData <- subset(data, chrom == chr)
    subData$regNum[1] <- offset
    if (nrow(subData) > 1) {
      temp <- (subData$start[2:nrow(subData)] - subData$end[1:(nrow(subData)-1)]) > gap
      dirCheck <- sign(subData[[dirCol]][2:nrow(subData)]) != sign(subData[[dirCol]][1:(nrow(subData)-1)])
      subData$regNum[2:nrow(subData)] <- offset + cumsum(temp | dirCheck)
    } 
    data[data$chrom == chr, ] <- subData
    offset <- max(data$regNum) + 1
  }
  out <- data.frame(
    chrom = aggregate(data$chrom, by=list(region=data$regNum), unique, simplify = TRUE)$x,
    start = aggregate(data$start, by=list(region=data$regNum), min, simplify = TRUE)$x,
    end = aggregate(data$end, by=list(region=data$regNum), max, simplify = TRUE)$x,
    count = aggregate(data$chrom, by=list(region=data$regNum), length, simplify = TRUE)$x,
    stringsAsFactors = FALSE
  )
  for (cn in otherCols) {
    out[[cn]] = aggregate(data[[cn]], by=list(region=data$regNum), summaryFunction, simplify = TRUE)$x
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title search nearest genomic regions given a set of genomic regions
#'@param data a data frame with at least four columns: chrom, start, end, and strand (encoded as 1 and -1)
#'@param otherRegions a data frame with the regions to search in, three columns: chrom, start, and end
#'@return a data frame with the upstream and downstream distances
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.closest.region <- function(data, otherRegions) {
  f.get.min.dist <- function(x, y) {
    startDist <- c(x[1]-y$end, y$start-x[1])
    endDist <- c(x[2]-y$end, y$start-x[2])
    if (is.list(startDist)) {
      startDist <- c()
    } else {
      startDist <- startDist[startDist>0]
    }
    if (is.list(endDist)) {
      endDist <- c()
    } else {
      endDist <- endDist[endDist>0]
    }
    minStart <- min(startDist, na.rm = TRUE)
    minEnd <-  min(endDist, na.rm = TRUE)
    out <- c(minStart, minEnd)
    if (x[3] == (-1)) { out <- rev(out) }
    return(out)
  }
  out <- list()
  for (chr in unique(data$chrom)) {
    subDat <- subset(data, chrom == chr)
    subReg <- subset(otherRegions, chrom == chr)
    if (nrow(subReg) == 0) {
      cat("no regions on chromosome", chr, "\n")
      out[[chr]] <- matrix(Inf, ncol = 2, nrow = nrow(subDat), dimnames = list(rownames(subDat), NULL))
      next
    }
    if (nrow(subDat) > 1) {
      out[[chr]] <- t(apply(subDat[,c("start", "end", "strand")], 1, function(x) f.get.min.dist(x, subReg)))
    } else {
      out[[chr]] <- f.get.min.dist(subDat[1, c("start", "end", "strand")], subReg)
    }
  }
  out <- do.call("rbind", out)
  colnames(out) <- c("upStream", "downStream")
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title convert a dataframe from wide format to long format
#'@param tab a dataframe in wide format
#'@param cols the columns which should be merged
#'@param newName the name of the merged column
#'@return a data frame in long format
#'@note aside the column <newName> with the merged data, there will be another column
#'called "factor_<newName>". The latter specifies the origin of the data in the <newName>
#'column.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.wide.to.long <- function(tab, cols, newName) {
  outList <- list()
  otherCols <- setdiff(colnames(tab), cols)
  if (length(otherCols) > 0) {
    for (cn in cols) {
      if (length(otherCols) == 1) { 
        temp <- data.frame(t = tab[,otherCols], stringsAsFactors = FALSE)
        colnames(temp) <- otherCols
      } else {
        temp <- as.data.frame(tab[,otherCols], stringsAsFactors = FALSE)
      }
      outList[[cn]] <- temp
      outList[[cn]][[newName]] <- tab[,cn]
      outList[[cn]][[paste(newName, 'factor', sep = '_')]] <- rep(cn, nrow(tab))
    }
  } else {
    for (cn in cols) {
      outList[[cn]] <- data.frame(A = tab[,cn], B = rep(cn, nrow(tab)))
      colnames(outList[[cn]]) <- c(newName, paste(newName, 'factor', sep = '_'))
    }
  }
  out <- do.call("rbind",outList)
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title convert a dataframe from long format to wide format
#'@param tab a dataframe in long format
#'@param futureRows the column containing the rownames of the wide table
#'@param futureCols the column containing the rownames of the wide table
#'@param futureEntries the column containing the values
#'@return a data frame in wide format
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.long.to.wide <- function(tab, futureRows, futureCols, futureEntries) {
  empty <- ifelse(is.numeric(tab[[futureEntries]]), 0, "")
  rn <- unique(as.character(tab[[futureRows]]))
  cn <- unique(as.character(tab[[futureCols]]))
  out <- as.data.frame(matrix(empty, nrow = length(rn), ncol = length(cn), dimnames = list(rn, cn)))
  temp <- split(tab[,c(futureCols, futureEntries)], tab[[futureRows]])
  for (r in names(temp)) {
    out[r, temp[[r]][,futureCols]] <- temp[[r]][,futureEntries]
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.genetodensity <- function(x,y) {
  require("MASS")
  est <- kde2d(x, y, n = 50) # needs MASS - calculates the two dimensional density
  if (sum(is.na(est$z)) > 0) {
    est <- kde2d(x, y, n = 50, h = c(bw.nrd0(x), bw.nrd0(y))) # needs MASS - calculates the two dimensional density
  }
  dvec <- as.vector(est$z) # note: as.vector takes columnwise, rows correspond to the value of x, columns to the value of y - the y index counts therefore 50 times
  xind <- findInterval(x, est$x) # find the intervals to which a gene belongs
  yind <- findInterval(y, est$y) # find the intervals to which a gene belongs
  vecind <- (yind-1)*50 + xind # this is a vector with indices for dvec. 
  dens <- dvec[vecind]
  names(dens) <- names(x)
  return(dens)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.genetodensitycolor <- function(x,y) {
  require("MASS")
  require("colorRamps")
  ## NOTE would not work if two spots have exactly the same density or?
  ## important for a nice picture is that one removes the grid points with no gene (the ones with the lowest density)
  grids <- 100
  out <- list()
  est <- kde2d(x, y, n = grids) # needs MASS - calculates the two dimensional density
  if (sum(is.na(est$z)) > 0) {
    est <- kde2d(x, y, n = grids, h = c(bw.nrd0(x), bw.nrd0(y))) # needs MASS - calculates the two dimensional density
  }
  dvec <- as.vector(est$z) # note: as.vector takes columnwise, rows correspond to the value of x, columns to the value of y - the y index counts therefore grids times
  xind <- findInterval(x, est$x) # find the intervals to which a gene belongs
  yind <- findInterval(y, est$y) # find the intervals to which a gene belongs
  vecind <- (yind-1)*grids + xind # this is a vector with indices for dvec. 
  dens <- dvec[vecind]
  dens_with_genes <- unique(dens)
  names(dens) <- names(x)
  out$dens <- dens
  out$denschar <- as.vector(dens, mode = "character")
  colgrad <- blue2red(length(dens_with_genes))
  #colgrad <- rainbow(length(dens_with_genes), start = 0, end = 1, alpha = 0.8)
  #colgrad <- heat.colors(length(dens_with_genes), alpha = 0.8)
  #colgrad <- terrain.colors(length(dens_with_genes), alpha = 0.8)
  #colgrad <- topo.colors(length(dens_with_genes), alpha = 0.8)
  #colgrad <- cm.colors(length(dens_with_genes), alpha = 0.8)
  sdens_with_genes <- sort(dens_with_genes, decreasing = FALSE)
  names(colgrad) <- sdens_with_genes
  out$cols <- colgrad
  genecols <- colgrad[as.vector(dens, mode = "character")]
  names(genecols) <- names(x)
  out$genecols <- genecols
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.yellowblueblack <- function(x) {
  #rg <- approx(c(0, 0.5, 1), c(1, 1/3, 0), n = x)$y
  #b <- approx(c(0, 0.5, 1), c(2/3, 2/3, 0), n = x)$y
  rg <- approx(c(0, 0.5, 1), c(1, 0, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 1, 0), n = x)$y
  return(rgb(rg, rg, b))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.yellowblackblue <- function(x) {
  #rg <- approx(c(0, 1, 0.5), c(1, 1/3, 0), n = x)$y
  #b <- approx(c(0, 1, 0.5), c(2/3, 2/3, 0), n = x)$y
  rg <- approx(c(0, 0.5, 1), c(1, 0, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 0, 1), n = x)$y
  return(rgb(rg, rg, b))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.yellowredblue <- function(x) {
  r <- approx(c(0, 0.5, 1), c(1, 1, 0), n = x)$y
  g <- approx(c(0, 0.5, 1), c(1, 0, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 0, 1), n = x)$y
  return(rgb(r, g, b))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.yellowblack <- function(x) {
  rg <- approx(c(0, 1), c(1, 0), n = x)$y
  b <- approx(c(0, 1), c(0, 0), n = x)$y
  return(rgb(rg, rg, b))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.yellowredblack <- function(x) {
  r <- approx(c(0, 0.5, 1), c(1, 1, 0), n = x)$y
  g <- approx(c(0, 0.5, 1), c(1, 0, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 0, 0), n = x)$y
  return(rgb(r, g, b))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.redwhiteblack <- function(x) {
  r <- approx(c(0, 0.5, 1), c(1, 1, 0), n = x)$y
  g <- approx(c(0, 0.5, 1), c(0, 1, 0), n = x)$y
  b <- approx(c(0, 0.5, 1), c(0, 1, 0), n = x)$y
  return(rgb(r, g, b))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.blackblueyellow <- function(x) {
  return(rev(f.yellowblueblack(x)))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.blueblackyellow <- function(x) {
  return(rev(f.yellowblackblue(x)))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.blueredyellow <- function(x) {
  return(rev(f.yellowredblue(x)))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.blackyellow <- function(x) {
  return(rev(f.yellowblack(x)))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.blackredyellow <- function(x) {
  return(rev(f.yellowredblack(x)))
}

#'@title undocumented internal function
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.blackwhitered <- function(x) {
  return(rev(f.redwhiteblack(x)))
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Draw an image() with the values printed as text
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.image.with.text <- function(x, y, z, xLabel, yLabel, mainLabel, useLog = FALSE, ...) {
  xChars <- as.character(x)
  yChars <- as.character(y)
  xPos <- 1:length(xChars)
  yPos <- 1:length(yChars)
  if (useLog) {
    toPlot <- log2(z+1)
  } else {
    toPlot <- z
  }
  image(xPos, yPos, toPlot, xlab = xLabel, ylab = yLabel, main = mainLabel, yaxt = "n", xaxt = "n", ...)
  axis(1, at = xPos, labels = xChars, outer = FALSE)
  axis(2, at = yPos, labels = yChars, outer = FALSE, las = 1)
  text(rep(xPos, each = length(yPos)), yPos, t(z), cex = 1)
  return(NULL)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Draw a histogram (frequency or density)
#'@param x values
#'@param xName x-axis label
#'@param useFreq plot frequencies instead of densities
#'@param doNotShowSummary omit printing mean/median/sd
#'@param useLog log2-transform the data
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.histogram <- function(x, y = c(), xName = "", useFreq = FALSE, doNotShowSummary = FALSE, useLog = FALSE, ...) 
{
  if (length(x) == length(y)) {
    toPlot <- list(x = x, y = y)
    class(toPlot) <- "density"
  } else {
    toPlot <- density(x)
  }
  if (useFreq) { 
    toPlot$y <- toPlot$y * length(x)
  }
  if (useLog) {
    if (useFreq) { toPlot$y <- log2(toPlot$y+1) }
    else { cat("log only possible for useFreq = TRUE\n") }
  }
  par(mar = c(5,4,1,1))
  plot(
    toPlot,
    #main = "",
    bty = "n",
    xaxs = "r",
    yaxs = "r",
    xlab = xName,
    ylab = "",
    las = 1,
    cex = 0.4,
    tck = 0.01,
    ...
  )
  if (!doNotShowSummary) {
    text(
      par("usr")[1] + par("usr")[2]/15,
      par("usr")[4] - par("usr")[4]/15,
      adj = c(0,1),
      labels = paste(c("median:", "mean:", "sd:"),
                     c(round(median(x), digits = 3),
                       round(mean(x), digits = 3),
                       round(sd(x), digits = 3)),
                     collapse = "\n", sep=" ")
    )
  }
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Draw a scatterplot
#'@param x x-values
#'@param y y-values
#'@param xName x-axis label
#'@param yName y-axis label
#'@param ... other arguments passed on to plot()
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.scatter <- function(x, y, xName = "", yName = "", cexFactor = 1, ...)
{
  par(mar = c(5,5,1,1))
  plot(
    y ~ x,
    #main = "",
    bty = "n",
    xaxs = "r",
    yaxs = "r",
    xlab = xName,
    ylab = yName,
    las = 1,
    cex = 0.2*cexFactor,
    tck = 0.01,
    pch = 19,
    ...
  )
  text(
    par("usr")[1] + par("usr")[2]/15,
    par("usr")[4] - par("usr")[4]/15,
    adj = c(0,1),
    labels = paste(c("corP:", "corS:", "n:"),
                   c(round(cor(x,y,method="pearson"), digits = 3),
                     round(cor(x,y,method="spearman"), digits = 3),
                     length(x)),
                   collapse = "\n", sep=" ")
  )
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Replace Inf/-Inf with the next highes/lowest value.
#'@param x a numeric vector
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.replace.Inf.by.next <- function(x) {
  if (sum(abs(x) == Inf, na.rm = TRUE) == 0) {return(x)}
  cat("replacing Inf/-Inf with the next highest/lowest value\n")
  notNA <- !is.na(x)
  minInf <- x == -Inf
  pluInf <- x == Inf
  nx <- x[notNA&(!minInf)&(!pluInf)]
  x[notNA&pluInf] <- max(nx)
  x[notNA&minInf] <- min(nx)
  return(x)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Draw a scatterplot with heat colors indicating the point density
#'@param x x-values
#'@param y y-values
#'@param xName x-axis label
#'@param yName y-axis label
#'@param tryToScale EXPERIMENTAL - replace outliers (10 times bigger than the median) with 10*median
#'@param ... other arguments passed on to plot()
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.density.scatter <- function(x, y, xName = "", yName = "", main = "", tryToScale = FALSE, ...) {
  x <- f.replace.Inf.by.next(x)
  y <- f.replace.Inf.by.next(y)
  if (tryToScale) {
    if (max(x) > 10*median(x)) {
      isExtreme <- x > 10*median(x)
      cat("limiting data to 10*median(x) - i.e. replacing", sum(isExtreme), "in percent:", mean(isExtreme)*100," extreme values\n")
      x[isExtreme] <- 10*median(x)
    }
    if (max(y) > 10*median(y)) {
      isExtreme <- y > 10*median(y)
      cat("limiting data to 10*median(y) - i.e. replacing", sum(isExtreme), "in percent:", mean(isExtreme)*100," extreme values\n")
      y[isExtreme] <- 10*median(y)
    }
  }
  #   xLimits <- range(x)
  #   if (max(x) > 10*median(x)) {
  #     cat("limiting plot to a quantile(x, 0.99)\n")
  #     xLimits[2] <- quantile(x, 0.99)
  #   }
  #   yLimits <- range(y)
  #   if (max(y) > 10*median(y)) {
  #     cat("limiting plot to a quantile(y, 0.99)\n")
  #     yLimits[2] <- quantile(y, 0.99)
  #   }
  colorGrad <- f.genetodensitycolor(x, y)
  # f.scatter(x, y, xName, yName, cexFactor = 1, col = colorGrad$genecols, main = main, xlim = xLimits, ylim = yLimits)
  f.scatter(x, y, xName, yName, cexFactor = 1, col = colorGrad$genecols, main = main, ...)
  return(NULL)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Draw scatter plots, sorted scatter plots and histograms for all columns in a data frame
#'@param data matrix or dataframe, not too many columns (there will be ncol(data)*ncol(data) plots in total)
#'@param rDir folder for the figure
#'@param rt threshold for the data (it will test >= threshold, except if the threshold is zero, then > 0)
#'@param prefix prefix for the name of the plot
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.allinone <- function(data, rDir, rt = 0, prefix = "", scaleToSample = FALSE) 
{
  hits.union <- rownames(data)[which(apply(f.comp(data,rt),1,sum) > 0)]
  data <- data[hits.union,]
  # get sorted data and make sure everything is a data frame
  s.data <- apply(data,2,sort)
  data <- as.data.frame(data)
  s.data <- as.data.frame(s.data)
  # some variables
  ns = ncol(data)
  sn = colnames(data)
  # use CairoSVG if you are on windows
  afname <- paste(prefix,'_allinone_',rt,'_',paste(sn,collapse='_'),sep='')
  if (nchar(afname) > 50) { afname <- paste(prefix,'_allinone_',rt,'_many',sep='') }
  afpath <- file.path(rDir,paste(afname, '.tiff', sep = ''))
  tiff(file=afpath, width = 600*ns, height = 600*ns, units = "px", pointsize = 20, compression = "lzw", bg = "white", type = GLOBAL_VARIABLE_TIFF_LIB)
  # plot everything
  layout(mat=matrix(c(1:(ns*ns)), nrow=ns, byrow = TRUE))
  cat(paste("plotting", ns, "rows\n"))
  cat(paste(c("|", rep(' ', ns),"|\n|"),sep='',collapse=''))
  for (i in c(1:ns))
  {
    for (k in c(1:ns))
    {
      if (scaleToSample) {
        xMin <- floor(min(data[,k]))
        xMax <- ceiling(max(data[,k]))
        yMin <- floor(min(data[,i]))
        yMax <- ceiling(max(data[,i]))
      } else {
        xMin <- floor(min(data))
        xMax <- ceiling(max(data))
        yMin <- floor(min(data))
        yMax <- ceiling(max(data))
      }
      if (i == k) 
      {
        # histogram in the diagonal
        f.histogram(data[f.comp(data[,k],rt),k], xName = sn[k], xlim = c(xMin, xMax), main = "")
      }
      else if (i < k)
      {	
        # scatterplot with densities on the top right
        if ((length(unique(data[,k][!is.na(data[,k])])) > 2) & (length(unique(data[,i][!is.na(data[,i])])) > 2)) {
          colorGrad <- f.genetodensitycolor(data[,k], data[,i])
        } else {
          colorGrad <- list(genecols = rep("black", nrow(data)))
        }
        f.scatter(data[,k], data[,i], sn[k], sn[i], cexFactor = 1.5, col = colorGrad$genecols, xlim = c(xMin, xMax), ylim = c(yMin, yMax), main = "")
        abline(0,1,col="black")
      }
      else if (i > k)
      {
        # sorted scatterplot on the bottom left
        f.scatter(s.data[,k], s.data[,i], sn[k], sn[i], cexFactor = 1.5, xlim = c(xMin, xMax), ylim = c(yMin, yMax), main = "")
        abline(0,1,col="darkorchid3")
      }
    }
    cat("#")
  }
  cat('|\n')
  dev.off()
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Draw a correlation matrix for multiple samples (correlation between multiple samples)
#'@param data matrix or dataframe with one sample per column.
#'@param rDir a directory where the figure is stored.
#'@param outfile name of the figure and the table (without file extension, _cor_mat.txt/.svg will be added).
#'@param corMethod method used to calculate correlation (pearson, spearman, kendall).
#'@param summaryFunction function used to summarize the individual 4C correlations (default is median).
#'@param useOnlyHighVar use only the rows with high variability (specified by highVarPercentile).
#'@param tryAutoColor try to get a nice color-gradient.
#'@param highVarPercentile if useOnlyHighVar is given, only the rows with variation greater than the highVarPercentile are used
#'(0.9 means that only the 10\% with the highest variation are used).
#'@param corUse see "use" in \code{\link{cor}}
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}
#'@export
f.generic.correlation.matrix <- function(data, rDir, outfile, corMethod = "pearson", useOnlyHighVar = FALSE, tryAutoColor = TRUE, highVarPercentile = 0.9, corUse = "everything", ...) {
  require("gplots")
  if (useOnlyHighVar) {
    varVec <- apply(data, 1, var)
    lowerBound <- quantile(varVec, highVarPercentile)
    data <- data[varVec>=lowerBound,]
  }
  sampleCorMat <- cor(data, method = corMethod, use = corUse)
  # some plotting parameters
  minCor <- min(sampleCorMat, na.rm = TRUE)
  maxCor <- max(sampleCorMat[sampleCorMat != 1], na.rm = TRUE)
  if (sum(is.na(c(minCor, maxCor))) > 0) { tryAutoColor <- FALSE }
  if (tryAutoColor) {
    if (minCor > 0.8) { minCor <- 0.8 } else { if (minCor > 0.7) { minCor <- 0.7 } else { if (minCor > 0.5) { minCor <- 0.5 }}}
    if (maxCor > 0.8) { maxCor <- 1 }   else { if (maxCor < 0.8) { maxCor <- 0.8 } else { if (maxCor < 0.7) { maxCor <- 0.7 } else { if (maxCor < 0.5) { maxCor <- 0.5 } }}}
  }
  numCols <- ceiling((maxCor-minCor)*50)
  if (GLOBAL_VARIABLE_USE_SVG) {
    svg(file.path(rDir, paste(outfile, '_cor_mat.svg', sep = '')), onefile = TRUE, bg = FALSE, antialias = "default", pointsize = 10, width = 12, height = 12)
  } else {
    tiff(file.path(rDir, paste(outfile, "_cor_mat.tiff", sep = '')), width = 1600, height = 1600, compression = "lzw", type = GLOBAL_VARIABLE_TIFF_LIB)
  }
  labelCex <- 1 + 1/log10(nrow(sampleCorMat))
  noteCex <- 4
  if (nrow(sampleCorMat) > 6) { noteCex <- 3 }
  if (nrow(sampleCorMat) > 12) { noteCex <- 2 }
  if (nrow(sampleCorMat) > 18) { noteCex <- 1 }
  if (nrow(sampleCorMat) > 30) { #no notes
    heatmap.2(sampleCorMat, col = f.blueredyellow(numCols), trace="none", scale = "none", margins = c(15,15), breaks = seq(minCor, maxCor, length.out = numCols+1), ...)
  } else {
    heatmap.2(sampleCorMat, col = f.blueredyellow(numCols), trace="none", scale = "none", margins = c(15,15), cellnote = round(sampleCorMat,2), notecol = "black", notecex = noteCex, breaks = seq(minCor, maxCor, length.out = numCols+1), cexRow = labelCex, cexCol = labelCex, ...)
  }
  dev.off()
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title combine a list of sets with a series of operations.
#'@param data a list with sets (vectors), names accordingly.
#'@param combination a vector with set operations (union, intersect, setdiff) followed by the result and the arguments for them.
#'Operations follow the order in the vector. Previous results can therefore be recycled. Setdiff is order-sensitive (in first but not the others).
#'@return a list with all results.
#'@note if you have a nested list, use \code{unlist(x, recursive = FALSE)}. The top-level names are then separated by a dot from the low level nodes.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.combine.a.list.of.sets <- function(data, combination) {
  opeIn <- c(which(combination %in% c("union", "intersect", "setdiff")), length(combination)+1)
  allResults <- combination[opeIn[1:(length(opeIn)-1)]+1]
  if (sum(allResults %in% names(data)) > 0) {
    cat("some result names are the same as the set names - stopping\n")
    return(NULL)
  }
  for (i in 1:(length(opeIn)-1)) {
    os <- opeIn[i]
    oe <- opeIn[i+1]-1
    ops <- combination[os]
    res <- combination[os+1]
    args <- combination[(os+2):oe]
    if (ops == "union") {
      data[[res]] <- Reduce(union, data[args])
    } else if (ops == "intersect") {
      data[[res]] <- Reduce(intersect, data[args])
    } else if (ops == "setdiff") {
      data[[res]] <- Reduce(setdiff, data[args])
    } else {
      cat("one should never end up here\n")
    }
  }
  out <- data[allResults]
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Formulate a simple contrast using two sets of regular expressions.
#'@param plusTermsRE a vector of regular expressions or strings for the terms which will be summed up.
#'@param minusTermsRE a vector of regular expressions or strings for the terms which will be substracted.
#'@param design a design matrix obtained with \code{\link{model.matrix}}.
#'@param invertPlus see note below
#'@param invertMinus see note below
#'@param fixed set to TRUE to disable regular expression mode of \code{\link{grep}}.
#'@return a string containing the contrast definition (to be used in \code{\link{makeContrast}}).
#'@note The function will use \code{\link{grep}} to retrieve all \code{colnames(design)} 
#'which fit to either the plusTermsRE or the minusTermsRE. The contrast is then formulated as
#'plusTerms/length(plusTerms) - minusTerms/length(minusTerms). If invertPlus/invertMinus are set
#'to TRUE, the terms which do NOT match the search criteria are extracted (e.g. if one would like
#'to test a certain set A versus all non-A, one can do that with:
#'\code{f.formulate.simple.contrast(setA, setA, design, invertMinus = TRUE)})
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.formulate.simple.contrast <- function(plusTermsRE, minusTermsRE, design, invertPlus = FALSE, invertMinus = FALSE, fixed = FALSE) {
  plusTerms <- sapply(plusTermsRE, function(x) grep(x, colnames(design), fixed = fixed, value = TRUE, invert = invertPlus))
  minusTerms <- sapply(minusTermsRE, function(x) grep(x, colnames(design), fixed = fixed, value = TRUE, invert = invertMinus))
  if (length(plusTermsRE) > 1) {
    if (invertPlus) {
      plusTerms <- Reduce(intersect, plusTerms)
    } else { 
      plusTerms <- Reduce(union, plusTerms)
    }
  } else {
    plusTerms <- unique(unlist(plusTerms))
  }
  if (length(minusTermsRE) > 1) {
    if (invertMinus) {
      minusTerms <- Reduce(intersect, minusTerms)
    } else { 
      minusTerms <- Reduce(union, minusTerms)
    }
  } else { 
    minusTerms <- unique(unlist(minusTerms))
  }
  out <- paste0(
    "(", paste(plusTerms, collapse = "+"), ")/", length(plusTerms), "-",
    "(", paste(minusTerms, collapse = "+"), ")/", length(minusTerms)
  )
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Normalize count data for plotting with edgeR.
#'@param data a data.frame or matrix object with counts.
#'@param logCPM return normalized logCPMs instead.
#'@return normalized counts (logged, optionally logCPM).
#'@note requires edgeR to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.normalize.counts.edgeR <- function(data, logCPM = FALSE) {
  require("edgeR")
  dge <- DGEList(counts=round(data))
  dge <- calcNormFactors(dge)
  if (logCPM) {
    out <- cpm(dge, prior.count=1, log=TRUE)
  } else {
    out <- log2(equalizeLibSizes(dge)$pseudo.counts+1)
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Normalize count data for plotting with limma (requires edgeR eventually)
#'@param data a data.frame or matrix object with counts.
#'@param design a design matrix made with model.matrix(formula(formulaString), data = sampleTab, contrasts.arg = NULL). See \code{\link{f.multi.level.limma}}.
#'@param method c("none", "TMM", "quantile").
#'@return normalized values
#'@note requires limma and edgeR (TMM) to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.normalize.counts.limma <- function(data, design, method = "none") {
  require("limma")
  if (!(method %in% c("none", "TMM", "quantile"))) { cat("method", method, "not available\n"); stop(); }
  if (method == "TMM") {
    require("edgeR")
    data <- DGEList(round(data))
    data <- calcNormFactors(data)
    v <- voom(data, design)
  } else {
    v <- voom(data, design, normalize.method = method)
  }
  return(v$E)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Normalize count data for plotting with DESeq2
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param formulaString a string holding the formula of the model. See also \code{\link{f.multi.level.DESeq}}.
#'to get blind dispersion estimates and normalization, one can supply ~1 as formula (TODO check this).
#'@return normalized values
#'@note requires DESeq2 to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.normalize.counts.DESeq <- function(data, sampleTab, formulaString, tryRescue = FALSE) {
  require("DESeq2")
  if (length(grep("0\\+|0 \\+", formulaString)) == 1) {
    cat("removing the no intercept (0+) term for DESeq to enable LFC shrinkage\n")
    formulaString <- gsub("0\\+|0 \\+", "", formulaString)
  }
  dds <- DESeqDataSetFromMatrix(countData = round(data[,rownames(sampleTab)]), colData = sampleTab, design = formula(formulaString))
  if (tryRescue) {
    dds <- tryCatch(DESeq(dds),
                    error = function(e) { DESeq(dds, fitType = "mean") },
                    finally = cat("###### NOTE: ran with tryRescue == TRUE\n"))
  } else {
    dds <- DESeq(dds)
  }
  out <- log2(DESeq2::counts(dds, normalized = TRUE)+1)
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Get all possible normalization given a certain experimental layout
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param formulaString a string holding the formula of the model.
#'@param design a design matrix made with model.matrix(formula(formulaString), data = sampleTab, contrasts.arg = NULL).
#'@return a list with normalized counts. 
#'Names are edgeR, DESeq_default, limma_none, limma_TMM, and limma_quantile.
#'@note requires edgeR, DESeq2, and limma to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@seealso \code{\link{f.normalize.counts.edgeR}}, \code{\link{f.normalize.counts.DESeq}}, and \code{\link{f.normalize.counts.limma}}.
#'@export
f.all.normalizations <- function(data, sampleTab, formulaString,  design) {
  data <- data[,rownames(sampleTab)]
  out <- list(
    edgeR = f.normalize.counts.edgeR(data),
    DESeq_default = f.normalize.counts.DESeq(data, sampleTab, formulaString),
    limma_none = f.normalize.counts.limma(data, design, "none"),
    limma_TMM = f.normalize.counts.limma(data, design, "TMM"),
    limma_quantile = f.normalize.counts.limma(data, design, "quantile")
  )
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Make some overview plots of RNA-Seq data.
#'@param data a data.frame or matrix object with counts (if required, normalize and or log beforehand)
#'@param rDir directory where things are stored.
#'@param filePrefix some prefix for the files.
#'@param skipScatters set to true if you want to skip the scatterplots.
#'@return NULL
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@seealso \code{\link{f.equalize.counts.edgeR}}...
#'@export
f.do.some.overview <- function(data, rDir, filePrefix, skipScatters = FALSE) {
  if (!skipScatters) { f.allinone(data, rDir, prefix = filePrefix) }
  #f.generic.correlation.matrix(data, rDir, paste0(filePrefix, "_spearman"), "spearman", FALSE)
  f.generic.correlation.matrix(data, rDir, paste0(filePrefix, "_spearman_onlyHighVar"), "spearman", TRUE)
  return(NULL)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title remove entries of a table with low values and variance.
#'@param data a data.frame or matrix object with counts.
#'@param minVal minimum value which should be seen at least minTimes per row.
#'@param minTimes see above.
#'@param lowerVarQuantileToRemove remove entries with var within this quantile.
#'@return filtered data.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.strip.data <- function(data, minVal = 5, minTimes = 1, lowerVarQuantileToRemove = 0, colsToStrip = colnames(data)) {
  before <- nrow(data)
  data <- as.data.frame(data[apply(data[,colsToStrip] >= minVal, 1, sum) >= minTimes,])
  inBetween <- nrow(data)
  if (lowerVarQuantileToRemove > 0) {
    vars <- apply(data[,colsToStrip], 1, var)
    lowerBound <- quantile(vars, lowerVarQuantileToRemove)
    data <- data[vars > lowerBound, ]
  }
  after <- nrow(data)
  cat("####### f.strip.data #######\n")
  cat("###", before, "entries at beginning\n")
  cat("###", inBetween, "entries after removal of low values\n")
  cat("###", after, "entries after variance filter\n")
  cat("####### f.strip.data #######\n")
  return(data)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'A class for storing results from any kind of differential expression method.
#'@field tool name of the tool used ("edgeR", "DESeq2", "limma").
#'@field method name of the method used in the tool (currently only for print purpose).
#'@field table data.frame with four columns: generic, logFC, pVal, and adjP.
#'@field isPairwise TRUE for pairwise, FALSE for contrast.
#'@field pairOrCont either the pair which was been compared, or the name of the contrast.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export c.DEGtab
#'@exportClass c.DEGtab
c.DEGtab <- setRefClass(
  "c.DEGtab",
  fields = list(
    tool = "character", # the DE tool (will be used to get the meaning of one of the columns)
    method = "character", # method used for the tool
    table = "data.frame", # table with generic, logFC, pVal, and adjP 
    isPairwise = "logical", # TRUE for pairwise, FALSE for contrast
    pairOrCont = "character", # either the name of the contrast or the two groups (equals to the "pair" in the DE functions)
    condA = "character", # either pair[1] or "down"
    condB = "character", # either pair[2] or "up"
    corner = "character" # either up in or name of the contrast
  ),
  methods = list(
    initialize = function(...) {
      "initialize the class"
      .self$condA <- "unInit"
      .self$condB <- "unInit"
      .self$corner <- "unInit"
      callSuper(...)
      .self$fill_conditions()
    },
    fill_conditions = function() {
      "fill in the conditions"
      if (.self$isPairwise) {
        .self$condA <- .self$pairOrCont[1]
        .self$condB <- .self$pairOrCont[2]
        .self$corner <- paste0("up in ", .self$condA, "/", .self$condB)
      } else {
        .self$condA <- "down"
        .self$condB <- "up"
        .self$corner <- paste(pairOrCont, ": down/up")
      }
    },
    get_table = function() { 
      "returns the table with 'generic' replaced by the actual variable name."
      toolToName <- c(edgeR = "logCPM", DESeq = "logBaseMean", limma = "aveExpr")
      out <- .self$table
      if (tool == "DESeq") {out$generic <- log2(out$generic+1)}
      colnames(out)[grep("generic", colnames(out))] <- toolToName[tool]
      return(out)
    },
    get_print_table = function() { 
      "like get_table, but rounds, sorts and formats the table as well."
      out <- .self$get_table()
      out <- out[with(out, order(pVal, 1/abs(logFC))),]
      out[,c(1,2)] <- round(out[,c(1,2)], 2)
      out$pVal <- format(out$pVal, digits = 3, scientific = TRUE)
      out$adjP <- format(out$adjP, digits = 3, scientific = TRUE)
      return(out)
    },
    get_significant_entries = function(rawPcut = 0.05, adjPcut = 0.05, LFCcut = ifelse(.self$tool == "DESeq", 0, 1)) {
      "returns a list of three gene sets (conditionA, conditionB, any)"
      maskA <- .self$table$logFC < 0 # up in pair[1] - or down in contrast 
      maskB <- .self$table$logFC > 0 # up in pair[2] - or up in contrast
      maskLFC <- abs(.self$table$logFC) >= LFCcut
      maskPVAL <- .self$table$pVal <= rawPcut
      maskADJP <- .self$table$adjP <= adjPcut
      maskLFC[is.na(maskLFC)] <- FALSE
      maskPVAL[is.na(maskPVAL)] <- FALSE
      maskADJP[is.na(maskADJP)] <- FALSE
      if (adjPcut == 1) { cat("ignoring adjusted p-values\n"); maskADJP[] <- TRUE; }
      maskAny <- (maskLFC & maskPVAL & maskADJP)
      out <- list()
      out$any <- rownames(.self$table)[maskAny]
      out[[.self$condA]] <- rownames(.self$table)[maskAny & maskA]
      out[[.self$condB]] <- rownames(.self$table)[maskAny & maskB]
      return(out)
    },
    get_stats = function(pValsToCheck = c(0.05, 0.01, 0.001), adjPsToCheck = c(0.1, 0.05, 0.01), logFCToCheck = c(0, 0.5, 1, 1.5, 2, 2.5, 3)) {
      "returns a 3D array with stats. The third dimension is c(conditionA, conditionB, 'simple')"
      tabRows <- c(paste("rawP_max_", pValsToCheck, sep = ''), paste("adjP_max_", adjPsToCheck, sep = ''))
      tabCols <- paste("logFC_min_", logFCToCheck, sep = '')
      tabLayer <- c(.self$condA, .self$condB, "simple")
      out <- array(0, c(length(tabRows), length(tabCols), length(tabLayer)), list(tabRows, tabCols, tabLayer))
      maskA <- table$logFC < 0 # up in pair[1] - or down in contrast 
      maskB <- table$logFC > 0 # up in pair[2] - or up in contras
      for (i in 1:length(logFCToCheck)) {
        maskLFC <- abs(.self$table$logFC) >= logFCToCheck[i]
        for (j in 1:length(pValsToCheck)) {
          k <- j + length(pValsToCheck)
          maskPVAL <- .self$table$pVal <= pValsToCheck[j]
          maskADJP <- .self$table$adjP <= adjPsToCheck[j]
          out[j, i, .self$condA] <- sum(maskLFC & maskPVAL & maskA, na.rm = TRUE)
          out[k, i, .self$condA] <- sum(maskLFC & maskADJP & maskA, na.rm = TRUE)
          out[j, i, .self$condB] <- sum(maskLFC & maskPVAL & maskB, na.rm = TRUE)
          out[k, i, .self$condB] <- sum(maskLFC & maskADJP & maskB, na.rm = TRUE)
        }
      }
      for (rn in tabRows) {
        out[rn,,"simple"] <- paste(out[rn,,.self$condA], out[rn,,.self$condB], sep = " / ")
      }
      return(out)
    },
    print_stats = function(pValsToCheck = c(0.05, 0.01, 0.001), adjPsToCheck = c(0.1, 0.05, 0.01), logFCToCheck = c(0, 0.5, 1, 1.5, 2, 2.5, 3)) {
      "prints the simple stats to the terminal"
      stats <- .self$get_stats(pValsToCheck, adjPsToCheck, logFCToCheck)
      cat(paste(c(.self$corner, colnames(stats)), collapse = "\t"), "\n")
      for (rn in rownames(stats)) {
        cat(paste(c(rn, stats[rn,,"simple"]), collapse = "\t"), "\n")
      }
    },
    print_stats_markdown = function(pValsToCheck = c(0.05, 0.01, 0.001), adjPsToCheck = c(0.1, 0.05, 0.01), logFCToCheck = c(0, 0.5, 1, 1.5, 2, 2.5, 3)) {
      "prints the simple stats to the terminal (in markdown format)"
      stats <- .self$get_stats(pValsToCheck, adjPsToCheck, logFCToCheck)
      cat("|", paste(c(.self$corner, colnames(stats)), collapse = " | "), "|\n")
      for (rn in rownames(stats)) {
        cat("|", paste(c(rn, stats[rn,,"simple"]), collapse = " | "), "|\n")
      }
    },
    add_data = function(x) {
      "adds new columns to the table"
      if (!is.data.frame(x)) {stop("expecting a data.frame object")}
      # check colnames
      colnames(x)[colnames(x) %in% colnames(.self$table)] <- paste0("added_", colnames(x)[colnames(x) %in% colnames(.self$table)])
      x <- x[intersect(rownames(x), rownames(.self$table)),]
      for (cn in colnames(x)) {
        .self$table[[cn]] <- NA
        .self$table[rownames(x), cn] <- x[,cn]
      }
      return(.self)
    }
  )
)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform pairwise comparison with edgeR.
#'@param data a data.frame or matrix object with counts.
#'@param groups the group labels for the samples (sorted accoring to the column names).
#'@param pair a vector with the two groups to be compared.
#'@param method the method used for dispersion estimates.
#'@param batch optionally a vector of same length as groups to specify a batch factor.
#'@return a \code{\link{c.DEGtab}} object. 
#'@note requires edgeR to be installed.
#'@references 
#'Robinson, M.D. and Oshlack, A. (2010)
#'A scaling normalization method for differential expression analysis of RNA-Seq data. \emph{Genome Biology} \bold{11}, R25.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.two.groups.edgeR <- function(data, groups, pair, method = "trended", batch = c()) {
  require("edgeR")
  # check args
  groups <- as.character(groups)
  if (!(method %in% c("common", "trended", "tagwise"))) { cat("undefined method:", method, "\n"); stop;}
  doBatch <- (length(batch) == length(groups))
  cat("running f.two.groups.edgeR with method =", method, "( doBatch is", doBatch, ")\n")
  # create DGEList - use only the samples that are relevant for the comparison
  mask <- groups %in% pair
  data <- data[,mask]
  groups <- groups[mask]
  if (doBatch) { batch <- batch[mask] }
  else { batch <- rep("B0", sum(mask)) }
  if (is.numeric(batch)) { batch <- paste0("B", batch) }
  groups <- factor(groups, levels = pair) # the first is the reference - so it is second - first
  batch <- factor(batch, levels = unique(batch)) # should be irrelevant though - except maybe if one would supply numbers
  # construct the DEGlist - add a batch factor if existing
  data <- new("DGEList", list(counts = round(data), samples = data.frame(GROUP = groups, BATCH = batch, lib.sizes = colSums(data), row.names = colnames(data))))
  formulaString <- ifelse(doBatch, "~BATCH+GROUP", "~GROUP")
  design <- model.matrix(formula(formulaString), data = data$samples, contrasts.arg = NULL)
  cat("testing", colnames(design)[ncol(design)], "vs its reference.\n")
  # calculate dispersion
  if (method == "common") { data <- estimateGLMCommonDisp(data, design) }
  if (method == "trended") { data <- estimateGLMTrendedDisp(data, design) }
  if (method == "tagwise") { 
    data <- estimateGLMCommonDisp(data, design)
    data <- estimateGLMTagwiseDisp(data, design)
  }
  # model fitting and adjusting p values
  fit <- glmFit(data, design)
  res <- glmLRT(fit, coef = ncol(design))$table # note - the default is the last column anyway
  res$adjustedPValue <- p.adjust(res$PValue, method = "BH")
  # modify/simplify the results to get a table valid for all testing functions
  deTab <- data.frame(
    generic = res$logCPM, 
    logFC = res$logFC,
    pVal = res$PValue,
    adjP = res$adjustedPValue,
    row.names = rownames(res)
  )
  out <- c.DEGtab(tool = "edgeR", method = method, table = deTab, isPairwise = TRUE, pairOrCont = as.character(pair))
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform pairwise comparison with DESeq2.
#'@param data a data.frame or matrix object with counts.
#'@param groups the group labels for the samples (sorted accoring to the column names).
#'@param pair a vector with the two groups to be compared.
#'@param method a dummy argument - not used.
#'@param batch optionally a vector of same length as groups to specify a batch factor.
#'@return a \code{\link{c.DEGtab}} object. 
#'@note requires DESeq2 to be installed.
#'@references 
#'Love, M.I. and Huber, W. and Anders, S. (2014)
#'Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. \emph{Genome Biology} \bold{15}, 550.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.two.groups.DESeq <- function(data, groups, pair, method = "default", batch = c()) {
  require("DESeq2")
  groups <- as.character(groups)
  doBatch <- (length(batch) == length(groups))
  cat("running f.two.groups.DESeq ( doBatch is", doBatch, ")\n")
  # create DESeq object - use only the samples that are relevant for the comparison
  mask <- groups %in% pair
  data <- data[,mask]
  groups <- groups[mask]
  if (doBatch) { batch <- batch[mask] }
  else { batch <- rep("B0", sum(mask)) }
  if (is.numeric(batch)) { batch <- paste0("B", batch) }
  groups <- factor(groups, levels = pair) # the first is the reference - so it is second - first
  batch <- factor(batch, levels = unique(batch)) # should be irrelevant though - except maybe if one would supply numbers
  sampleTab <- data.frame(GROUP = groups, BATCH = batch, row.names = colnames(data))
  # construct the DEseqData - add a batch factor if existing
  formulaString <- ifelse(doBatch, "~BATCH+GROUP", "~GROUP")
  dds <- DESeqDataSetFromMatrix(countData = round(data), colData = sampleTab, design = formula(formulaString))
  # do the fitting and testing
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds, contrast = c("GROUP", rev(pair)))) # order in pair needs to be switched
  # modify/simplify the results to get a table valid for all testing functions
  deTab <- data.frame(
    generic = res$baseMean,
    logFC = res$log2FoldChange,
    pVal = res$pvalue,
    adjP = res$padj,
    row.names = rownames(res)
  )
  out <- c.DEGtab(tool = "DESeq", method = method, table = deTab, isPairwise = TRUE, pairOrCont = as.character(pair))
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform pairwise comparison with limma.
#'@param data a data.frame or matrix object with counts.
#'@param groups the group labels for the samples (sorted accoring to the column names).
#'@param pair a vector with the two groups to be compared
#'@param method method for normalization (none, TMM, or quantile)
#'@param batch optionally a vector of same length as groups to specify a batch factor.
#'@return a \code{\link{c.DEGtab}} object. 
#'@note requires limma to be installed.
#'@references 
#'Ritchie, M.E. and Phipson, B. and Wu, D. and Hu, Y. and Law, C.W. and Shi, W. and Smyth, G. (2015)
#'limma powers differential expression analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research} \bold{43}, e47.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.two.groups.limma <- function(data, groups, pair, method = "TMM", batch = c()) {
  require("limma")
  # check args
  groups <- as.character(groups)
  if (!(method %in% c("TMM", "quantile", "none"))) { cat("undefined method:", method, "\n"); stop;}
  doBatch <- (length(batch) == length(groups))
  cat("running f.two.groups.limma with method =", method, "( doBatch is", doBatch, ")\n")
  # create DGEList - use only the samples that are relevant for the comparison
  mask <- groups %in% pair
  data <- data[,mask]
  groups <- groups[mask]
  if (doBatch) { batch <- batch[mask] }
  else { batch <- rep("B0", sum(mask)) }
  if (is.numeric(batch)) { batch <- paste0("B", batch) }
  groups <- factor(groups, levels = pair) # the first is the reference - so it is second - first
  batch <- factor(batch, levels = unique(batch)) # should be irrelevant though - except maybe if one would supply numbers
  sampleTab <- data.frame(GROUP = groups, BATCH = batch, row.names = colnames(data))
  # round and get the design
  data <- round(data)
  formulaString <- ifelse(doBatch, "~BATCH+GROUP", "~GROUP")
  design <- model.matrix(formula(formulaString), sampleTab, contrasts.arg = NULL)
  cat("testing", colnames(design)[ncol(design)], "vs its reference.\n")
  # apply normalizations and voom transformation
  if (method == "TMM") {
    d <- DGEList(counts = data)
    d <- calcNormFactors(d)
    v <- voom(d, design)
  } else {
    v <- voom(data, design, normalize.method = method)
  }
  # fit and test
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = ncol(design), number = Inf, adjust="BH")
  # modify/simplify the results to get a table valid for all testing functions
  deTab <- data.frame(
    generic = res$AveExpr,
    logFC = res$logFC,
    pVal = res$P.Value,
    adjP = res$adj.P.Val,
    row.names = rownames(res)
  )
  out <- c.DEGtab(tool = "limma", method = method, table = deTab, isPairwise = TRUE, pairOrCont = as.character(pair))
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform all pairwise comparisons in an experiment with several groups.
#'@param data a data.frame or matrix object with counts.
#'@param samples vector with the sample names - same length and sorting as groups.
#'@param groups vector with the groups the samples belong to.
#'@param DEFUN function used for the pairwise comparison (\code{\link{f.two.groups.edgeR}}, \code{\link{f.two.groups.DESeq}} or \code{\link{f.two.groups.limma}}).
#'@param method is passed on to the method argument of the DEFUN.
#'@param batch optionally a vector of same length as groups to specify a batch factor.
#'@return a list with all the result from DEFUN. The names equal to GroupA_vs_GroupB (for the formulas, this means groupB-groupA). 
#'@note requires edgeR, DESeq2, and limma to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@seealso \code{\link{f.two.groups.edgeR}}, \code{\link{f.two.groups.DESeq}}, and \code{\link{f.two.groups.limma}}.
#'@export
f.all.pairwise.comparisons <- function(data, samples, groups, DEFUN, method = "none", batch = c()) {
  groups <- as.character(groups)
  doBatch <- (length(batch) == length(groups))
  temp <- table(groups)
  okGroups <- names(temp)[temp > 1]
  samples <- samples[groups %in% okGroups]
  groups <- groups[groups %in% okGroups]
  if (doBatch) { batch <- batch[groups %in% okGroups] }
  uniqueGroups <- unique(groups)
  out <- list()
  for (i in 1:(length(uniqueGroups)-1)) {
    groupA <- uniqueGroups[i]
    for (j in (i+1):length(uniqueGroups)){
      groupB <- uniqueGroups[j]
      samplesToUse <- c(samples[groups %in% groupA], samples[groups %in% groupB])
      groupsToUse <- c(groups[groups %in% groupA], groups[groups %in% groupB])
      if (doBatch) { batchToUse <- c(batch[groups %in% groupA], batch[groups %in% groupB]) }
      else { batchToUse <- c() }
      out[[paste(groupA, "_vs_", groupB, sep = '')]] <- DEFUN(data[,samplesToUse], groupsToUse, c(groupA, groupB), method, batchToUse)
    }
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform multiple series of pairwise comparisons.
#'@param data a data.frame or matrix object with counts.
#'@param samples vector with the sample names - same length and sorting as groups.
#'@param groups vector with the groups the samples belong to.
#'@param batch optionally a vector of same length as groups to specify a batch factor.
#'@return a list with the results of \code{\link{f.all.pairwise.comparisons}}.
#'Names are edgeR_common, edgeR_tagwise, edgeR_trended, DESeq_default, limma_none, limma_TMM, and limma_quantile.
#'@note requires edgeR, DESeq2, and limma to be installed.
#'This function will give a list with the method first and then all the pairs. For the other way around, check \code{\link{f.multiple.two.group.comparisons.pairwise}}
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@seealso \code{\link{f.two.groups.edgeR}}, \code{\link{f.two.groups.DESeq}}, and \code{\link{f.two.groups.limma}}.
#'@export
f.multiple.all.pairwise.comparisons <- function(data, samples, groups, batch = c()) {
  groups <- as.character(groups)
  if (length(unique(samples)) != length(unique(groups))) {
    out <- list(
      edgeR_common = f.all.pairwise.comparisons(data, samples, groups, f.two.groups.edgeR, "common", batch),
      edgeR_tagwise = f.all.pairwise.comparisons(data, samples, groups, f.two.groups.edgeR, "tagwise", batch),
      edgeR_trended = f.all.pairwise.comparisons(data, samples, groups, f.two.groups.edgeR, "trended", batch),
      DESeq_default = f.all.pairwise.comparisons(data, samples, groups, f.two.groups.DESeq, "default", batch),
      limma_none = f.all.pairwise.comparisons(data, samples, groups, f.two.groups.limma, "none", batch),
      limma_TMM = f.all.pairwise.comparisons(data, samples, groups, f.two.groups.limma, "TMM", batch),
      limma_quantile = f.all.pairwise.comparisons(data, samples, groups, f.two.groups.limma, "quantile", batch)
    )
  } else {
    if (length(batch) == length(groups)) { cat("unicates are already poor enough - if there's a batch effect on top of it, there is not solution\n"); }
    out <- list(DESeq_default = f.all.pairwise.comparisons(data, samples, groups, f.two.groups.DESeq))
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform a series of two group comparisons.
#'@param data a data.frame or matrix object with counts.
#'@param groups the group labels for the samples (sorted accoring to the column names).
#'@param pair a vector with the two groups to be compared.
#'@param batch optionally a vector of same length as groups to specify a batch factor.
#'@return a list with results from the individual test methods. 
#'Names are edgeR_common, edgeR_tagwise, edgeR_trended, DESeq_default, limma_none, limma_TMM, and limma_quantile.
#'@note requires edgeR, DESeq2, and limma to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@seealso \code{\link{f.two.groups.edgeR}}, \code{\link{f.two.groups.DESeq}}, and \code{\link{f.two.groups.limma}}.
#'@export
f.multiple.two.group.comparisons <- function(data, groups, pair, batch = c()) {
  groups <- as.character(groups)
  if (sum(groups %in% pair) > 2) {
    out <- list(
      edgeR_common = f.two.groups.edgeR(data, groups, pair, "common", batch),
      edgeR_tagwise = f.two.groups.edgeR(data, groups, pair, "tagwise", batch),
      edgeR_trended = f.two.groups.edgeR(data, groups, pair, "trended", batch),
      DESeq_default = f.two.groups.DESeq(data, groups, pair, "default", batch),
      limma_none = f.two.groups.limma(data, groups, pair, "none", batch),
      limma_TMM = f.two.groups.limma(data, groups, pair, "TMM", batch),
      limma_quantile = f.two.groups.limma(data, groups, pair, "quantile", batch)
    )
  } else {
    if (length(unique(batch)) == length(unique(groups))) { cat("unicates are already poor enough - if there's a batch effect on top of it, there is not solution\n"); }
    out <- list(DESeq_default = f.two.groups.DESeq(data, groups, pair))
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform a pairwise series of two group comparisons.
#'@param data a data.frame or matrix object with counts.
#'@param samples vector with the sample names - same length and sorting as groups.
#'@param groups vector with the groups the samples belong to.
#'@param batch optionally a vector of same length as groups to specify a batch factor.
#'@return a list with the results of \code{\link{f.multiple.two.group.comparisons}}. The names equal to GroupA_vs_GroupB (for the formulas, this means groupB-groupA). 
#'@note requires edgeR, DESeq2, and limma to be installed. 
#'This function will give a list with the pair first and then all the methods. For the other way around, check \code{\link{f.multiple.all.pairwise.comparisons}}
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@seealso \code{\link{f.two.groups.edgeR}}, \code{\link{f.two.groups.DESeq}}, and \code{\link{f.two.groups.limma}}.
#'@export
f.multiple.two.group.comparisons.pairwise <- function(data, samples, groups, batch = c()) {
  groups <- as.character(groups)
  doBatch <- (length(batch) == length(groups))
  temp <- table(groups)
  okGroups <- names(temp)[temp > 1]
  samples <- samples[groups %in% okGroups]
  groups <- groups[groups %in% okGroups]
  if (doBatch) { batch <- batch[groups %in% okGroups] }
  uniqueGroups <- unique(groups)
  out <- list()
  for (i in 1:(length(uniqueGroups)-1)) {
    groupA <- uniqueGroups[i]
    for (j in (i+1):length(uniqueGroups)){
      groupB <- uniqueGroups[j]
      samplesToUse <- c(samples[groups %in% groupA], samples[groups %in% groupB])
      groupsToUse <- c(groups[groups %in% groupA], groups[groups %in% groupB])
      if (doBatch) { batchToUse <- c(batch[groups %in% groupA], batch[groups %in% groupB]) }
      else { batchToUse <- c() }
      out[[paste(groupA, "_vs_", groupB, sep = '')]] <- f.multiple.two.group.comparisons(data[,samplesToUse], groupsToUse, c(groupA, groupB), batchToUse)
    }
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Combine two or more columns of a sample table for DE testing with multiple groups.
#'@param sampleTab a table with samples as rownames, variables as colnames, and their levels as entires.
#'@param colsToUse the columns which should be combined.
#'@return the sample table with an additional column
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.combine.factors.of.sample.table <- function(sampleTab, colsToUse = colnames(sampleTab)) {
  newFac <- paste(paste(colsToUse, collapse = 'x'), '__', sep = '')
  cat("new factor name:", newFac, "\n")
  sampleTab[[newFac]] <- apply(sampleTab[,colsToUse], 1, function(x) paste(x, collapse = '_'))
  return(sampleTab)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title get the names resulting from a DESeq2 model - just a helper
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'@param formulaString a string holding the formula of the model - results .
#'@return a vector with the resultnames given the sampleTable and the formula.
#'@note all variables will be converted to factors!
#'useful for \code{\link{f.multiple.groups.DESeq}}, \code{\link{f.multiple.groups.edgeR}}, and \code{\link{f.multiple.groups.limma}}
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.result.names.and.model.matrix <- function(sampleTab, formulaString) {
  require("DESeq2")
  dds <- makeExampleDESeqDataSet(n = 100, m = nrow(sampleTab))
  for (fn in colnames(sampleTab)) {
    dds[[fn]] <- as.factor(sampleTab[[fn]])
  }
  if (length(grep("0\\+|0 \\+", formulaString)) == 1) {
    cat("removing the no intercept (0+) term for DESeq as it seems to be unnecessary\n")
    formulaString <- gsub("0\\+|0 \\+", "", formulaString)
  }
  design(dds) <- formula(formulaString)
  dds <- DESeq(dds)
  out <- list(
    resNames = resultsNames(dds),
    modelMatrix = attr(dds, "modelMatrix"),
    modelMatrixType = attr(dds, "modelMatrixType")
  )
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform multi-level comparisons with DESeq2. This function is NOT FOR MULTIPLE FACTORS.
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param formulaString a string holding the formula of the model.
#'@param contrastsToTest a contrast matrix (see examples).
#'@param design normally a design matrix - only here for compatibility.
#'@param method also for compatibility.
#'@return list with \code{\link{c.DEGtab}} objects, named after the name of the contrast.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.multi.level.DESeq <- function(data, sampleTab, formulaString, contrastsToTest, design = "none", method = "default", tryRescue = FALSE) {
  require("DESeq2")
  cat("running multi-level comparison with DESeq2\n")
  if (length(grep("0\\+|0 \\+", formulaString)) == 1) {
    cat("removing the no intercept (0+) term for DESeq to enable LFC shrinkage\n")
    formulaString <- gsub("0\\+|0 \\+", "", formulaString)
  }
  contrastsToTest <- as.data.frame(contrastsToTest)
  dds <- DESeqDataSetFromMatrix(countData = round(data[,rownames(sampleTab)]), colData = sampleTab, design = formula(formulaString))
  if (tryRescue) {
    dds <- tryCatch(DESeq(dds),
                    error = function(e) { DESeq(dds, fitType = "mean") },
                    finally = cat("###### NOTE: ran with tryRescue == TRUE\n"))
  } else {
    dds <- DESeq(dds)
  }
  conRes <- lapply(contrastsToTest, function(x) as.data.frame(results(dds, contrast = c(0,x)))) # the zero is for the intercept added in DESeq2
  conRes <- lapply(conRes, function(x) data.frame(generic = x$baseMean, logFC = x$log2FoldChange, pVal = x$pvalue, adjP = x$padj, row.names = rownames(x)))
  out <- list()
  for (cont in names(conRes)) {
    out[[cont]] <- c.DEGtab(tool = "DESeq", method = method, table = conRes[[cont]], isPairwise = FALSE, pairOrCont = as.character(cont))
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform multi-level comparisons with edgeR. This function is NOT FOR MULTIPLE FACTORS.
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param formulaString a string holding the formula of the model.
#'@param contrastsToTest a contrast matrix (see examples).
#'@param design a design matrix made with model.matrix(formula(formulaString), data = sampleTab, contrasts.arg = NULL).
#'@param method the method used for dispersion estimates (common, trended, tagwise).
#'@return list with \code{\link{c.DEGtab}} objects, named after the name of the contrast.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.multi.level.edgeR <- function(data, sampleTab, formulaString, contrastsToTest, design, method = "trended") {
  require("edgeR")
  cat("running multi-level comparison with edgeR\n")
  # check method argument
  if (!(method %in% c("common", "trended", "tagwise"))) { cat("undefined method:", method, "\n"); stop;}
  # prepare data - one needs the column with the factor
  contrastsToTest <- as.data.frame(contrastsToTest)
  multiLevelFactor <- paste0(unlist(strsplit(colnames(design)[ncol(design)], "__", fixed = TRUE))[1], "__")
  cat("identified", multiLevelFactor, "as the column of interest\n")
  data <- data[,rownames(sampleTab)]
  temp <- list(
    counts = round(data),
    samples = data.frame(sampleTab[[multiLevelFactor]], lib.sizes = colSums(data), row.names = rownames(sampleTab))
  )
  colnames(temp$samples) <- c(multiLevelFactor, "lib.sizes")
  data <- new("DGEList", temp)
  # calculate dispersion
  if (method == "common") { data <- estimateGLMCommonDisp(data, design) }
  if (method == "trended") { data <- estimateGLMTrendedDisp(data, design) }
  if (method == "tagwise") { 
    data <- estimateGLMCommonDisp(data, design)
    data <- estimateGLMTagwiseDisp(data, design)
  }
  # model fitting
  fit <- glmFit(data, design)
  conRes <- lapply(contrastsToTest, function(x) as.data.frame(glmLRT(fit, contrast = x)$table))
  conRes <- lapply(conRes, function(x) data.frame(generic = x$logCPM, logFC = x$logFC, pVal = x$PValue, adjP = p.adjust(x$PValue, "BH"), row.names = rownames(x)))
  out <- list()
  for (cont in names(conRes)) {
    out[[cont]] <- c.DEGtab(tool = "edgeR", method = method, table = conRes[[cont]], isPairwise = FALSE, pairOrCont = as.character(cont))
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform multi-level comparisons with limma. This function is NOT FOR MULTIPLE FACTORS.
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param formulaString a string holding the formula of the model.
#'@param contrastsToTest a contrast matrix (see examples).
#'@param design a design matrix made with model.matrix(formula(formulaString), data = sampleTab, contrasts.arg = NULL).
#'@param method the method used for normalization (none, TMM, quantile).
#'@return list with \code{\link{c.DEGtab}} objects, named after the name of the contrast.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.multi.level.limma <- function(data, sampleTab, formulaString, contrastsToTest, design, method = "TMM") {
  require("limma")
  cat("running multi-level comparison with limma\n")
  # check method argument
  if (!(method %in% c("none", "TMM", "quantile"))) { cat("undefined method:", method, "\n"); stop;}
  # make sure data is sorted correctly
  data <- data[,rownames(sampleTab)]
  # apply normalizations and voom transformation
  if (method == "TMM") {
    data <- DGEList(round(data))
    data <- calcNormFactors(data)
    v <- voom(data, design)
  } else {
    v <- voom(data, design, normalize.method = method)
  }
  # fit and test
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrastsToTest)
  fit <- eBayes(fit)
  # extract results
  conRes <- list()
  for (cont in colnames(contrastsToTest)) {
    conRes[[cont]] <- topTable(fit, coef = cont, number = Inf, adjust = "BH")
  }
  # and simplify
  conRes <- lapply(conRes, function(x) data.frame(generic = x$AveExpr, logFC = x$logFC, pVal = x$P.Value, adjP = x$adj.P.Val, row.names = rownames(x)))
  out <- list()
  for (cont in names(conRes)) {
    out[[cont]] <- c.DEGtab(tool = "limma", method = method, table = conRes[[cont]], isPairwise = FALSE, pairOrCont = as.character(cont))
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform a series of multi level tests.
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param formulaString a string holding the formula of the model.
#'@param contrastsToTest a contrast matrix (see examples).
#'@param design a design matrix made with model.matrix(formula(formulaString), data = sampleTab, contrasts.arg = NULL).
#'@return a list with results from the individual test methods. 
#'Names are edgeR_common, edgeR_tagwise, edgeR_trended, DESeq_default, limma_none, limma_TMM, and limma_quantile.
#'@note requires edgeR, DESeq2, and limma to be installed.
#'@references 
#'Robinson, M.D. and Oshlack, A. (2010)
#'A scaling normalization method for differential expression analysis of RNA-Seq data. \emph{Genome Biology} \bold{11}, R25.
#'Love, M.I. and Huber, W. and Anders, S. (2014)
#'Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. \emph{Genome Biology} \bold{15}, 550.
#'Ritchie, M.E. and Phipson, B. and Wu, D. and Hu, Y. and Law, C.W. and Shi, W. and Smyth, G. (2015)
#'limma powers differential expression analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research} \bold{43}, e47.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@seealso \code{\link{f.multi.level.edgeR}}, \code{\link{f.multi.level.DESeq}}, and \code{\link{f.multi.level.limma}}.
#'@export
f.multiple.multi.level.comparisons <- function(data, sampleTab, formulaString, contrastsToTest, design) {
  out <- list(
    edgeR_common = f.multi.level.edgeR(data, sampleTab, formulaString, contrastsToTest, design, "common"),
    edgeR_tagwise = f.multi.level.edgeR(data, sampleTab, formulaString, contrastsToTest, design, "tagwise"),
    edgeR_trended = f.multi.level.edgeR(data, sampleTab, formulaString, contrastsToTest, design, "trended"),
    DESeq_default = f.multi.level.DESeq(data, sampleTab, formulaString, contrastsToTest),
    limma_none = f.multi.level.limma(data, sampleTab, formulaString, contrastsToTest, design, "none"),
    limma_TMM = f.multi.level.limma(data, sampleTab, formulaString, contrastsToTest, design, "TMM"),
    limma_quantile = f.multi.level.limma(data, sampleTab, formulaString, contrastsToTest, design, "quantile")
  )
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title A helper function for the two-by-two factor/level analysis functions
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param factorA the name of the first factor in the model (a colname in sampleTab)
#'@param factorB the name of the second factor in the model (a colname in sampleTab)
#'@param batchCol optionally the name of the column holding the info on the batch.
#'@return a list with the formulaString, the design, and the defaultContrasts.
#'@seealso \code{\link{f.two.by.two.edgeR}}, \code{\link{f.two.by.two.DESeq}}, and \code{\link{f.two.by.two.limma}}.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.prepare.two.by.two.analysis <- function(sampleTab, factorA, factorB, batchCol = "") {
  # check args
  if (!(is.factor(sampleTab[[factorA]]) & is.factor(sampleTab[[factorB]]))) { cat("convert the variables in the sampleTab into factors first.\n"); stop(); }
  if (!((length(levels(sampleTab[[factorA]])) == 2) & (length(levels(sampleTab[[factorB]])) == 2))) { cat("this is only available for 2x2 factor designs.\n"); stop(); }
  doBatch <- nchar(batchCol) > 0
  if (doBatch & !(is.factor(sampleTab[[batchCol]]))) { cat("convert the batch variable into a factor first.\n"); stop(); }
  A_first <- levels(sampleTab[[factorA]])[1]
  A_second <- levels(sampleTab[[factorA]])[2]
  B_first <- levels(sampleTab[[factorB]])[1]
  B_second <- levels(sampleTab[[factorB]])[2]
  formulaString <- ifelse(doBatch, paste0("~", batchCol, "+", factorA, "*", factorB), paste0("~", factorA, "*", factorB))
  design <- model.matrix(formula(formulaString), data = sampleTab)
  numNonsense <- ifelse(doBatch, length(levels(sampleTab[[batchCol]])), 1)
  toAdd <- rep(0, numNonsense)
  defaultContrasts <- data.frame(
    A_within_B_first = c(toAdd,1,0,0),
    A_within_B_second = c(toAdd,1,0,1),
    B_within_A_first = c(toAdd,0,1,0),
    B_within_A_second = c(toAdd,0,1,1),
    A_cross_B = c(toAdd,0,0,1),
    row.names = colnames(design)
  )
  colnames(defaultContrasts) <- c(
    paste0(A_second, "_vs_", A_first, "_in_", B_first),
    paste0(A_second, "_vs_", A_first, "_in_", B_second),
    paste0(B_second, "_vs_", B_first, "_in_", A_first),
    paste0(B_second, "_vs_", B_first, "_in_", A_second),
    paste0(factorA, "_x_", factorB)
  )
  defaultContrasts
  out <- list(formulaString = formulaString, design = design, defaultContrasts = defaultContrasts)
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform a two-by-two factor/level analysis with DESeq
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param factorA the name of the first factor in the model (a colname in sampleTab)
#'@param factorB the name of the second factor in the model (a colname in sampleTab)
#'@param method for compatibility only.
#'@param batchCol optionally the name of the column holding the info on the batch.
#'@return list with \code{\link{c.DEGtab}} objects, named after the name of the contrast.
#'@note this analysis is based on standard design matrices. The sampleTab should be prepared in a way that the levels of the factors are clear.
#'@seealso \code{\link{f.prepare.two.by.two.analysis}}
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.two.by.two.DESeq <- function(data, sampleTab, factorA, factorB, method = "default", batchCol = "") {
  require("DESeq2")
  helper <- f.prepare.two.by.two.analysis(sampleTab, factorA, factorB, batchCol)
  data <- round(data[,rownames(sampleTab)])
  dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleTab, design = formula(helper$formulaString))
  dds <- DESeq(dds)
  conRes <- lapply(helper$defaultContrasts, function(x) as.data.frame(results(dds, contrast = x)))
  conRes <- lapply(conRes, function(x) data.frame(generic = x$baseMean, logFC = x$log2FoldChange, pVal = x$pvalue, adjP = x$padj, row.names = rownames(x)))
  out <- list()
  for (cont in names(conRes)) {
    out[[cont]] <- c.DEGtab(tool = "DESeq", method = method, table = conRes[[cont]], isPairwise = FALSE, pairOrCont = as.character(cont))
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform a two-by-two factor/level analysis with EdgeR
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param factorA the name of the first factor in the model (a colname in sampleTab)
#'@param factorB the name of the second factor in the model (a colname in sampleTab)
#'@param method used for dispersion estimates: common, trended, or tagwise.
#'@param batchCol optionally the name of the column holding the info on the batch.
#'@return list with \code{\link{c.DEGtab}} objects, named after the name of the contrast.
#'@note this analysis is based on standard design matrices. The sampleTab should be prepared in a way that the levels of the factors are clear.
#'@seealso \code{\link{f.prepare.two.by.two.analysis}}
#'@examples \dontrun{
#'rDir <- "/Users/marc/Projects/Raffaella/1941_RNA-Seq"
#'data <-  read.table(file.path(rDir, "Rcount_comb.txt"), sep = '\t')
#'sampleTab <- read.csv(file.path(rDir, "sampleTab.csv"), stringsAsFactors = FALSE, row.names = 1)
#'sampleTab$BA <- paste0("B", rep(1:3, 4))
#'sampleTab$BA <- factor(sampleTab$BA, levels = c("B1", "B2", "B3"))
#'sampleTab$COL <- factor(sampleTab$COL, levels = c("adherent", "sphere"))
#'sampleTab$INH <- factor(sampleTab$INH, levels = c("DMSO", "ICR"))
#'res <- f.two.by.two.edgeR(data, sampleTab, "COL", "INH", "trended")
#'resBatch <- f.two.by.two.edgeR(data, sampleTab, "COL", "INH", "trended", "BA")
#'}
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.two.by.two.edgeR <- function(data, sampleTab, factorA, factorB, method = "trended", batchCol = "") {
  require("edgeR")
  if (!(method %in% c("common", "trended", "tagwise"))) { cat("undefined method:", method, "\n"); stop;}
  helper <- f.prepare.two.by.two.analysis(sampleTab, factorA, factorB, batchCol)
  data <- round(data[,rownames(sampleTab)])
  sampleTab$lib.sizes <- apply(data, 2, sum)
  data <- new("DGEList", list(counts = data, samples = sampleTab))
  # calculate dispersion
  if (method == "common") { data <- estimateGLMCommonDisp(data, helper$design) }
  if (method == "trended") { data <- estimateGLMTrendedDisp(data, helper$design) }
  if (method == "tagwise") { 
    data <- estimateGLMCommonDisp(data, helper$design)
    data <- estimateGLMTagwiseDisp(data, helper$design)
  }
  fit <- glmFit(data, helper$design)
  conRes <- lapply(helper$defaultContrasts, function(x) as.data.frame(glmLRT(fit, contrast = x)$table))
  conRes <- lapply(conRes, function(x) data.frame(generic = x$logCPM, logFC = x$logFC, pVal = x$PValue, adjP = p.adjust(x$PValue, "BH"), row.names = rownames(x)))
  out <- list()
  for (cont in names(conRes)) {
    out[[cont]] <- c.DEGtab(tool = "edgeR", method = method, table = conRes[[cont]], isPairwise = FALSE, pairOrCont = as.character(cont))
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform a two-by-two factor/level analysis with limma
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param factorA the name of the first factor in the model (a colname in sampleTab)
#'@param factorB the name of the second factor in the model (a colname in sampleTab)
#'@param method used for dispersion estimates: none, TMM, or quantile.
#'@param batchCol optionally the name of the column holding the info on the batch.
#'@return list with \code{\link{c.DEGtab}} objects, named after the name of the contrast.
#'@note this analysis is based on standard design matrices. The sampleTab should be prepared in a way that the levels of the factors are clear.
#'@seealso \code{\link{f.prepare.two.by.two.analysis}}
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.two.by.two.limma <- function(data, sampleTab, factorA, factorB, method = "TMM", batchCol = "") {
  require("limma")
  if (!(method %in% c("none", "TMM", "quantile"))) { cat("undefined method:", method, "\n"); stop;}
  helper <- f.prepare.two.by.two.analysis(sampleTab, factorA, factorB, batchCol)
  data <- round(data[,rownames(sampleTab)])
  if (method == "TMM") {
    data <- DGEList(data)
    data <- calcNormFactors(data)
    v <- voom(data, helper$design)
  } else {
    v <- voom(data, helper$design, normalize.method = method)
  }
  # fit and test
  fit <- lmFit(v, helper$design)
  fit <- contrasts.fit(fit, helper$defaultContrasts)
  fit <- eBayes(fit)
  # extract results
  conRes <- list()
  for (cont in colnames(helper$defaultContrasts)) {
    conRes[[cont]] <- topTable(fit, coef = cont, number = Inf, adjust = "BH")
  }
  # and simplify
  conRes <- lapply(conRes, function(x) data.frame(generic = x$AveExpr, logFC = x$logFC, pVal = x$P.Value, adjP = x$adj.P.Val, row.names = rownames(x)))
  out <- list()
  for (cont in names(conRes)) {
    out[[cont]] <- c.DEGtab(tool = "limma", method = method, table = conRes[[cont]], isPairwise = FALSE, pairOrCont = as.character(cont))
  }
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Perform a series of two-by-two factor/level tests.
#'@param data a data.frame or matrix object with counts.
#'@param sampleTab a data.frame with grouping variables:
#'rownames are samples, colnames are variables for the model, and entries are factor levels.
#'only the variable in the model is relevant - and this should only be one. 
#'@param factorA the name of the first factor in the model (a colname in sampleTab)
#'@param factorB the name of the second factor in the model (a colname in sampleTab)
#'@param batchCol optionally the name of the column holding the info on the batch.
#'@return a list with results from the individual test methods. 
#'Names are edgeR_common, edgeR_tagwise, edgeR_trended, DESeq_default, limma_none, limma_TMM, and limma_quantile.
#'@note this analysis is based on standard design matrices. The sampleTab should be prepared in a way that the levels of the factors are clear.
#'Requires edgeR, DESeq2, and limma to be installed.
#'@references 
#'Robinson, M.D. and Oshlack, A. (2010)
#'A scaling normalization method for differential expression analysis of RNA-Seq data. \emph{Genome Biology} \bold{11}, R25.
#'Love, M.I. and Huber, W. and Anders, S. (2014)
#'Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. \emph{Genome Biology} \bold{15}, 550.
#'Ritchie, M.E. and Phipson, B. and Wu, D. and Hu, Y. and Law, C.W. and Shi, W. and Smyth, G. (2015)
#'limma powers differential expression analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research} \bold{43}, e47.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@seealso \code{\link{f.multi.level.edgeR}}, \code{\link{f.multi.level.DESeq}}, and \code{\link{f.multi.level.limma}}.
#'@export
f.multiple.two.by.two.comparisons <- function(data, sampleTab, factorA, factorB, batchCol = "") {
  data <- data[,rownames(sampleTab)]
  out <- list(
    edgeR_common = f.two.by.two.edgeR(data, sampleTab, factorA, factorB, "common", batchCol),
    edgeR_tagwise = f.two.by.two.edgeR(data, sampleTab, factorA, factorB, "tagwise", batchCol),
    edgeR_trended = f.two.by.two.edgeR(data, sampleTab, factorA, factorB, "trended", batchCol),
    DESeq_default = f.two.by.two.DESeq(data, sampleTab, factorA, factorB, "default", batchCol),
    limma_none = f.two.by.two.limma(data, sampleTab, factorA, factorB, "none", batchCol),
    limma_TMM = f.two.by.two.limma(data, sampleTab, factorA, factorB, "TMM", batchCol),
    limma_quantile = f.two.by.two.limma(data, sampleTab, factorA, factorB, "quantile", batchCol)
  )
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Write out an overview, the individual results, and a table with gene descriptions to a workbook.
#'@param DEGtabList a list with \code{\link{c.DEGtab}} objects. Names will be used as sheet names.
#'@param rDir directory where the workbook will be saved.
#'@param filePrefix a prefix for the name of the workbook ("_DE_results.xlsx" will be added).
#'@param mart the mart to use (e.g. useMart("ensembl", dataset = "mmusculus_gene_ensembl")).
#'@param addGeneCol set to true if the rownames in data were of type ensembl_transcript_id.
#'@param onlySigRawP set to true if only entires with a raw pValue < 0.05 should be saved.
#'@note requires XLConnect to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.write.DEGtabs.to.workbook <- function(DEGtabList, rDir, filePrefix, mart = NA, addGeneCol = FALSE, onlySigRawP = FALSE) {
  require("XLConnect")
  # load workbook 
  wb <- loadWorkbook(file.path(rDir, paste0(filePrefix, "_DE_results.xlsx")), create = TRUE)
  # check what the rownames are
  contentName <- ifelse(addGeneCol, "transcripts", "genes")
  if (is.na(mart)) {contentName <- "feature"}
  # write the overview tables into one sheet
  cat("writing overview\n")
  sheetName <- "overview"
  createSheet(wb, name = sheetName)
  placeOnRow <- 1
  for (resName in names(DEGtabList)) {
    DEGtab <- DEGtabList[[resName]]
    temp <- DEGtab$get_stats()[,,"simple"]
    placeRef <- paste("", sheetName, resName, sep = '_')
    createName(wb, name = placeRef, formula = paste0(sheetName, "!$A$", placeOnRow))
    cornerField <- ifelse(resName == paste(DEGtab$tool, DEGtab$method, sep = '_'), paste0(resName, ": ", DEGtab$corner), DEGtab$corner) # if it is a method collection, this should be added
    writeNamedRegion(wb, temp, name = placeRef, rownames = cornerField)
    placeOnRow <- placeOnRow + nrow(temp) + 2
  }
  for (i in 1:(ncol(temp)+1)) {setColumnWidth(wb, sheet = sheetName, column = i)}
  # write the detailed test results
  cat("writing detailed test results\n")
  allEntries <- c()
  placeOnRow <- 1
  for (resName in names(DEGtabList)) {
    DEGtab <- DEGtabList[[resName]]
    temp <- DEGtab$get_print_table()
    if (onlySigRawP) {
      toKeep <- DEGtab$get_significant_entries(0.05, 1, 0)
      temp <- temp[toKeep$any,]
    }
    allEntries <- union(allEntries, rownames(temp))
    sheetName <- paste0("_", resName)
    if (nchar(sheetName) > 31) { sheetName <- substr(sheetName, 1, 31); cat("limited a name to 31 characters\n"); }
    createSheet(wb, name = sheetName)
    placeRef <- paste("", sheetName, contentName, sep = '_')
    createName(wb, name = placeRef, formula = paste0(sheetName, "!$A$", placeOnRow))
    writeNamedRegion(wb, temp, name = placeRef, rownames = contentName)
    for (i in 1:(ncol(temp)+1)) {setColumnWidth(wb, sheet = sheetName, column = i)}
  }
  # write a sheet with the gene annotations
  if (!is.na(mart)) {
    cat("searching and writing gene descriptions\n")
    entryType <- ifelse(addGeneCol, "ensembl_transcript_id", "ensembl_gene_id")
    geneDescriptions <- getGene(allEntries, type = entryType, mart = mart)
    symbolName <- grep("_symbol", colnames(geneDescriptions), value = TRUE)
    if (addGeneCol) {
      colsToTake <- c("ensembl_transcript_id", "ensembl_gene_id", symbolName, "description")
    } else {
      colsToTake <- c("ensembl_gene_id", symbolName, "description")
    }
    geneDescriptions <- geneDescriptions[,colsToTake]
    sheetName <- "gene_description"
    createSheet(wb, name = sheetName)
    placeOnRow <- 1
    placeRef <- paste("", sheetName, contentName, sep = '_')
    createName(wb, name = placeRef, formula = paste0(sheetName, "!$A$", placeOnRow))
    writeNamedRegion(wb, geneDescriptions, name = placeRef, rownames = NULL)
    for (i in 1:(ncol(geneDescriptions)-1)) {setColumnWidth(wb, sheet = sheetName, column = i)} # no scaling for descriptions
  }
  # save, clear, and gc
  saveWorkbook(wb)
  rm(wb)
  gc()
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Write out set collections to a workbook.
#'@param data any table with some data.
#'@param setCollection the list returned by \code{\link{f.combine.a.list.of.sets}}
#'@param rDir directory where things are stored.
#'@param filePrefix some prefix for the files.
#'@param mart the mart to use (e.g. useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
#'@param addGeneCol set to true if the rownames in data were of type ensembl_transcript_id.
#'@note requires XLConnect to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.write.set.collection.to.workbook <- function(data, setCollection, rDir, filePrefix, mart = NA, addGeneCol = FALSE) {
  require("XLConnect")
  # load workbook 
  wb <- loadWorkbook(file.path(rDir, paste0(filePrefix, "_SET_collection.xlsx")), create = TRUE)
  # check what the rownames are
  contentName <- ifelse(addGeneCol, "transcripts", "genes")
  if (is.na(mart)) {contentName <- "feature"}
  # retrieve annotation
  cat("downloading annotation\n")
  if (!is.na(mart)) {
    entryType <- ifelse(addGeneCol, "ensembl_transcript_id", "ensembl_gene_id")
    geneDescriptions <- getGene(rownames(data), type = entryType, mart = mart)
    symbolName <- grep("_symbol", colnames(geneDescriptions), value = TRUE)
    if (addGeneCol) {
      colsToTake <- c("ensembl_transcript_id", "ensembl_gene_id", symbolName, "description")
    } else {
      colsToTake <- c("ensembl_gene_id", symbolName, "description")
    }
    geneDescriptions <- geneDescriptions[,colsToTake]
  } else {
    geneDescriptions <- NA
  }
  # write gene data and description for all the sets
  cat("writing gene sets\n")
  placeOnRow <- 1
  for (geneSetName in names(setCollection)) {
    if (!(length(setCollection[[geneSetName]]) > 1)) {
      cat("skipping", geneSetName, "because it is too small\n")
      next
    }
    if (!is.na(geneDescriptions)) {
      descMask <- geneDescriptions[,1] %in% setCollection[[geneSetName]]
      geneNames <- geneDescriptions[,1][descMask]
      temp <- cbind(data[geneNames,], geneDescriptions[descMask,])
    } else {
      temp <- data[setCollection[[geneSetName]],]
    }
    sheetName <- paste0("_", geneSetName)
    if (nchar(sheetName) > 31) { sheetName <- substr(sheetName, 1, 31); cat("limited a name to 31 characters\n"); }
    createSheet(wb, name = sheetName)
    placeRef <- paste("", sheetName, contentName, sep = '_') # it's not allower to have a number at beginning...
    createName(wb, name = placeRef, formula = paste0(sheetName, "!$A$", placeOnRow))
    writeNamedRegion(wb, temp, name = placeRef, rownames = contentName)
    for (i in 1:ncol(temp)) {setColumnWidth(wb, sheet=sheetName, column=i)}
  }
  # save, clear, and gc
  saveWorkbook(wb)
  rm(wb)
  gc()
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Write out set collections to a workbook.
#'@param data any table with data
#'@param setCollection a list of sets - can be for example the list returned by \code{\link{f.combine.a.list.of.sets}}
#'or \code{\link{f.extract.significant.entries.pairwise}}. Nested lists should be unlisted first.
#'@param rDir directory where things are stored.
#'@param filePrefix some prefix for the files.
#'@param mart the mart to use (e.g. useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
#'@param addGeneCol set to true if the rownames in data were of type ensembl_transcript_id.
#'@note requires XLConnect to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.write.multicolumn.sets.to.workbook <- function(data, setCollection, rDir, filePrefix, mart = NA, addGeneCol = FALSE) {
  require("XLConnect")
  # check what the rownames are
  contentName <- ifelse(addGeneCol, "transcripts", "genes")
  if (is.na(mart)) {contentName <- "feature"}
  # combine table and sets
  cat("combining tables and sets - removing genes not present in any set:\n")
  data <- as.data.frame(data)
  placeOnRow <- 1
  for (geneSetName in names(setCollection)) {
    data[[geneSetName]] <- 0
    if (length(setCollection[[geneSetName]]) > 0) {
      data[setCollection[[geneSetName]], geneSetName] <- 1
    }
  }
  before <- nrow(data)
  data <- data[apply(data[,names(setCollection)], 1, sum) > 0,]
  cat(before - nrow(data), "- so there are still", nrow(data), "left\n")
  # retrieve annotation
  if (!is.na(mart)) {
    cat("downloading annotation\n")
    entryType <- ifelse(addGeneCol, "ensembl_transcript_id", "ensembl_gene_id")
    geneDescriptions <- getGene(rownames(data), type = entryType, mart = mart)
    symbolName <- grep("_symbol", colnames(geneDescriptions), value = TRUE)
    if (addGeneCol) {
      colsToTake <- c("ensembl_transcript_id", "ensembl_gene_id", symbolName, "description")
    } else {
      colsToTake <- c("ensembl_gene_id", symbolName, "description")
    }
    geneDescriptions <- geneDescriptions[,colsToTake]
    data <- cbind(data[geneDescriptions[,1],], geneDescriptions)
  }
  # load workbook and write everything
  wb <- loadWorkbook(file.path(rDir, paste0(filePrefix, "_SET_multicolumn.xlsx")), create = TRUE)
  cat("writing to workbook\n")
  sheetName <- "data"
  createSheet(wb, name = sheetName)
  placeOnRow <- 1
  placeRef <- paste("", sheetName, contentName, sep = '_')
  createName(wb, name = placeRef, formula = paste0(sheetName, "!$A$", placeOnRow))
  writeNamedRegion(wb, data, name = placeRef, rownames = contentName)
  for (i in 1:ncol(data)) {setColumnWidth(wb, sheet=sheetName, column=i)}
  # save, clear, and gc
  saveWorkbook(wb)
  rm(wb)
  gc()
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Combine individual DE tabs and write all into one sheet.
#'@param data any table with additional data.
#'@param DEGtabList a list with \code{\link{c.DEGtab}} objects.
#'@param rDir directory where things are stored.
#'@param filePrefix some prefix for the files.
#'@param pVal threshold for the pValue
#'@param adjP threshold for the adjusted pValue
#'@param logFC threshold for the logFC
#'@param mart the mart to use (e.g. useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
#'@param addGeneCol set to true if the rownames in data were of type ensembl_transcript_id.
#'@note requires XLConnect to be installed.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.combine.DEGtabs.and.write.to.workbook <- function(data, DEGtabList, rDir, filePrefix, pVal = 1, adjP = 0.1, logFC = 0, mart = NA, addGeneCol = FALSE) {
  require("XLConnect")
  # check what the rownames are
  contentName <- ifelse(addGeneCol, "transcripts", "genes")
  if (is.na(mart)) {contentName <- "feature"}
  # combine table and sets
  cat("combining tables and sets - removing genes not present in any set:\n")
  data <- as.data.frame(data)
  placeOnRow <- 1
  for (compName in names(DEGtabList)) {
    temp <- DEGtabList[[compName]]$get_table()
    newCols <- paste0(compName, '.', colnames(temp))
    data[,newCols] <- NA
    sigEntries <- DEGtabList[[compName]]$get_significant_entries(pVal, adjP, logFC)$any
    temp <- temp[sigEntries,]
    if (nrow(temp) > 0) { data[rownames(temp), newCols] <- temp }
  }
  before <- nrow(data)
  checkNA <- paste0(names(DEGtabList), '.', "logFC")
  data <- data[apply(!is.na(data[,checkNA]), 1, sum) > 0,]
  cat(before - nrow(data), "- so there are still", nrow(data), "left\n")
  # retrieve annotation
  if (!is.na(mart)) {
    cat("downloading annotation\n")
    entryType <- ifelse(addGeneCol, "ensembl_transcript_id", "ensembl_gene_id")
    geneDescriptions <- getGene(rownames(data), type = entryType, mart = mart)
    symbolName <- grep("_symbol", colnames(geneDescriptions), value = TRUE)
    if (addGeneCol) {
      colsToTake <- c("ensembl_transcript_id", "ensembl_gene_id", symbolName, "description")
    } else {
      colsToTake <- c("ensembl_gene_id", symbolName, "description")
    }
    geneDescriptions <- geneDescriptions[,colsToTake]
    data <- cbind(data[geneDescriptions[,1],], geneDescriptions)
  }
  # load workbook and write everything
  wb <- loadWorkbook(file.path(rDir, paste0(filePrefix, "_DEresults_multicolumn.xlsx")), create = TRUE)
  cat("writing to workbook\n")
  sheetName <- "data"
  createSheet(wb, name = sheetName)
  placeOnRow <- 1
  placeRef <- paste("", sheetName, contentName, sep = '_')
  createName(wb, name = placeRef, formula = paste0(sheetName, "!$A$", placeOnRow))
  writeNamedRegion(wb, data, name = placeRef, rownames = contentName)
  for (i in 1:ncol(data)) {setColumnWidth(wb, sheet=sheetName, column=i)}
  # save, clear, and gc
  saveWorkbook(wb)
  rm(wb)
  gc()
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Get chromosomes which should be kept.
#'@return a character vector with valid chromosomes
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.internal.get.relevant.chromosomes <- function(chromList) {
  out <- c("MT", "M", "C", "X", "Y", as.character(1:100))
  out <- c(out, paste("Chr", out, sep = ''), paste("chr", out, sep =''))
  chromList <- unique(chromList)
  alsoOK <- grep("[[:alpha:]]{2}[[:digit:]]{6}\\.[[:digit:]]{1}", chromList, value = TRUE)
  out <- c(out, alsoOK)
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Retrieve promoter regions.
#'@param mart the mart to use (e.g. useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
#'@param rDir directory where the file is stored
#'@param us upstream bp
#'@param ds downstream bp
#'@return the promoter regions
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.promoter.regions <- function(mart, rDir, us = 2e3, ds = 5e2) {
  require("biomaRt")
  toRetrieve <- c("ensembl_gene_id",
                  "gene_biotype",
                  "ensembl_transcript_id",
                  "transcript_biotype",
                  "chromosome_name",
                  "strand", 
                  "transcription_start_site")
  data <- getBM(attributes = toRetrieve, mart = mart)
  out <- data.frame(GeneID = data$ensembl_gene_id,
                    Chr = data$chromosome_name,
                    Start = data$transcription_start_site - data$strand*us,
                    End = data$transcription_start_site + data$strand*ds, 
                    Strand = ifelse(data$strand == 1, "+", "-"), 
                    GeneType = data$gene_biotype,
                    TransType = data$transcript_biotype,
                    stringsAsFactors = FALSE, row.names = data$ensembl_transcript_id)
  out[out$Strand == "-", c("Start", "End")] <- out[out$Strand == "-", c("End", "Start")]
  # remove weird chromosomes
  toKeep <- f.internal.get.relevant.chromosomes(unique(out$Chr))
  before <- nrow(out)
  out <- out[out$Chr %in% toKeep,]
  cat("removed", before-nrow(out), "entries which were not on regular chromosomes\n")
  # save to outfile
  outFile <- file.path(rDir, paste("promoter_", us, "_", ds, ".txt", sep = ''))
  write.table(out, outFile, sep = '\t', quote = FALSE)
  cat("saved promoters in", outFile, "\n")
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Retrieve TES regions.
#'@param mart the mart to use (e.g. useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
#'@param rDir directory where the file is stored
#'@param us upstream bp
#'@param ds downstream bp
#'@return the TES regions
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.TES.regions <- function(mart, rDir, us = 2e3, ds = 2e3) {
  require("biomaRt")
  toRetrieve <- c("ensembl_gene_id",
                  "gene_biotype",
                  "ensembl_transcript_id",
                  "transcript_biotype",
                  "chromosome_name",
                  "strand", 
                  "transcript_end")
  data <- getBM(attributes = toRetrieve, mart = mart)
  out <- data.frame(GeneID = data$ensembl_gene_id,
                    Chr = data$chromosome_name,
                    Start = data$transcript_end - data$strand*us,
                    End = data$transcript_end + data$strand*ds, 
                    Strand = ifelse(data$strand == 1, "+", "-"), 
                    GeneType = data$gene_biotype,
                    TransType = data$transcript_biotype,
                    stringsAsFactors = FALSE, row.names = data$ensembl_transcript_id)
  out[out$Strand == "-", c("Start", "End")] <- out[out$Strand == "-", c("End", "Start")]
  # remove weird chromosomes
  toKeep <- f.internal.get.relevant.chromosomes(unique(out$Chr))
  before <- nrow(out)
  out <- out[out$Chr %in% toKeep,]
  cat("removed", before-nrow(out), "entries which were not on regular chromosomes\n")
  # save to outfile
  outFile <- file.path(rDir, paste("TES_", us, "_", ds, ".txt", sep = ''))
  write.table(out, outFile, sep = '\t', quote = FALSE)
  cat("saved TES in", outFile, "\n")
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Retrieve gene regions (from gene start to gene end).
#'@param mart the mart to use (e.g. useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
#'@param rDir directory where the file is stored
#'@return the gene regions
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.gene.regions <- function(mart, rDir) {
  require("biomaRt")
  toRetrieve <- c("ensembl_gene_id",
                  "gene_biotype",
                  "chromosome_name",
                  "start_position",
                  "end_position",
                  "strand")
  data <- getBM(attributes = toRetrieve, mart = mart)
  out <- data.frame(GeneID = data$ensembl_gene_id,
                    Chr = data$chromosome_name,
                    Start = data$start_position,
                    End = data$end_position, 
                    Strand = ifelse(data$strand == 1, "+", "-"), 
                    GeneType = data$gene_biotype,
                    stringsAsFactors = FALSE, row.names = data$ensembl_gene_id)
  # remove weird chromosomes
  toKeep <- f.internal.get.relevant.chromosomes(unique(out$Chr))
  before <- nrow(out)
  out <- out[out$Chr %in% toKeep,]
  cat("removed", before-nrow(out), "entries which were not on regular chromosomes\n")
  # save to outfile
  outFile <- file.path(rDir, paste("genes.txt", sep = ''))
  write.table(out, outFile, sep = '\t', quote = FALSE)
  cat("saved genes in", outFile, "\n")
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Retrieve transcript regions (from transcript start to transcript end).
#'@param mart the mart to use (e.g. useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
#'@param rDir directory where the file is stored
#'@return the transcript regions
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.transcript.regions <- function(mart, rDir) {
  require("biomaRt")
  toRetrieve <- c("ensembl_gene_id",
                  "gene_biotype",
                  "ensembl_transcript_id",
                  "transcript_biotype",
                  "transcript_start",
                  "transcript_end",
                  "chromosome_name",
                  "strand")
  data <- getBM(attributes = toRetrieve, mart = mart)
  out <- data.frame(GeneID = data$ensembl_gene_id,
                    Chr = data$chromosome_name,
                    Start = data$transcript_start,
                    End = data$transcript_end, 
                    Strand = ifelse(data$strand == 1, "+", "-"), 
                    GeneType = data$gene_biotype,
                    TransType = data$transcript_biotype,
                    stringsAsFactors = FALSE, row.names = data$ensembl_transcript_id)
  # remove weird chromosomes
  toKeep <- f.internal.get.relevant.chromosomes(unique(out$Chr))
  before <- nrow(out)
  out <- out[out$Chr %in% toKeep,]
  cat("removed", before-nrow(out), "entries which were not on regular chromosomes\n")
  # save to outfile
  outFile <- file.path(rDir, paste("transcripts.txt", sep = ''))
  write.table(out, outFile, sep = '\t', quote = FALSE)
  cat("saved transcripts in", outFile, "\n")
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Retrieve exon regions (from exon start to exon end).
#'@param mart the mart to use (e.g. useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
#'@param rDir directory where the file is stored
#'@return the transcript regions
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.get.exon.regions <- function(mart, rDir) {
  require("biomaRt")
  toRetrieve <- c("ensembl_gene_id",
                  "gene_biotype",
                  "ensembl_transcript_id",
                  "transcript_biotype",
                  "ensembl_exon_id",
                  "exon_chrom_start",
                  "exon_chrom_end",
                  "rank",
                  "chromosome_name",
                  "strand")
  data <- getBM(attributes = toRetrieve, mart = mart)
  out <- data.frame(GeneID = data$ensembl_gene_id,
                    Chr = data$chromosome_name,
                    Start = data$exon_chrom_start,
                    End = data$exon_chrom_end, 
                    Strand = ifelse(data$strand == 1, "+", "-"), 
                    GeneType = data$gene_biotype,
                    TransType = data$transcript_biotype,
                    TransID = data$ensembl_transcript_id,
                    ExonID = data$ensembl_exon_id,
                    Rank = data$rank,
                    stringsAsFactors = FALSE, row.names = NULL)
  # remove weird chromosomes
  toKeep <- f.internal.get.relevant.chromosomes(unique(out$Chr))
  before <- nrow(out)
  out <- out[out$Chr %in% toKeep,]
  cat("removed", before-nrow(out), "entries which were not on regular chromosomes\n")
  # rownames only here (exon_id is not unique)
  rownames(out) <- paste("generic", 1:nrow(out), sep = '_')
  # save to outfile
  outFile <- file.path(rDir, paste("exons.txt", sep = ''))
  write.table(out, outFile, sep = '\t', quote = FALSE)
  cat("saved exons in", outFile, "\n")
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
### functions related to the workflow
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Create Rcount project files
#'@param sampleName name of the sample
#'@param annotationXML path to the Rcount annotation file (.xml)
#'@param readsInfile path to the input bam-file (.bam)
#'@param readsOutfile path to the output bam-file (.bam)
#'@param countTableFile path to the output count table (.txt)
#'@param OUTprojectXML path to the Rcount project file (.xml)
#'@param useMulti: true or false (character!)
#'@param useStrand: false, sense or antisense
#'@param minReads: see note
#'@param maxDist: see note
#'@param minBelowMaxDist: see note
#'@return NULL - creates a Rcount project file.
#'@note Genes must have at least "minReads" reads in total and "minBelowMaxDist"_" reads
#'within the first "maxDist" bps at the 3' end (see reference for details).
#'The parameters may be chosen according to the total number of reads and
#'the genome/transcriptome size.
#'@references 
#'Schmid, M.W. and Grossniklaus, U. (2015)
#'Rcount: simple and flexible RNA-Seq read counting.\emph{Bioinformatics} \bold{31}: 436-437.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.write.Rcount.project.file <- function(sampleName, annotationXML, readsInfile, readsOutfile, countTableFile, OUTprojectXML, 
                                        useMulti = "true", useStrand = "false", minReads = 5, maxDist = 250, minBelowMaxDist = 1) {
  lines <- c('<?xml version="1.0" encoding="UTF-8"?>',
             '<p502project version="1.0">',
             paste0('<name entry="', sampleName,'">'),
             '<files entry="">',
             paste0('<dataBaseInfile entry="', annotationXML, '"/>'),
             '<dataBaseOutfile entry=""/>',
             paste0('<readsInfile entry="', readsInfile, '"/>'),
             paste0('<readsOutfile entry="', readsOutfile, '"/>'),
             paste0('<countTableFile entry="', countTableFile, '"/>'),
             '</files>',
             '<parameters entry="">',
             paste0('<multi entry="', useMulti, '"/>'),
             paste0('<stranded entry="', useStrand, '"/>'),
             paste0('<minReads entry="', minReads, '"/>'),
             paste0('<maxDist entry="', maxDist, '"/>'),
             paste0('<minBelowMaxDist entry="', minBelowMaxDist, '"/>'),
             '</parameters>',
             '<region entry="">',
             '<useRegion entry="false"/>',
             '<regionStartName entry=""/>',
             '<regionStart entry="0"/>',
             '<regionEndName entry=""/>',
             '<regionEnd entry="0"/>',
             '</region>',
             '<settings entry="">',
             '<indexStepSize entry="10000"/>',
             '<bufferSizeBAM entry="200000"/>',
             '<bufferSizeMAP entry="200000"/>',
             '<bufferSizeOUT entry="200000"/>',
             '</settings>',
             '</name>',
             '</p502project>')
  output <- file(OUTprojectXML, open = "w")
  for (outLine in lines) {
    writeLines(outLine, output)
  }
  close(output)
  return(NULL)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Sort several bam files in parallel
#'@param bamFiles a vector with the bam files to be sorted (with or without the directory).
#'@param bamDirIn path to a folder where the aligned reads (bamFiles) are stored.
#'@param bamOutDir path to a folder where the sorted bam files will be stored.
#'@param numCores number of parallel cores/threads.
#'@param maxMemPerCore the maximal amount of RAM used by one thread (in Mb).
#'@return a list with two vectors (success, fail). Files listed in out$fail were not sorted successfully.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.sort.bam.parallel <- function(bamFiles, bamDirIn, bamDirOut, numCores = 4, maxMemPerCore = 3600) {
  require("Rsamtools")
  require("parallel")
  
  # Create a list for mclapply
  toSort <- list()
  for (curFile in bamFiles) {
    curFile <- basename(curFile) # remove directory if present
    curName <- gsub("\\.bam$", "", curFile) # remove extension
    infile <- file.path(bamDirIn, curFile)
    outfileNoExt <- file.path(bamDirOut, paste0(curName, "_srt"))
    outfile <- paste0(outfileNoExt, ".bam")
    toSort[[curName]] <- c(infile, outfileNoExt, outfile)
  }
  
  # Internal function for sorting and indexing
  sortAndIndexBam <- function(entry, mmpc) {
    sortBam(entry[1], entry[2], maxMemory = mmpc)
    indexBam(entry[3])
  }
  
  # Sort the reads in the bam files according to 
  # chromosome and position, and create an index.
  cat("Sorting bam files...\n")
  cat("Note that warning messages containing\n")
  cat("[bam_sort_core] merging from X files\n")
  cat("can be ignored.\n")
  mclapply(toSort, function(x) sortAndIndexBam(x, maxMemPerCore), mc.cores = numCores)
  cat("Note that warning messages containing\n")
  cat("[bam_sort_core] merging from X files\n")
  cat("can be ignored.\n")
  
  # Check if all files are present
  out <- list(success = c(), fail = c())
  for (curEntry in toSort) {
    if (file.exists(curEntry[3]) & file.exists(paste0(curEntry[3], ".bai"))) {
      out$success <- c(out$success, curEntry[3])
    } else {
      out$fail <- c(out$fail, curEntry[1])
    }
  }

  # Return the list with successfully sorted and failed files.
  return(out)
}

#'@title Downloads samples from SRA
#'@param readsDir path to a folder where the raw read files will be stored.
#'@param metaDB path to the SRA database with the metadata ("SRAmetadb.sqlite"); will be downloaded if the file does not exist.
#'IMPORTANT: this file is sometimes not downloaded correctly or per se corrupt. If you encounter errors like "table does not exist",
#'it is a strong indication that something is wrong with the file. Try downloading it again and if the error persists, contact the 
#'developers of the SRAdb package.
#'@param mySamplesFile path to a file with SRA-Identifiers (one ID per line, nothing else). 
#'All short read data linked to these SRA-IDs will be downloaded.
#'@param autoAnnotation path to a *.csv file which will contain the automatically retrieved annotation from the SRA database.
#'This is the file which has to be edited by the user and which will guide the entire workflow. IMPORTANT:
#'Fill in all columns with entries "toBeProvided", check/modify the parameters and finally erase the columns starting with "SRA".
#'@param colsToAdd a vector with column names which should be added to the autoAnnotion file (e.g. c("GENOTYPE", "TREATMENT")).
#'These columns will need to be filled up by the user and are then used during the analysis of differential expression.
#'@return TRUE if there is no error.
#'@note This is a top level wrapper function - see github.com/MWSchmid/RNAseq_protocol for details.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.wf.download.from.SRA <- function(readsDir, metaDB, mySamplesFile, autoAnnotation, colsToAdd) {
  require("SRAdb")

  # Check for a local copy of the database and
  # download it if it does not exist.
  if (!file.exists(metaDB)) {
    cat("Could not fine a local copy of the SRA database.\n")
    cat("Downloading it now...\n")
    getSRAdbFile(dirname(metaDB))
  }
  
  # Open a connection to the database
  cat("Connecting to the SRA database...\n")
  sra_con <- dbConnect(SQLite(), metaDB)
  
  # Here we assume that you have a precompiled list of 
  # samples/experiments/studies which are of interest to you.
  # Read in the list of experiments you are interested in.
  cat("Reading the file with the SRA-IDs of interest...\n")
  samples <- scan(mySamplesFile, what = "character")
  
  # You may have a collection of experiments and individual runs. 
  # However, only individual runs can be downloaded. Note that 
  # the conversion of a certain SRA ID type to another requires
  # the different query types not to be mixed.
  cat("Unifying SRA-IDs (converting to SRRxxx)...\n")
  sampleTypes <- unique(substr(samples, 1, 3))
  splitSamples <- lapply(sampleTypes, function(x)
    grep(paste0("^", x), samples, value = TRUE))
  conversions <- lapply(splitSamples, function(x)
    sraConvert(x, sra_con = sra_con ))
  allSRAids <- do.call("rbind", conversions)
  
  # Retrieve the remaining metadata to determine
  # the sequencing platform, the read length, etc.
  cat("Retrieving information about the SRA-IDs...\n")
  SRAdata <- lapply(allSRAids$run, function(x)
    getSRA(search_terms = x, out_types = "sra", sra_con))
  SRAdata <- do.call("rbind", SRAdata)
  rownames(SRAdata) <- SRAdata$run # for data access
  SRAdata$approxRL <- SRAdata$bases/SRAdata$spots

  # Search for strand-specific samples with regular
  # expressions. Note that this gives no guarantee
  # for strand specificity - in any case, this
  # should be verified manually.
  cat("Trying to guess whether samples are strand-specific...\n")
  posExp <- "strand[[:space:]{0,1}|-]specific"
  negExp <- paste0("no[[:alpha:]{0,1}][[:space:]{0,1}|-]", posExp)
  strSpec <- apply(SRAdata, 1, function(x)
    length(grep(posExp, x)) > 0)
  notStSp <- apply(SRAdata, 1, function(x)
    length(grep(negExp, x)) > 0)
  strSpec[notStSp] <- FALSE
  sum(strSpec)  # the number of strand specific libraries
  SRAdata$strSpec <- strSpec

  # Check if there are FASTQ files available for the samples.
  # For the samples without FASTQ file, you will download the
  # SRA archives and convert them into FASTQ (this is a lot 
  # slower than downloading FASTQ files).
  cat("Checking for FASTQ and SRA availability...\n")
  fastqInfo <- getFASTQinfo(allSRAids$run, sra_con, 'ftp')
  samplesWithFastq <- unique(with(fastqInfo, run[!is.na(fastq_ID)]))
  noFastqSamples <- unique(with(fastqInfo, run[is.na(fastq_ID)]))
  sraInfo <- getSRAinfo(noFastqSamples, sra_con)
  noFastqButSRA <- unique(with(sraInfo, run[!is.na(ftp)]))
  noFiles <- unique(with(sraInfo, run[is.na(ftp)]))
  if (length(noFiles) > 0) {
    cat(length(noFiles), "samples do not have a FASTQ or SRA file:\n")
    cat(paste(noFiles, collapse = '\n'), '\n')
    cat("Visit http://www.ncbi.nlm.nih.gov/ and check the SRA-IDs listed above.\n")
    cat("Note: These samples will be removed from the sample table.\n")
    # SRAdata <- subset(SRAdata, !(run %in% noFiles)) # it's done automatically below
  }
  
  # Download the FASTQ files.
  # Note that an error like "no such table: fastq" is likely
  # due to a corrupt SRAdb file. Delete SRAmetadb.sqlite
  # and download it again. If the error persists, contact
  # the SRAdb developers.
  if (length(samplesWithFastq) > 0) {
    cat("Downloading FASTQ files...\n")
    getFASTQfile(samplesWithFastq, sra_con, readsDir, 'ftp')
  }

  # If necessary, download the SRA archives of the remaining
  # samples and convert them to FASTQ. Note that the conversion
  # requires the "SRA toolkit" to be installed:
  # www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
  if (length(noFastqButSRA) > 0) {
    cat("Downloading SRA files...\n")
    getSRAfile(noFastqButSRA, sra_con, readsDir, 'sra' )
    commands <- sapply(noFastqButSRA, function(x)
      paste0("fastq-dump --gzip --split-files -O ",
             readsDir, " ", file.path(readsDir, x), ".sra"))
    sapply(commands, system)
  }
  
  # Get a vector of all FASTQ file names and
  # remove the reverse reads of the PE samples.
  cat("Removing reverse reads of the paired-end samples...\n")
  allFastq <- list.files(readsDir, ".fastq.gz")
  revFastq <- list.files(readsDir, "_2.fastq.gz")
  onlyForwardFastq <- setdiff(allFastq, revFastq)
  sampleIDs <- gsub("_1|.fastq.gz", "", onlyForwardFastq)
  
  # Samples with paired-end sequencing have twice the
  # coverage - therefore adjust the approximated
  # read length of these samples.
  cat("Adjusting the read-length of the paired-end samples...\n")
  PES <- gsub("_2|.fastq.gz", "", revFastq)
  SRAdata[PES, "approxRL"] <- SRAdata[PES, "approxRL"]/2
  
  # Create a table with information for the alignment
  # and the preprocessing. Store the sampleID, read length,
  # platform type, FASTQ file names, sample/library names, 
  # strand specificity and additional attributes in a csv file. 
  # Even though the sample/library names and the attributes
  # are not complete for each sample, they may help renaming
  # all samples.
  cat("Generating table with the automatically retrieved annotation...\n")
  sortedSRAdata <- SRAdata[sampleIDs,] # this ensures that the sample has a fastq file
  tabForProcessing <- data.frame(
    sampleID   = sampleIDs,
    readLen    = round(sortedSRAdata$approxRL, 0),
    numReads   = sortedSRAdata$spots,
    platform   = sortedSRAdata$platform,
    fastqFile  = onlyForwardFastq,
    phredOffset = 33,
    useMulti   = "true",
    useStrand  = as.logical(sortedSRAdata$strSpec),
    minReads   = 5,
    maxDist    = 250,
    minBelowMaxDist = 1, 
    processedReads = 0,
    alignedReads = 0,
    processedReadsFC = 0,
    assignedReadsFC = 0,
    stringsAsFactors = FALSE
  )
  for (toAdd in colsToAdd) { tabForProcessing[[toAdd]] <- "toBeProvided" }
  tabForProcessing$mySampleName <-"toBeProvided"
  tabForProcessing$SRAsamName <- sortedSRAdata$sample_name
  tabForProcessing$SRAlibName <- sortedSRAdata$library_name
  tabForProcessing$SRAattrib <- sortedSRAdata$sample_attribute

  # Optional: Verify that the PHRED offset is 33 and not 64.
  # This is not strictly necessary as the files from SRA are
  # generally having an offset of 33 (they otherwise convert
  # the encoding). However, there are still some cases with
  # an offset of 64.
  cat("Verifying whether the PHRED-offsets are correct...\n")
  rownames(tabForProcessing) <- tabForProcessing$sampleID
  for (cSam in rownames(tabForProcessing)) {
    fastqFile <- file.path(readsDir,
                           tabForProcessing[cSam, "fastqFile"])
    phredOffset  <- f.check.phred.offset(fastqFile)
    tabForProcessing[cSam, "phredOffset"] <- phredOffset
  }
  
  # Save the table.
  cat("Writing the sample annotation table...\n")
  write.csv(tabForProcessing, autoAnnotation, row.names = FALSE)
  
  # Return TRUE for success.
  return(TRUE)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Align samples with Rsubread
#'@param myAnnotation path to the *.csv with the curated annotation (based on the file created by \code{\link{f.wf.download.from.SRA}})
#'@param readsDir path to a folder where the raw read files are stored.
#'@param bamDir path to a folder where the aligned reads (bam files) will be stored.
#'@param nucIndex path to the nucleotide index for Rsubread.
#'@param nucIndex path to the color index for Rsubread (if existing).
#'@param readType type of experiment the reads come from - either "rna" for RNAseq, or "dna" for genomic sequencing, ChIPseq, etc.
#'@param ... additional arguments passed on to \code{\link{align}}. IMPORTANT!
#'For example:
#'nthreads = 4, unique = FALSE, nBestLocations = 10
#'Note that the following \code{\link{align}} parameters are set internally and cannot be specified again:
#'index, readFile1, readFile2, output_file, type, phredOffset 
#'@return TRUE if there is no error.
#'@note This is a top level wrapper function - see github.com/MWSchmid/RNAseq_protocol for details.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.wf.align.reads.with.Rsubread <- function(myAnnotation, readsDir, bamDir, nucIndex, colIndex = "", readType = c("rna", "dna"), ...) {
  require("Rsubread")
  # Import the sample annotation.
  cat("Importing sample annotation...\n")
  sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
  rownames(sampleTab) <- sampleTab$sampleID
  
  # Align the reads - allow up to 10 alignments per read.
  # Note that we redirect the messages from the aligner
  # into a character vector to extract the number of
  # mapped reads. However, the drawback is that we 
  # will not see the progress messages from the aligner
  # in real-time.
  cat("Aligning reads to the genome...\n")
  for (cSam in rownames(sampleTab)) {
    cat("[ALIGN] processing", cSam, "...\n")
    isSOLiD <- (sampleTab[cSam, "platform"] == "ABI_SOLID")
    pOff <- sampleTab[cSam, "phredOffset"]
    cIdx <- ifelse(isSOLiD, colIndex, nucIndex)
    inFile <- file.path(readsDir, sampleTab[cSam, "fastqFile"])
    outFile <- file.path(bamDir, paste0(cSam, ".bam"))
    if ("fastqFileReverse" %in% colnames(sampleTab)) {
      revReads <- file.path(readsDir, sampleTab[cSam, "fastqFileReverse"])
    } else {
      revReads <- NULL
    }
    if (file.exists(outFile) & (sampleTab[cSam,"processedReads"]>0)) {next} # TODO
    stdout <- capture.output(
      align(cIdx, readfile1 = inFile, readfile2 = revReads,
            output_file = outFile, type = readType,
            phredOffset = pOff, ...)
    )
    # extract the number of processed and mapped reads
    readPattern <- "(Processed|Mapped).*reads"
    readMatches <- regexpr(readPattern, stdout)
    readCountsStr <- regmatches(stdout, readMatches)
    numberPattern <- "[,0123456789]+"
    numberMatches <- regexpr(numberPattern, readCountsStr)
    readCounts <- regmatches(readCountsStr, numberMatches)
    sampleTab[cSam, "processedReads"] <- gsub(",", "", readCounts[1])
    sampleTab[cSam, "alignedReads"] <- gsub(",", "", readCounts[2])
    write.csv(sampleTab, myAnnotation, row.names = FALSE) # TODO
  }
  # calculate alignment percentages
  sampleTab$alPerc <- 100*with(sampleTab, as.numeric(alignedReads)/as.numeric(processedReads))
  
  # Save the table again to store the number of processed
  # and mapped reads as well. "processedReads" should be
  # equal to "numReads", and the number of "alignedReads"
  # may give you an indication whether the samples are
  # fine (note that alignment rates around 40 % are not
  # necessarily a problem - it depends on the protocol).
  cat("Updating sample annotation table (adding number of aligned reads).\n")
  write.csv(sampleTab, myAnnotation, row.names = FALSE)
  
  # Draw a histogram with the alignment percentages.
  cat("Drawing a histogram with the alignment rates.\n")
  counter <- 0
  histoGramFig <- file.path(dirname(), paste0("alignmentRates_", counter, ".png"))
  while (file.exists(histoGramFig)) {
    counter <- counter + 1
    histoGramFig <- file.path(dirname(), paste0("alignmentRates_", counter, ".png"))
  }
  cat("see:", histoGramFig, "\n")
  png(histoGramFig)
  hist(sampleTab$alPerc, xlab = "alignment rate (in %)",
       ylab = "no. samples", main = "", las = 1)
  dev.off()

  # Display the samples with a low alignment rate.
  cat("The following samples have an alignment rate below 20 %.\n")
  colsOfInt <- c("readLen", "platform", "mySampleName", "alPerc")
  print(subset(sampleTab, alPerc < 20)[,colsOfInt])
  
  # Return TRUE for success.
  return(TRUE)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Count reads per gene with featureCounts
#'@param myAnnotation path to the *.csv with the curated annotation (based on the file created by \code{\link{f.wf.download.from.SRA}})
#'@param bamDir path to a folder where the aligned reads (bam files) are stored.
#'@param gAnnFile path to the file with the genome annotation (GTF/GFF/SAF format).
#'@param countTabsFC path to a folder where the individual count tables (one per sample) will be stored. Note that this folder
#'should not contain any other files (at least no other *.csv files).
#'@param countTabFC path to a *.csv table which will contain the expression values for all samples.
#'@param minOverlapFrac miminal fraction of a read overlapping with a genomic feature to be assigned (e.g., if 0.7 with 100 bp reads,
#'at least 70 bps need to be within a genomic feature - if not, the read is not assigned).
#'@param ... additional arguments passed on to \code{\link{featureCounts}}. IMPORTANT!
#'For example: 
#'nthreads = 4, allowMultiOverlap = TRUE, largestOverlap = TRUE
#'Note that the following \code{\link{featureCounts}} parameters are set internally and cannot be specified again:
#'files, annot.ext, isGTFAnnotationFile, strandSpecific, minOverlap
#'@return TRUE if there is no error.
#'@note This is a top level wrapper function - see github.com/MWSchmid/RNAseq_protocol for details.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.wf.count.reads.with.featureCounts <- function(myAnnotation, bamDir, gAnnFile, countTabsFC, countTabFC,
                                                minOverlapFrac = .7, ...) {
  require("Rsubread")

  # Import the sample annotation.
  cat("Importing sample annotation...\n")
  sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
  rownames(sampleTab) <- sampleTab$sampleID
  
  # Count the number of (unique!) reads with featureCounts().
  # I strongly recommend to read the documentation of this
  # function (type ?featureCounts).
  # Note that we redirect the messages from featureCounts
  # into a character vector to extract the number of
  # assigned reads. However, the drawback is that we 
  # will not see the progress messages in real-time.
  cat("Assigning reads to genomic features...\n")
  sT <- c("false" = 0, "sense" = 1, "antisense" = 2)
  for (cSam in rownames(sampleTab)) {
    cat("[ASSIGN] processing", cSam, "...\n")
    bamFile <- file.path(bamDir, paste0(cSam, ".bam"))
    outFile <- file.path(countTabsFC, paste0(cSam, ".txt"))
    SSmode <- tolower(sampleTab[cSam, "useStrand"])
    SSmode <- sT[SSmode]
    minBPinExon <- ceiling(minOverlapFrac*sampleTab[cSam, "readLen"])
    isGTForGFF <- grepl("\\.gtf$|\\.gff$", gAnnFile)
    stdout <- capture.output(
      countRes <- featureCounts(bamFile,
                                annot.ext = gAnnFile,       # reference annotation (.gtf|.gff|.saf)
                                isGTFAnnotationFile = isGTForGFF, # TRUE if ref. anno. is ".gtf" of ".gff"
                                strandSpecific = SSmode,    # set strand-specificity
                                minOverlap = minBPinExon,
                                ...)
    )
    # Save the count table for the current sample.
    colnames(countRes$counts) <- cSam
    write.table(countRes$counts, outFile, sep = '\t', quote = FALSE)
    # extract the number of processed and assigned reads
    readPattern <- "(Total|Successfully).*reads"
    readCountsStr <- grep(readPattern, stdout, value = TRUE)
    numberPattern <- "[0123456789]+"
    numberMatches <- regexpr(numberPattern, readCountsStr)
    readCounts <- regmatches(readCountsStr, numberMatches)
    sampleTab[cSam, "processedReadsFC"] <- readCounts[1]
    sampleTab[cSam, "assignedReadsFC"] <- readCounts[2]
  }
  # calculate assignment percentages
  sampleTab$asPercFC <- 100*with(sampleTab, as.numeric(assignedReadsFC)/as.numeric(processedReadsFC))
  
  # Save the sample annotation table again to store the 
  # number of processed and successfully assigned reads.
  # The number of "processedReadsFC" should equal to
  # "alignedReads", and the number of "assignedReadsFC"
  # may give you an indication whether a sample has
  # too many reads in intergenic regions.
  cat("Updating sample annotation table (adding number of assigned reads).\n")
  write.csv(sampleTab, myAnnotation, row.names = FALSE)
  
  # Read-in all the tables and merge the files which
  # correspond to only one sample (some samples were
  # sequenced on multiple lanes). Windows users must
  # switch to the "native Windows R" to run the code
  # below.
  cat("Merging split samples and combining the individual count tables...\n")
  counts_FC <- f.read.featureCounts(countTabsFC)
  counts_FC <- f.summarize.columns(
    counts_FC,
    data.frame(sample = rownames(sampleTab),
               group = sampleTab$mySampleName),
    sum)
  write.csv(counts_FC, countTabFC)
  
  # Draw a histogram with the assigment percentages.
  cat("Drawing a histogram with the assignment rates.\n")
  counter <- 0
  histoGramFig <- file.path(dirname(), paste0("assignmentRatesFC_", counter, ".png"))
  while (file.exists(histoGramFig)) {
    counter <- counter + 1
    histoGramFig <- file.path(dirname(), paste0("assignmentRatesFC_", counter, ".png"))
  }
  cat("see:", histoGramFig, "\n")
  png(histoGramFig)
  hist(sampleTab$asPerc, xlab = "assignment rate (in %)",
       ylab = "no. samples", main = "", las = 1)
  dev.off()
  
  # Display the samples with a low assignment rate.
  cat("The following samples have an assignment rate below 50 %.\n")
  colsOfInt <- c("readLen", "platform", "mySampleName", "asPerc")
  print(subset(sampleTab, asPerc < 50)[,colsOfInt])
  
  # Return TRUE for success.
  return(TRUE)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Sort bam files for Rcount
#'@param myAnnotation path to the *.csv with the curated annotation (based on the file created by \code{\link{f.wf.download.from.SRA}})
#'@param bamDirIn path to a folder where the aligned reads (bam files) are stored.
#'@param bamDirOut path to a folder where the sorted bam files will be stored.
#'@param numCores number of parallel cores/threads.
#'@param maxMemPerCore the maximal amount of RAM used by one thread (in Mb).
#'@return TRUE if there is no error.
#'@note This is a top level wrapper function - see github.com/MWSchmid/RNAseq_protocol for details.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.wf.sort.bam.files <- function(myAnnotation, bamDirIn, bamDirOut = bamDirIn, numCores = 4, maxMem = 3600) {
  require("Rsamtools")
  
  # Import the sample annotation.
  cat("Importing sample annotation...\n")
  sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
  rownames(sampleTab) <- sampleTab$sampleID
  
  # Sort the reads in the bam files according to 
  # chromosome and position, and create an index.
  bamFiles <- paste0(rownames(sampleTab), ".bam")
  sorted <- f.sort.bam.parallel(bamFiles, bamDirIn, bamDirOut,
                                numCores, maxMemPerCore)

  # Check if all files were sorted and indexed successfully.
  if (length(sorted$fail) > 0) {
    cat("some files were not sorted or indexed.\n")
    cat(paste0(sorted$fail, collapse = '\n'), '\n')
    out <- FALSE
  } else {
    out <- TRUE
  }
  
  return(out)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Set weights of reads with multiple alignments using Rcount
#'@param myAnnotation path to the *.csv with the curated annotation (based on the file created by \code{\link{f.wf.download.from.SRA}})
#'@param RcountMultir path to the Rcount-multireads executable
#'@param bamDirIn path to a folder where the sorted bam files are stored.
#'@param bamDirOut path to a folder where the bam files with the weighted reads will be stored.
#'@return TRUE if there is no error.
#'@note This is a top level wrapper function - see github.com/MWSchmid/RNAseq_protocol for details.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.wf.weight.multireads.with.Rcount <- function(myAnnotation, RcountMultir, bamDirIn, bamDirOut = bamDirIn) {
  # Import the sample annotation.
  cat("Importing sample annotation...\n")
  sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
  rownames(sampleTab) <- sampleTab$sampleID
  
  # Distribute multireads based on the number of 
  # unique reads with Rcount-multireads
  cat("Setting weights of multireads...")
  for (cSam in rownames(sampleTab)) {
    cat("[WEIGHT] processing", cSam, "...\n")
    inFile <- file.path(bamDirIn, paste0(cSam, "_srt.bam"))
    outFile <- file.path(bamDirOut, paste0(cSam, "_wts.bam"))
    calcWeights <- tolower(sampleTab[cSam, "useMulti"])
    calcWeights <- ifelse(calcWeights == "true", "y", "n")
    alloDist <- sampleTab[cSam, "readLen"]
    command <- paste0(RcountMultir, ' -c "', inFile, '","',
                      outFile, '",', calcWeights, ',', alloDist)
    system(command)
  }
  
  # Return TRUE for success.
  return(TRUE)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Create Rcount project files.
#'@param myAnnotation path to the *.csv with the curated annotation (based on the file created by \code{\link{f.wf.download.from.SRA}})
#'@param RcountPrePro path to the "prepareXMLfileForRcount.py" script.
#'@param RcountProDir path to a folder where the Rcount project files will be stored.
#'@param RcountAnnota path to the genome annotation in the Rcount xml format.
#'@param countTabsRC path to a folder where the individual count tables (one per sample) will be stored.
#'@param bamDirIn path to a folder where the bam files with the weighted reads are stored.
#'@param bamDirOut path to a folder where the bam files with the mapping information will be stored.
#'@return TRUE if there is no error.
#'@note This is a top level wrapper function - see github.com/MWSchmid/RNAseq_protocol for details.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.wf.create.Rcount.projects <- function(myAnnotation, RcountPrePro, RcountProDir, RcountAnnota, countTabsRC, bamDirIn, bamDirOut = bamDirIn) {
  # Import the sample annotation.
  cat("Importing sample annotation...\n")
  sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
  rownames(sampleTab) <- sampleTab$sampleID
  
  # Create Rcount project files
  cat("Creating Rcount project files...\n")
  for (cSam in rownames(sampleTab)) {
    cat("[CREATE PROJECT] processing", cSam, "...\n")
    proFile <- file.path(RcountProDir, paste0(cSam, ".xml"))
    inBam <- file.path(bamDirIn, paste0(cSam, "_wts.bam"))
    outBam <- file.path(bamDirOut, paste0(cSam, "_mpd.bam"))
    countTab <- file.path(countTabsRC, paste0(cSam, ".txt"))
    useMulti <- tolower(sampleTab[cSam, "useMulti"])
    useStrand <- tolower(sampleTab[cSam, "useStrand"])
    minReads <- sampleTab[cSam, "minReads"]
    maxDist <- sampleTab[cSam, "maxDist"]
    minBelowMaxDist <- sampleTab[cSam, "minBelowMaxDist"]
    command <- paste("python", RcountPrePro, 
                     RcountAnnota, inBam, outBam, 
                     countTab, cSam, useMulti, 
                     useStrand, minReads, maxDist, 
                     minBelowMaxDist, proFile)
    system(command)
  }

  # Return TRUE for success.
  return(TRUE)
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#'@title Count reads per gene with Rcount
#'@param myAnnotation path to the *.csv with the curated annotation (based on the file created by \code{\link{f.wf.download.from.SRA}})
#'@param RcountDistri path to the Rcount-distribute executable
#'@param RcountProDir path to a folder where the Rcount project files will be stored.
#'@param countTabsRC path to a folder where the individual count tables (one per sample) will be stored. Note that this folder
#'should not contain any other files (at least no other *.txt files).
#'@param countTabRC path to a *.csv table which will contain the expression values for all samples.
#'@return TRUE if there is no error.
#'@note This is a top level wrapper function - see github.com/MWSchmid/RNAseq_protocol for details.
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.wf.count.reads.with.Rcount <- function(myAnnotation, RcountDistri, RcountProDir, countTabsRC, countTabRC) {
  # Import the sample annotation.
  cat("Importing sample annotation...\n")
  sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
  rownames(sampleTab) <- sampleTab$sampleID
  
  # Or run it via the command line:
  cat("Counting reads per gene...\n")
  proFiles <- sapply(rownames(sampleTab), function(x) 
    file.path(RcountProDir, paste0(x, ".xml")))
  proFiles <- paste(proFiles, collapse = ',')
  command <- paste0(RcountDistri, " -c ", proFiles)
  system(command)
  
  # Read-in all the tables and merge the files which
  # correspond to only one sample (some samples were
  # sequenced on multiple lanes). Windows users must
  # switch to the "native Windows R" to run the code
  # below.
  cat("Merging split samples and combining the individual count tables...\n")
  counts_RC <- f.read.Rcount(countTabsRC)
  counts_RC <- f.summarize.columns(
    counts_RC,
    data.frame(sample = rownames(sampleTab),
               group = sampleTab$mySampleName),
    sum)
  write.csv(counts_RC, countTabRC)
  
  # Return TRUE for success.
  return(TRUE)
}


