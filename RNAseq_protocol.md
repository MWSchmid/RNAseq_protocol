# RNA-Seq data analysis protocol: combining in-house and publicly available data

This page serves as an up-to-date mirror of the workflow published in [not yet published]().

Please note that I could not include the full text due to licence restrictions.

If you encounter problems, please report them with the "Issues" tool (on top of this page).

## Material

### Installing R and dependencies on Linux (Ubuntu-like)
Visit [cran.r-project.org](https://cran.r-project.org/) and follow the instructions on the web site to get the newest stable version of R (the default version in the package manager is frequently outdated). In brief, open a terminal and type:
```SH
## add the repository to your software sources
# sudo add-apt-repository "deb http://<MIRR>/bin/linux/ubuntu <VERS>/"
# <VERS>: Ubuntu version (code name; e.g. "trusty" for Ubuntu 14.04)
# <MIRR>: A mirror listed on https://cran.r-project.org/mirrors.html 
sudo add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/"

## add the authentication key
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

## update the software package list
sudo apt-get update

## install R
sudo apt-get install r-base r-base-dev r-cran-rjava

## install all other dependencies
sudo apt-get install libxml2-dev libssl-dev libcurl3-openssl-dev qt5-default
```

### Installing R and dependencies on MacOS
Install XCode (AppStore) and XQuartz ([www.xquartz.org](https://www.xquartz.org/)). Then download the R-installer from [cran.r-project.org](https://cran.r-project.org/) and install R according to the installation instructions.

### Installing R and dependencies on Windows
Rsubread and many other NGS tools do not run natively on Windows. An option is to use Cygwin, a Unix-like environment and command-line interface for Windows. Download the most recent installer from ([cygwin.com](https://cygwin.com/)) and start the setup. Go through the dialog until you can select the packages you would like to install. Select the following packages in addition to what is per default selected (for the libraries ending on "-devel", double check if the version without "-devel" is selected as well).

+ Archive:
    * libbz2-devel: BZip file de/compressor
+ Devel:
    * gcc-core: GNU Compiler Collection (C, OpenMP)
    * gcc-g++: GNU Compiler Collection (C++)
    * gcc-fortran: GNU Compiler Collection (Fortran)
    * make: The GNU version of the 'make' utility
+ Libs:
    * libcurl-devel: Multi-protocol file transfer library (development)
    * libiconv-devel: Unicode iconv() implementation
    * libicu-devel: IBM Internationalization Component for Unicode
    * libintl-devel: GNU Internationalization runtime library
    * liblzma-devel: LZMA de/compressor library (development)
    * libpcre-devel: Perl Compatible Regular Expressions library development
    * libtirpc-devel: A port of Sun's Transport-Independent RPC library
    * libxml2-devel: GNOME XML library (development)
    * zlib-devel: Gzip de/compression library (development)
+ Net:
    * openssl-devel: A general purpose cryptography toolkit with TLS implementation (development)
+ Science:
    * R: R Statistical computing language

After completing the Cygwin installation, start Cygwin, type "R" and press enter to start an R console. At least for the alignment of the short reads (i.e. usage of Rsubread), you need to use this console. Note that packages installed in the Cygwin R console are not available to any other native R installation (i.e. R installed with the installer available on [cran.r-project.org](https://cran.r-project.org/)). Aside the R within Cygwin, I recommend to install the "native Windows R" as well.

### Other software and packages (all platforms)
1. Bioconductor packages and other dependencies. Start R and type:
    ```R
    source("http://bioconductor.org/biocLite.R")
    biocLite()
    biocLite(c("SRAdb", "Rsamtools", "Rsubread", 
               "biomaRt", "DESeq2", "edgeR",
               "limma", "XLConnect", "gplots",
               "colorRamps", "ShortRead"))
    # Note for Windows users: 
    # I recommend using the "Cygwin R" only for 
    # the data pre-processing. Only SRAdb, Rsamtools
    # and Rsubread are required for this part.
    # All the other libraries should be installed 
    # on the "native Windows R".
    # Important for Windows users:
    # Current versions of Rsubread may fail to compile
    # on Cygwin and multithreading may be limited as well.
    # It is unclear if and when this will be fixed.
    # You can download version 1.16.1 from 
    # bioconductor.org/packages/3.0/bioc/
    # src/contrib/Rsubread_1.16.1.tar.gz
    # and install it in R with:
    install.packages("/path/to/Rsubread_1.16.1.tar.gz")
    ```
2. Download [RNAseqWrapper](RNAseqWrapper_0.99.0.tar.gz?raw=true) and install it in R with:
    ```R
    install.packages("/path/to/RNAseqWrapper_0.99.0.tar.gz",
                     repos = NULL, type = "source")
    # Note for Windows users: 
    # Install this package only in the "native Windows R".
    ```
3. Download the archive [workFlowData.zip](workFlowData.zip?raw=true) and unpack it.
4. The workflow includes featureCounts to count the number of reads per gene. However, if you would like to specifically address the problem of reads aligning with multiple locations in the genome (multireads) or reads aligning with positions where two or more genes overlap (ambiguous reads), you can use Rcount. Download the archive matching your operating system from [github.com/MWSchmid/Rcount](https://github.com/MWSchmid/Rcount) and unpack it. 

### Reference sequence and annotation

```R
# Specify the local file paths and the links.
gSeqFileGZ <- "/path/to/genome/sequence.fasta(.gz)"
gAnnFileGZ <- "/path/to/genome/annotation.gtf|gff|bed(.gz)"
gSeqLink <- "http:|ftp://link/to/genome/sequence.fasta"
gAnnLink <- "http:|ftp://link/to/genome/annotation.gtf"

# Download the file.
download.file(gSeqLink, gSeqFileGZ, "auto")
download.file(gAnnLink, gAnnFileGZ, "auto")

# If the files were compressed, unpack them either manually
# or with a system command (Linux, MacOS and Cygwin). Note
# that this will fail on a native Windows R installation.
system(paste0("gunzip ", gSeqFileGZ))
system(paste0("gunzip ", gAnnFileGZ))
```


## Methods

### Download publicly available data

```R
library("SRAdb")

# Specify several file paths and directories in which
# you would like to work in.
readsDir <- "/path/to/the/raw/read/files"
metaDB <- "/path/to/metadata/database/SRAmetadb.sqlite"
mySamplesFile <- "/path/to/the/file/with/all/mySamples.txt"
autoAnnotation <- "/path/to/autoGenAnnotation.csv"

# Check for a local copy of the database and
# download it if it does not exist.
if (!file.exists(metaDB)) { getSRAdbFile(dirname(metaDB)) }

# Windows users: if the download fails while using the
# Cygwin-R, try:
if (!file.exists(metaDB)) { 
  getSRAdbFile(dirname(metaDB), method = "libcurl")
}

# Open a connection to the database
sra_con <- dbConnect(SQLite(), metaDB)

# Here we assume that you have a precompiled list of 
# samples/experiments/studies which are of interest to you.
# Read in the list of experiments you are interested in.
samples <- scan(mySamplesFile, what = "character")

# You may have a collection of experiments and individual runs. 
# However, only individual runs can be downloaded. Note that 
# the conversion of a certain SRA ID type to another requires
# the different query types not to be mixed.
sampleTypes <- unique(substr(samples, 1, 3))
splitSamples <- lapply(sampleTypes, function(x)
  grep(paste0("^", x), samples, value = TRUE))
conversions <- lapply(splitSamples, function(x)
  sraConvert(x, sra_con = sra_con ))
allSRAids <- do.call("rbind", conversions)

# Retrieve the remaining metadata to determine
# the sequencing platform, the read length, etc.
SRAdata <- lapply(allSRAids$run, function(x)
  getSRA(search_terms = x, out_types = "sra", sra_con))
SRAdata <- do.call("rbind", SRAdata)
rownames(SRAdata) <- SRAdata$run # for data access
SRAdata$approxRL <- SRAdata$bases/SRAdata$spots
summary(SRAdata$spots)    # the total number of reads
summary(SRAdata$approxRL) # the average read lengths
table(SRAdata$platform)   # the sequencing platform
SRAdata$sample_name[1:10] # the names of the first 10 samples
# The sample names are not always too informative
# and I recommend renaming them afterwards.

# Search for strand-specific samples with regular
# expressions. Note that this gives no guarantee
# for strand-specificity - in any case, this
# should be verified manually.
posExp <- "strand[[:space:]{0,1}|-]specific"
negExp <- paste0("no[[:alpha:]{0,1}][[:space:]{0,1}|-]", posExp)
strSpec <- apply(SRAdata, 1, function(x)
  length(grep(posExp, x)) > 0)
notStSp <- apply(SRAdata, 1, function(x)
  length(grep(negExp, x)) > 0)
strSpec[notStSp] <- FALSE
sum(strSpec)  # the number of strand-specific libraries
SRAdata$strSpec <- strSpec

# Check if there are FASTQ files available for the samples.
# For the samples without FASTQ file, you will download the
# SRA archives and convert them into FASTQ (this is a lot 
# slower than downloading FASTQ files).
fastqInfo <- getFASTQinfo(allSRAids$run, sra_con, 'ftp')

# Check if the paths were retrieved correctly. If there is
# no column called "fastq_ID" (e.g., FASTQ_FILES with all
# entries being "1"), the database may be broken.
head(fastqInfo)

samplesWithFastq <- unique(with(fastqInfo, run[!is.na(fastq_ID)]))
noFastqSamples <- unique(with(fastqInfo, run[is.na(fastq_ID)]))
if (length(noFastqSamples) > 0) {
  sraInfo <- getSRAinfo(noFastqSamples, sra_con)
  noFastqButSRA <- unique(with(sraInfo, run[!is.na(ftp)]))
  noFiles <- unique(with(sraInfo, run[is.na(ftp)]))
} else {
  noFiles <- c()
}
length(noFiles) # number of samples without files
# if there are samples without FASTQ or SRA file, 
# go to http://www.ncbi.nlm.nih.gov/, under "Databases"
# select "SRA" and check the status of the IDs for which
# there are no files available.

# Download the FASTQ files.
# Note that an error like "no such table: fastq" is likely
# due to a corrupt SRAdb file. Delete SRAmetadb.sqlite
# and download it again. If the error persists, contact
# the SRAdb developers.
getFASTQfile(samplesWithFastq, sra_con, readsDir, 'ftp')

# If necessary, download the SRA archives of the remaining
# samples and convert them to FASTQ. Note that the conversion
# requires the "SRA toolkit" to be installed:
# www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
getSRAfile(noFastqButSRA, sra_con, readsDir, 'sra' )
commands <- sapply(noFastqButSRA, function(x)
  paste0("fastq-dump --gzip --split-files -O ",
         readsDir, " ", file.path(readsDir, x), ".sra"))
sapply(commands, system)

# Get a vector of all FASTQ file names and
# remove the reverse reads of the PE samples.
allFastq <- list.files(readsDir, ".fastq.gz")
revFastq <- list.files(readsDir, "_2.fastq.gz")
onlyForwardFastq <- setdiff(allFastq, revFastq)
sampleIDs <- gsub("_1|.fastq.gz", "", onlyForwardFastq)

# Samples with paired-end sequencing have twice the
# coverage - therefore adjust the approximated
# read length of these samples.
PES <- gsub("_2|.fastq.gz", "", revFastq)
SRAdata[PES, "approxRL"] <- SRAdata[PES, "approxRL"]/2

# Create a table with information for the alignment
# and the preprocessing. Store the sampleID, read length,
# platform type, FASTQ file names, sample/library names, 
# strand-specificity and additional attributes in a csv file. 
# Even though the sample/library names and the attributes
# are not complete for each sample, they may help renaming
# all samples.
sortedSRAdata <- SRAdata[sampleIDs,] # this ensures that the
tabForProcessing <- data.frame(      # sample has a FASTQ file
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
  processedAlignmentsFC = 0,
  assignedAlignmentsFC = 0,
  TISSUE     = "toBeProvided",   # these columns depend on the
  CELLTYPE   = "toBeProvided",   # analysis. In another case
  STAGE      = "toBeProvided",   # it may be GENOTYPE, TREATMENT 
  BLOCK      = "toBeProvided",   # or PERTURBATION.
  mySampleName = "toBeProvided",
  SRAsamSRS  = sortedSRAdata$sample,
  SRAexpSRX  = sortedSRAdata$experiment,
  SRAsamName = sortedSRAdata$sample_name,
  SRAlibName = sortedSRAdata$library_name,
  SRAattrib  = sortedSRAdata$sample_attribute,
  stringsAsFactors = FALSE
)

# Optional: Verify that the PHRED offset is 33 and not 64.
# This is not strictly necessary as the files from SRA are
# generally having an offset of 33 (they otherwise convert
# the encoding). However, there are still some cases with
# an offset of 64.
library("RNAseqWrapper")
rownames(tabForProcessing) <- tabForProcessing$sampleID
for (cSam in rownames(tabForProcessing)) {
  fastqFile <- file.path(readsDir,
                         tabForProcessing[cSam, "fastqFile"])
  phredOffset  <- f.check.phred.offset(fastqFile)
  tabForProcessing[cSam, "phredOffset"] <- phredOffset
}

# Save the table.
write.csv(tabForProcessing, autoAnnotation, row.names = FALSE)

# Optional: Generate quality control reports for all the
# fastq files using the "ShortRead" package. By default,
# the reports are based on 1 million randomly chosen reads.
library("ShortRead")
set.seed(123)
reportDir <- "/path/to/the/reports"
qaSummary <- qa(readsDir, "fastq.gz", type = "fastq")
report_html(qaSummary, dest = reportDir)
```

### Manual curation of the sample annotation

see [not yet published]() for detailed instructions

### Alignment of short reads with the reference genome

```R
library("Rsubread")

# Specify several file paths and directories in which
# you would like to work in. 
gSeqFile <- "/path/to/genome/sequence.fasta"
nucIndex <- "/base/name/of/the/genome/nucIndex" 
colIndex <- "/base/name/of/the/genome/colIndex" 
myAnnotation <- "/path/to/the/curatedAnnotation.csv"
readsDir <- "/path/to/the/raw/read/files"
bamDir <- "/path/to/aligned/reads/files"

# Build the indices Rsubread. Check the documentation
# for the details on the arguments (?buildIndex).
# The index for the Illumina (nucleotide-encoded) reads:
buildindex(nucIndex, gSeqFile,
           gappedIndex = TRUE,  # FALSE on high-RAM machines
           memory = 8000,       # max 80 % of all your RAM
           colorspace = FALSE)

# The index for the SOLiD (color-encoded) reads:
buildindex(colIndex, gSeqFile,
           gappedIndex = TRUE,
           memory = 8000,
           colorspace = TRUE)

# Import the sample annotation.
sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
rownames(sampleTab) <- sampleTab$sampleID

# Align the reads - allow up to 10 alignments per read.
# Note that we redirect the messages from the aligner
# into a character vector to extract the number of
# mapped reads. However, the drawback is that we 
# will not see the progress messages from the aligner
# in real-time.
for (cSam in rownames(sampleTab)) {
  cat("processing", cSam, "...\n")
  isSOLiD <- (sampleTab[cSam, "platform"] == "ABI_SOLID")
  pOff <- sampleTab[cSam, "phredOffset"]
  cIdx <- ifelse(isSOLiD, colIndex, nucIndex)
  inFile <- file.path(readsDir, sampleTab[cSam, "fastqFile"])
  outFile <- file.path(bamDir, paste0(cSam, ".bam"))
  if (file.exists(outFile)&(sampleTab[cSam,"processedReads"]>0)){next}
  stdout <- capture.output( # delete this line to disable redirecting
    subjunc(cIdx, inFile, output_file = outFile,
            nthreads = 4,        # adjust the number of CPUs 
            unique = FALSE,      # set to TRUE to avoid multireads
            nBestLocations = 10, # max number of multireads
            phredOffset = pOff   # PHRED offset is either 33 or 64
            )
  )                         # delete this line to disable redirecting
  # extract the number of processed and mapped reads
  readPattern <- "(Processed|Mapped).*reads"
  readMatches <- regexpr(readPattern, stdout)
  readCountsStr <- regmatches(stdout, readMatches)
  numberPattern <- "[,0123456789]+"
  numberMatches <- regexpr(numberPattern, readCountsStr)
  readCounts <- regmatches(readCountsStr, numberMatches)
  sampleTab[cSam, "processedReads"] <- gsub(",", "", readCounts[1])
  sampleTab[cSam, "alignedReads"] <- gsub(",", "", readCounts[2])
  write.csv(sampleTab, myAnnotation, row.names = FALSE)
}

# Save the table again to store the number of processed
# and mapped reads as well. "processedReads" should be
# equal to "numReads", and the number of "alignedReads"
# may give you an indication whether the samples are
# fine (note that alignment rates around 40 % are not
# necessarily a problem - it depends on the protocol).
write.csv(sampleTab, myAnnotation, row.names = FALSE)

# Optional: Draw a histogram with the alignment percentages.
sampleTab$alPerc <- 100*with(sampleTab,
  as.numeric(alignedReads)/as.numeric(processedReads))
png("/path/to/the/figure.png")
hist(sampleTab$alPerc, xlab = "alignment rate (in %)",
     ylab = "no. samples", main = "", las = 1)
dev.off()

# Optional: Display the samples with a low alignment rate.
colsOfInt <- c("readLen", "platform", "mySampleName", "alPerc")
subset(sampleTab, alPerc < 20)[,colsOfInt] # e.g., less than 20%
```

### Counting unique reads per gene with featureCount

```R
library("Rsubread")

# Specify several file paths and directories in which
# you would like to work in. 
myAnnotation <- "/path/to/the/curatedAnnotation.csv"
gAnnFile <- "/path/to/genome/annotation.gtf|gff|saf"
bamDir <- "/path/to/aligned/reads/files"
countTabsFC <- "/path/to/featureCounts/tables"
countTabFC <- "/path/to/the/final/featureCounts_table.csv"

# NOTE: the directory with the count tables from featureCounts
# (variable countTabsFC) should not contain anything else
# aside the count tables (at least no other *.csv files).

# Import the sample annotation.
sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
rownames(sampleTab) <- sampleTab$sampleID

# Count the number of (unique!) reads with featureCounts().
# I strongly recommend to read the documentation of this
# function (type ?featureCounts).
# Note that we redirect the messages from featureCounts
# into a character vector to extract the number of
# assigned reads. However, the drawback is that we 
# will not see the progress messages in real-time.
sT <- c("false" = 0, "sense" = 1, "antisense" = 2)
for (cSam in rownames(sampleTab)) {
  cat("processing", cSam, "...\n")
  bamFile <- file.path(bamDir, paste0(cSam, ".bam"))
  outFile <- file.path(countTabsFC, paste0(cSam, ".txt"))
  SSmode <- tolower(sampleTab[cSam, "useStrand"])
  SSmode <- sT[SSmode]
  minBPinExon <- ceiling(.7*sampleTab[cSam, "readLen"])
  stdout <- capture.output( # delete this line to disable redirecting
    countRes <- featureCounts(bamFile,
        annot.ext = gAnnFile,       # reference annotation
        isGTFAnnotationFile = TRUE, # TRUE if ref. anno. is GTF or GFF
        GTF.featureType = "exon",   # map alignments to exons
        GTF.attrType = "gene_id",   # use geneIDs to group exons
        useMetaFeatures = TRUE,     # yes, summarize counts per gene
        allowMultiOverlap = TRUE,   # count ambiguous reads
        nthreads = 4,               # adjust the number of CPUs
        strandSpecific = SSmode,    # set strand-specificity
        # Windows (Rsubread 1.16.1): minReadOverlap instead of minOverlap
        minOverlap = minBPinExon,   # min num of BPs within exons
        # Windows (Rsubread 1.16.1): remove largestOverlap argument
        largestOverlap = TRUE)      # see documentation
  )                         # delete this line to disable redirecting
  # Save the count table for the current sample.
  colnames(countRes$counts) <- cSam
  write.table(countRes$counts, outFile, sep = '\t', quote = FALSE)
  # extract the number of processed and assigned reads
  readPattern <- "(Total|Successfully).*reads"
  readCountsStr <- grep(readPattern, stdout, value = TRUE)
  numberPattern <- "[0123456789]+"
  numberMatches <- regexpr(numberPattern, readCountsStr)
  readCounts <- regmatches(readCountsStr, numberMatches)
  sampleTab[cSam, "processedAlignmentsFC"] <- readCounts[1]
  sampleTab[cSam, "assignedAlignmentsFC"] <- readCounts[2]
}

# Save the sample annotation table again to store the 
# number of processed and successfully assigned reads.
# The number of "processedAlignmentsFC" should equal to
# "alignedReads", and the number of "assignedAlignmentsFC"
# may give you an indication whether a sample has
# too many reads in intergenic regions.
write.csv(sampleTab, myAnnotation, row.names = FALSE)

# Read-in all the tables and merge the files which
# correspond to only one sample (some samples were
# sequenced on multiple lanes). Windows users must
# switch to the "native Windows R" to run the code
# below.
library("RNAseqWrapper")
counts_FC <- f.read.featureCounts(countTabsFC)
counts_FC <- f.summarize.columns(
  counts_FC,
  data.frame(sample = rownames(sampleTab),
             group = sampleTab$mySampleName),
  sum)
write.csv(counts_FC, countTabFC)

# Optional: Draw a histogram with the assignment percentages.
sampleTab$asPerc <- 100*with(sampleTab,
  as.numeric(assignedAlignmentsFC)/as.numeric(processedAlignmentsFC))
png("/path/to/the/figure.png")
hist(sampleTab$asPerc, xlab = "assignment rate (in %)",
     ylab = "no. samples", main = "", las = 1)
dev.off()

# Optional: Display the samples with a low assignmend rate.
colsOfInt <- c("readLen", "platform", "mySampleName", "asPerc")
subset(sampleTab, asPerc < 50)[,colsOfInt] # e.g., less than 50 %
```

### Counting weighted alignments with Rcount

```R
library("Rsamtools")
library("RNAseqWrapper")

# Specify several file paths and directories in which
# you would like to work in. 
myAnnotation <- "/path/to/the/curatedAnnotation.csv"
bamDir <- "/path/to/aligned/reads/files"
countTabsRC <- "/path/to/Rcount/tables"
countTabRC <- "/path/to/the/final/Rcount_table.csv"
RcountFormat <- "/path/to/Rcount-format(.exe)"
RcountMultir <- "/path/to/Rcount-multireads(.exe)"
RcountDistri <- "/path/to/Rcount-distribute(.exe)"
RcountProDir <- "/path/to/directory/with/Rcount/project/files"

# NOTE: the directory with the count tables from Rcount
# (variable countTabsRC) should not contain anything else
# aside the count tables (at least no other *.txt files).

# Create a database for Rcount with Rcount-format.
# See the Rcount user guide for details.
system(RcountFormat, wait = FALSE)

# Set the path for the annotation you just created.
RcountAnnota <- "/path/to/Rcount/annotation.xml" 

# Import the sample annotation.
sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
rownames(sampleTab) <- sampleTab$sampleID

# Sort the reads in the bam files according to 
# chromosome and position, and create an index.
for (cSam in rownames(sampleTab)) {
  cat("processing", cSam, "...\n")
  inFile <- file.path(bamDir, paste0(cSam, ".bam"))
  outFile <- file.path(bamDir, paste0(cSam, "_srt"))
  sortBam(inFile, outFile, maxMemory = 8192)  # adjust RAM
  indexBam(paste0(outFile, ".bam"))
}
# Note that warning messages containing 
# [bam_sort_core] merging from X files
# can be ignored.

# Distribute multireads based on the number of 
# unique reads with Rcount-multireads
for (cSam in rownames(sampleTab)) {
  cat("processing", cSam, "...\n")
  inFile <- file.path(bamDir, paste0(cSam, "_srt.bam"))
  outFile <- file.path(bamDir, paste0(cSam, "_wts.bam"))
  calcWeights <- tolower(sampleTab[cSam, "useMulti"])
  calcWeights <- ifelse(calcWeights == "true", "y", "n")
  alloDist <- sampleTab[cSam, "readLen"]
  command <- paste0(RcountMultir, ' -c "', inFile, '","',
                    outFile, '",', calcWeights, ',', alloDist)
  system(command)
}

# Create Rcount project files
for (cSam in rownames(sampleTab)) {
  cat("processing", cSam, "...\n")
  f.write.Rcount.project.file(
    sampleName = cSam,
    annotationXML = RcountAnnota,
    readsInfile = file.path(bamDir, paste0(cSam, "_wts.bam")),
    readsOutfile = file.path(bamDir, paste0(cSam, "_mpd.bam")),
    countTableFile = file.path(countTabsRC, paste0(cSam, ".txt")),
    OUTprojectXML = file.path(RcountProDir, paste0(cSam, ".xml")),
    useMulti = tolower(sampleTab[cSam, "useMulti"]),
    useStrand = tolower(sampleTab[cSam, "useStrand"]),
    minReads = sampleTab[cSam, "minReads"],
    maxDist = sampleTab[cSam, "maxDist"],
    minBelowMaxDist = sampleTab[cSam, "minBelowMaxDist"])
}

# Distribute and count alignments. Note that this
# may be done with the GUI instead of the command
# line prompt. The advantage is that the read 
# statistics are added to the project files and
# browsable within the GUI.

# GUI - start it with a double click or the system
# command below, open all the projects and press
# "run all".
system(RcountDistri, wait = FALSE)

# Or run it via the command line:
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
counts_RC <- f.read.Rcount(countTabsRC)
counts_RC <- f.summarize.columns(
  counts_RC,
  data.frame(sample = rownames(sampleTab),
             group = sampleTab$mySampleName),
  sum)
write.csv(counts_RC, countTabRC)
```

### Identification of differentially expressed genes

see [not yet published]() for detailed instructions

```R
# Writing the MS Excel workbooks needs quite some RAM.
# On machines with less than 8GB of RAM, you can either
# set the ensembl variable below to NA (ensembl <- NA)
# or write csv files instead of workbooks.
options(java.parameters = "-Xmx8g")
library("biomaRt")
library("edgeR")
library("DESeq2")
library("limma")
library("XLConnect")
library("RNAseqWrapper")

# Specify several file paths and directories in which
# you would like to work in.
myAnnotation <- "/path/to/the/curatedAnnotation.csv"
countTabFC <- "/path/to/featureCounts_table.csv"
countTabRC <- "/path/to/Rcount_table.csv"
rDir <- "/directory/where/results/are/stored"

# Import the sample annotation.
sampleTab <- read.csv(myAnnotation, stringsAsFactors = FALSE)
rownames(sampleTab) <- sampleTab$sampleID

# Currently there are still some samples with
# several runs in the annotation. We need to remove
# them latest now.
colsToKeep <- c("TISSUE", "CELLTYPE", "STAGE", "BLOCK",
                "mySampleName")
sampleTab <- unique(sampleTab[,colsToKeep])
rownames(sampleTab) <- sampleTab$mySampleName

# Choose a biomart database (only if you were using
# the corresponding reference genome and annotation).
# To display available marts and datasets for animals
# and plants (there may also be other hosts):
listMarts(host = "www.ensembl.org")     # animals
listMarts(host = "plants.ensembl.org")  # plants

# Connect to a database and check whether there
# is a dataset for Arabidopsis available:
ensembl <- useMart("plants_mart", host = "plants.ensembl.org")
ensemblDatasets <- listDatasets(ensembl)
ensemblDatasets[grep("Arabidopsis", ensemblDatasets$description),]

# Finally connect to the database for A. thaliana:
ensembl <- useDataset("athaliana_eg_gene", mart = ensembl)

# Set the class of the biomart DB to "ensembl".
# Note that this is a simple error-workaround
# which is only necessary for some hosts.
biomaRt:::martBM(ensembl) <- "ensembl"

# Create a single factor which summarizes all the
# information given in the annotation table. Note
# that TISSUE, CELLTYPE, STAGE and BLOCK can be
# deleted/replaced/extended as required. See the
# paragraph on manual curation of the sample 
# annotation above.
sampleTab$SFAC__ <- with(sampleTab, 
                         paste(TISSUE, CELLTYPE,
                               STAGE, BLOCK,
                               sep = '_'))

# Convert the new variable to a factor.
sampleTab$SFAC__ <- factor(sampleTab$SFAC__)

# Specify a generic model and the design matrix.
formulaString <- "~0+SFAC__"
design <- model.matrix(formula(formulaString),
                       data = sampleTab,
                       contrasts.arg = NULL)


# Specify the contrasts (comparisons) you are 
# interested in. Instead of writing the entire
# formula by hand, one can define the formulas
# using regular expressions and the function
# f.formulate.simple.contrast. See ?regex for
# more details on regular expressions. However,
# you can also specify the names directly without
# using regular expressions (in this case you need
# to use fixed=TRUE in f.formulate.simple.contrast).
# See ?f.formulate.simple.contrast for more details.
# The function will specify the formula that it equals
# to the samples in the first argument minus the samples
# in the second argument (each scaled to the number of
# samples). In the example contrast "femVSmal" below,
# a positive log fold-change therefore means that the
# gene is expressed at a higher level in the female
# tissues compared to the male tissues.
# Define groups of samples based on their names
# and/or name patterns. The available names are
# given by the design matrix:
colnames(design)

# Based on the names, you can define an expression
# which specifies all the "green" tissues of a plant:
RE_green <- "__flo|__inflo|__leaf|__SAM|__sdl|__sil"

# Define some more groups:
RE_root <- "__root"         # all root tissues
RE_emb  <- "__emb"          # all embryonic tissues
RE_femGam <- "__femGam"     # all cells of the female gametophyte
RE_malGam <- "__malGam"     # all samples from the male gametophyte

# And a group which includes any gametophytic tissue/cell type.
# Note that there are several ways to achieve the same:
RE_anyGam <- c("__malGam", "__femGam") # two exact words
RE_anyGam <- "__[[:alpha:]]{3}Gam"     # __<three characters>Gam
RE_anyGam <- "__femGam|__malGam"       # __femGam OR __malGam
myCont <- makeContrasts(
  f.formulate.simple.contrast(         # rootVSgreen
    RE_root, RE_green, design),
  f.formulate.simple.contrast(         # femVSmal
    RE_femGam, RE_malGam, design),
  f.formulate.simple.contrast(         # gamVSspo
    RE_anyGam, RE_anyGam, design,
    invertMinus = TRUE),
  f.formulate.simple.contrast(         # gamVSemb
    RE_anyGam, RE_emb, design),
  levels = design)

# Give short names to each contrast. 
colnames(myCont) <- c("rootVSgreen", "femVSmal",
                      "gamVSspo", "gamVSemb")

# Load the featureCounts and Rcount data sets.
myDataList <- list(
  feaCou = read.csv(countTabFC, row.names = 1),
  Rcount = read.csv(countTabRC, row.names = 1)
)

# Remove genes with less than five reads in all 
# samples.
myDataList <- lapply(myDataList, function(x)
  f.strip.data(x, minVal = 5, minTimes = 1))

# Draw a sample correlation matrix. Given the size
# the data set, skip drawing all scatter plots at once.
for (dT in names(myDataList)) {
  logData <- log2(myDataList[[dT]]+1)
  f.do.some.overview(logData, rDir, dT,
                     skipScatters = TRUE)
}

# Normalize the data and save all tables.
myNormData <- lapply(myDataList, function(x) 
  f.all.normalizations(x, sampleTab, formulaString, design))
for (dT in names(myNormData)) {
  for (nM in names(myNormData[[dT]])) {
    write.csv(myNormData[[dT]][[nM]],
      file.path(rDir, paste0(dT, '_', nM, "_normalized.csv")))
  }
}

# Calculate the average expression level across a 
# specific condition (average across replicates).
byTab <- data.frame(sample = rownames(sampleTab),
                    group = sampleTab$SFAC__,
                    stringsAsFactors = FALSE)
meanTabs <- list()
for (dT in names(myDataList)) {
  meanTabs[[dT]] <- lapply(myNormData[[dT]], function(x)
    f.summarize.columns(x, byTab, mean))
}

# Test for differential expression.
results <- lapply(myDataList, function(x)
  f.multiple.multi.level.comparisons(x, sampleTab,
                                     formulaString,
                                     myCont, design))

# Write the results into MS Excel workbooks.
for (dT in names(results)) {
  for (dM in names(results[[dT]])) {
    f.write.DEGtabs.to.workbook(results[[dT]][[dM]], rDir,
                                paste0(dT, '_', dM), ensembl)
  }
}
# Note that warning messages containing 
# is.na() applied to non-(list or vector)...
# can be ignored.
```
