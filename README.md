# RNAseq_protocol
From short reads to differential expression

I'm currently working on a workflow for RNA-Seq data analysis. It covers all steps from retrieval of sequences from SRA down to the differential expression analysis. It will be purely based on R and it looks like it will run on MacOSX, Linux and Windows. 

If (legally) possible, I will deposit it here somewhen in April 2016 (well - whenever I get the answer from the official side...).

# RNAseqWrapper

I wrote a collection of wrapper functions which will be used in the workflow as well. Independent of the workflow, I added several generic examples for RNAseq data analysis.

If you are on an Ubuntu-like system and if you would like to have the newest version of R:
```SH
### install the newest version of R
## add the repository to your software sources
# sudo add-apt-repository "deb http://<MIRR>/bin/linux/ubuntu <VERS>/"
# <VERS>: Ubuntu version (code name; e.g. "trusty" for Ubuntu 14.04)
# <MIRR>: A mirror listed on https://cran.r-project.org/mirrors.html 
sudo add-apt-repository "deb https://stat.ethz.ch/CRAN/bin/linux/ubuntu trusty/"

## add the authentication key
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

## update the software package list
sudo apt-get update

## install R
sudo apt-get install r-base r-base-dev r-cran-rjava
```

RNAseqWrapper imports several other packages. Install the following dependencies within R:
```R
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("biomaRt", "DESeq2", "edgeR", "limma", "XLConnect", "gplots", "colorRamps", "SRAdb"))
```

Download [RNAseqWrapper](RNAseqWrapper_0.99.0.tar.gz?raw=true) and install it:
```R
install.packages("/path/to/RNAseqWrapper_0.99.0.tar.gz", repos=NULL)
# you may need to specify the type:
install.packages("/path/to/RNAseqWrapper_0.99.0.tar.gz", repos=NULL, type="source")
```

# RNAseq data analysis examples

Please note that the examples are based on real data, which however is not yet publicly available. Thus the examples are tested, but I can't supply the test data. Please contact me if you encounter a problem.

## Pairwise comparisons

In the most simple case you may have only two groups of samples (e.g. wild-type vs mutant, treated vs mock control or tissue A vs tissue B):

[two group comparison](examples/TGNB.md)

However, you might have processed the samples on different days (hopefully on each day from both conditions...). In this cases it's sometimes beneficial to include this information as a batch effect:

[two group comparison with batch effect](examples/TGWB.md)

or a "real-world example": [paired data from the cancer genome atlas](examples/TCGA_pairedSamples.md)

Eventually you have more than two groups for which you would like to do all pairwise comparisons:

[multiple two group comparisons](examples/MTGNB.md) [REQUIRES EACH GROUP TO BE REPLICATED]

However, you might have processed these samples as well in batches:

[multiple two group comparisons with batch effect](examples/MTGWB.md) [REQUIRES EACH GROUP TO BE REPLICATED]

## Two-by-two crossed factorial design

You may have two factors in your experiment, e.g. genetic background (wildType vs mutant) and drug treatment (mock vs drug). For cases where one (or both) of the factors have more than two levels (e.g. many different mutants), I recommend to use the approach below. However, in case of two two-level-factors we may also used a crossed factorial design (with/without batch):

[two-by-two factorial design](examples/TBTNB.md)

[two-by-two factorial design with batch effect](examples/TBTWB.md)

## Single/compound factor with multiple levels

Finally, your experiment may have several factors with two or more levels. Multifactorial models tend to be quickly quite complex and non-intuitive. Alternatively one can combine all factors and their levels into one single factor with multiple levels. The comparisons of interest can then be done using linear contrasts:

[single/compound factor with several levels](examples/MLNB.md) [SOME GROUPS MAY BE UNREPLICATED]

(there is no simple one-fits-all solution for the single/compound factor with several levels plus batch effect - well, one can write it out, but it's really lengthy then. And finally - using a model with batch effect is not the best choice in all cases. It's mainly beneficial if the batch effect is strong. Otherwise it is safe (and sometimes better) not to use the batch factor)

## (anticipated) FAQs

1. I don't have any replicates, what can I do now?
    * edgeR and limma will not work, use the individual DESeq2 functions (the single/compound factor should work in any case as long as there are some groups with replicates). But - just because DESeq2 may run, it does not mean that you should do it.
    * think about a possible pseudoreplication. It's not very clean, but you sometimes you have the option to treat two different samples as replicates.
    * if it was your doing, don't do it again ;)
2. My count table does not have counts per genes but counts per transcripts - is there anything I need to consider?
    * Yes - in the function `f.write.DEGtabs.to.workbook()` you should set `addGeneCol=TRUE`.

## About c.DEGtab

The results of the DE analysis are stored as c.DEGtab object. This object is implemented as a reference class and has the following functions associated with it:

```R
DETAB # a c.DEGtab object

# add more data to the tables (this would be saved in the Excel workbooks as well):
DETAB$add_data(x) # x must be a data.frame - merging is done based on rownames()

# convert the DETAB to a regular data.frame:
deData <- DETAB$get_table()

# similar, but round, sort and apply some formatting:
deData <- DETAB$get_print_table()

# extract a list with three gene sets (upregulated in any, conditionA or conditionB):
sigGenes <- DETAB$get_significant_entries(rawPcut = 1,    # raw P-value cutoff
                                          adjPcut = 0.05, # FDR cutoff
                                          LFCcut = 0)     # LFC threshold

# get some simple information on the number of genes upregulated in one or the other
# condition. The function returns a 3D array with the third dimension being conditionA,
# conditionB and "simple". 
stats <- DETAB$get_stats()

# or print them directly:
DETAB$print_stats()

# (markdown is possible as well)
DETAB$print_stats_markdown()
```





