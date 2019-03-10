#!/usr/bin/env Rscript

# PODNAME: combineData.R
# ABSTRACT: combineData.R

if(!require(getopt)){install.packages("getopt", repos = "http://cran.us.r-project.org")}
library("getopt")

if(!require(devtools)){install.packages("devtools", repos = "http://cran.us.r-project.org")}
library("devtools")
if(!require(NNLSexperiment)){install_github("sbastkowski/NNLSexperiment")}
library("NNLSexperiment")

options(warn=-1)

options(width=80)

opt = getopt(matrix( c('help', 'h', 0, "logical",
                       'otu_table', 'u', 1, "character",
                       'mapping', 'm', 1, "character",
                       'starter_communities', 's', 1 , "character",
                       'target_communities', 't', 1, "character",
                       'pool', 'p', 1, "character",
                       'outfile', 'o', 1, "character",
                       'level', 'l', 1, "character"), ncol=4, byrow=TRUE ) );


if(!is.null(opt$help) || is.null(opt$otu_table) || is.null(opt$mapping) || is.null(opt$starter_communities) || is.null(opt$target_communities))
{
  cat(paste("Usage: run_example.R [-h] --otu_table otu.biom --mapping mapping.txt [--output output_prefix]\n\n"));
  writeLines(c(strwrap("Reads in OTU table from biom file and calculates contribution of starter microbial communities to target communities."),
               "\n\nRequired Arguments:\n",
               strwrap("--otu_table : OTU table in biom format."),
               strwrap("--mapping : mapping filename in txt format."),
               strwrap("--starter_communities : names of samples that refer to starter communities seperated by a space."),
               strwrap("--target_communities : names of samples that refer to target communities seperated by a space."),
               "\nOptional Arguments:\n",
               strwrap("--output : prefix for output files. (Default: output"),
               strwrap("--pool : samples to be pooled (average is used). A string specifying a pool as follows: poolname1=sample1, sample2, sample3; poolname2=sample4, sample5, sample6, sample7; ... (Default:  none)"),
               strwrap("--level : level to what to collapse the OTU table. (Default: Family)"),"\n"))
  q(status=1);
}

if ( is.null(opt$outfile ) ) {
  writeLines(c(strwrap("No out file name supplied. Output will be saved as output_solutionWeights.txt"),"\n"))
  opt$outfile = "output"
}

if ( is.null(opt$level ) ) {
  writeLines(c(strwrap("No taxanomic level supplied. The OTUs are collapsed to Family level."),"\n"))
  opt$level = "Family"
}

if ( is.null(opt$lpool ) ) {
  writeLines(c(strwrap("No pooling of samples."),"\n"))
}

weightSolutions=run_nnls_analysis(opt$otu_table, opt$mapping, opt$starter_communities, opt$target_communities, opt$level, opt$pool)
write.table(weightSolutions, sep="\t", paste0(opt$outfile, "_solutionWeights.txt"))
ggplot2::ggsave(file = paste0(opt$outfile, "_plots.pdf") , make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")


