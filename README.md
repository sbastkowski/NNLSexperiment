# NNLSexperiment

A set of tools to analyse and plot the output from a BMG FLUOStar OMEGA plate reader.  


## Contents
  * [Introduction](#introduction)
  * [Installation of R package from Github](#installation-of-r-package-from-github)
  * [Script Usage](#script-usage)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)

## Introduction
NNLSexperiment is an R package to calculate and plot the contribution of "starter" microbial communities to a "target" (e.g. mixed) community using a non-negative least square approach. This requires several replications of the experiment. For more details on the method and what kind of experiment it can be used for please see https://doi.org/10.1016/j.cub.2017.09.056. 

## Installation of R package from Github

install.packages("devtools", repos = "http://cran.us.r-project.org")

library("devtools")

install_github("sbastkowski/NNLSexperiment")

library("NNLSexperiment")


For running the Scripts no package installation is necessary, but it requires R (at least version 3.5.2)

## Script Usage

This script reads in OTU table from biom file, a mapping file as well as optional parameters and calculates contribution of starter microbial communities to target communities, which will be outputted in a csv file containing the contributions and a pdf with barplots illustrating the contributions.

### Command-line usage instructions:

./run_example.R [-h] --otu_table otu.biom --mapping mapping.txt [--output output_prefix]



  writeLines(c(strwrap("Reads in OTU table from biom file and calculates contribution of starter microbial communities to target communities."),
               "\n\nRequired Arguments:\n",strwrap("--mapping : mapping filename in txt format."),
               strwrap("--otu_table : OTU table in biom format."),
               "\nOptional Arguments:\n",
               strwrap("--output : prefix for output files."),
               strwrap("--starter_communities : names of samples that refer to starter communities seperated by a space."),
               strwrap("--target_communities : names of samples that refer to target communities seperated by a space."),
               strwrap("--level : level to what to collapse the OTU table. Default: Family"),"\n"))
  q(status=1);

### Required:

--otu_table : OTU table in biom format.


--mapping : mapping filename in txt format.


### Optional:


--output : prefix for output files.

--starter_communities : names of samples that refer to starter communities seperated by a space.

--target_communities : names of samples that refer to target communities seperated by a space.

--level : level to what to collapse the OTU table. Default: Family


Please see examples of input format in Data folder. 

A help menu can be accessed by running the script with -h.

## License
GKinect is free software, licensed under [GPLv3](https://github.com/sbastkowski/GKinect/blob/master/software_license).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sbastkowski/NNLSexperiment/issues) or email sarah.bastkowski@quadram.ac.uk
