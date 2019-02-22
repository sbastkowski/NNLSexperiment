# NNLSexperiment

A set of tools to reads in OTU table from biom file and calculate contribution of starter microbial communities to target communities.  


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

./run_example.R [-h] --otu_table otu.biom --mapping mapping.txt --starter_communities "x y z" --target_communities "a b c d" [--output output_prefix] [--level taxanomix_level] [--pool]


### Required:

--otu_table /-u : OTU table in biom format.

--mapping / -m : mapping filename in txt format.

--starter_communities / -s : names of samples that refer to starter communities seperated by a space. Note, that if you pool samples, you will need to use the poolnames.

--target_communities / -t : names of samples that refer to target communities seperated by a space. Note, that if you pool samples, you will need to use the poolnames.


### Optional:


--output / -o : prefix for output files. (Default: output)

--level / -l : level to what to collapse the OTU table. (Default: Family)
      Options: "Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species"

--pool  / -p : samples to be pooled (average is used). A string specifying a pool as follows: 
      "poolname1=sample1, sample2, sample3; poolname2=sample4, sample5, sample6, sample7; ..."


Please see examples of input format in Data folder. 

A help menu can be accessed by running the script with -h.

## License
NNLSexperiment is free software, licensed under [GPLv3](https://github.com/sbastkowski/NNLSexperiment/blob/master/software_license).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sbastkowski/NNLSexperiment/issues) or email sarah.bastkowski@quadram.ac.uk
