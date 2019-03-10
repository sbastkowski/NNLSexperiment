FROM debian:testing

# Install required debian packages
RUN apt-get update -qq && apt-get install -y sudo git r-base libssl-dev libcurl4-openssl-dev curl r-cran-httr

# Install required R packages
RUN Rscript -e "install.packages(\"getopt\", repos = \"http://cran.us.r-project.org\")"
RUN Rscript -e "install.packages(\"devtools\", repos = \"http://cran.us.r-project.org\")"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" -e "biocLite(c('phyloseq'))"
RUN Rscript -e "library(\"devtools\")" -e "install_github(\"sbastkowski/NNLSexperiment\")"

# Get repo (add to root)
RUN git clone https://github.com/sbastkowski/NNLSexperiment.git

# Set /work as working directory... this is where we should mount external file systems
WORKDIR /work

# Default app to run
ENTRYPOINT ["Rscript", "/NNLSexperiment/Scripts/run_example.R"]

