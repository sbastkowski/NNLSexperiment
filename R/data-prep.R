read_data <- function(biomfile, mappingfile, pool=NULL, cellcounts=NULL) {

  data = phyloseq::import_biom(biomfile)
  samples = phyloseq::import_qiime_sample_data(mappingfile)
  myTaxTable <- phyloseq::tax_table(data)
  colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
  otu = phyloseq::otu_table(data)
  tax = phyloseq::tax_table(myTaxTable)
  myPhyloSeq <- phyloseq::phyloseq(otu,tax)

  if(!missing(cellcounts)){
    myPhyloSeq=norm_by_cellcount(myPhyloSeq, cellcountvalues = cellcounts)
  }

  if(!is.null(pool)){
    pool <- pool_to_df(myPhyloSeq, pool)
    OTU_pooled=pool_samples(myPhyloSeq, pool)
    myPhyloSeq <- phyloseq::merge_phyloseq(OTU_pooled, tax)
  }

  return(myPhyloSeq)
}

#Normalise by cellcount

norm_by_cellcount <- function(phyloseqObj, cellcountvalues) {

  pooledOTU=rep(0,nrow(phyloseq::otu_table(phyloseqObj)))

  for(i in 1:length(cellcountvalues)){
    for(j in 1:nrow(phyloseqObj@otu_table@.Data)){
      phyloseqObj@otu_table@.Data[j,i]=phyloseqObj@otu_table@.Data[j,i]*cellcountvalues[i]
    }
  }

  return(phyloseqObj)
}

# Input: OTU table and matrix that indicates which samples should be pooled,
# row is one pool and columns (number of samples) indicates which Samples are in this pool
# if entry is 0 its not in pool, if sample is 1 it is in pool.
# single samples are treated as a pool as well, ie one row with zeros except one 1
pool_samples<-function(phyloseqObj, pool) {


  pooledOTU <- rep(0,nrow(phyloseq::otu_table(phyloseqObj)))
  poolMatrix <- data.matrix(pool)
  for(i in 1:nrow(poolMatrix)){
    C <- matrix(rep(0, nrow(phyloseq::otu_table(phyloseqObj))))
    for(j in 1:ncol(poolMatrix)) {
      C <- C + poolMatrix[i,j]*(phyloseq::otu_table(phyloseqObj)[,j])
    }
    pooledOTU <- cbind.data.frame(pooledOTU,C/sum(poolMatrix[i,]))

  }
  row.names(pooledOTU)=row.names(phyloseq::tax_table(phyloseqObj))
  pooledOTU=pooledOTU[,-1]
  colnames(pooledOTU) <- row.names(pool)

  return(phyloseq::otu_table(pooledOTU, taxa_are_rows = TRUE))
}


# Check if this works
pool_to_df <- function(myPhyloseqObj, pool){

  pools <- strsplit(pool, ";")
  samples <- colnames(myPhyloseqObj@otu_table@.Data)
  pooldf <- data.frame(matrix(rep(0, length = length(pools[[1]]) * length(samples)), ncol = length(samples), nrow = length(pools[[1]])))
  colnames(pooldf) <- samples
  poolnames <- vector(mode = "character", length = length(pools[[1]]))
  for (i in 1:length(pools[[1]])) {

    poolSeperated <- strsplit(pools[[1]][i], "=")
    poolnames[i] <- trimws(poolSeperated[[1]][1])
    poolIndividual <- strsplit(poolSeperated[[1]][2], ",")
    for (j in 1:length(poolIndividual[[1]])) {

      sampleInPool <- trimws(poolIndividual[[1]][j])
      pooldf[i, sampleInPool] <- 1

    }
  }
  row.names(pooldf) <- poolnames
  return(pooldf)

}

make_pool<-function(size){

  pool=matrix(rep(0,300), ncol=30, nrow=10)
  for(i in 1:nrow(pool)){
    for(j in 1:size){
      pool[i,((i-1)*size + j)]=1
    }
  }
  return(pool)
}
