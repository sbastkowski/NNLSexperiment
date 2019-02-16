read_data <- function(biomfile, mappingfile, pool=NULL, cellcounts=NULL) {

  data = import_biom(biomfile)
  samples = phyloseq::import_qiime_sample_data(mappingfile)
  myTaxTable <- phyloseq::tax_table(data)
  colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
  otu = phyloseq::otu_table(data)
  tax = phyloseq::tax_table(myTaxTable)
  myPhyloSeq <- phyloseq::phyloseq(otu,tax)

  if(!missing(cellcounts)){
    myPhyloSeq=norm_by_cellcount(myPhyloSeq, cellcountvalues = cellcounts)
  }

  if(!missing(pool)){
    OTU_pooled=pool_samples(myPhyloSeq, pool)
    myPhyloSeq <- phyloseq::merge_phyloseq(OTU_pooled, tax)
  }

  return(myPhyloSeq)
}

#Normalise by cellcount

norm_by_cellcount <- function(phyloseqObj, cellcountvalues) {

  pooledOTU=rep(0,nrow(otu_table(phyloseqObj)))

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

  pooledOTU=rep(0,nrow(otu_table(phyloseqObj)))

  for(i in 1:nrow(pool)){
    C=matrix(rep(0, nrow(otu_table(phyloseqObj))))
    for(j in 1:ncol(pool)) {
      C=C+ pool[i,j]*otu_table(phyloseqObj)[,j]
    }
    pooledOTU=cbind.data.frame(pooledOTU,C/sum(pool[i,]))

  }
  row.names(pooledOTU)=row.names(tax_table(phyloseqObj))
  pooledOTU=pooledOTU[,-1]

  return(otu_table(pooledOTU, taxa_are_rows = TRUE))
}

# Little helper to create pooling matrix
# Need to specify how to present the pooling or add function to pull this out of mapping file
make_pool<-function(size){

  pool=matrix(rep(0,300), ncol=30, nrow=10)
  for(i in 1:nrow(pool)){
    for(j in 1:size){
      pool[i,((i-1)*size + j)]=1
    }
  }
  return(pool)
}

