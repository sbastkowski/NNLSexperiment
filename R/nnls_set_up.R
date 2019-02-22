set_up_A <- function(myPhyloSeq, level, starter) {

  starters <- names_to_index(myPhyloSeq, starter)
  bacteria_level <- phyloseq::tax_glom(myPhyloSeq, taxrank=level)
  bacteria_level_df <- as.data.frame(phyloseq::get_taxa(phyloseq::otu_table(bacteria_level)))
  bacteria_level_starterComm_df  <-  bacteria_level_df[,starters]
  bacteria_family_matrix_A <- data.frame(as.matrix(bacteria_level_starterComm_df))
  colnames(bacteria_family_matrix_A)=colnames(myPhyloSeq@otu_table@.Data)[starters]

  return(bacteria_family_matrix_A)
}

set_up_b <- function(myPhyloSeq, level, target) {

  targets <- names_to_index(myPhyloSeq, target)
  my_bvectors <- list()
  bacteria_level <- phyloseq::tax_glom(myPhyloSeq, taxrank=level)
  bacteria_level_df <- as.data.frame(phyloseq::get_taxa(phyloseq::otu_table(bacteria_level)))
  for (i in 1:length(targets))
  {
    my_bvectors[[i]] <-bacteria_level_df[,targets[i]]
  }

  names(my_bvectors)=colnames(myPhyloSeq@otu_table@.Data)[targets]
  return(my_bvectors)

}

run_nnls<-function(matrixA, vectorb){

  my_weightsList <- list()

  for (i in 1:length(vectorb))
  {

    weightObject <-  nnls::nnls(as.matrix(matrixA),vectorb[[i]])

    my_weightsList[[i]] <- data.frame(weightObject$x, colnames=colnames(matrixA))
    names(my_weightsList)[i]=names(vectorb)[i]
  }

  return(my_weightsList)
}

names_to_index <- function(myPhyloSeq, starter) {

  starterSamplesNames <- strsplit(starter, " ")
  indices <- vector(mode = "integer", length = length(starterSamplesNames[[1]]))
  for( i in 1: length(starterSamplesNames[[1]])) {
    indices[i] <- which(colnames(myPhyloSeq@otu_table@.Data) == starterSamplesNames[[1]][i])
  }

  return(indices)
}
