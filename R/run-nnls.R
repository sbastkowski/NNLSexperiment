set_up_A <- function(myPhyloSeq, level, singles) {

  bacteria_level <- phyloseq::tax_glom(myPhyloSeq, taxrank=level)
  bacteria_level_df <- as.data.frame(phyloseq::get_taxa(otu_table(bacteria_level)))
  bacteria_level_singleComm_df  <-  bacteria_level_df[,singles]
  bacteria_family_matrix_A <- data.frame(as.matrix(bacteria_level_singleComm_df))
  colnames(bacteria_family_matrix_A)=colnames(myPhyloSeq@otu_table@.Data)[singles]

  return(bacteria_family_matrix_A)
}

set_up_b <- function(myPhyloSeq, level, mixed) {

  my_bvectors <- list()
  bacteria_level <- tax_glom(myPhyloSeq, taxrank=level)
  bacteria_level_df <- as.data.frame(get_taxa(otu_table(bacteria_level)))
  for (i in 1:length(mixed))
  {
    my_bvectors[[i]] <-bacteria_level_df[,mixed[i]]
  }

  names(my_bvectors)=colnames(myPhyloSeq@otu_table@.Data)[mixed]
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


