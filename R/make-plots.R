weight_to_percent<-function(weights){

  percentages=list()
  for (i in 1:length(weights))
  {
    percent=vector(mode="integer", length = length(weights[[i]]$weightObject.x))
    for (j in 1:length(weights[[i]]$weightObject.x)){
      percent[j] <-  weights[[i]]$weightObject.x[j]/sum(weights[[i]]$weightObject.x)
    }
    percentages[[i]]=percent
  }

  return(percentages)
}

make_barcharts <- function(weightSolutions, percent=TRUE){


  if (percent==TRUE){
    p=weight_to_percent(weights = weightSolutions)
    trellis.objects.list = list()

    for (i in 1:length(weightSolutions))
    {
      trellis.objects.list[[i]] <- lattice::barchart(as.table(p[[i]]),
                                            main=names(weightSolutions)[i],
                                            horizontal=FALSE, col="steelblue", ylab="Weight", xlab = "Samples", aspect=1,
                                            scales=list(x=list(rot=70, labels=weightSolutions[[i]]$colnames, cex=1.1)))
    }

    multiPageGrobs <-  gridExtra::marrangeGrob(grobs = trellis.objects.list, nrow=4, ncol=2)

  }else{

    trellis.objects.list = list()

    for (i in 1:length(weightSolutions))
    {
      trellis.objects.list[[i]] <- lattice::barchart(as.table(weightSolutions[[i]]$weightObject.x),
                                            main=names(weightSolutions)[i],
                                            horizontal=FALSE, col="steelblue", ylab="Weight", xlab = "Samples", aspect=1,
                                            scales=list(x=list(rot=70, labels=weightSolutions[[i]]$colnames, cex=1.1)))
    }

    multiPageGrobs <- gridExtra::marrangeGrob(grobs = trellis.objects.list, nrow=4, ncol=2)

  }

  return(multiPageGrobs)
}
