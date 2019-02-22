run_nnls_analysis <- function( biomfile, mapping, starter, target, level = "Family", pool = NULL){

  myDataObj=read_data(biomfile, mapping, pool)
  A=set_up_A(myDataObj, level, starter)
  b=set_up_b(myDataObj, level, target)
  weightSolutions=run_nnls(A, b)

  return(weightSolutions)

}
