
#####################################################################
#Test data-prep
#####################################################################

test_that("read_data",{
  expectedOTU=read.csv("Testdata/Expected_read_data_OTUtable.csv", header=TRUE, row.names = 1)
  expectedTax=read.csv("Testdata/Expected_read_data_Taxtable.csv", header=TRUE, row.names = 1)
  to_test=read_data("Testdata/otu_table_complete.txt.biom", "Testdata/mappingFile.txt")
  to_test_OTU=as.data.frame(to_test@otu_table@.Data)
  to_test_Tax=as.data.frame(to_test@tax_table@.Data)
  expect_equal( to_test_OTU, expectedOTU)
  expect_equal( dim(to_test_Tax), dim(expectedTax))
})

test_that("pool_to_df",{
  pool="p1=S01,S02,S03;p2=S04,S05,S06; p3=S07,S08,S09; t1=S10,S11,S12; t2=S13,S14,S15; t3=S16,S17,S18; t4=S19,S20,S21; t5=S22,S23,S24;t6=S25,S26,S27;t7=S28,S29,S30;"
  OTU=phyloseq::import_biom("Testdata/otu_table_complete.txt.biom")
  to_test=pool_to_df(OTU,pool)
  expected=read.csv("Testdata/test_pool_df.csv", header=TRUE, row.names = 1)
  expect_equal(expected, to_test)
})

test_that("pool_samples",{
  pool="p1=S01,S02,S03;p2=S04,S05,S06; p3=S07,S08,S09; t1=S10,S11,S12; t2=S13,S14,S15; t3=S16,S17,S18; t4=S19,S20,S21; t5=S22,S23,S24;t6=S25,S26,S27;t7=S28,S29,S30;"
  OTU=phyloseq::import_biom("Testdata/otu_table_complete.txt.biom")
  pool_df=pool_to_df(OTU,pool)
  to_test=as.data.frame(pool_samples(OTU, pool_df))
  expected=read.csv("Testdata/expectedOTU_pooled.csv", header=TRUE, row.names = 1)
  expect_equal(expected, to_test)
})

#####################################################################
#Test setting up the problem and running solver
#####################################################################

test_that("set_up_A",{

  level="Family"
  starter= "S1 S2 S3"
  OTU=phyloseq::import_biom("Testdata/SimOTU.biom")
  taxmat=read.csv("Testdata/SimTax.csv",row.names = 1, stringsAsFactors = FALSE)
  taxmat=as.matrix(taxmat)
  rownames(taxmat) <- rownames(OTU)
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  TAX = phyloseq::tax_table(taxmat)
  myObj <-phyloseq::phyloseq(OTU,TAX)
  to_test=set_up_A(myObj, level, starter)
  expected=read.csv("Testdata/expected_A.csv",header=TRUE,row.names = 1, stringsAsFactors = FALSE)
  expect_equal(expected, to_test)
})

test_that("set_up_b",{
  level="Family"
  target= "S4 S5 S6 S7 S8 S9 S10"
  OTU=phyloseq::import_biom("Testdata/SimOTU.biom")
  taxmat=read.csv("Testdata/SimTax.csv",row.names = 1, stringsAsFactors = FALSE)
  taxmat=as.matrix(taxmat)
  rownames(taxmat) <- rownames(OTU)
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  TAX = phyloseq::tax_table(taxmat)
  myObj <- phyloseq::phyloseq(OTU,TAX)
  to_test=as.data.frame(set_up_b(myObj, "Family", target))
  expected=read.csv("Testdata/expected_b.csv",header=TRUE,row.names = 1, stringsAsFactors = FALSE)
  expect_equal(expected, to_test)
          })

test_that("run_nnls",{
  level="Family"
  starter= "S1 S2 S3"
  target= "S4 S5 S6 S7 S8 S9 S10"
  OTU=phyloseq::import_biom("Testdata/SimOTU.biom")
  taxmat=read.csv("Testdata/SimTax.csv",row.names = 1, stringsAsFactors = FALSE)
  taxmat=as.matrix(taxmat)
  rownames(taxmat) <- rownames(OTU)
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  TAX = phyloseq::tax_table(taxmat)
  myObj <- phyloseq::phyloseq(OTU,TAX)
  A=set_up_A(myObj, level, starter)
  b=set_up_b(myObj, "Family", target)
  to_test=as.data.frame(run_nnls(A, b))
  expected=read.csv("Testdata/expected_run_nnls_solutions.csv",header=TRUE,row.names = 1)
  expect_equal(expected, to_test)

})

#####################################################################
#Test ploting
#####################################################################

test_that("weight_to_percent",{
  level="Family"
  starter= "S1 S2 S3"
  target= "S4 S5 S6 S7 S8 S9 S10"
  OTU=phyloseq::import_biom("Testdata/SimOTU.biom")
  taxmat=read.csv("Testdata/SimTax.csv",row.names = 1, stringsAsFactors = FALSE)
  taxmat=as.matrix(taxmat)
  rownames(taxmat) <- rownames(OTU)
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  TAX = phyloseq::tax_table(taxmat)
  myObj <- phyloseq::phyloseq(OTU,TAX)
  A=set_up_A(myObj, level, starter)
  b=set_up_b(myObj, "Family", target)
  weights=run_nnls(A, b)
  to_test=as.data.frame(weight_to_percent(weights))
  expected=read.csv("Testdata/expected_weights_to_percentage.csv", header=TRUE, row.names = 1)
  expect_equal(expected, to_test)
})


test_that("make_barcharts",{
  level="Family"
  starter= "S1 S2 S3"
  target= "S4 S5 S6 S7 S8 S9 S10"
  OTU=phyloseq::import_biom("Testdata/SimOTU.biom")
  taxmat=read.csv("Testdata/SimTax.csv",row.names = 1, stringsAsFactors = FALSE)
  taxmat=as.matrix(taxmat)
  rownames(taxmat) <- rownames(OTU)
  colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  TAX = phyloseq::tax_table(taxmat)
  myObj <- phyloseq::phyloseq(OTU,TAX)
  A=set_up_A(myObj, level, starter)
  b=set_up_b(myObj, "Family", target)
  weights=run_nnls(A, b)
  to_test=make_barcharts(weights, percent=TRUE)
  expect_true(length(to_test[[1]]$grobs)==8)
})

#####################################################################
#Test running the whole thing
#####################################################################


test_that("Test run_nnls_analysis",{


  level="Family"
  biom="Testdata/otu_table_complete.txt.biom"
  mapping="Testdata/mappingFile.txt"
  pool="p1=S01,S02,S03;p2=S04,S05,S06; p3=S07,S08,S09; t1=S10,S11,S12; t2=S13,S14,S15; t3=S16,S17,S18; t4=S19,S20,S21; t5=S22,S23,S24;t6=S25,S26,S27;t7=S28,S29,S30;"
  starter = "p1 p2 p3"
  target = "t1 t2 t3 t4 t5 t6 t7"
  to_test=as.data.frame(run_nnls_analysis( biom, mapping, starter, target, level = "Family", pool))

  expected=read.csv("Testdata/expected_run_full_nnls_solutions.csv",header=TRUE,row.names = 1)
  expect_equal(expected, to_test)
})
