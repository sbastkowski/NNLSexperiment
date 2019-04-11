
library("phyloseq")

# Simulate OTU tables that show little overlap

# n is number of OTUs
# k is number of samples

# i-j is range for counts
# We need to define overlap o, maybe in %

# First simple example:
# n=100, k=3, i=5, j=1000, o=0%
# all samples are weighted equally

realdata = import_biom("Data/otu_table_complete.txt.biom")
myTaxTable <- tax_table(realdata)
colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
OTU = otu_table(realdata)
OTU=OTU@.Data
# Estimate parameter from real data
muS=mean(OTU)
thetaS=MASS::theta.ml(OTU, mu=muS)

starter= "S1 S2 S3"
target= "S4 S5 S6"

# Create datastructure and sample for single communities, no overlap

SimOTU=matrix(0L, nrow = 90, ncol = 6)
SimOTU[1:30,1]=MASS::rnegbin(n=30, mu=muS, theta=thetaS)
SimOTU[31:60,2]=MASS::rnegbin(n=30, mu=muS, theta=thetaS)
SimOTU[61:90,3]=MASS::rnegbin(n=30, mu=muS, theta=thetaS)

# Add mixed communities
SimOTU[,4]=round(1/3*SimOTU[,1]+1/3*SimOTU[,2]+1/3*SimOTU[,3])
SimOTU[,5]=round(1/10*SimOTU[,1]+1/10*SimOTU[,2]+8/10*SimOTU[,3])
SimOTU[,6]=round(1/100*SimOTU[,1]+1/100*SimOTU[,2]+98/100*SimOTU[,3])

rownames(SimOTU) <- paste0("OTU", 1:nrow(SimOTU))
colnames(SimOTU) <- paste0("S", 1:ncol(SimOTU))

#Get taxonomy

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(SimOTU), ncol = 7)
rownames(taxmat) <- rownames(SimOTU)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")


SimOTU = otu_table(SimOTU, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

myObj <- phyloseq(SimOTU,TAX)


# Perform analysis

A=set_up_A(myObj, "Family", starter)
b=set_up_b(myObj, "Family", target)
weightSolutions=run_nnls(A, b)
write.table(weightSolutions, sep="\t", paste0("Test1", "_solutionWeights.txt"))
ggplot2::ggsave(file = paste0("Test1", "_plots.pdf") , make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")

#################################################################

#Overlapping OTUs 10 per pair

#################################################################

SimOTU=matrix(0L, nrow = 90, ncol = 6)
SimOTU[1:30,1]=MASS::rnegbin(n=30, mu=muS, theta=thetaS)
SimOTU[21:30,2]=SimOTU[21:30,1]
SimOTU[31:50,2]=MASS::rnegbin(n=30, mu=muS, theta=thetaS)
SimOTU[41:50,3]=SimOTU[41:50,2]
SimOTU[51:70,3]=MASS::rnegbin(n=30, mu=muS, theta=thetaS)

# Add mixed communities
SimOTU[,4]=round(1/3*SimOTU[,1]+1/3*SimOTU[,2]+1/3*SimOTU[,3])
SimOTU[,5]=round(1/10*SimOTU[,1]+1/10*SimOTU[,2]+8/10*SimOTU[,3])
SimOTU[,6]=round(1/100*SimOTU[,1]+1/100*SimOTU[,2]+98/100*SimOTU[,3])

rownames(SimOTU) <- paste0("OTU", 1:nrow(SimOTU))
colnames(SimOTU) <- paste0("Sample", 1:ncol(SimOTU))

#Get taxonomy

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(SimOTU), ncol = 7)
rownames(taxmat) <- rownames(SimOTU)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")


library("phyloseq")
OTU = otu_table(SimOTU, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

myObj <- phyloseq(OTU,TAX)


A=set_up_A(myObj, "Family", starter)
b=set_up_b(myObj, "Family", target)
weightSolutions=run_nnls(A, b)
write.table(weightSolutions, sep="\t", "Test2_refactored_solutionWeights.txt")
ggplot2::ggsave(file = "Test2.pdf", make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")

#################################################################
# Overlap of 10 between all samples
#################################################################

SimOTU=matrix(0L, nrow = 90, ncol = 6)
SimOTU[1:30,1]=sample(5:1000, 30, replace=TRUE)
SimOTU[21:30,2]=SimOTU[21:30,1]
SimOTU[31:50,2]=sample(5:1000, 20, replace=TRUE)
SimOTU[21:30,3]=SimOTU[21:30,2]
SimOTU[51:70,3]=sample(5:1000, 20, replace=TRUE)

# Add mixed communities
SimOTU[,4]=round(1/3*SimOTU[,1]+1/3*SimOTU[,2]+1/3*SimOTU[,3])
SimOTU[,5]=round(1/10*SimOTU[,1]+1/10*SimOTU[,2]+8/10*SimOTU[,3])
SimOTU[,6]=round(1/100*SimOTU[,1]+1/100*SimOTU[,2]+98/100*SimOTU[,3])

rownames(SimOTU) <- paste0("OTU", 1:nrow(SimOTU))
colnames(SimOTU) <- paste0("Sample", 1:ncol(SimOTU))

#Get taxonomy

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(SimOTU), ncol = 7)
rownames(taxmat) <- rownames(SimOTU)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")


OTU = otu_table(SimOTU, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

myObj <- phyloseq(OTU,TAX)

#Perform analysis

A=set_up_A(myObj, "Family", starter)
b=set_up_b(myObj, "Family", target)
weightSolutions=run_nnls(A, b)
write.table(weightSolutions, sep="\t", "Test3_refactored_solutionWeights.txt")
ggplot2::ggsave(file = "Test3.pdf", make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")

#################################################################

# Overlap of 20 between all and 10 between each pairs

#################################################################

SimOTU=matrix(0L, nrow = 60, ncol = 6)
SimOTU[1:30,1]=sample(5:1000, 30, replace=TRUE)
SimOTU[41:60,1]=sample(5:1000, 20, replace=TRUE)
SimOTU[1:30,2]=SimOTU[1:30,1]
SimOTU[31:40,2]=sample(5:1000, 10, replace=TRUE)
SimOTU[51:60,2]=sample(5:1000, 10, replace=TRUE)
SimOTU[1:20,3]=SimOTU[1:20,1]
SimOTU[31:40,3]=SimOTU[31:40,2]
SimOTU[41:50,3]=SimOTU[41:50,1]
SimOTU[51:60,3]=sample(5:1000, 10, replace=TRUE)

# Add mixed communities
SimOTU[,4]=round(1/3*SimOTU[,1]+1/3*SimOTU[,2]+1/3*SimOTU[,3])
SimOTU[,5]=round(1/10*SimOTU[,1]+1/10*SimOTU[,2]+8/10*SimOTU[,3])
SimOTU[,6]=round(1/100*SimOTU[,1]+1/100*SimOTU[,2]+98/100*SimOTU[,3])

rownames(SimOTU) <- paste0("OTU", 1:nrow(SimOTU))
colnames(SimOTU) <- paste0("Sample", 1:ncol(SimOTU))

#Get taxonomy

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(SimOTU), ncol = 7)
rownames(taxmat) <- rownames(SimOTU)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

OTU = otu_table(SimOTU, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

myObj <- phyloseq(OTU,TAX)

# Set up system of lin equations

A=set_up_A(myObj, "Family", starter)
b=set_up_b(myObj, "Family", target)
weightSolutions=run_nnls(A, b)
write.table(weightSolutions, sep="\t", "Test5_refactored_solutionWeights.txt")
ggplot2::ggsave(file = "Test5.pdf", make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")

#################################################################

# Overlap of 20 between all and 10 between each pairs different abundances

#################################################################

SimOTU=matrix(0L, nrow = 60, ncol = 6)
SimOTU[1:30,1]=sample(5:1000, 30, replace=TRUE)
SimOTU[41:60,1]=sample(5:1000, 20, replace=TRUE)
SimOTU[1:40,2]=sample(5:1000, 40, replace=TRUE)
SimOTU[51:60,2]=sample(5:1000, 10, replace=TRUE)
SimOTU[1:20,3]=sample(5:1000, 20, replace=TRUE)
SimOTU[31:60,3]=sample(5:1000, 30, replace=TRUE)

# Add mixed communities
SimOTU[,4]=round(1/3*SimOTU[,1]+1/3*SimOTU[,2]+1/3*SimOTU[,3])
SimOTU[,5]=round(1/10*SimOTU[,1]+1/10*SimOTU[,2]+8/10*SimOTU[,3])
SimOTU[,6]=round(1/100*SimOTU[,1]+1/100*SimOTU[,2]+98/100*SimOTU[,3])

rownames(SimOTU) <- paste0("OTU", 1:nrow(SimOTU))
colnames(SimOTU) <- paste0("Sample", 1:ncol(SimOTU))

#Get taxonomy

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(SimOTU), ncol = 7)
rownames(taxmat) <- rownames(SimOTU)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

OTU = otu_table(SimOTU, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

myObj <- phyloseq(OTU,TAX)

# Perform analysis

A=set_up_A(myObj, "Family", starter)
b=set_up_b(myObj, "Family", target)
weightSolutions=run_nnls(A, b)
write.table(weightSolutions, sep="\t", "Test5_2_refactored_solutionWeights.txt")
ggplot2::ggsave(file = "Test5_2.pdf", make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")

#################################################################

# all the same, sanity check

#################################################################

SimOTU=matrix(0L, nrow = 60, ncol = 6)
SimOTU[1:60,1]=sample(5:1000, 60, replace=TRUE)
SimOTU[1:60,2]=SimOTU[1:60,1]
SimOTU[1:60,3]=SimOTU[1:60,2]


# Add mixed communities
SimOTU[,4]=round(1/3*SimOTU[,1]+1/3*SimOTU[,2]+1/3*SimOTU[,3])
SimOTU[,5]=round(1/10*SimOTU[,1]+1/10*SimOTU[,2]+8/10*SimOTU[,3])
SimOTU[,6]=round(1/100*SimOTU[,1]+1/100*SimOTU[,2]+98/100*SimOTU[,3])

rownames(SimOTU) <- paste0("OTU", 1:nrow(SimOTU))
colnames(SimOTU) <- paste0("Sample", 1:ncol(SimOTU))

#Get taxonomy

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(SimOTU), ncol = 7)
rownames(taxmat) <- rownames(SimOTU)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")


OTU = otu_table(SimOTU, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
myObj <- phyloseq(OTU,TAX)

# Perform analysis

A=set_up_A(myObj, "Family", starter)
b=set_up_b(myObj, "Family", target)
weightSolutions=run_nnls(A, b)
write.table(weightSolutions, sep="\t", "Test6_refactored_solutionWeights.txt")
ggplot2::ggsave(file = "Test6.pdf", make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")

#################################################################

# Using real data

#################################################################


pool="p1=S01,S02,S03;p2=S04,S05,S06; p3=S07,S08,S09; t1=S10; t2=S10; t3=S10; t4=S10; t5=S10;t6=S10;t7=S10"
starter="p1 p2 p3"
target="t1 t2 t3 t4 t5 t6 t7"

data = realdata
myTaxTable <- tax_table(data)
colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
TAX = tax_table(myTaxTable)
OTU = otu_table(data)[,1:10]

myObj <- phyloseq(OTU,TAX)
pool=pool_to_df(myObj,pool)
pooledOTU=pool_samples(myObj, pool)
OTU = otu_table(pooledOTU)

#Replace with simulated counts
OTU[,4] = round(1/3*OTU[,1]+1/3*OTU[,2]+1/3*OTU[,3])
OTU[,7] = round(1/10*OTU[,1]+1/10*OTU[,2]+8/10*OTU[,3])
OTU[,10] = round(1/100*OTU[,1]+1/100*OTU[,2]+98/100*OTU[,3])
OTU[,6] = round(1/10*OTU[,1]+8/10*OTU[,2]+1/10*OTU[,3])
OTU[,9] = round(1/100*OTU[,1]+98/100*OTU[,2]+1/100*OTU[,3])
OTU[,5] = round(8/10*OTU[,1]+1/10*OTU[,2]+1/10*OTU[,3])
OTU[,8] = round(98/100*OTU[,1]+1/100*OTU[,2]+1/100*OTU[,3])

TAX = tax_table(myTaxTable)
myObj <- phyloseq(OTU,TAX)

# Perform analysis

A=set_up_A(myObj, "Family", starter)
b=set_up_b(myObj, "Family", target)
weightSolutions=run_nnls(A, b)
write.table(weightSolutions, sep="\t", "Test7_refactored_solutionWeights.txt")
ggplot2::ggsave(file = "Test7.pdf", make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")


#################################################################

# Consider presence and absence

#################################################################

pool="p1=S01,S02,S03;p2=S04,S05,S06; p3=S07,S08,S09; t1=S10,S11,S12; t2=S13,S14,S15; t3=S16,S17,S18; t4=S19,S20,S21; t5=S22,S23,S24;t6=S25,S26,S27;t7=S28,S29,S30;"
starter="p1 p2 p3"
target="t1 t2 t3 t4 t5 t6 t7"

presence_absence<-function(pooledOTU) {


    for(i in 1:nrow(pooledOTU)){
      for(j in 1:ncol(pooledOTU)){
        if(pooledOTU[i,j]!=0){
          pooledOTU[i,j]=1
        }
      }
    }
  row.names(pooledOTU)=row.names(tax_table(phyloseqObj))


  return(otu_table(pooledOTU, taxa_are_rows = TRUE))
}

data = realdata
myTaxTable <- tax_table(data)
colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
OTU = otu_table(data)

TAX = tax_table(myTaxTable)
myPhyloSeq <- phyloseq(OTU,TAX)

pool=pool_to_df(myPhyloSeq,pool)
pooledOTU=pool_samples(myPhyloSeq, pool)
OTU_pooled=presence_absence(pooledOTU)
myObj <- merge_phyloseq(OTU_pooled, TAX)

# Set up system of lin equations

A=set_up_A(myObj, "Family", starter)
b=set_up_b(myObj, "Family", target)
weightSolutions=run_nnls(A, b)
write.table(weightSolutions, sep="\t", "refactored_solutionWeightsPresenceAbsence.txt")
ggplot2::ggsave(file = "PresenceAbsence.pdf", make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")

#################################################################

#Presence ansence on sim data

#################################################################
pool="P1=S01,S02,S03; P2=S04,S05,S06; P3=S07,S08,S09; P4=S10; P5=S11; P6=S12; P7=S13; P8=S14; P9=S15; P10=S16"
starter="P1 P2 P3"
target="P4 P5 P6 P7 P8 P9 P10"
data = realdata
myTaxTable <- tax_table(data)
colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
TAX = tax_table(myTaxTable)
OTU = otu_table(data)[,1:16]

myObj <- phyloseq(OTU,TAX)
pool=pool_to_df(myObj,pool)
pooledOTU=pool_samples(myObj, pool)
OTU = otu_table(pooledOTU)

OTU[,4] = round(1/3*OTU[,1]+1/3*OTU[,2]+1/3*OTU[,3])
OTU[,7] = round(1/10*OTU[,1]+1/10*OTU[,2]+8/10*OTU[,3])
OTU[,10] = round(1/100*OTU[,1]+1/100*OTU[,2]+98/100*OTU[,3])
OTU[,6] = round(1/10*OTU[,1]+8/10*OTU[,2]+1/10*OTU[,3])
OTU[,9] = round(1/100*OTU[,1]+98/100*OTU[,2]+1/100*OTU[,3])
OTU[,5] = round(8/10*OTU[,1]+1/10*OTU[,2]+1/10*OTU[,3])
OTU[,8] = round(98/100*OTU[,1]+1/100*OTU[,2]+1/100*OTU[,3])

OTU_pooled=presence_absence(OTU)
myObj <- merge_phyloseq(OTU_pooled, TAX)

# Set up system of lin equations

A=set_up_A(myObj, "Family", starter)
b=set_up_b(myObj, "Family", target)
weightSolutions=run_nnls(A, b)
write.table(weightSolutions, sep="\t", "SIMrefactored_solutionWeightsPresenceAbsence.txt")
ggplot2::ggsave(file = "SIMPresenceAbsence.pdf", make_barcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")
#####################################################################

#Assess the difference between the seeds (simple difference)

pool="P1=S01,S02,S03; P2=S04,S05,S06; P3=S07,S08,S09; P4=S10; P5=S11; P6=S12"
starter="P1 P2 P3"
target="P4 P5 P6"

data = realdata
myTaxTable <- tax_table(data)
colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
TAX = tax_table(myTaxTable)
OTU = otu_table(data)[,1:12]

myObj <- phyloseq(OTU,TAX)
pool=pool_to_df(myObj,pool)
pooledOTU=pool_samples(myObj, pool)
OTU = otu_table(pooledOTU)
dist12=0
dist13=0
dist23=0
for(i in 1:nrow(myObj@otu_table)){
  dist12=dist12+abs(OTU@.Data[i,1]-OTU@.Data[i,2])
  dist13=dist13+abs(OTU@.Data[i,1]-OTU@.Data[i,3])
  dist23=dist23+abs(OTU@.Data[i,2]-OTU@.Data[i,3])
}
#Assess the difference between the seeds using presence absence(simple difference)
pool="P1=S01,S02,S03; P2=S04,S05,S06; P3=S07,S08,S09; P4=S10; P5=S11; P6=S12; P7=S13; P8=S14; P9=S15; P10=S16"
starter="P1 P2 P3"
target="P4 P5 P6 P7 P8 P9 P10"
data = realdata
myTaxTable <- tax_table(data)
colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
TAX = tax_table(myTaxTable)
OTU = otu_table(data)[,1:16]

myObj <- phyloseq(OTU,TAX)
pool=pool_to_df(myObj,pool)
pooledOTU=pool_samples(myObj, pool)
OTU = otu_table(pooledOTU)
OTU=presence_absence(OTU)

dist12p=0
dist13p=0
dist23p=0
for(i in 1:nrow(myObj@otu_table)){
  dist12p=dist12p+abs(OTU@.Data[i,1]-OTU@.Data[i,2])
  dist13p=dist13p+abs(OTU@.Data[i,1]-OTU@.Data[i,3])
  dist23p=dist23p+abs(OTU@.Data[i,2]-OTU@.Data[i,3])
}



# Counting number of species

#pool="P1=S01,S02,S03; P2=S04,S05,S06; P3=S07,S08,S09; P4=S10; P5=S11; P6=S12; P7=S13; P8=S14; P9=S15; P10=S16"

data = realdata
myTaxTable <- tax_table(data)
colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
TAX = tax_table(myTaxTable)
OTU = otu_table(data)

SpeciesSeed1=sum(OTU@.Data[,1])
SpeciesSeed2=sum(OTU@.Data[,2])
SpeciesSeed3=sum(OTU@.Data[,3])

SpeciesSeed4=sum(OTU@.Data[,4])
SpeciesSeed5=sum(OTU@.Data[,5])
SpeciesSeed6=sum(OTU@.Data[,6])

SpeciesSeed7=sum(OTU@.Data[,7])
SpeciesSeed8=sum(OTU@.Data[,8])
SpeciesSeed9=sum(OTU@.Data[,9])

SeedA=SpeciesSeed1+SpeciesSeed2+SpeciesSeed3
SeedB=SpeciesSeed4+SpeciesSeed5+SpeciesSeed6
SeedC=SpeciesSeed7+SpeciesSeed8+SpeciesSeed9

OTU=presence_absence(OTU)

SpeciesSeed1=sum(OTU@.Data[,1])
SpeciesSeed2=sum(OTU@.Data[,2])
SpeciesSeed3=sum(OTU@.Data[,3])

SpeciesSeed4=sum(OTU@.Data[,4])
SpeciesSeed5=sum(OTU@.Data[,5])
SpeciesSeed6=sum(OTU@.Data[,6])

SpeciesSeed7=sum(OTU@.Data[,7])
SpeciesSeed8=sum(OTU@.Data[,8])
SpeciesSeed9=sum(OTU@.Data[,9])

SeedA=SpeciesSeed1+SpeciesSeed2+SpeciesSeed3
SeedB=SpeciesSeed4+SpeciesSeed5+SpeciesSeed6
SeedC=SpeciesSeed7+SpeciesSeed8+SpeciesSeed9

#Find average abundance per seed

subA_1=subset(OTU@.Data[,1], OTU@.Data[,1]!=0)
subA_2=subset(OTU@.Data[,2], OTU@.Data[,2]!=0)
subA_3=subset(OTU@.Data[,3], OTU@.Data[,3]!=0)

mean(c(subA_1,subA_2,subA_3))

subB_1=subset(OTU@.Data[,4], OTU@.Data[,4]!=0)
subB_2=subset(OTU@.Data[,5], OTU@.Data[,5]!=0)
subB_3=subset(OTU@.Data[,6], OTU@.Data[,6]!=0)

mean(c(subB_1,subB_2,subB_3))

subC_1=subset(OTU@.Data[,7], OTU@.Data[,7]!=0)
subC_2=subset(OTU@.Data[,8], OTU@.Data[,8]!=0)
subC_3=subset(OTU@.Data[,9], OTU@.Data[,9]!=0)

mean(c(subC_1,subC_2,subC_3))

