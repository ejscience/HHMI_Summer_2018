#Load the data from Dan
load("~/R/HHMI EXROP 18 R Folder/Lohman_RNA-Seq_Workspace/Heritability_Input.RData")
#I have here the 2 data frames: (1)rnaDesign.filt (2)counts.mean.filt

#Load up rethinking package
library(rethinking)

#adding the counts of the gene of interest to the dataframe
d <- counts.mean.filt["ENSGACG00000017764",] #take only the gene I want here it is Ndufs8
dd<- as.list(d) #make it into a list of values
MHC.Dataframe <- rnaDesign.filt #make new dataframe that will have all the information
MHC.Dataframe$counts <- dd #add the list to the new dataframe which is in the same order as the counts data

#adding the theoretical maximum reads variable
n_1 <- colSums(counts.mean.filt) #get the sum for each of the columns
MHC.Dataframe$n_max <- n_1 #place these values into the new variable n_max

#create dummy variable for Infection
MHC.Dataframe$inf <- ifelse(MHC.Dataframe$Population == "Infected", 1, 0)

#create an ID system for the fish
MHC.Dataframe$id <- 1:nrow(MHC.Dataframe)

#making a dataframe of just the data I need
MHC.Dataframe.small <- as.data.frame(MHC.Dataframe[,c("Batch","counts","n_max","id","inf")])
