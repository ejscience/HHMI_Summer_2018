#Load the data from Dan
load("~/R/HHMI EXROP 18 R Folder/Lohman_RNA-Seq_Workspace/Heritability_Input.RData")
#I have here the 2 data frames: (1)rnaDesign.filt (2)counts.mean.filt

#Load up rethinking package
library(rethinking)

#adding the counts of the gene of interest to the dataframe
d <- counts.mean.filt["ENSGACG00000000336",] #take only the gene I want here it is MHC II
dd<- as.list(d) #make it into a list of values
MHC.Dataframe <- rnaDesign.filt #make new dataframe that will have all the information
MHC.Dataframe$counts <- dd #add the list to the new dataframe which is in the same order as the counts data

#adding the theoretical maximum reads variable
n_1 <- colSums(counts.mean.filt) #get the sum for each of the columns
MHC.Dataframe$n_max <- n_1 #place these values into the new variable n_max

#create dummy variable for Roberts
MHC.Dataframe$Rob <- ifelse(MHC.Dataframe$Population == "Roberts", 1, 0)

#create an ID system for the fish
MHC.Dataframe$id <- 1:nrow(MHC.Dataframe)

##### Trying to make an ID for batch...
MHC.Dataframe$batch_id <- as.numeric(MHC.Dataframe$Batch)

#making a dataframe of just the data I need
MHC.Dataframe.small <- as.data.frame(MHC.Dataframe[,c("Batch","counts","n_max","id","Rob")])


#the model:

#made the list first to write out all my information

flist <- alist(
  #likelihood
  counts ~ dbetabinom(n_max, pbar, theta),
  
  #linear model
  logit(pbar) <- A + b_Rob * Rob,
  A <- a + a_batch[Batch] + a_ind[id],
  
  #adaptive prior
  a_batch[Batch] ~ dnorm(0, sigma_batch), 
  a_ind[id] ~ dnorm(0,sigma_id),
  
  #fixed priors
  a ~ dnorm(0,10),
  b_Rob ~ dnorm(0,10),
  theta ~ dexp(1),
  sigma_batch ~ dcauchy(0,1),
  sigma_id ~ dcauchy(0,1)
)

#then put the stan model after so that I can read the details of it easier
mMHC1 <- map2stan(flist, data=MHC.Dataframe.small, iter = 2000, chains = 1)
