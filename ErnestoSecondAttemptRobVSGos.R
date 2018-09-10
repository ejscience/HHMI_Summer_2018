#Newest attempt to make the data frames to analyze ----------------------

#Load the data from Dan
load("~/R/HHMI EXROP 18 R Folder/Lohman_RNA-Seq_Workspace/Heritability_Input.RData")
#I have here the 2 data frames: (1)rnaDesign.filt (2)counts.mean.filt

#Load up rethinking package
library(rethinking)

#make the count data in the same design as the rnaDesign.filt
transformed_count <- as.data.frame(t(counts.mean.filt))

#make new dataframe that will have all the information I will add
MHC.Dataframe <- rnaDesign.filt

#take only the gene I want here it is MHC II
MHC.Dataframe$MHC_counts <- as.integer(transformed_count$ENSGACG00000000336)


#finding the theoretical count maximum
n_1 <- rowSums(transformed_count) #get the sum for each of the rows for the transformed_count
MHC.Dataframe$n_max <- as.integer(n_1) #place these values into the new variable n_max

#create dummy variable for Roberts
MHC.Dataframe$Rob <- as.integer(ifelse(MHC.Dataframe$Population == "Roberts", 1, 0))

#create an ID system for the fish
MHC.Dataframe$fish_id <- as.integer(1:nrow(MHC.Dataframe))

#Trying to make an ID (dummy variable) for batch...
MHC.Dataframe$batch_id <- as.integer(MHC.Dataframe$Batch)

#check the structure of the dataframe
str(MHC.Dataframe)




#### ***FIXED AUG 7th!!! BetaBinom Model w/o fish_id and w/o hyperpriors ------------

#making a dataframe #making a dataframe of just the data I need
MHC.Dataframe.small <- as.data.frame(MHC.Dataframe[,c("batch_id","MHC_counts","Rob","n_max")])

#made the list first to write out all my information

flist <- alist(
  #likelihood
  MHC_counts ~ dbetabinom(n_max, pbar, theta),
  
  #linear model
  logit(pbar) <- A + b_Rob * Rob,
  A <- a + a_batch[batch_id],
  
  #adaptive prior
  a_batch[batch_id] ~ dnorm(avg_batch, sigma_batch), 
  
  #fixed priors
  a ~ dnorm(0,10),
  b_Rob ~ dnorm(0,10),
  theta ~ dcauchy(0,1),
  avg_batch ~ dnorm(0,10),
  sigma_batch ~ dcauchy(0,1)
)

#then put the stan model after so that I can read the details of it easier
mMHC1 <- map2stan(flist, data=MHC.Dataframe.small, iter = 2000, chains = 1)


#### Poisson Model -------------------

#making a dataframe of just the data I need
MHC.Dataframe.small <- as.data.frame(MHC.Dataframe[,c("batch_id","MHC_counts","fish_id","Rob")])

#Check the structure again
str(MHC.Dataframe.small)

flist <- alist(
  #likelihood
  counts ~ dpois( lambda ),
  
  #linear model
  log(lambda) <- A + bR * Rob,
  A <- a + a_batch[batch_id] + a_ind[fish_id],
  
  #adaptive prior
  a_batch[batch_id] ~ dnorm(a_bt, sigma_batch), 
  a_ind[fish_id] ~ dnorm(a_fs,sigma_id),
  
  #fixed priors
  c(a, a_bt, a_fs) ~ dnorm(0,10),
  bR ~ dnorm(0,10),
  sigma_batch ~ dcauchy(0,1),
  sigma_id ~ dcauchy(0,1)
)

#then put the stan model after so that I can read the details of it easier
mMHC1 <- map2stan(flist, data=MHC.Dataframe.small, iter = 2000, chains = 1)




