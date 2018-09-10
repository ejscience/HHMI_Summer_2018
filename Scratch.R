##First Attempts RobVSGos Scratch-----------
#Load the data from Dan
load("~/R/HHMI EXROP 18 R Folder/Lohman_RNA-Seq_Workspace/Heritability_Input.RData")

#Load up rethinking package
library(rethinking)

head(counts.mean.filt)
head(rnaDesign.filt)

load("~/R/HHMI EXROP 18 R Folder/Lohman_RNA-Seq_Workspace/Heritability_Input.RData")
d <- counts.mean.filt["ENSGACG00000000336",] #take only the gene I want
dd<- as.list(d) #make it into a list of values
MHC.Dataframe <- rnaDesign.filt #make new dataframe that will have all the information
MHC.Dataframe$counts <- dd #add the list to the new dataframe which is in the same order as the counts data

str(MHC.Dataframe)
range(MHC.Dataframe$counts)

#get the theoretical maximum. 
n_1 <- colSums(counts.mean.filt) #get the sum for each of the columns
MHC.Dataframe$n_max <- n_1 #place these values into the new variable n_max


#We are going to make a

alist(
  #likelihood
  counts ~ dbetabinom(n, pbar, theta),
  
  
  #linear model
  
  #adaptive prior
  
  #fixed priors
  theta ~ dexp(1)
)


transformed.counts <- as.data.frame(t(counts.mean.filt))





##Rethinking stuff as reference for how to make a model--------------

## R code 13.22
library(rethinking)
data(chimpanzees)
d <- chimpanzees
d$recipient <- NULL
d$block_id <- d$block

m13.6 <- map2stan(
  alist(
    # likeliood
    pulled_left ~ dbinom(1,p),
    
    # linear models
    logit(p) <- A + (BP + BPC*condition)*prosoc_left,
    A <- a + a_actor[actor] + a_block[block_id],
    BP <- bp + bp_actor[actor] + bp_block[block_id],
    BPC <- bpc + bpc_actor[actor] + bpc_block[block_id],
    
    # adaptive priors
    c(a_actor,bp_actor,bpc_actor)[actor] ~
      dmvnorm2(0,sigma_actor,Rho_actor),
    c(a_block,bp_block,bpc_block)[block_id] ~
      dmvnorm2(0,sigma_block,Rho_block),
    
    # fixed priors
    c(a,bp,bpc) ~ dnorm(0,1),
    sigma_actor ~ dcauchy(0,2),
    sigma_block ~ dcauchy(0,2),
    Rho_actor ~ dlkjcorr(4),
    Rho_block ~ dlkjcorr(4)
  ) , data=d , iter=5000 , warmup=1000 , chains=3 , cores=3 )






#MAP shit I am playing with messing around with the model--------

m.MHC <- map2stan(
  alist(
    #likelihood
    counts ~ dbetabinom(n_max, pbar, theta),
    
    #linear model
    logit(pbar) <- A + b_Rob * Rob,
    A <- a + a_batch[batch] + a_ind[id],
    
    

    
    #adaptive prior
    a_batch[batch] ~ dnorm(0, sigma_batch), 
    a_ind[id] ~ dnorm(0,sigma_id),
    
    #fixed priors
    a ~ dnorm(0,10),
    b_Rob ~ dnorm(0,10),
    theta ~ dexp(1),
    sigma_batch ~ dcauchy(0,1),
    sigma_id ~ dcauchy(0,1)
)
)




#Trying to make an ID for batch---------------
#adding the counts of the gene of interest to the dataframe
d <- counts.mean.filt["ENSGACG00000000336",] #take only the gene I want here it is MHC II
dd<- as.list(d) #make it into a list of values
MHC.Dataframe <- rnaDesign.filt #make new dataframe that will have all the information
MHC.Dataframe$counts <- dd #add the list to the new dataframe which is in the same order as the counts data

#adding the theoretical maximum reads variable
n_1 <- as.list(colSums(counts.mean.filt)) #get the sum for each of the columns
MHC.Dataframe$n_max <- n_1 #place these values into the new variable n_max

#create dummy variable for Roberts
MHC.Dataframe$Rob <- ifelse(MHC.Dataframe$Population == "Roberts", 1, 0)

#create an ID system for the fish
MHC.Dataframe$id <- 1:nrow(MHC.Dataframe)

##### Trying to make an ID for batch...
MHC.Dataframe$batch_id <- as.numeric(MHC.Dataframe$Batch)

#making a dataframe of just the data I need
MHC.Dataframe.small <- as.data.frame(MHC.Dataframe[,c("batch_id","counts","n_max","id","Rob")])


#the model:

#made the list first to write out all my information

flist <- alist(
  #likelihood
  counts ~ dbetabinom(n_max, pbar, theta),
  
  #linear model
  logit(pbar) <- A + b_Rob * Rob,
  A <- a + a_batch[batch_id] + a_ind[id],
  
  #adaptive prior
  a_batch[batch_id] ~ dnorm(0, sigma_batch), 
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








### NOTES ON ATTEMPT
#I tried making the batch_id no soln
#I tried making the n_max a list then putting it into the MHC II
#I tried making switching the n_max and pbar





#### OKAY IMMA TRY POISSON...----------------


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
MHC.Dataframe$fish_id <- 1:nrow(MHC.Dataframe)

##### Trying to make an ID (dummy variable) for batch...
MHC.Dataframe$batch_id <- as.numeric(MHC.Dataframe$Batch)

#making a dataframe of just the data I need
MHC.Dataframe.small <- as.data.frame(MHC.Dataframe[,c("batch_id","counts","fish_id","Rob")])



flist <- alist(
  #likelihood
  counts ~ dpois(m_lambda),
  
  #linear model
  log(m_lambda) <- A + bR * Rob,
  A <- a + a_batch[batch_id] + a_ind[fish_id],
  
  #adaptive prior
  a_batch[batch_id] ~ dnorm(0, sigma_batch), 
  a_ind[fish_id] ~ dnorm(0,sigma_id),
  
  #fixed priors
  a ~ dnorm(0,10),
  bR ~ dnorm(0,10),
  sigma_batch ~ dcauchy(0,1),
  sigma_id ~ dcauchy(0,1)
)
mMHC1 <- map2stan(flist, data=MHC.Dataframe.small, iter = 2000, chains = 1)



##### I'm going to try to remove the a_ind[id]---------------

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

##### Trying to make an ID for batch...
MHC.Dataframe$batch_id <- as.numeric(MHC.Dataframe$Batch)

#making a dataframe of just the data I need
MHC.Dataframe.small <- as.data.frame(MHC.Dataframe[,c("batch_id","counts","n_max","Rob")])


#the model:

#made the list first to write out all my information

flist <- alist(
  #likelihood
  counts ~ dbetabinom(n_max, pbar, theta),
  
  #linear model
  logit(pbar) <- a + a_batch[batch_id] + b_Rob * Rob,
  
  #adaptive prior
  a_batch[batch_id] ~ dnorm(0, sigma_batch), 
  
  #fixed priors
  a ~ dnorm(0,10),
  b_Rob ~ dnorm(0,10),
  theta ~ dexp(1),
  sigma_batch ~ dcauchy(0,1),
)

#then put the stan model after so that I can read the details of it easier
mMHC1 <- map2stan(flist, data=MHC.Dataframe.small)

## Now I am going to try without hyperpriors----------

flist <- alist(
  #likelihood
  counts ~ dpois( lambda ),
  
  #linear model
  log(lambda) <- A + bR * Rob,
  A <- a + a_batch[batch_id] + a_ind[fish_id],
  
  #adaptive prior
  a_batch[batch_id] ~ dnorm(0, sigma_batch), 
  a_ind[fish_id] ~ dnorm(0,sigma_id),
  
  #fixed priors
  a ~ dnorm(0,10),
  bR ~ dnorm(0,10),
  sigma_batch ~ dcauchy(0,1),
  sigma_id ~ dcauchy(0,1)
)

#then put the stan model after so that I can read the details of it easier
mMHC1 <- map2stan(flist, data=MHC.Dataframe.small, iter = 2000, chains = 1)




#### BetaBinom Model w/o fish_id ------------

#making a dataframe of just the data I need
MHC.Dataframe.small <- as.data.frame(MHC.Dataframe[,c("batch_id","MHC_counts","Rob","n_max")])

#made the list first to write out all my information

flist <- alist(
  #likelihood
  MHC_counts ~ dbetabinom(n_max, pbar, theta),
  
  #linear model
  logit(pbar) <- A + b_Rob * Rob,
  A <- a + a_batch[batch_id],
  
  #adaptive prior
  a_batch[batch_id] ~ dnorm(0, 10), 
  
  #fixed priors
  a ~ dnorm(0,10),
  b_Rob ~ dnorm(0,10),
  theta ~ dcauchy(0,1)
)
mMHC1 <- map2stan(flist, data=MHC.Dataframe.small, iter = 2000, chains = 1)


##Trying to do the hyperpriors -------------



###Working with length of the warmup-----------
{
mW1KI10K <- map2stan(flist, data = work.frame, warmup = 1e3, iter= )
mW2KI10K <- map2stan(flist, data = work.frame, warmup = 2e3, iter= )
mW4KI10K <- map2stan(flist, data = work.frame, warmup = 4e3, iter= )

mWI
mWI
mWI
}

i = 50


results[ i , 1] <- mean(post$b_RobInf)

## "sd_beta" 
results[ i , 2] <- sd(post$b_RobInf)

# "89_hpdi_hi"
results[ i , 3] <- HPDI(post$b_RobInf)[2]

# "89_hpdi_lo" 
results[ i , 4] <- HPDI(post$b_RobInf)[1]

# "log2_mean_beta"
results[ i , 5] <- 2^(mean(post$b_RobInf))

# "log2_sd_beta"
results[ i , 6] <- 2^(sd(post$b_RobInf))

# "log2_89_hpdi_hi"
results[ i , 7] <- 2^(HPDI(post$b_RobInf)[2])

# "log2_89_hpdi_lo"
results[ i , 8] <- 2^(HPDI(post$b_RobInf)[1])
####log2 change values---------
log2value <- 2^post$b_RobInf
mean(log2value)
