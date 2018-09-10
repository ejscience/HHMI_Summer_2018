#------------
#Load the data from Dan
load("~/R/HHMI EXROP 18 R Folder/Lohman_RNA-Seq_Workspace/Heritability_Input.RData")
#I have here the 2 data frames: (1)rnaDesign.filt (2)counts.mean.filt

#Load up rethinking package
library(rethinking)

###Making the SMALL forloop---------
small.counts.mean.filt <- counts.mean.filt[1:50,] # took rows 1 to 50 to try out the forloop and transformed it to get the genes

t.small.counts.mean.filt <- as.data.frame(t(small.counts.mean.filt))
n_genes <- nrow(small.counts.mean.filt)
gene_names <- as.list(small.counts.mean.filt[,1])


results.summary <- matrix(nrow = n_genes , ncol = 8, dimnames = list(, c("mean_beta", "sd_beta","89_hpdi_hi", "89_hpdi_lo", "log2_mean_beta", "log2_sd_beta", "log2_89_hpdi_hi", "log2_89_hpdi_lo")))



results.summary <- matrix(nrow = n_genes , ncol = 8, dimnames = list(, c("mean_beta", "sd_beta","89_hpdi_hi", "89_hpdi_lo", "log2_mean_beta", "log2_sd_beta", "log2_89_hpdi_hi", "log2_89_hpdi_lo")))
#This forloop wil be for the actual analysis
for(i in 1:n_genes) {
  
  
  
  ####Model
  
  work.frame <- as.data.frame(test.frame[,c("batch_id","counts","Rob","n_max", "family_id", "infected_id")])
  
  flist <- alist(
    #likelihood
    counts ~ dbetabinom(n_max, pbar, theta),
    
    #linear model
    logit(pbar) <- A + b_rob * Rob + b_Inf*infected_id + b_RobInf*Rob*infected_id,
    A <- a + a_batch[batch_id] + a_family[family_id],
    
    #adaptive prior
    a_batch[batch_id] ~ dnorm(0, sigma_batch), 
    a_family[family_id] ~ dnorm(0, sigma_family),
    
    #hyperpriors
    sigma_batch ~ dcauchy(0,2),
    sigma_family ~ dcauchy(0,2),
    
    #fixed priors
    a ~ dnorm(0,10),
    b_rob ~ dnorm(0,10),
    b_Inf ~ dnorm(0,10),
    b_RobInf ~ dnorm(0,10),
    theta ~ dexp(1)
  )
  
  #then put the stan model after so that I can read the details of it easier
  model.out <- map2stan(flist, data=work.frame, iter = 2000, chains = 1)
  
  #### Output data
  post <- extract.samples(mInt1) #take the samples here
  
  
  
  
  
  
}

###Making the LARGE forloop--------

design.data <- rnaDesign.filt
gene.data <- counts.mean.filt[1:50,] #small data set
t.gene.data<- as.data.frame(t(gene.data)) #transformed full data set

results.matrix <- x

#creating the dummy variables

{
  #create dummy variable for Roberts
  design.data$Rob <- as.integer(ifelse(design.data$Population == "Roberts", 1, 0))
  
  #create dummy variable for Infected
  design.data$infected_id <- as.integer(ifelse(design.data$Infected == "YES", 1, 0))
  
  #create dummy variable family ID
  design.data$family_id <- as.integer(design.data$Family)
  
  #Trying to make an dummy variable for batch...
  design.data$batch_id <- as.integer(design.data$Batch)
  
  #create an ID system for the fish
  design.data$fish_id <- as.integer(1:nrow(design.data))
  
  #finding the theoretical count maximum
  n_1 <- rowSums(t.gene.data) #get the sum for each of the rows for the transformed_count
  design.data$n_max <- as.integer(n_1) #place these values into the new variable n_max
  
}

n_genes <- n_genes <- nrow(gene.data)
work.frame <- work.frame <- as.data.frame(design.data[,c("batch_id","Rob", "family_id", "infected_id", "n_max")])

for(i in 1:n_genes) {
  
  ####Model
  
  model.frame <- work.frame #this has "batch_id","Rob", "family_id", "infected_id"
  
  #I need "counts"
  
  model.frame$counts <- as.integer(gene.data[i,])
  
  
  flist <- alist(
    #likelihood
    counts ~ dbetabinom(n_max, pbar, theta),
    
    #linear model
    logit(pbar) <- A + b_rob * Rob + b_Inf*infected_id + b_RobInf*Rob*infected_id,
    A <- a + a_batch[batch_id] + a_family[family_id],
    
    #adaptive prior
    a_batch[batch_id] ~ dnorm(0, sigma_batch), 
    a_family[family_id] ~ dnorm(0, sigma_family),
    
    #hyperpriors
    sigma_batch ~ dcauchy(0,2),
    sigma_family ~ dcauchy(0,2),
    
    #fixed priors
    a ~ dnorm(0,10),
    b_rob ~ dnorm(0,10),
    b_Inf ~ dnorm(0,10),
    b_RobInf ~ dnorm(0,10),
    theta ~ dexp(1)
  )
  
  #then put the stan model after so that I can read the details of it easier
  model.out <- map2stan(flist, data=model.frame, iter = 2000, chains = 1)
  
  #### Output data
  post <- extract.samples(model.out) #take the samples here