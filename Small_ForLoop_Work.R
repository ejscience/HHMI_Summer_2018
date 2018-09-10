load("~/R/HHMI EXROP 18 R Folder/Lohman_RNA-Seq_Workspace/Heritability_Input.RData")
#I have here the 2 data frames: (1)rnaDesign.filt (2)counts.mean.filt

#Load up rethinking package
library(rethinking)

#getting my working data set
design.data <- rnaDesign.filt
gene.data <- counts.mean.filt[1:50,] #small data set
t.gene.data<- as.data.frame(t(gene.data)) #transformed full data set

n_genes <- n_genes <- nrow(gene.data)

results.summary <- matrix(nrow = n_genes , ncol = 8)

#creating ALL the dummy variables

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


work.frame <- as.data.frame(design.data[,c("batch_id","Rob", "family_id", "infected_id", "n_max")])

for(i in 1:n_genes) {
  
  ####Model
  print(i)
  
  model.frame <- work.frame #this has "batch_id","Rob", "family_id", "infected_id"
  
  #I need "counts"
  
  model.frame$counts <- as.integer(gene.data[ i ,])
  
  
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
  model.out <- map2stan(flist, data= model.frame, iter = 1500, chains = 1, cores = 3)
  
  #### Output data
  post <- extract.samples(model.out) #take the samples here
  
  ## "mean beta"
  results.summary[ i , 1] <- mean(post$b_RobInf)
  
  ## "sd_beta" 
  results.summary[ i , 2] <- sd(post$b_RobInf)
  
  # "89_hpdi_hi"
  results.summary[ i , 3] <- HPDI(post$b_RobInf)[2]
  
  # "89_hpdi_lo" 
  results.summary[ i , 4] <- HPDI(post$b_RobInf)[1]
  
  # "log2_mean_beta"
  results.summary[ i , 5] <- 2^(mean(post$b_RobInf))
  
  # "log2_sd_beta"
  results.summary[ i , 6] <- 2^(sd(post$b_RobInf))
  
  # "log2_89_hpdi_hi"
  results.summary[ i , 7] <- 2^(HPDI(post$b_RobInf)[2])
  
  # "log2_89_hpdi_lo"
  results.summary[ i , 8] <- 2^(HPDI(post$b_RobInf)[1])
  
  
  
}

colnames(results.summary) <- c("mean_beta", "sd_beta","hpdi89_hi", "hpdi89_lo", "log2_mean_beta", "log2_sd_beta", "log2_hpdi89_hi", "log2_hpdi89_lo")

write.csv(results.summary, file = "results_summary_genes1-50.csv")


###------ DOING THE WHOLE THING...

design.data <- rnaDesign.filt
gene.data <- counts.mean.filt
t.gene.data<- as.data.frame(t(gene.data)) #transformed full data set

n_genes <- n_genes <- nrow(gene.data)

results.summary.full <- matrix(nrow = n_genes , ncol = 8)

for(i in 1236:n_genes) {
  
  ####Model
  print(i)
  
  model.frame <- work.frame #this has "batch_id","Rob", "family_id", "infected_id"
  
  #I need "counts"
  
  model.frame$counts <- as.integer(gene.data[ i ,])
  
  
  flist <- alist(
    #likelihood
    counts ~ dbetabinom(n_max, pbar, theta),
    
    #linear model
    logit(pbar) <- A + b_rob * Rob + b_Inf*infected_id + b_RobInf*Rob*infected_id,
    A <- a + a_batch[batch_id] + a_family[family_id],
    
    #adaptive prior
    a_batch[batch_id] ~ dnorm(0.01, sigma_batch), 
    a_family[family_id] ~ dnorm(0.01, sigma_family),
    
    #hyperpriors
    sigma_batch ~ dcauchy(0.01,2),
    sigma_family ~ dcauchy(0.01,2),
    
    #fixed priors
    a ~ dnorm(0.01,10),
    b_rob ~ dnorm(0.01,10),
    b_Inf ~ dnorm(0.01,10),
    b_RobInf ~ dnorm(0.01,10),
    theta ~ dexp(1)
  )
  
  #then put the stan model after so that I can read the details of it easier
  model.out <- map2stan(flist, data= model.frame, iter = 1500, chains = 1, cores = 1)
  
  #### Output data
  post <- extract.samples(model.out) #take the samples here
  
  ## "mean beta"
  results.summary.full[ i , 1] <- mean(post$b_RobInf)
  
  ## "sd_beta" 
  results.summary.full[ i , 2] <- sd(post$b_RobInf)
  
  # "89_hpdi_hi"
  results.summary.full[ i , 3] <- HPDI(post$b_RobInf)[2]
  
  # "89_hpdi_lo" 
  results.summary.full[ i , 4] <- HPDI(post$b_RobInf)[1]
  
  # "log2_mean_beta"
  results.summary.full[ i , 5] <- 2^(mean(post$b_RobInf))
  
  # "log2_sd_beta"
  results.summary.full[ i , 6] <- 2^(sd(post$b_RobInf))
  
  # "log2_89_hpdi_hi"
  results.summary.full[ i , 7] <- 2^(HPDI(post$b_RobInf)[2])
  
  # "log2_89_hpdi_lo"
  results.summary.full[ i , 8] <- 2^(HPDI(post$b_RobInf)[1])
  
  
  
}

colnames(results.summary.full) <- c("mean_beta", "sd_beta","hpdi89_hi", "hpdi89_lo", "log2_mean_beta", "log2_sd_beta", "log2_hpdi89_hi", "log2_hpdi89_lo")

write.csv(results.summary.full, file = "results_summary_genes50-9078.csv")
