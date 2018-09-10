###Preparation work||| checked on 8-13 and it's all good ---------------
#Load the data from Dan
load("~/R/HHMI EXROP 18 R Folder/Lohman_RNA-Seq_Workspace/Heritability_Input.RData")
#I have here the 2 data frames: (1)rnaDesign.filt (2)counts.mean.filt

#Load up rethinking package
library(rethinking)

#make the count data in the same design as the rnaDesign.filt
transformed_count <- as.data.frame(t(counts.mean.filt))

#make new dataframe that will have all the information I will add
test.frame <- rnaDesign.filt

#take only the gene I want here it is CXCL19
test.frame$counts <- as.integer(transformed_count$ENSGACG00000019078)

#finding the theoretical count maximum
n_1 <- rowSums(transformed_count) #get the sum for each of the rows for the transformed_count
test.frame$n_max <- as.integer(n_1) #place these values into the new variable n_max

#create dummy variable for Roberts
test.frame$Rob <- as.integer(ifelse(test.frame$Population == "Roberts", 1, 0))

#create dummy variable for Infected
test.frame$infected_id <- as.integer(ifelse(test.frame$Infected == "YES", 1, 0))

#create dummy variable for male
test.frame$male <- as.integer(ifelse(test.frame$Sex == "M", 1, 0))

#create dummy variable family ID
test.frame$family_id <- as.integer(test.frame$Family)

#Trying to make an dummy variable for batch...
test.frame$batch_id <- as.integer(test.frame$Batch)

#create an ID system for the fish
test.frame$fish_id <- as.integer(1:nrow(test.frame))

#check the structure of the dataframe
str(test.frame)


####Model----------------

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
mInt1 <- map2stan(flist, data=work.frame, iter = 2000, chains = 1)


#### Cleaning the Model and Looking Around ------------

#using 3 chains to check whether or not they converge to a similar distribution
mInt2 <- map2stan(flist, data=work.frame, iter = 1500, chains = 3, cores = 3) #chains=3, distributed across 3 cores to alleviate time

#view the tabular output
precis(mInt2)
#The n_eff is a little low. Which may mean that the chain is just inefficient... 
  #The Rhat is generally about 1 so most of the chains did converge

#view the trace plots for the chains...
plot(mInt2)
#here they all seem to converge to a similar point. So the chains seem to be doing okay

#trying more chains and longer iter
mInt3 <- map2stan(flist, data=work.frame, iter = 4000, chains = 4, cores = 4) #chains=3, distributed across 3 cores to alleviate time
#my n_eff obv was longer with an iter of 4000 but the proportion of the n_eff to the iter was higher
  ## before for theta we had 750/1500 and then it went to 2859/4000 (.5-.714)



##POST-CHECK
post.mInt1 <- extract.samples(mInt1)
str(post.mInt1)
dens(post.mInt1$b_RobInf)
HPDI(post.mInt1$b_RobInf)

###read out of the densities for the interaction HPDI
#|0.89      0.89| 
#  -0.5770093  0.8578638
# readout for a random binom with these same st.dev and mean
#|0.89      0.89| 
#  -0.5677317  0.8736436 


#### checking the 3 different models with amounts of cores playing a role within each other
model.compare <- compare(mInt1, mInt2, mInt3) #building the WAIC model comparison machinery
model.compare #view the weights of each item

#        WAIC pWAIC dWAIC weight   SE  dSE
#mInt1 1636.5   0.5   0.0   0.52 8.66   NA
#mInt3 1636.8   0.5   0.3   0.44 8.66 0.06
#mInt2 1641.9   0.6   5.5   0.03 8.64 0.21

## model weights are the same basically for models 1 and 3 so I really don't need that many cores running this data

plot(model.compare, SE = TRUE, dSE= TRUE)
coeftab(mInt1, mInt2, mInt3)

#                mInt1   mInt2   mInt3 
#a              -4.54   -4.44   -4.55
#b_rob          -0.07   -0.01   -0.07
#b_Inf          -0.13   -0.17   -0.13
#b_RobInf        0.12    0.06    0.14
#theta          22.18   20.07   22.30
#nobs              98      98      98

### all the coefficients are basically the same so there is not much loss for the information using all that many cores

# Messing around with priors-------------------


# grabbing the posterior check information to be used ----------------
post <- extract.samples(mInt1) #making a group of 10K samples to be used as our posterior distribution
mean(post$b_RobInf) #output is 0.116941 which is the log 2 of 1.08 which is very close to the -1.17 fold change we saw
sd(post$b_RobInf) #.4346966
dens(post$b_RobInf, col="green")
?ab

?rbetabinom

dnorm()
?dnorm
dens(rnorm(1e4, mean = .1, sd = .4), col= "red")
ablines(dens(post$b_RobInf, col="green"))

dens(rnorm(1e4, mean = .1, sd = .4), col= "red", )
par(new=TRUE)
dens(post$b_RobInf, col="green", norm.comp = TRUE, show.HPDI = .50
     )


#### Dan Phone Call Information 8/11/2018-----------------------

### Model will be 

counts ~ dbetabinom(n_max, pbar, theta)
logit(pbar) <- A + br*rob + bi*inf + b_ri * rob * inf 
A <- a + a_batch[batch] + a_family[family]
a_batch[batch] ~ dnorm(0, sigma_batch)
a_family[family] ~ dnorm(0, sigma_family)
c(sigma_batch, sigma_family) ~ dcauchy(0,1)
a ~ dnorm(0,10)
c(br,bi,b_ri) ~ dnorm(0,10)
theta ~ dexp(1?)

###Then use a forloop

for(i in n:n_genes) {
  #grab the gene of interest in some meaningful way...
  #then take the gene and push it though your analysis!
  #make a matrix of the size of your output 
  results.output <- matrix(nrow = n_genes,ncol = X)
  #X is the number of output variables I want (the mu's and sigma's of the posterior, point estimate and such)
}

### Model as a list 



flist <- alist(
  #likelihood
  counts ~ dbetabinom(n_max, pbar, theta)
  
  #linear model
  logit(pbar) <- A + b_rob * rob + b_Inf*infected + b_RobInf*Rob*infected,
  A <- a + a_batch[batch_id] + a_family[family_id],
  
  #adaptive prior
  a_batch[batch_id] ~ dnorm(0, sigma_batch), 
  a_family[family] ~ dnorm(0, sigma_family),
  
  #hyperpriors
  sigma_batch ~ dcauchy(0,1),
  sigma_family ~ dcauchy(0,1),
  
  #fixed priors
  a ~ dnorm(0,10),
  b_Rob ~ dnorm(0,10),
  b_Inf ~ dnorm(0,10),
  b_RobInf ~ dnorm(0,10),
  theta ~ dexp(1)
)

