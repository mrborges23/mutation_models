# Infering population sizes with boundary and recurrent mutations 

In this tutorial, we will learn how to use Bayesian estimators for the effective population size with boundary and recurrent mutation models. 
These models differ in that the BM model only permits mutations in monomorphic or fixed states, neglecting their contribution to frequency shifts between polymorphic states. 
The RM model is thus more appropriate for species with higher levels of diversity, such as the generality of prokaryotes. 

If you want to know more about the mechanics of these models, our should read the first section of:
* Borges (2024) Evalutating the rarity of mutations in population genetics inferences biorXiv

To run these estimators, you will need two (or three) things:
* the empirical sampled frequency spectrum or the allelic counts
* the mutation rate
* if you want to consider possible errors in the mutation rate, its standard deviation is also necessary


## The allelic counts 

Imagine that you inspect the allelic content (alleles 0 and 1) of a sample of 10 individuals at a particular genomic position. You can see that the frequency of a certain allele may vary between 0 and $M$, these extremes representing when none or all the sampled individuals have the said allele. If you further sample $S$ genomic positions and count the number of positions in which allele 1 appears in $0$, $1$, $...$ and $M$ of the sampled individuals, you have obtained the empirical sampled site frequency spectrum. 

A more concrete example follows. Imagine that you sample 12 genomic positions in 5 individuals with the following (not very realistic) sequences 

```
Individual 1 |  A A A A A C C C C C C A
Individual 2 |  A A A A C C C C C C A A
Individual 3 |  A A A C C C C C C A A A
Individual 4 |  A A C C C C C C A A A A
Individual 5 |  A C C C C C C A A A A A
```

Taking A=1 and C=0, we can further sum the values in each column and we get ```(5,4,3,2,1,0,0,1,2,3,4,5)```, which is the frequency of the 1 allele per genomic position. All that the sampled site frequency spectrum does is summarize this vector: we can see that the allele was not observed in any of the individuals in 2 positions, the same for when it appered in only one individual, and so forth. The sampled site frequency spectrum is thus  ````(2, 2, 2, 2, 2, 2)````. We note that this is a rather unlikely sampled frequency spectrum; usually, the monomorphic counts (first and last elements) dominate the counts.

## Running the estimators

In order to run the estimators, we first need to load into our R session the script with all the necessary functions. We called this file ````mutation_models_functions.R```, and you can find it in this repository. 
Open a new script file and write:
```{r}
source("mutation_models_functions.R")
```

Now that all the functions are available in our environment, let us run the estimators. 
These employ an MCMC simulator to estimate the effective population size. You can:
* estimate the effective population size with boundary mutations and a fixed mutation rate: ```mcmc_boundary```
* estimate the effective population size with recurrent mutations and a fixed mutation rate: ```mcmc_recurrent```
* estimate the effective population size with boundary mutations and an uncertain mutation rate: ```mcmc_boundary_umu```
* estimate the effective population size with recurrent mutations and an uncertain mutation rate: ```mcmc_recurrent_umu```

Let us now see how these functions work. We will use data from our case study, the allele counts of the third codon position of the bacterium *Pseudomonas fluorescens*.
First, we load the allelic counts in R. You can write down the vector or use the function ```scan()```, which allows you to read in a vector.
```{r}
# allelic counts
counts <- c(182747, 32606, 23601, 37990, 27474, 20903, 29517, 47825, 31377, 44219, 203259)
```
We then set the mutation rate and its standard deviation. The standard deviation is only needed if we want to incorporate uncertainty in the mutation rate during the inferences.
```{r}
# mutation rate
mu  <- 0.00000000009300
sd  <- 0.00000000000660
```
We are now able to run the estimators:
```{r}
# seed
set.seed(1)

# boundary mutation model with fixed mu
mcmc_boundary1 <- mcmc_boundary(I,counts,mu)
mcmc_boundary2 <- mcmc_boundary(I,counts,mu)

# recurrent mutation model with fixed mu
mcmc_recurrent1 <- mcmc_recurrent(I,counts,mu)
mcmc_recurrent2 <- mcmc_recurrent(I,counts,mu)

# boundary mutation model with uncertain mu
mcmc_boundary_u1 <- mcmc_boundary_umu(I,counts,mu,sd/2)
mcmc_boundary_u2 <- mcmc_boundary_umu(I,counts,mu,sd/2)

# recurrent mutation model with uncertain mu
mcmc_recurrent_u1 <- mcmc_recurrent_umu(I,counts,mu,sd/2)
mcmc_recurrent_u2 <- mcmc_recurrent_umu(I,counts,mu,sd/2)
```
Before summarizing the estimated population size, let us first check whether the MCMC chains converged and mixed properly. 
You can use the following code to produce trace plots.

```{r}
par(mfrow=c(2,2))
plot(mcmc_boundary1,main="BM model: mu fixed",xlab="iterations",ylab="N")
points(mcmc_boundary2,col=2)

plot(mcmc_recurrent1,main="RM model: mu fixed",xlab="iterations",ylab="N")
points(mcmc_recurrent2,col=2)

plot(mcmc_boundary_u1,main="BM model: mu uncertain",xlab="iterations",ylab="N")
points(mcmc_boundary_u2,col=2)

plot(mcmc_recurrent_u1,main="RM model:mu uncertain",xlab="iterations",ylab="N")
points(mcmc_recurrent_u2,col=2)
```

<p align="center">
<img src="https://github.com/mrborges23/mutation_models/blob/main/mcmc_effective_population_size.png" alt="drawing" width="700"/>
</p>

As it seems that the chains have converged properly, we are in conditions to estimate the effective population size.
We will summarize the estimates for each model using the median and the 99% credible interval. 
```
quantile(c(mcmc_boundary1,mcmc_boundary2), prob=c(0.01,.5,0.99) )
quantile(c(mcmc_recurrent1,mcmc_recurrent2), prob=c(0.01,.5,0.99) )
quantile(c(mcmc_boundary_u1,mcmc_boundary_u2), prob=c(0.01,.5,0.99) )
quantile(c(mcmc_recurrent_u1,mcmc_recurrent_u2), prob=c(0.01,.5,0.99) )
```

| quantile        |      1%  |    5%    |     99%  |
|-----------------|----------|----------|----------|
| BM mu fixed     | 10.02115 | 10.03159 | 10.04137 | 
| BM mu fixed     | 9.416650 | 9.418904 | 9.421308 |
| BM mu uncertain | 9.998025 | 10.03176 | 10.06866 | 
| BM mu uncertain | 9.384474 | 9.419235 | 9.455736 |

As expected, the medians of each models are similar, but because in some of these models we accounts for uncertatinties in the mutation rate, the credible intervals are larger.
