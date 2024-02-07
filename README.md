# Infering population sizes with boundary and recurrent mutations 

In this tutorial we learn how to use Bayesian estimators of the effective population size with a boundary and recurrent mutation models. 
These models differ in the fact the BM model only permits mutations in monomorphic or fixed states, neglecting their contribution to frequency shifts between polymorphic states. 
The RM model is thus more appropriate for species with higher levels of diversity like the generality of procariotes. 
We note that in order to assess the diversity levels of a species Watersons theta can be computed. A theta higher tahn 0.05 is indication that a recurrent model should be used.

If you want to know more about the mechanics of these models, our should read the first section of:
* Borges (2024) Evalutating the rarity of mutations in population genetics inferences biorXiv

To run these estimator you will need two (or three) things:
* the sampled frequency spectrum or the allelic counts
* the mutation rate
* if you want to consider possible errors in the mutation rate, its standard deviation is also necessary


## The allelic counts 
Image that you inspect the allelic content (alleles 0 and 1) of a sampled of 10 inidividuals in a particular genomic position. You can see that the frequency of a certain allele may vary between 0 and $M$, these extremes representing when none or all the sampled individuals have the said allele. If you further sample $S$ genomic positions and count the number of positions in which the allele 1 appears in $0$, $1$, $...$ and $M$ of sampled individuals, you have obtained the empirical sampled site frequency spectrum. An visual example follows.

Imagine that you sample 12 genomic positions five individuals with the following (very realistic) sequences 

```
Individual 1 |  A A A A A C C C C C C A
Individual 2 |  A A A A C C C C C C A A
Individual 3 |  A A A C C C C C C A A A
Individual 4 |  A A C C C C C C A A A A
Individual 5 |  A C C C C C C A A A A A
```

Taking A=1 and C=0, we can further sum the columns and we get ```(5,4,3,2,1,0,0,1,2,3,4,5)``` which is the frequency of the 1 allele per genomic position. All that the sampled SFS does is summarizing this vector:
Allele frequency 0 1 2 3 4 5
sampled SFS      2 2 2 2 2 2 

We note that this is a rather unlikely sampled frequency spectrum; usually the monomorphic counts dominate.

# Running the estimators

In order to run the estimators we first to load into our R session the script with all the necessary functions. We called this file ````mutation_models_functions.R``` and you can find it in this repository. 
Open a new file script and write:
```{R}
source("mutation_models_functions.R")
```
 
