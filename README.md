# Infering population sizes with boundary and recurrent mutations 

In this tutorial we learn how to use Bayesian estimators of the effective population size with a boundary and recurrent mutation models. 
These models differ in the fact the BM model only permits mutations in monomorphic or fixed states, neglecting their contribution to frequency shifts between polymorphic states. 
The RM model is thus more appropriate for species with higher levels of diversity like the generality of procariotes. 
We note that in order to assess the diversity levels of a species Watersons theta can be computed. A theta higher tahn 0.05 is indication that a recurrent model should be used.

If you want to know more about the mechanics of these models, our should read the first section of:
* Borges (2024) Evalutating the rarity of mutations in population genetics inferences biorXiv

I order to start 


## The allelic counts 
Image that you inspect the allelic content (alleles 0 and 1) of a sampled of 10 inidividuals in a particular genomic position. You can see that the frequency of a certain allele may vary between 0 and $M$, these extremes representing when none or all the sampled individuals have it. 
If we further sampled S genomic positions and count the number of positions in which the allele 1 appears 0, 1, ..., M times among the sampled individuals, you have the empirical sampled frequency spectrum. An visuall example follows. 

Imagine that I sample 12 genomic positions in a four individuals with the following (very realistic) sequences 

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

 
