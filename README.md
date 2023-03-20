
# NGSremix

NGSremix is a software tool for estimating relatedness coefficients from admixed populations with low depth sequencing data. When working with low depth sequencing data it takes as input genotype likelihoods in [beagle format](http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#Beagle_format). NGSremix also integrates [relateAdmix](https://github.com/aalbrechtsen/relateAdmix), the predecessor to estimate relatedness form admixed populations, that assumed genotypes are known or are accuratelly called. When you have called genotypes you can use them as input in [binary plink format](https://www.cog-genomics.org/plink/1.9/formats#bed).

NGSremix uses known admixture proportions and populaiton allele frequencies to account for the confounding effects of admixture in relatedness inference. To do so it also requires as input a file with admixture proportions (Q) and a file with ancestral population allele frequencies (F). With low depth data Q and F can be estiamted with [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) that takes also as input genotype likelihoods in beagle format. When having called genotypes, Q and F can be estimated with [ADMIXTURE](https://genome.cshlp.org/content/19/9/1655.long) (it would work with the esitamtes from other software such as STRUCTURE). Importantly, the admixture proporitons need to be accurate and reflect actual population structure, otherwise relatedness estimates might not be accurate.


# Install

First clone the code
```bash
git clone https://github.com/KHanghoj/NGSremix.git
```

go to the scr folder that contains the C++ files and type 

```
cd NGSremix/src
make -f CPP_Makefile
```


#  Arguments

To get the arguments available for `NGSremix` type the following:

```bash
./NGSremix

Arguments:
	-plink     [str] name of the binary plink file (excluding the .bed)
	-beagle    [str] name of the gzipped beagle file
	-fname     [str] Ancestral population frequencies
	-qname     [str] Admixture proportions
	-o         [str] name of the output file

Paired and Parental Ancestry Estimation
	-bothanc 1 [int] Estimates paired and parental ancestries for each individual. Note Relatedness will NOT be estimated. [Default -bothanc 0]

Setup:
	-P         [int] Number of threads
	-seed      [uint]
	-tol       [float] Lower tolerance for excluding admixture/paired ancestry estimates [Default: 0.001]
	-select    [Comma separated (1,2,3) and/or range with dash (1-3), 1-based indexes]
	-notcool 1 [int] Disables paired ancestry. [Default -notcool 0] 
Legacy:
	-F 1       if you want to estimate inbreeding
	-autosomeMax 22	 autosome ends with this chromsome


```

# Examples

We include two small datasets as test examples to run the method.

## Genotype likelihoods

The example for genotype likelihoods has 6 individuals and 10000 sites. The example is a subset of individuals shown in figure 3 of the [NGSremix paper](https://academic.oup.com/g3journal/article/11/8/jkab174/6279082#304747454). 


The ancestry proportions and allele frequencies are already estimated with [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix) using 3 source populations with a larger dataset that includes unrelated individuals. The dataset is not included here, but we assume it is called `beagle.gz` and include the NGSadmix commnad for guidance:

```bash

NGSadmix -likes beagle.gz -K 3 -outfiles 3

```

This will generate a file with admixture proporitons Q called `3.qopt` and a file with population allele frequencies F called `3.fopt.gz`. The script `subsample.py` was used to keep only 6 related indivdiuals and subsample to 10000 sites. We can then use NGSremix to estiamte relatedness:


```bash

../src/NGSremix -beagle test.beagle.gz -qname test.3.qopt -fname test.3.fopt.gz -seed 1 -o res -P 1

```


## Called genotypes
The example for called genotypes consists of 50 individuals that are admixed from 2 source populations. 

The example requires that ADMIXTURE software and NGSremix are installed

```
cd data

# First run Admixture using a plink ".bed" as input to produce population specific allele 
# frequencies (smallPlink.2.P)  and individual ancestry proportions (smallPlink.2.Q).
# (note other programs can be used instead of Admixture, e.g. Structure and FRAPPE)
admixture smallPlink.bed 2 

# Then run NGSremix with plink bed, bim and fam files plus the Admixture output files as input
../src/NGSremix -plink smallPlink -f smallPlink.2.P -q smallPlink.2.Q -P 20

```
