
# NGSremix
Estimating relatedness coefficients from admixed populations. It can take called genotypes in PLINK format and genotype likelihoods in a beagle format

The method is implemented in an R package and as a commandline based C++ program embeded in the R package. 

# Install

First clone the code
```bash
git clone https://github.com/KHanghoj/NGSremix.git
```

go to the scr folder that contains the C++ files and type 

```
cd NGSremix/scr
make -f CPP_Makefile
```


#  Arguments

To get the arguments available for `NGSremix` type the following:

```bash
./NGSremix
```

# Examples

## Genotype likelihoods
The example for genotype likelihoods of 6 individuals and 10000 sites. The ancestry proportions and allele frequencies are already estimated with `ngsadmix` using 3 source populations. The example is a subset of individuals shown in figure 3 of the `NGSremix` paper. 

```bash

cd dataNGS

./run.sh

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
 



**ADD EXAMPLE FOR GENOTYPE LIKELIHOOD DATA ANGSD - NGSADMIX - NGSremix**
