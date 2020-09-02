[![Build Status](https://travis-ci.org/aalbrechtsen/relateAdmix.svg?branch=master)](https://travis-ci.org/aalbrechtsen/relateAdmix)

# ngsRelateAdmix
Estimating relatedness coefficients from admixed populations. It can take called genotypes in PLINK format and genotype likelihoods in a beagle format

The method is implemented in an R package and as a commandline based C++ program embeded in the R package. 

# Install
### R package using devtools

If you have the devtools packages (https://github.com/hadley/devtools) installed in R then you can install the package i R directly from github

```
library(devtools)
install_github("KHanghoj/ngsRelateAdmix")
```

### To compile the C++ version
download the code

```
git clone https://github.com/KHanghoj/ngsRelateAdmix.git
```

go to the scr folder that contains the C++ files 
type 

```
cd ngsRelateAdmix/scr
make -f CPP_Makefile
```
### R package without devtools

If you do not have the devtools package (and dont want to install it) then you will have to build the R package 

first download the code (you need to have a clean version without the compiled c++ code)
```
git clone https://github.com/KHanghoj/ngsRelateAdmix.git
```

```
R CMD build ngsRelateAdmix
```

# getting started

 * The manual can be found on the wiki [http://www.popgen.dk/software/index.php/RelateAdmix]

### in R
```
library(ngsRelateAdmix)
example(relate))
```

### C++ from commandline - examples
After installing the program you can try running it on the example data set in the data folder, which consists of 50 individuals that are admixed from 2 source populations.

If you are in the src folder where you installed relateAdmix and you have the software ADMIXTURE installed this can be done as follows: 
```
cd ../data

# First run Admixture using a plink ".bed" as input to produce population specific allele 
# frequencies (smallPlink.2.P)  and individual ancestry proportions (smallPlink.2.Q).
# (note other programs can be used instead of Admixture, e.g. Structure and FRAPPE)
admixture smallPlink.bed 2 

# Then run relateAdmix with plink bed, bim and fam files plus the Admixture output files as input
../src/relateAdmix -plink smallPlink -f smallPlink.2.P -q smallPlink.2.Q -P 20

```
 

**ADD EXAMPLE FOR GENOTYPE LIKELIHOOD DATA**

