# relateAdmix
Estimating relatedness coefficients from admixed populations

The method is implemented in an R package and as a commandline based C++ program embeded in the R package. 

# Install
### R package using devtools

If you have the devtools packages (https://github.com/hadley/devtools) installed in R then you can install the package i R directly from github

```
library(devtools)
install_github("aalbrechtsen/relateAdmix")
```

### To compile the C++ version
download the code

```
git clone https://github.com/aalbrechtsen/relateAdmix.git
```

go to the scr folder that contains the C++ files 
type 

```
cd relateAdmix/scr
mv CPP_Makefile Makefile
make
```
### R package without devtools

If you do not have the devtools package (and dont want to install it) then you will have to build the R package 

first download the code (you need to have a clean version without the compiled c++ code)
```
git clone https://github.com/aalbrechtsen/relateAdmix.git
```

```
R CMD build relateAdmix
```

or  build and install

```
R CMD build relate
R CMD INSTALL Relate_<add version number>.tar.gz
```


# getting started

 * The manual can be found on the wiki [http://www.popgen.dk/software/index.php/RelateAdmix]

### in R
```
library(relateAdmix)
example(relate))
```

### C++ from commandline
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
 

