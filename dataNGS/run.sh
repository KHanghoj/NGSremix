../src/NGSremix -beagle test.beagle.gz -qname test.3.qopt -fname test.3.fopt.gz -seed 1 -o res -P 1
md5sum -c truth/res.md
