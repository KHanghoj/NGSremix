../src/NGSremix -beagle test.beagle.gz -qname test.3.qopt -fname test.3.fopt.gz -seed 1 -o res_all -P 1

awk 'BEGIN{OFS="\t"}
     NR==1 {print $1,$2,$3,$4,$5; next}
     {printf "%s\t%s\t%.3f\t%.3f\t%.3f\n", $1, $2, $3, $4, $5}' \
     res_all > res

md5sum -c truth/res.md
