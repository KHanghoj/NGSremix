#!/usr/bin/env python3
import numpy as np
import pandas as pd

np.random.seed(seed=1)
orig = pd.read_csv("res.cool.txt", sep="\t").query("k1>0.8").head(3)
k = orig[["ind1", "ind2"]].values.flatten()
k.sort()
header = ["marker", "allele1", "allele2"]

rows_zero = []

for x in k:
    header.append(f"Ind{x-1}")
    rows_zero.append(x-1)
    for y in range(1,3):
        header.append(f"Ind{x-1}.{y}")

fs = pd.read_csv("3.fopt.gz", sep=" ", header=None)

s = np.random.choice(fs.shape[0], size=10000, replace=False)
s.sort()

fs.loc[s].to_csv("test.3.fopt.gz", index=False, header=False, sep=" ")

pd.read_csv("beagle.gz", sep="\t", usecols=header).loc[s].to_csv("test.beagle.gz", index=False, header=True, sep="\t")


pd.read_csv("3.qopt", sep=" ", header=None).loc[rows_zero].to_csv("test.3.qopt", index=False, header=False, sep=" ")

pd.read_csv("names.list", header=None).loc[rows_zero].to_csv("test.names.list", index=False, header=False)
