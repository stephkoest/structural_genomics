#!/bin/python
"""
Python script that takes a3m alignment file with single line protein sequence.
Outputs a identity filtered a3m file without redunda occurrance of query seqs.
"""

import sys

lines = []
fastadict = {}
i = 0
query = ""

with open(sys.argv[1]) as f:
    lines = f.readlines()
    cllines = [s.strip() for s in lines]
    query = cllines[i]
    print(query)
    print(cllines[2])
    for line in cllines:
        if line.startswith(">"):
            if not cllines[i + 1].startswith(">"):
                if not query in cllines[i]:
                    fastadict[cllines[i]] = [cllines[i].split()[2], cllines[i+1]]
        i+=1

odfastadict = dict(sorted(fastadict.items(), key=lambda item: item[1][0]))

sorted_dict = dict(sorted(fastadict.items(), key=lambda kv: float(kv[1][0]), reverse=True))

for k in sorted_dict.items():
    print(k[0])
    print(k[1][1])
