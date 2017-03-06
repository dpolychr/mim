MIM: Maximal Inexact Matches
===

<b>Description</b>: Given two sequences x and y, MIM finds all maximal inexact matches between
x and y.

<b>Installation</b>: To compile MIM, please follow the instructions given in file INSTALL.
```
 MIM <options>
 Standard (Mandatory):
  -a, --alphabet               <str>      'DNA' for nucleotide  sequences  or 'PROT' for protein  sequences.
  -i, --input-file             <str>      MultiFASTA input filename.
  -o, --output-file            <str>      Output filename with maximal inexact matches.
  -l, --min-seq-length         <int>      Minimum length of match.
  -k, --max-error-size         <int>      Maximum error size between matches.
 Optional:
  -M, --longest-inc-matches    <int>      Choose 1 to return all longest increasing maximal inexact matches
                                          or 0 to return all maximal inexact matches. Default: 0.
  -c, --min-cluster-size       <int>      Minimum number of MIM in each cluster when M=1. Default: 5.
  -r, --rev-compliment         <int>      Choose 1 to compute reverse compliment matches or 0 otherwise. Default: 0.
 Number of threads: 
  -T, --threads                <int>      Number of threads to use. Default: 1.
```

<b>License</b>: GNU GPLv3 License; Copyright (C) 2017 Lorraine A.K. Ayad, Chang Liu, and Solon P. Pissis.

