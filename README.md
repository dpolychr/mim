MIM: Maximal Inexact Matches
===

<b>Description</b>: Given two sequences x and y, MIM finds all maximal inexact matches between
x and y.

<b>Installation</b>: To compile MIM, please follow the instructions given in file INSTALL.
```
 MIM <options>
 Standard (Mandatory):
  -a, --alphabet              <str>     'DNA' for nucleotide  sequences  or 'PROT' for protein  sequences.
  -i, --input-file            <str>     MultiFASTA input filename.
  -o, --output-file           <str>     Output filename with maximal inexact matches.
  -q, --q-size                <int>     q-gram length.
  -l, --min-seq-length        <int>     Minimum length of match.
  -k, --min-error-size        <dbl>     Minimum error size between matches.
 Number of threads:
  -T, --threads               <int>     Number of threads to use. Default: 1.
```

<b>License</b>: GNU GPLv3 License; Copyright (C) 2017 Lorraine A.K. Ayad, Chang Liu, and Solon P. Pissis.

