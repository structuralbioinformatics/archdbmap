ArchDBmap
---------

This **python2 only** library creates the executable ``archdbmap`` that can map, through sequence similarity, Smotifs from the ArchDB14 database into the
queries (being PDB or FASTA inputs).

Install with:

```
pip install https://github.com/structuralbioinformatics/archdbmap/archive/v0.2.zip
```

There are 3 files needed to properly execute this script: the FASTA and FASTA.index version of the ArchDB14 and access to the ArchDB14 SQL.
The sequence files and the SQL if you don't have SQL access to the actual ArchDB14) can be downloaded at http://sbi.upf.edu/archdb/download

It also requires a **blast** installation.

Check on the available options of the script with:

```
archdbmap -h
```

If of use, please cite:

Bonet, J., Planas-Iglesias, J., Garcia-Garcia, J., Marin-Lopez, M.A., Fernandez-Fuentes, N., Oliva, B. ArchDB 2014: structural classification of loops in proteins. Nucleic Acids Res. http://doi.org/10.1093/nar/gkt1189
