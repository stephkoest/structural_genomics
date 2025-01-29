# structural_genomics

Collection of scripts used for analyses in "Structure-based inference of eukaryotic complexity in Asgard archaea."


##

### Necessary conda environments

```bash
conda install -n mmseqs2_v14 -c bioconda mmseqs2::14.7e284
```
### other software

mmseqs2 version for colabfold search: fd1837b600c57278bcfb2ac1ac7f024e458c0606

## Genetic search including a user-provided database

### mmseqs_search.sh
Adapted from the ColabFold genetic search shell script to include
an external database.
