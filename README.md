PhenoRank
=====================
Reducing study bias in gene prioritisation through simulation


PhenoRank reference
----------
[Cornish, A.J., David, A. and Sternberg M.J.E. PhenoRank: reducing study bias in gene prioritisation through simulation. Bioinformatics (2018).][2]


Package requirements
----------
- A Linux or Mac operating system.
- At least 1Gb free disk space.
- Python version 2.7.
- Python packages: numpy (tested with version 1.11.0), pandas (tested with version 0.19.0), scipy (tested with version 0.18.1), sklearn (tested with version 0.0).


Package installation
----------
First download the PhenoRank repository.

The `PhenoRank` python package should then be built by navigating to the main directory and entering the following two commands:

```
python setup.py build
python setup.py install
```

PhenoRank can then be run using `run_PhenoRank.py`, as described below. An implementation of the [PRINCE][1] algorithm can also be run using `run_PRINCE.py`.


Running PhenoRank
----------
PhenoRank is run using the `run_PhenoRank.py` script contained in the main directory. To run PhenoRank, it is necessary to specify either 1) an OMIM ID for the query disease or 2) a set of phenotype terms describing the query disease.

Running PhenoRank using an OMIM term specifying the query disease:

`python run_PhenoRank.py -o results.tsv -d OMIM:606070`

Running PhenoRank using a set of phenotype terms describing the query disease:

`python run_PhenoRank.py -o results.tsv -p 'HP:0007354;HP:0002460;HP:0001739'`

PhenoRank can also be run with additional parameters, as described below.


PhenoRank parameters
----------
These following parameters are accepted by `run_PhenoRank.py`:


##### `-o --file_output`
Required. Name of the file to which the results are written.


##### `-d --omim_id`
OMIM term specifying the query disease. Either this or --phenotype_ids should be specified.


##### `-p --phenotype_ids`
Set of phenotype terms describing the query disease. Should be semi-colon separated. Either this or --omim_id should be specified.


##### `-n --nsim`
Number of simulated diseases to use. Default is *1,000*.


##### `-r --r`
Restart probability to use in the RWR algorithm. Default is *0.1*.


##### `-i --ni`
Number of iterations of the RWR algorithm to complete. Default is *20*.


##### `-g --gene_mask`
Can be used in benchmarking PhenoRank. If an Ensembl gene ID is specified here, then the association between this gene and the query disease is removed from the data used by PhenoRank before gene prioritisation. Can only be used when the query disease is specified using an OMIM term.


PhenoRank results
----------
The results file generated by PhenoRank contains these columns:

- GENE: Ensembl gene identifier.
- SCORE UNRANKED UNPROP: Phenotypic-relevance score of each gene (before propagation across the PPI network) for the query disease.
- SCORE UNRANKED PROP: Phenotypic-relevance score of each gene (after propagation across the PPI network) for the query disease.  
- PVALUE: P-value for each gene. Generated by comparing the phenotypic-relevance scores (after propagation across the PPI network) for the query disease against the phenotypic-relevance scores for each simulated disease.


Running PRINCE
----------
PRINCE is run using the `run_PRINCE.py` script contained in the main directory. To run PRINCE, it is necessary to specify an OMIM ID for the query disease.

Running PRINCE using an OMIM term specifying the query disease:

`python run_PRINCE.py -o results.tsv -d OMIM:606070`

PRINCE can also be run with additional parameters, as described below.


PRINCE parameters
----------
These following parameters are accepted by `run_PRINCE.py`:


##### `-o --file_output`
Required. Name of the file to which the results are written.


##### `-d --omim_id`
OMIM term specifying the query disease. Either this or --phenotype_ids should be specified.


##### `-a --a`
Alpaca parameter value to use. Default is *0.5*.


##### `-n --n`
Number of iterations to complete. Default is *20*.


##### `-c --c`
C paramter value to use. Default is *-15*.


##### `-g --gene_mask`
Can be used in benchmarking PRINCE. If an Ensembl gene ID is specified here, then the association between this gene and the query disease is removed from the data used by PRINCE before gene prioritisation. Can only be used when the query disease is specified using an OMIM term.


PRINCE results
----------
The results file generated by PRINCE contains these columns:

- GENE: Ensembl gene identifier.
- SCORE: Gene score 
- Y: Phenotypic-relevance score of each gene (before propagation across the PPI network).

[1]: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000641
[2]: https://doi.org/10.1093/bioinformatics/bty028
