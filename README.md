# exscan

## Introduction

`exscan` is bioinformatics pipeline for the exploration and annotation of exons
and their specific features. It takes a nucleotide `fasta` file to extract all
of the open reading frames (ORFs), uses profile hidden Markov models (hmm) to detect exons displaying specific domains and finally filters the results based on specific criteria.

1. Translate all ORFs of a nucleotide sequence ([`seqkit2`](< https://doi.org/10.1002/imt2.191>)).
2. Query each ORFs against a profile HMM database ([`hmmscan`](<http://hmmer.org/>))
3. Perform different operations on each query results. (([`python`](<https://www.python.org/>), [`biopython`](<https://biopython.org/>)), [`jq`](<https://jqlang.github.io/jq/>), [bedtools](<https://bedtools.readthedocs.io/en/latest/>)).

## Usage

You can run the pipeline using:

```bash
nextflow run main.nf
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Pipeline output

#TODO
To see the results of an example test run with a full size dataset refer to

For more details about the output files and reports, please refer to the
[output documentation](docs/output.md).

## Credits

`exscan` was originally written by [Joan LLuis Pons Ramon](<mail>)

We thank the following people for their extensive assistance in the development of this pipeline:

#TODO
EU funding?

