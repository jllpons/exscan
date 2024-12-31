# exscan.nf - An extended sequence scanner

> [!IMPORTANT]
> This project is still under development.
> After first release, a license and Zenodo DOI will be added.

## Introduction

**`exscan.nf`** is a bioinformatics pipeline designed to scan **DNA or protein** sequences for key features of interest, with a focus on aiding in genomic annotation and comparative genomics studies.

The pipeline is implemented using [Nextflow](https://www.nextflow.io), and
performs the following steps:

1. **Translation of ORFs (if input is DNA)**  
   Uses [`seqkit2`](<https://doi.org/10.1002/imt2.191>) to translate sequences into all possible open reading frames (ORFs).

2. **Profile HMM Search**  
   Queries each translated ORF or raw protein sequence against a profile HMM database using [`hmmscan`](<http://hmmer.org/>).

3. **Perform different operations on each query result**. Operations include: 
    - Filtering the results by e-value, score, and coverage.
    - Selecting the best hit for each ORF, sequence fragment, or full protein
      sequence.
    - Comparing hits with a GFF file to retrain the features intersecting with the hits.
    - Generating FASTA files contaning aligned sequences for each hit, corresponding ORFs, or the original sequence.

All operations are handled via [`python`](<https://www.python.org/>), [`biopython`](<https://biopython.org/>), [`jq`](<https://jqlang.github.io/jq/>), and [`bedtools`](<https://bedtools.readthedocs.io/en/latest/>).

## Usage

You can run the pipeline using:

```bash
nextflow run main.nf --fasta sequences.fasta --hmmdb hmmdb.hmm
```

Or alternatively:

```bash
nextflow run main.nf -params-file param_files/params.yaml
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Pipeline output

#TODO
To see the results of an example test run with a full size dataset refer to

For more details about the output files and reports, please refer to the
[output documentation](docs/output.md).

## Credits

`exscan.nf` was originally written by [Joan LLuis Pons Ramon](<mail>) at the
[Station Biologique de Roscoff](<https://www.sb-roscoff.fr/en/team-algal-genetics>).

We thank the following people for their extensive assistance in the development of this pipeline:
- #TODO

This work was supported by the HORIZON–MSCA-2022-DN program of the European Commission under the Grant Agreement No 101120280.

## License

Still to be decided.

## References

1. Wei Shen, Botond Sipos, and Liuyang Zhao, “SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing,” Imeta 3, no. 3 (June 2024): e191, https://doi.org/10.1002/imt2.191.
2. Sean R. Eddy, “Accelerated Profile HMM Searches,” ed. William R. Pearson, Plos Computational Biology 7, no. 10 (October 2011): e1002195, https://doi.org/10.1371/journal.pcbi.1002195.
3. Peter J. A. Cock et al., “Biopython: Freely Available Python Tools for Computational Molecular Biology and Bioinformatics,” Bioinformatics 25, no. 11 (June 2009): 1422–23, https://doi.org/10.1093/bioinformatics/btp163.
4. Aaron R. Quinlan and Ira M. Hall, “BEDTools: A Flexible Suite of Utilities for Comparing Genomic Features,” Bioinformatics 26, no. 6 (March 2010): 841–42, https://doi.org/10.1093/bioinformatics/btq033.
