# Design Documentation

## Table of Contents

1. [Overview](#overview)
2. [Component Breakdown](#component-breakdown)
    1. [Processes](#processes)
    2. [Subworkflows](#subworkflows)
    3. [Workflows](#workflows)

## Overview

This document outlines the structure and design choices of the Nextflow pipeline, explaining the organization of components and the rationale behind them. Many of the ideas come from the [nf-core documentation](<https://nf-co.re/docs/usage/getting_started/terminology>).

## Component Breakdown

### Processes

Processes are the most basic building blocks of the pipeline and are designed to be as atomic as possible. They are located in `./modules/local/`.

### Subworkflows

Subworkflows are chains of multiple processes that provide specific functionality. They are used for:

1. Reusability: Facilitating the reuse of common process sequences.
2. Readability: Simplifying the workflow by replacing long lists of processes with concise subworkflows.

Currently, only `domtblop.py` operations are implemented as subworkflows. All subworkflows are located in `./subworkflows/`.


### Workflows

Workflows are end-to-end pipelines that utilize processes and subworkflows to achieve specific goals.

There are two workflows in this project:

1. `./main.nf`: The entry point of the pipeline. It sets up and validates the environment before running the actual pipeline and handles termination.
2. `./workflows/exscan.nf`: Contains the main pipeline where processes and subworkflows are connected to accomplish the pipeline's objectives.

## Parameters

#TODO

All parameters are defined in `./nextflow.config`.

Required parameters are defined with `null` as defalut.
More specific parameters are defined with default values (e.g. params.group =
3750)

Parameters can be provided either as:
- CLI arguments `$ nextflow run main.nf --<parameter> <value>`
- In a `params.yaml` or `params.json` file `$ nextflow run main.nf -params-file
  <params.yaml>`
- In a `params.csv` file `$ nextflow run main.nf --input <params.csv>`.

Notice that if you can split each chr into one fasta, and map them to each
gff/bam file using the `--input <params.csv>` option. This way, each chr will be
processed in parallel.

>[WARNING]
> The `-c` opti #TODO

## Configuration

#TODO

