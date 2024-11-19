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

