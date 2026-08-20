# MarkerGeneAnalysis16sDADA2

A Nextflow pipeline that processes 16S marker gene sequencing runs into amplicon sequence variants (ASVs) and OTU-level taxonomic assignments using DADA2.

## Overview

This pipeline takes a list of SRA run accessions, downloads the corresponding read data, and runs them through the [DADA2](https://benjjneb.github.io/dada2/) workflow: quality filtering, error-model learning, denoising into amplicon sequence variants, and taxonomic assignment against a reference training set. It is used within VEuPathDB's microbiome/marker-gene data workflows to turn raw 16S sequencing runs into per-run ASV/OTU tables with taxonomic and species-level assignments. Each run accession is downloaded with `fasterq-dump`, quality-filtered, used to build a per-run error model, denoised into a feature table, and finally merged and classified against the supplied training and species-assignment reference FASTAs.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Docker (enabled by default via `nextflow.config`)

The pipeline runs on the `veupathdb/markergeneanalysis16sdada2` container image (which bundles R, DADA2, and the pipeline's R scripts), except for the download step, which uses `veupathdb/bowtiemapping` (providing the SRA toolkit's `fasterq-dump`).

## Usage

```
nextflow run VEuPathDB/MarkerGeneAnalysis16sDADA2 \
  -r main \
  --studyIdFile /path/to/SRAIDS.tsv \
  --trainingSet /path/to/trainingSet.fa \
  --speciesAssignment /path/to/speciesAssignment.fa \
  --outputDir /path/to/output \
  -resume
```

The pipeline has a single named workflow entry point, `markerGeneAnalysis`, which is invoked automatically by the top-level `workflow` block in `main.nf`. For each run accession listed in `studyIdFile`, it downloads the reads, filters them, builds an error model, denoises reads into ASVs, and merges/assigns taxonomy.

## Key Parameters

| Parameter | Description |
| --- | --- |
| `params.studyIdFile` | Path to a TSV file with a `run_accession` column listing the SRA run accessions to process. |
| `params.platform` | Sequencing platform (e.g. `illumina`); passed to the filtering and error-model steps. |
| `params.isPaired` | Whether the reads are paired-end. |
| `params.trimLeft` / `params.trimLeftR` | Bases to trim from the start of the forward/reverse reads during filtering. |
| `params.truncLen` / `params.truncLenR` | Length to truncate the forward/reverse reads to during filtering. |
| `params.maxLen` | Maximum read length allowed during filtering. |
| `params.mergeTechReps` | Whether to merge technical replicates when building the ASV feature table. |
| `params.nValue` | Number of reads dereplicated at a time when building the DADA2 error model (`derepFastq`); lowering this reduces memory usage on large runs. |
| `params.trainingSet` | Path to the reference FASTA used for taxonomic assignment (`assignTaxonomy`). |
| `params.speciesAssignment` | Path to the reference FASTA used for species-level assignment (`addSpecies`). |
| `params.outputDir` | Directory where the final per-run output files are published. |

## Output

For each run accession, three files named `<run_accession>_output`, `<run_accession>_output.bootstraps`, and `<run_accession>_output.full` are published to `params.outputDir`, containing the merged ASV/OTU table with taxonomic assignments, assignment bootstrap confidence values, and the full (unfiltered) assignment output, respectively.
