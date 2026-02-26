# GATK SNP calling pipeline

Small Nextflow pipeline containing the GATK scripts of josieparis that perform read QC, alignment (BWA), marking duplicates and SNP calling with GATK HaplotypeCaller and FreeBayes, then compare and refine calls.

**Repository layout**
- `main.nf` — Nextflow workflow entrypoint.
- `modules/` — pipeline process modules (e.g. `gatk.nf`).
- `data/` — example inputs (FASTQ, reference, metadata).
- `work/`, `logs/`, `snp/`, `test1/`, `test2/` — runtime and results directories produced by runs.

**Requirements**
- Java 11+ and Nextflow (tested with Nextflow DSL2).
- BWA, Samtools, Picard, GATK, FreeBayes available either on PATH or via Singularity/containers.
- Singularity images not present in the repository (`*.sif`) but Def files for build if present.
- Sufficient CPU and RAM for alignment and variant calling.

Recommended: use the provided Singularity images (e.g. `gatk_snp.sif`, `samtools.sif`) or configure `process.container` in Nextflow.

Quick checklist before running:
- Prepare `--data` path with paired FASTQ files named like `SAMPLE_1.fastq` and `SAMPLE_2.fastq`.
- Provide `--metadata` file (sample metadata used by read-group and merging steps).
- Provide `--reference` pointing to the reference FASTA (with accompanying index files `.fai`, dict, etc.).

Quick start

1. Run locally (example):

```bash
nextflow run main.nf \
  --data data \
  --metadata data/metadata_2006.tsv \
  --reference data/IHEM_04380_LY1996_0339_1_reference_assembly.fasta \
  --outdir GATK_results \
  --cpus 8
```


Key pipeline parameters (set via CLI)
- `--data` (required): directory containing paired FASTQ files matching `*_{1,2}.fastq`.
- `--metadata` (required): metadata TSV used by read-group/merge steps.
- `--reference` (required): path to reference FASTA (must have `.fai` and dict).
- `--outdir` (default `GATK_results`): pipeline output directory.
- `--cpus` (default `4`): total CPUs; pipeline computes per-task allocation.
- `--merged` (optional): if set, triggers sample merging workflow.
- QC params: `--galore_quality`, `--galore_length`.

Outputs
- Final VCFs and intermediate GVCFs appear under `snp/` and `GATK_results` (depending on run config).
- Per-process logs are written to `work/` and `logs/` folders; see the run log for details (Nextflow `.log`).

Troubleshooting
- Missing reference indices: create `.fai` with `samtools faidx` and dict with Picard `CreateSequenceDictionary`.
- Container issues: ensure Singularity is installed and images are accessible.
- If Nextflow fails on resource assignment, adjust `--cpus` or per-process `cpus` settings in `main.nf`.

Notes
- This pipeline assumes paired-end Illumina reads and uses GATK HaplotypeCaller plus FreeBayes for cross-caller comparison and refinement.
- See `main.nf` for process names and included modules in `modules/gatk.nf`.

Contact
- For questions or issues, open a GitHub issue or contact the repository owner.

License
- No license specified; add a `LICENSE` file if you intend to publish this repository.
