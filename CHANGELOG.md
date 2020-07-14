# nf-core/mag: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## v1.1.0 - TBD

### `Added`

- Add host read removal with Bowtie 2 and according custom section to MultiQC
- Add separate MultiQC section for FastQC after preprocessing
- Add social preview image
- Compress assembly files
- Add MetaBAT2 RNG seed parameter `--metabat_rng_seed` and set the default to 1 which ensures reproducible binning results
- Add parameters `--megahit_fix_cpu_1`, `--spades_fix_cpus` and `--spadeshybrid_fix_cpus` to ensure reproducible results from assembly tools

### `Fixed`

- Fix links in README
- Fix MetaBAT2 binning discards unbinned contigs [#27](https://github.com/nf-core/mag/issues/27)
- Fix missing MultiQC when `--skip_quast` or `--skip_busco` was specified
- Fix missing channels when `--keep_phix` is specified
- Added missing parameters to summary
- Updated links to minikraken db
- Fixed kraken2 dp preparation: allow different names for compressed archive file and contained folder as for some minikraken dbs
- Fixed channel joining for multiple samples causing MetaBAT2 error [#32](https://github.com/nf-core/mag/issues/32)
- Update MetaBAT2 from v2.13 to v2.15
- Fix number of threads used by MetaBAT2 program `jgi_summarize_bam_contig_depths`
- No more ignoring errors in SPAdes assembly
- No more ignoring of BUSCO errors

### `Deprecated`

- Change depreciated parameters: singleEnd -> single_end, igenomesIgnore -> igenomes_ignore

## v1.0.0 - 2019/12/20

Initial release of nf-core/mag, created with the [nf-core](http://nf-co.re/) template.

### `Added`

- short and long reads QC (fastp, porechop, filtlong, fastqc)
- Lambda and PhiX detection and filtering (bowtie2, nanolyse)
- Taxonomic classification of reads (centrifuge, kraken2)
- Short read and hybrid assembly (megahit, metaspades)
- metagenome binning (metabat2)
- QC of bins (busco, quast)
- annotation (cat/bat)
