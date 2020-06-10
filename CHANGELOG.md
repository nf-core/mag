# nf-core/mag: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## v1.1.0 - TBD

### `Added`

- Add host read removal with Bowtie 2
- Add separate MultiQC section for FastQC after preprocessing
- Add social preview image

### `Fixed`

- Fix links in README
- Fix MetaBAT2 binning discards unbinned contigs [#27](https://github.com/nf-core/mag/issues/27)
- Fix missing MultiQC when `--skip_quast` or `--skip_busco` was specified
- Fix missing channels when `--keep_phix` is specified

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
