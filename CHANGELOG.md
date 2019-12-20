# nf-core/mag: Changelog

## v1.0.0 - 2019/12/20

Initial release of nf-core/mag, created with the [nf-core](http://nf-co.re/) template.
As this release the pipeline will have the following functionailities:

- short and long reads QC (fastp, porechop, filtlong, fastqc)
- Lambda and PhiX detection and filtering (bowtie2, nanolyse)
- Taxonomic classification of reads (centrifuge, kraken2)
- Short read and hybrid assembly (megahit, metaspades)
- metagenome binning (metabat2)
- QC of bins (busco, quast)
- annotation (cat/bat)
