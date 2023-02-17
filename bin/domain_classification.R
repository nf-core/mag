#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)

parser <- OptionParser()
parser <- add_option(parser, c("-t", "--classification_file"),
            action = "store",
            type = "character",
            metavar = "character",
            help = "The out.txt tsv file of per-contig classifications from Tiara.")
parser <- add_option(parser, c("-s", "--contig_to_bin"),
            action = "store",
            type = "character",
            metavar = "character",
            help = "A tsv file with two columns, bin and contig, listing the contig membership for each bin.")
parser <- add_option(parser, c("-j", "--join_prokaryotes"),
            action = "store_true",
            type = "logical",
            default = TRUE,
            metavar = "logical",
            help = "Use an general prokaryote classification instead of separating Archaea and Bacteria.")
parser <- add_option(parser, c("-a", "--assembler"),
            action = "store",
            type = "character",
            metavar = "character",
            help = "Assembler used to assemble the contigs. 'MEGAHIT' or 'SPAdes' only.")
parser <- add_option(parser, c("-o", "--output_prefix"),
            action = "store",
            type = "character",
            metavar = "character",
            help = "Prefix for the output classification table name.")
args <- parse_args(parser)

## optparse doesn't have a required flag so exit if we don't get given a file
if(is.null(args$classification_file)) {
    stop("Tiara classification file not provided.")
}
if(is.null(args$contig_to_bin)) {
    stop("Contig to bin file not provided.")
}
if(is.null(args$assembler)) {
    stop("Assembler not provided.")
}
if(!(args$assembler %in% c("MEGAHIT", "SPAdes"))) {
    stop("Invalid assembler provided.")
}

classifications <- read_tsv(args$classification_file, na = c("NA", "n/a"))

## DASTOOL_FASTATOCONTIG2BIN drops everything after a space in contig names, which MEGAHIT uses
## Modify Tiara classifications contig IDs accordingly so we can merge
if(args$assembler == "MEGAHIT"){
    classifications$sequence_id <- str_split(classifications$sequence_id, " ", simplify = TRUE)[,1]
}

contig_to_bin <- read_tsv(args$contig_to_bin, col_names = c("sequence_id", "BinID"))

## combination of left_join and filter collectively eliminate unclassified contigs
classifications <- classifications |>
    left_join(contig_to_bin) |>
    filter(!is.na(BinID)) |>
    select(sequence_id,
            BinID,
            Archaea = arc,
            Bacteria = bac,
            Eukarya = euk,
            Organelle = org,
            Unknown = unk1)

if(args$join_prokaryotes) {
    classifications <- classifications |>
        mutate(Prokarya = Archaea + Bacteria) |>
        select(sequence_id, BinID, Prokarya, Eukarya, Organelle, Unknown)
    n_classifications <- 4
} else {
    n_classifications <- 5
}

find_classification <- function(softmax_values, join_prokaryotes = TRUE) {
    if(join_prokaryotes) {
        classifications <- c("prokarya", "eukarya", "organelle", "unknown")
    } else {
        classifications <- c("archaea", "bacteria", "eukarya", "organelle", "unknown")
    }
    return(classifications[which.max(softmax_values)])
}

## Identify the columns to softmax
prob_columns <- 2:(2 + n_classifications - 1)

## Calculate softmax probabilites based on summed bin probabilities for each category
softmax_probabilities <- classifications |>
    group_by(BinID) |>
    summarise(across(all_of(prob_columns), sum), .groups = "drop") |>
    rowwise() |>
    mutate(denominator = sum(exp(c_across(all_of(prob_columns))))) |>
    mutate(across(all_of(prob_columns), \(x) exp(x)/denominator),
            classification = find_classification(c_across(all_of(prob_columns)),
                                join_prokaryotes = args$join_prokaryotes)) |>
    select(-denominator)

## Keep just the classifications so we can loop over more easily
softmax_classifications <- select(softmax_probabilities, BinID, classification)

write_tsv(softmax_probabilities, paste0(args$output_prefix, ".binclassification.tsv"))
write_tsv(softmax_classifications, "bin2classification.tsv")
