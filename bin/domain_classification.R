#!/usr/bin/env Rscript

## Written by Jim Downie and released under the MIT license.
## See git repository (https://github.com/nf-core/mag) for full license text.

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

find_classification <- function(probabilities, join_prokaryotes = TRUE) {
    if(join_prokaryotes) {
        classifications <- c("prokarya", "eukarya", "organelle", "unknown")
    } else {
        classifications <- c("archaea", "bacteria", "eukarya", "organelle", "unknown")
    }
    return(classifications[which.max(probabilities)])
}

classify_bins <- function(tiara, contig2bin, join_prokaryotes){
    ## Some assemblers produce contigs with spaces in the name
    ## Depending on the binner, everything after the first space is sometimes dropped
    ## Make sure that we drop everything after a possible space before doing anything else to allow merging
    tiara$sequence_id <- word(tiara$sequence_id)
    contig2bin$sequence_id <- word(contig2bin$sequence_id)

    if(join_prokaryotes) {
        n_classifications <- 4
    } else {
        n_classifications <- 5
    }

    ## combination of left_join and filter collectively eliminate unclassified contigs
    tiara <- tiara |>
        left_join(contig2bin) |>
        filter(!is.na(BinID)) |>
        select(sequence_id,
                BinID,
                Archaea = arc,
                Bacteria = bac,
                Eukarya = euk,
                Organelle = org,
                Unknown = unk1)

    if(join_prokaryotes) {
        tiara <- tiara |>
            mutate(Prokarya = Archaea + Bacteria) |>
            select(sequence_id, BinID, Prokarya, Eukarya, Organelle, Unknown)
    }

    ## Identify the columns to softmax
    prob_columns <- 2:(2 + n_classifications - 1)

    ## Calculate softmax probabilites based on summed bin probabilities for each category
    softmax_probabilities <- tiara |>
        group_by(BinID) |>
        summarise(across(all_of(prob_columns), sum), .groups = "drop") |>
        rowwise() |>
        mutate(denominator = sum(exp(c_across(all_of(prob_columns))))) |>
        mutate(across(all_of(prob_columns), \(x) exp(x)/denominator),
                classification = find_classification(c_across(all_of(prob_columns)),
                                                join_prokaryotes = join_prokaryotes)) |>
        select(-denominator)

    ## A bin may have no classified contigs if all contigs are below the minimum
    ## Tiara length threshold
    all_bins <- unique(contig2bin$BinID)
    unclassified_bins <- all_bins[!(all_bins %in% softmax_probabilities$BinID)]

    ## Assign these as unclassified
    if(length(unclassified_bins) > 0) {
        if(join_prokaryotes == TRUE){
            unclassified_bins_tbl <- tibble(
                BinID = unclassified_bins,
                Prokarya = NA,
                Eukarya = NA,
                Organelle = NA,
                Unknown = NA,
                classification = "unknown"
            )
        } else {
            unclassified_bins_tbl <- tibble(
                BinID = unclassified_bins,
                Bacteria = NA,
                Archaea = NA,
                Eukarya = NA,
                Organelle = NA,
                Unknown = NA,
                classification = "unknown"
            )
        }
        softmax_probabilities <- bind_rows(softmax_probabilities, unclassified_bins_tbl)
    }

    return(softmax_probabilities)
}

classifications <- read_tsv(args$classification_file, na = c("NA", "n/a"))
contig_to_bin <- read_tsv(args$contig_to_bin, col_names = c("sequence_id", "BinID"))

results <- classify_bins(
  tiara = classifications,
  contig2bin = contig_to_bin,
  join_prokaryotes = args$join_prokaryotes
)

## Keep just the classifications so we can loop over more easily
results_basic <- select(results, BinID, classification)

## write outputs
write_tsv(results, paste0(args$output_prefix, ".binclassification.tsv"))
write_tsv(results_basic, "bin2classification.tsv", col_names = FALSE)

## write out package versions
packageVersion("tidyverse") |> as.character() |> writeLines("tidyverse_version.txt")
