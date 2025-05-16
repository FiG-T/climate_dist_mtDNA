# Looking to process extracted GenBank information (from python script)

# Author: FG-T

# Required libraries
library(dplyr)
library(stringi)
library(phylotools)

# Getting country list: 
country_info <- feather::read_feather(
  "data/countries_subcontinents.feather")

# List of countries and subcontinents with associated coordinate info. 
# Cannot be pulled directly from GitHub (not currently possible with feather files)

# create one long regular expression of all countries
country_list <- stringr::str_flatten(
  string = country_info$region,
  collapse = "|"
) 

# add in "old" country names
country_list <- stringr::str_c(
  country_list, 
  "Czechia", "United Kingdom", "Great Britain", "United States" , "Turkije", 
  "Holland",
  sep = "|"
)

## Import files
# Change this path to match your output file.
mito_extracted_1_10 <- read.delim("data/genbank/genbank_output_1-10.tsv")
mito_extracted_10_30 <- read.delim("data/genbank/genbank_output_10-30.tsv")
mito_extracted_30_60 <- read.delim("data/genbank/genbank_output_30-60.tsv")
mito_extracted_60_64 <- read.delim("data/genbank/genbank_output_60-64.tsv")


# create function to processes extracted gb files

extract_processor <- function(
          input
          ) {
  
  output <- paste(
    input, "_processed", 
    sep = ""
  )
  
  output <- input %>%
    filter(
      Organism == "Homo sapiens" # exclude ancient DNA for now
    ) %>%
    mutate(
      Note = stringr::str_extract(
        string  = Note, 
        pattern = country_list
      )
      #. Extract names from Notes column
    ) %>%
    mutate(
      Country = stringr::str_extract(
        string  = Country, 
        pattern = country_list
      )
      #. Extract names from country columns
      #. This will remove regions/ additional geographical info
    ) %>%
    tidyr::unite(
      col = country,
      c(Country, Note), 
      sep = "", 
      na.rm = TRUE
    ) %>% 
    dplyr::select(
      Accession, country, Sequence
    ) %>%
    filter(
      country != ""
    ) %>%
    tidyr::unite(
      col = id,
      c(Accession, country)
    ) %>%
    mutate(Sequence = stringr::str_to_upper(Sequence))
}

# repeat for each of the separate files
mito_extracted_1_10_processed <- extract_processor(mito_extracted_1_10)
mito_extracted_10_30_processed <- extract_processor(mito_extracted_10_30)
mito_extracted_30_60_processed <- extract_processor(mito_extracted_30_60)
mito_extracted_60_64_processed <- extract_processor(mito_extracted_60_64)

# combine into a single file: 
mito_processed <- rbind(
  mito_extracted_1_10_processed,
  mito_extracted_10_30_processed,
  mito_extracted_30_60_processed,
  mito_extracted_60_64_processed
)

mito_processed <- dplyr::distinct(
  mito_processed
)

# convert to fasta
names(mito_processed) <- c("seq.name", "seq.text")

mito_processed$seq.name <- gsub(
  pattern = "UK", 
  x = mito_processed$seq.name, 
  replacement = "United Kingdom"
)

mito_processed$seq.name <- gsub(
  pattern = " ", 
  x = mito_processed$seq.name, 
  replacement = "_"
)

library(stringr)
# list problem letters in alignment 
problem_letter_list <- stringr::str_c(
  "B|D|H|K|M|N|R|S|V|W|Y"
)

# detect and remove samples with problem letters in the sequences 
mito_processed <- mito_processed %>%
  mutate(
    problem.letter = stringr::str_detect(
      string = seq.text, 
      pattern = problem_letter_list
    )
  ) %>%
  filter(
    problem.letter == "FALSE"
  ) %>%
  dplyr::select(
    seq.name, seq.text
  )

  
# make upper case
mito_processed <- mito_processed %>% 
  rowwise() %>%
  mutate(seq.text = stringr::str_to_upper(seq.text))

# write into a fasta file:
phylotools::dat2fasta(
    mito_processed,
    "data/fasta/complete_output_09.2023.fasta"
  )
