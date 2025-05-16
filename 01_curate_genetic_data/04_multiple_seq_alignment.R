# ---------------- Performing multiple sequence alignments ----------------

# This scripts contains the test code to align sequences using the MAFFT package (see sections below).  

# MAFFT was found to perform best (fewest gaps, fastest alignment), and is thus advised for future use. Following the initial alignment problematic sequences were filtered out, with the alignment subsequently being rerun. 

## Libraries required:    -----

library(stringr)
library(dplyr)
library(phylotools)
library(msa)
library(ape)
library(seqinr)
library(DECIPHER)
library(Biostrings)

## --------------------------------------------------------------------------
##  Checking the MAFFT alignment -------------------------------------------- 

# Alignment completed using MAFFT of the command line: 
#  To run in terminal:  
{'mafft --thread n <input_file.fasta> <output_file.fasta>'}

# Or to use a reference genome (this will help improve efficiency and accuracy)
{'mafft --6merpair --thread n --keeplength input_file.fasta ref.fasta > output.fasta '}

# run on 05.09.2023
{'mafft --6merpair --thread 6  --addfragments  
  data/fasta/complete_output_09.2023.fasta
  data/fasta/human_mtDNA_reference_seq.fasta 
  > data/fasta/aln_mafft_incRef_2_09_2023.fasta'}

# read MAFFT alignment into R
mafft_algn <- ape::read.dna(
  file = "data/fasta/aln_mafft_incRef_2_09_2023.fasta",
  format = "fasta"
)
# Check alignment
ape::checkAlignment(mafft_algn)

##  Filtering the MAFFT alignment ------------------------------------------ 

# The resulting alignment has a number of sites where there is a single base 
# called (with all other sequences thus having a gap at this point).  To improve
# the alignment, these problem sequences should be identified and removed. 

# This filtering process should be repeated until there is stability in the 
# sequences.

# round 1 input: "data/fasta/fasta/aln_mafft_incRef_2_09_2023.fasta"
# round 2 input: "data/fasta/aln_mafft_incRef_filtered_09_2023.fasta"
# round 3 imput: "data/fasta/aln_mafft_incRef_filtered_2_09_2023.fasta"

mafft_algn <- Biostrings::readDNAStringSet(
  file = "data/fasta/aln_mafft_incRef_2_09_2023.fasta"
)

# Convert to a matrix 
algn_matrix <-  as.matrix(mafft_algn)

# Calculate the number of sequences and alignment length
num_sequences <- nrow(algn_matrix)
alignment_length <- ncol(algn_matrix)

# Create a matrix to store gap information for each position
gap_counts <- matrix(
  0,
  nrow = alignment_length, 
  ncol = 1)

# Loop through each position in the alignment:

for (position in 1:alignment_length) {
  # Count the number of gaps at the current position
  gap_counts[position, 1] <- sum(
    algn_matrix[, position] == "-"
    )
}

mac <- 10
threshold <- 1-(mac/num_sequences) # % cutoff
threshold <- threshold * num_sequences # number of  gaps limit

# list values where there are few
gap_points <- which(
  gap_counts >= threshold, 
  #arr.ind = TRUE
)

outlier_sequences <- rownames(algn_matrix)[which(
  rowSums(
    algn_matrix[, c(gap_points)] == '-'
  ) != 2 # this number should equal the number of positions in gap_points
)]

outlier_sequences <- algn_matrix[, c(gap_points)]
outlier_sequences <- rowSums(outlier_sequences == '-')
outlier_sequences <- which(outlier_sequences != 2) # this should too ^

outlier_sequences <- as.vector(names(outlier_sequences))

## If you are NOT on the first round, skip to below...
##
# reading in the original fasta file...(if on the first round...)
og_fasta <- seqinr::read.fasta(
  file = "data/fasta/complete_output_09.2023.fasta"
)

# filter out the outlier sequences defined above:
filtered_fasta <- og_fasta[-which(names(og_fasta) %in% outlier_sequences)]

# Write the filtered FASTA file
seqinr::write.fasta(
  sequences = filtered_fasta, 
  names = names(filtered_fasta), 
  file = "data/fasta/complete_output_filtered_09.2023.fasta"
)  # then re-align this file... (using MAFFT)

#  run on 14.10.2023 -- using the filtered fasta files (and standard parameters)
{'mafft --6merpair --thread 6  --addfragments  
  data/fasta/complete_output_filtered_09.2023.fasta
  data/fasta/human_mtDNA_reference_seq.fasta 
  > data/fasta/aln_mafft_incRef_filtered_09_2023.fasta'}

## Skip to here ....

# reading in the latest fasta file...
og_fasta <- seqinr::read.fasta(
  file = "data/fasta/complete_output_filtered_2_09.2023.fasta"
)

# filter out the outlier sequences defined above:
filtered_fasta <- og_fasta[-which(names(og_fasta) %in% outlier_sequences)]

# Write the filtered FASTA file
seqinr::write.fasta(
  sequences = filtered_fasta, 
  names = names(filtered_fasta), 
  file = "data/fasta/complete_output_filtered_3_09.2023.fasta"
)  # then re-align this file... (using MAFFT)

