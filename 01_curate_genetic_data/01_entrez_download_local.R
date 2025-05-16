# ----------------  Download NCBI files onto local computer ----------------

#  This pulls data directly from NCBI using the "rentrez" package and downloads
#  all the genbank files (in groups of a specified size). 

# set cran mirror to RStudio maintained global server
options(repos = c(CRAN = "https://cloud.r-project.org/"))


#install.packages("rentrez") 
library(rentrez)

#  NCBI  API default is 3 requests per second - this may result in large requests 
#  taking a long time to process.  This can be improved to 10 requests per 
#  second once you register for a (personal) API key (available once you sign up 
#  for an NCBI account). 

set_entrez_key("[your-key-here]")
# set the key for each session
# use your own personal key

# ---------------------------------------------------------------------------------------

# to complete the search term:
mito_search <- entrez_search(
  db = "nucleotide",
  term = "(00000015400[SLEN] : 00000016600[SLEN]) AND Homo[Organism] AND mitochondrion
  [FILT] AND (15400[SLEN] : 17000[SLEN])",
  use_history = TRUE
) 

# This search term defines the length (ensuring only full mtDNA sequences arereturned), the organism, as well as the organelle. Rather than downloading this search result (which will impose an storage-related cutoff), the output is stored directly on the NCBI server (\`use_history = TRUE\`). 

#To download the sequences the resulting search item can be used to fetch the associated GenBank documents.

# check the search has successfully been conducted
mito_search$web_history

# This loops through items in the search term (for the specified range) and appends all the GenBank files into a singular text file. This was done in multiple chunks for practicality reasons (it's easier and quicker to spot if there was an error in smaller groups.

# 4 chuncks used: 
# 1 - 10,000
# 10,0001 - 30,000 
# 30,001 - 60,000
# 60,001 - 64,300

# Using a web history within an "entrez_fetch" loop: 
for( seq_start in seq(
  60001, # start of chunk
  64300, # end of chunk
  100
  )){
  recs <- entrez_fetch(
    db ="nucleotide", 
    web_history = mito_search$web_history,
    #id = mito_search$ids,
    rettype ="gb", 
    retmax = 100,  # download 100 files at a time
    retstart=seq_start)
  Sys.sleep(0.1)          # to ensure NCBI is not overloaded.
  cat(recs, file="/data/genbank/mito_gb_60-64_.txt", append=TRUE)
  cat(seq_start + 99, "GenBank files downloaded\r") # append onto document
}

# This output file can be processed using the python script to extract the country information etc.

# The end result of this process is a group of .txt files that contain a very large number of concatenated GenBank files. The relevant information is extracted from these in the following stages.


