## Reading in MT loci data 

# taking the csv from the MITOMAP website [accessed Nov.2023], save into filed
mt_loci_pos <- read.csv(
  file = "data/GenomeLoci_MITOMAP_Foswiki_08.2024.csv"
)

# add an extra column
mt_loci_pos$classification <- mt_loci_pos$Shorthand

# tRNA is tricky as the shorthand names are very generic (and thus are not easily exactly matched with a regular expression) - this matches to the 'description' column. 
mt_loci_pos$classification[
  grepl(
    pattern = "tRNA.*",
    mt_loci_pos$Description
  )
] <- "tRNA"

# create a tibble with MITOMAP detailed names and desired replacement general categories 
mt_loci_names <- tidyr::tibble(
  regex_used = c(
    "12S|16S", "ND.*","Cytb" ,"CO.*", "ATP.*", "NC.*", ".*TF.", ".*HVS.*",
    "CR:OH.*" , ".*CS.*", "CR:mt.H", ".*PH.*", "SHLP.", ".*Control.*"
  ), 
  replace_with = c(
    "rRNA", "C_I", "C_III", "C_IV", "C_V", "NC", "mtTF1", "HVS",
    "OH", "CSB", "mtH", "mtPH", "SHLP", "C_R"
  )
)

# create a vector of mt regions (this will be used to define levels)
mt_classification_levels <- c(
  "C_I", "C_III", "C_IV", "C_V", "ATT", "C_R", "HVS", "CSB", "OH", "mtPH", "mtH",
  "mtTF1","rRNA", "Humanin","SHLP", "tRNA", "NC"
)

# for each row replace the complex with the associated general category defined above: 
for (i in 1:nrow(mt_loci_names)) {
  mt_loci_pos$classification <- stringi::stri_replace_all(
    str = mt_loci_pos$classification,  # the column to overwrite
    regex = mt_loci_names$regex_used[i],  # the MITOMAP name 
    replacement = mt_loci_names$replace_with[i] # new general/clean category
  )
}

# duplicate the CR and ATT rows (to control for the fact that they loop around)
mt_loci_origin <- mt_loci_pos[grepl(
  pattern = "C_R|ATT|.*7S.*", 
  mt_loci_pos$classification
),]
mt_loci_origin$Starting <- 1
#mt_loci_origin$aln_start <- 1

# replace ending with end of 'loop'
mt_loci_pos$Ending[mt_loci_pos$Ending %in% c(499, 576, 191)] <- 16569

mt_loci_pos <- dplyr::bind_rows(
  mt_loci_pos, 
  mt_loci_origin
)

# calaculate the length of each region
mt_loci_pos$length <- (mt_loci_pos$Ending-mt_loci_pos$Starting) + 1

# calculate the total length per classification
mt_loci_length <- mt_loci_pos %>%
  group_by(classification) %>%
  reframe(classification_length = sum(length))

# add in classification length
mt_loci_pos <- dplyr::left_join(
  x = mt_loci_pos, 
  y = mt_loci_length, 
  by = "classification"
)

# reorder columns
mt_loci_pos <- mt_loci_pos %>%
  select(
    Map.Locus, Starting, Ending, length, Shorthand, 
    Description, classification, classification_length
  )
# convert positions to numeric
mt_loci_pos$Starting <- as.numeric(mt_loci_pos$Starting)
mt_loci_pos$Ending <- as.numeric(mt_loci_pos$Ending)

feather::write_feather(
  x = mt_loci_pos, 
  path = "data/feather/mtDNA_loci_all_positions_classifications.feather"
)
