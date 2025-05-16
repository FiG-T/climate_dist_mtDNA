# Climate associated natural selection in the human mitochondrial genome. 

Finley Grover-Thomas $^1$  , Lucy van Dorp $^1$, Francois Balloux $^1$ , Aida Moran Andrés $^1$,  Maria Florencia Camus $^{1*}$

$^1$ : Department of Genetics, Evolution and Environment, University College London, London, UK 

$*$ Corresponding author: 
f.camus@ucl.ac.uk; ORCID: https://orcid.org/0000-0003-0626-6865 

This repository contains a copy of the scripts required to perform Generalised Linear Mixed Models (GLMMs) and TreeWAS analysis, in addition to initial data processing and curation. 

## Summary 
Each of the 6 subdirectories contains a set of scripts to perform one suite of actions. Broadly speaking, these actions are...  

* **01_curate_genetic_data** 
  1. Download full mitochondrial DNA sequences from NCBI
  2. Extract sequence and location data from these files
  3. Perform a multiple sequence alignment 
  4. Convert this to vcf and genotype matrix formats
  5. Extract Principal Components from the genotype matrix
  6. Create a tree using FastTree2 on the multiple sequence alignment. 

* **02_curate_metadata** 
  1. Extract climate variables from WORLDCLIM and coordinates from the worldmap. 
  2. Extract paleoclimate variables using the *pastclim* package. 
  3. Extract measures of plant diversity and numbers of outbreaks. 
  4. Calculate mean values per country for each of the extracted variables. 
  5. Download mtDNA positions from the MITOMAP database

* **03_LatPrecip_GLMM_analysis**
  1. Merge genotype and environmental data. 
  2. Perform Generalised Linear Mixed Model analysis to compare models with *latitude* and *annual precipitation* to models with shared ancestry alone. 
  3. Cluster highly correlated SNPs together and plot the results
  4. Repeat the above analyis without sequences from Australia, New Zealand, and the Americas. 

* **04_PaleoBio_GLMM_analysis** 
  1. Merge genotype and paleobio data 
  2. Perform GLMM analysis with the additional paleobio climatic variables
  3. Plot these results
  4. Perform additional GLMMs checks using environmental PC1

* **05_TreeWAS_analysis** 
  1. Prepare data for TreeWAS analysis 
  2. Run TreeWAS and plot results

* **06_supplementary_figures**
  1. Create plots comparing allele frequencies and paleobio climate variables
  2. Assess correlations between environmental variables 
  3. Plot trees showing mutation distributions
  4. Plot PCAs

## Comments

It is advised that the scripts are run in the order listed.  Most scripts should call in at the start all required variables, however if this does not work please run the previous script. 

The majority of the code is written in R, however please note that some code chunks (particuarly within the 01_curate_genetic_data directory) are written in python or bash. This is indicated at the start of the chunk. In such cases please ensure that the code is able to call the required programs (e.g.: Fasttree2 or MAFFT), these may need to first be downloaded locally. 
