# Colour Palettes 

# This script contains the colour palettes used in R to help keep 
# schemes consistent across scripts. 


# -----  mtDNA colour palette -----
# Colours per mtDNA loci categorisation

mtDNA_palette <- c(
  C_I = "orchid", 
  C_III = "darkorchid2", 
  C_IV = "darkorchid4",
  C_V = "deepskyblue3",
  ATT = "darkorange",
  C_R = "snow4",
  HVS = "firebrick2", CSB = "firebrick3", 
  OH = "darkorange2", mtPH = "darkorange3", mtH = "darkorange4",
  mtTF1= "orange", 
  rRNA = "limegreen", 
  Humanin = "chartreuse3", SHLP = "chartreuse2", 
  tRNA = "chartreuse4", 
  NC = "grey60"
)

# ----- Environmental covariate palette -----
env_var_palette <- c(
  lat = "black",
  precip_yr = "grey30",
  tmp_max = "grey70", tmp_min = "grey50", 
  mean_T_yr_paleo = "navy",
  T_seasonality_paleo = "maroon",
  T_wetQ_paleo = "maroon1", T_dryQ_paleo = "maroon3",
  mean_T_warmQ_paleo = "hotpink3", mean_T_coldQ_paleo = "hotpink4",
  precip_yr_paleo = "turquoise4", 
  precip_dryQ_paleo = "turquoise2", precip_coldQ_paleo = "turquoise",
  npp_paleo = "forestgreen", 
  lai_paleo = "darkolivegreen" , endemic_area = "darkolivegreen3",
  outbreaks_area = "red4"
)

# ----- Continental Colours ----- 

FiGT_continent_palette <- c(
  Africa = "deepskyblue4",
  Asia = "firebrick2", 
  Europe ="orchid3", 
  `North America` ="chartreuse4",
  Oceania = "turquoise",
  `South America` = "chartreuse3"
  
)

# ----- Older env palettes ----- 
meta_var_palette <- c(
  # stage 1 variable
  lat = "firebrick3",
  precip_yr = "deepskyblue3", 
  
  # Thermo & diet var
  T_seasonality_paleo = "darkorchid",
  # Thermo variable
  T_wetQ_paleo = "hotpink2", 
  mean_T_coldQ_paleo = "hotpink4",
  # Diet variables climate
  T_dryQ_paleo = "orchid", 
  precip_dryQ_paleo = "darkorchid4", 
  # Diet plant measures
  lai_paleo = "chartreuse3", 
  endemic_area = "chartreuse4",
  npp_paleo = "palegreen4",
  # Disease varible
  outbreaks_area = "darkorange2"
)

meta_var_palette_2 <- c(
  T_seasonality_paleo = "black",
  # Thermo variable
  T_wetQ_paleo = "turquoise4", 
  mean_T_coldQ_paleo = "turquoise3",
  precip_dryQ_paleo = "turquoise", 
  # Diet plant measures
  endemic_area = "yellow3",
  npp_paleo = "yellow"
)