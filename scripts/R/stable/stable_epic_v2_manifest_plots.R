
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                Swifthoof:: 
#               Plot EPIC v2 Data for Validation/Comercial Sake
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# General Package Install:
#
# packages_vec <- c( "doParallel", "tidyverse", "optparse", "plyr",
#                    "glue", "matrixStats", "scales", "stringr", "Rcpp", "R6",
#                    "ggplot2", "GGally" )
# install.packages( pkgs = packages_vec, 
#                   dependencies = TRUE )
#
# Bioconductor Install::
#

rm(list=ls(all=TRUE))

# Options, Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("doParallel", quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("sesame", quietly = TRUE) ) )

# test_idats_dir <- "/Users/bretbarnes/Documents/data/idats/idats_AKE-MVP-Failed-v1"
# open_ses_dat <- sesame::openSesame( sesame::searchIDATprefixes( test_idats_dir, recursive = TRUE ), platform = "EPIC" )
# open_ses_ssets <- sesame::openSesame( sesame::searchIDATprefixes( test_idats_dir, recursive = TRUE ), platform = "EPIC" )
# open_ses_dat %>% as.data.frame() %>% tibble::as_tibble( rownames = "Probe_ID" ) %>% dplyr::summarise( dplyr::across( -Probe_ID , ~ sum(is.na(.x)) /865918 ) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Set Run Environment:: RStudio/Command-Line
#                            Source All Functions
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par  <- list()
args <- commandArgs(trailingOnly = FALSE)

par$src_path <- NULL
par$run_mode <- args[1]
par$date_str <- Sys.Date() %>% as.character()
par$prgm_dir <- 'swifthoof'
par$prgm_tag <- 'stable_epic_v2_manifest_plots'
par$verbose  <- 3
local_paths  <- c( 
  "/Users/bbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
  "/Users/bretbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
  "/illumina/scratch/darkmatter/tools/Workhorse-Unstained/scripts/R",
  "/repo/Workhorse-Unstained/scripts/R" )

if ( par$run_mode == "RStudio" ) {
  for (path in local_paths) if ( dir.exists(path) ) par$src_path <- path
} else {
  par$exe_path <- 
    base::normalizePath( base::substring( args[grep("--file=", args)], 8 ) )
  par$scr_path <- base::dirname( par$exe_path )
  par$src_path <- base::dirname( par$scr_path )
  par$run_mode <- 'Command_Line'
}
stopifnot( length(par$src_path) > 0 )
stopifnot( dir.exists(par$src_path) )

par$fun_path <- file.path( par$src_path, "functions" )
stopifnot( dir.exists(par$fun_path) )

src_ret <- base::lapply( base::list.files( 
  par$fun_path, pattern = "\\.R$", full.names = TRUE ), source )
par <- source_functions( pars = par, vb = par$verbose )
par <- params_check( pars = par, args = args, 
                     prgm_aux_check = FALSE, vb = par$verbose )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Get Program Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par$run_name <- "EPIC_v2"
par$run_name <- "GSIBIOINFO-597"
par$run_name <- "Embark_v3"
par$run_name <- "GSIBIOINFO-638"

# Multiple copies of v0 (bk and bk24)
par$version <- 0
par$version <- 1
# par$version <- 2
# par$version <- 3

opt <- NULL
opt <- swifthoof_sesame_options( pars = par, args = args, vb = par$verbose )
vb  <- opt$verbose
vt  <- 0
tc  <- 0

p0  <- vb > vt + 0
p1  <- vb > vt + 1
p2  <- vb > vt + 2
p4  <- vb > vt + 4
p8  <- vb > vt + 8

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c( 'run_mode',
               'src_path', 'scr_path', 'exe_path', 'prgm_dir', 'prgm_tag' )
opt_reqs <- c( 'out_path', # 'ref_path', 'ref_file', 'ref_build', 'ref_species',
               'Rscript', 'verbose' )

opt$rcpp <- 2
opt$rcpp <- 0
#
# TBD:: Return Rcpp value instead of boolean::
#
prgm_dat <- program_init( name = par$prgm_tag,
                          opts = opt, opt_reqs = opt_reqs,
                          pars = par, par_reqs = par_reqs, 
                          rcpp = opt$rcpp,
                          vb=vb, vt=3, tc=0 )

opt <- prgm_dat$opt
par <- prgm_dat$par
opt_tib <- prgm_dat$opt_tib
par_tib <- prgm_dat$par_tib

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Initialize Run Objects
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tt <- timeTracker$new()

pmssg <- glue::glue("[{par$prgm_tag}]:")
pwarn <- glue::glue("{RET}[{par$prgm_tag}]: Warning:")
perrs <- glue::glue("{RET}[{par$prgm_tag}]: ERROR:")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#               Pre-processing:: Initialize Parameter Vectors
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

gst_manifest_vec <- NULL
gst_manifest_vec <- stringr::str_split( 
  string = opt$genomestudio, 
  pattern = ",", simplify = FALSE ) %>% 
  base::unlist() %>% as.vector()

ses_manifest_vec <- NULL
ses_manifest_vec <- stringr::str_split( 
  string = opt$sesame, 
  pattern = ",", simplify = FALSE ) %>% 
  base::unlist() %>% as.vector()

manifest_key_vec <- NULL
manifest_key_vec <- stringr::str_split( 
  string = opt$manifest_name, 
  pattern = ",", simplify = FALSE ) %>% 
  base::unlist() %>% as.vector()

manifest_key_cnt <- manifest_key_vec %>% length()
manfiest_max_cnt <- max( length(gst_manifest_vec),length(ses_manifest_vec) )
if ( manfiest_max_cnt != manifest_key_cnt )
  cat(glue::glue("{pmssg} {BRK}{RET}",
                 "{pmssg} {S15}   FAILED Manifest Counts NOT equal!{RET}",
                 "{pmssg} {S25} max({manfiest_max_cnt}) != key({manifest_key_cnt}) {RET}",
                 "{pmssg} {BRK}{RET}"))
stopifnot( manfiest_max_cnt == manifest_key_cnt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Pre-processing:: Load Manifests/Controls
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Load Controls to be added::
man_ctls_tib <- NULL
man_ctls_tib <- readr::read_rds( file = opt$controls )
# man_ctls_tib <- readr::read_rds( file.path( opt$top_path, "data/manifests/methylation/bgz/epic_ctls.rds" ) )

# Load both EPIC v1/v2 Analytical Manifests::
#   TBD:: Replace with global manifest and manifest identification in c++
manifest_tibs <- base::list()
if ( opt$run_name %>% stringr::str_starts("GSIBIOINFO-638") ) {
  for ( ii in c(1:manfiest_max_cnt) ) {
    manifest_tibs[[manifest_key_vec[ii]]] <- load_epicv2_manifests( 
      sesame_csv = ses_manifest_vec[ii], 
      genome_csv = gst_manifest_vec[ii],
      name = manifest_key_vec[ii],
      ctls = man_ctls_tib,
      out_dir = opt$out_path, 
      run_tag = opt$run_name, 
      reload = opt$reload, 
      reload_min = 0, 
      parallel = opt$parallel,
      vb=vb,vt=vt+2,tc=tc )
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Load Sample Sheets
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# A tibble: 327 × 19 (all ACA sample sheets)
# A tibble: 168 × 21 (unique)
sample_sheet_tib <- NULL

if ( opt$run_name %>% stringr::str_starts("GSIBIOINFO-638") ) {
  sample_sheet_tib <- build_epicv2_sample_sheet( 
    sheet_str = opt$sample_sheets, 
    idat_path = opt$idat_path, 
    sentrix_unq = TRUE,
    out_dir = opt$out_path, 
    run_tag = opt$run_name, 
    reload = opt$reload, 
    reload_min = 10, 
    parallel = opt$parallel, 
    vb=vb,vt=vt+4,tc=tc )
} else if ( opt$run_name %>% stringr::str_starts("Embark_v3") ) {
  idat_tib <- sesame::searchIDATprefixes( 
    dir.name = file.path( opt$idat_path ), 
    recursive = TRUE ) %>% 
    as.data.frame() %>% 
    tibble::as_tibble( rownames = "Sentrix_Name" ) %>% 
    purrr::set_names( c("Sentrix_Name","Prefix") )
  
  sheet_rm_vec <- c("206712840012_R01C01","206712840015_R01C01","206712840015_R04C02","206712840015_R06C02")
  
  sample_sheet_tib <- readr::read_csv( file = opt$sample_sheets[1], show_col_types = FALSE ) %>%
    dplyr::inner_join( idat_tib, by=c("Sentrix_Name") ) %>% 
    dplyr::filter( !Sentrix_Name %in% sheet_rm_vec )
  
  unq_man_key <- sample_sheet_tib$Manifest_Key %>% unique()
  manifest_tibs[[unq_man_key]] <- NULL
  manifest_tibs[[unq_man_key]] <- readr::read_csv( 
    file = opt$sesame[1], show_col_types = FALSE ) %>%
    dplyr::bind_rows( man_ctls_tib )
}
# sample_sheet_tib %>% dplyr::group_by(Sample_Group,Sample_Base) %>% dplyr::summarise( Count=n() ) %>% print(n=1000)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Reduce Sample Sheets
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Sample Count Filtering:: May Want to Remove this based on Experiment...
#  - NOTE:: OK for EPIC v2 ("GSIBIOINFO-638")
#
opt$min_samp_cnt <- 3
sample_sheet_tib <- sample_sheet_tib %>%
  dplyr::add_count( Sample_Name, name="Sample_Rep" ) %>% 
  dplyr::filter( Sample_Rep >= opt$min_samp_cnt )

# Sample Sheets by Sentrix_Name::
sentrix_sheets <- NULL
sentrix_sheets <- sample_sheet_tib %>% split(.$Sentrix_Name)

# Sample Sheets by Manifest::
manifest_sheets <- NULL
manifest_sheets <- sample_sheet_tib %>% split(.$Manifest_Key )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Clean EPIC v2 Manifest with v1 Intersection::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

epi1_man_tmp <- NULL
epi1_man_tmp <- manifest_tibs[["EPIC_v1"]] %>%
  dplyr::select( Probe_ID ) %>%
  dplyr::mutate( v1 = "v1" )

epic_join_tib <- NULL
epic_join_tib <- manifest_tibs[["EPIC_v2"]] %>% 
  dplyr::filter( !Probe_ID %>% stringr::str_starts("^ctl_") ) %>% 
  dplyr::select( Probe_ID:AlleleB_ProbeSeq ) %>%
  dplyr::select( -mask ) %>%
  dplyr::left_join( epi1_man_tmp, by=c("Probe_ID") ) %>%
  dplyr::mutate( Loci_ID = Probe_ID %>% stringr::str_remove("_.*$") ) %>%
  dplyr::select( Probe_ID,Loci_ID, dplyr::everything() )

#
# Wanding Masks::
#
wan_rm1_tib <- NULL
wan_rm1_txt <- file.path( opt$top_path, "Projects/EPIC_v2/masks/Wanding_rm_list_33268.txt" )
wan_rm1_tib <- readr::read_tsv( file = wan_rm1_txt, show_col_types = FALSE )

# sesameData::sesameDataList()
ses_epic_add_dat <- NULL
ses_epic_add_tib <- NULL
ses_epic_add_dat <- sesameData::sesameDataGet( title = "EPIC.address" )
ses_epic_add_tib <- ses_epic_add_dat$ordering %>% tibble::as_tibble()

ses_epic_add_sum <- ses_epic_add_tib %>% 
  dplyr::group_by( mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

#
# Incorrect Matching::
#
# ses_epic_add_tib %>% dplyr::filter( U %in% wan_rm1_tib$AddressA_ID | M %in% wan_rm1_tib$AddressA_ID )
#
# wan_rm1_tib %>% 
#   dplyr::anti_join( ses_epic_add_tib, by=c("AddressA_ID"="U") ) %>% 
#   dplyr::anti_join( ses_epic_add_tib, by=c("AddressA_ID"="M") )

# A tibble: 866,553 × 5
aca_mask_man_tib <- NULL
aca_mask_man_tib <- epic_join_tib %>% 
  dplyr::filter( U %in% wan_rm1_tib$AddressA_ID | M %in% wan_rm1_tib$AddressA_ID )

epic_v1_mask_tib <- NULL
epic_v1_mask_tib <- manifest_tibs[["EPIC_v1"]] %>% 
  dplyr::select( Probe_ID, U ) %>%
  dplyr::inner_join( ses_epic_add_tib %>% dplyr::rename( Loci_ID = Probe_ID),
                     by=c("U") ) %>%
  dplyr::select( Probe_ID,Loci_ID, U,M, col,mask )

epic_v1_mask_sum <- epic_v1_mask_tib %>% 
  dplyr::filter( Probe_ID %in% aca_mask_man_tib$Probe_ID ) %>% 
  dplyr::group_by( mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
epic_v1_mask_sum %>% print()

#
# Get most recent Sesame mask_recommended::
#
sdf_rec_masked_tib <- NULL
sdf_rec_masked_tib <- sesameDataGet('EPIC.1.SigDF') %>%
  sesame::qualityMask( ) %>%
  tibble::as_tibble()

sdf_rec_masked_sum <- NULL
sdf_rec_masked_sum <- sdf_rec_masked_tib %>% 
  dplyr::group_by( mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
sdf_rec_masked_sum %>% print()

sdf_aca_masked_sum <- NULL
sdf_aca_masked_sum <- sdf_rec_masked_tib %>% 
  dplyr::filter( Probe_ID %in% aca_mask_man_tib$Loci_ID ) %>% 
  dplyr::group_by( mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
sdf_aca_masked_sum %>% print()

# All Filter Types::
maks_lists <- sesame::listAvailableMasks("EPIC")[,c("mask_name","mask_recommended","mask_group")]

ses_rec_masked_tib <- NULL
ses_rec_masked_tib <- sdf_rec_masked_tib %>% 
  dplyr::select( Probe_ID,mask ) %>% 
  dplyr::rename( Loci_ID = Probe_ID, 
                 Sesame_Mask = mask ) %>%
  dplyr::right_join( epic_join_tib, by=c("Loci_ID")  ) %>%
  dplyr::mutate( 
    Sesame_Mask = Sesame_Mask %>% as.character(),
    dplyr::across( Sesame_Mask, ~tidyr::replace_na(.x, "v2" ) ) )

ses_rec_masked_sum <- ses_rec_masked_tib %>% 
  dplyr::group_by( Sesame_Mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
ses_rec_masked_sum %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                             Pre-processing
#                                   END
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                   BEG
#                 Start Previously Processed Data Loading
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$dat_path <- file.path( opt$top_path, "scratch/swifthoof_sesame_test_mask/GSIBIOINFO-638-v0" )

uhm_top_v1_tib <- NULL
uhm_top_v1_csv <- file.path( opt$dat_path, "mask/iDNPBV.negs-1.poob-1/iDNPBV.negs-1.poob-1.Lightning.Auto.EPIC.1.EPIC_v1.MeTritration.Epigen.uhm.r2.mask.dat.csv.gz" )
uhm_top_v1_tib <- readr::read_csv( file = uhm_top_v1_csv, show_col_types = FALSE ) %>% 
  dplyr::mutate( uhm_v1_mask = mask )
uhm_top_v1_sum <- NULL
uhm_top_v1_sum <- uhm_top_v1_tib %>% dplyr::group_by( uhm_v1_mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

uhm_top_v2_tib <- NULL
uhm_top_v2_csv <- file.path( opt$dat_path, "mask/iDNPBV.negs-1.poob-1/iDNPBV.negs-1.poob-1.Lightning.Auto.EPIC.2.EPIC_v2.MeTritration.Epigen.uhm.r2.mask.dat.csv.gz" )
uhm_top_v2_tib <- readr::read_csv( file = uhm_top_v2_csv, show_col_types = FALSE ) %>% 
  dplyr::mutate( uhm_v2_mask = mask )
uhm_top_v2_sum <- NULL
uhm_top_v2_sum <- uhm_top_v2_tib %>% dplyr::group_by( uhm_v2_mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

epic_uhm_man_tib <- NULL
epic_uhm_man_tib <- epic_join_tib %>%
  dplyr::right_join( dplyr::select( uhm_top_v2_tib, Probe_ID,uhm_v2_mask ), 
                     by=c("Probe_ID") ) %>%
  dplyr::left_join( dplyr::select( uhm_top_v1_tib, Probe_ID,uhm_v1_mask ), 
                    by=c("Probe_ID") )

epic_uhm_v2_sum <- NULL
epic_uhm_v2_sum <- epic_uhm_man_tib %>% 
  dplyr::group_by( v1,uhm_v2_mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
# A tibble: 4 × 3
# v1    uhm_v2_mask  Count
# <chr> <lgl>        <int>
# 1 v1    FALSE       660977
# 2 v1    TRUE         93357
# 3 NA    FALSE       173392
# 4 NA    TRUE         42078
epic_uhm_v21_sum <- NULL
epic_uhm_v21_sum <- epic_uhm_man_tib %>% 
  dplyr::group_by( v1,uhm_v2_mask,uhm_v1_mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
# A tibble: 6 × 4
# v1    uhm_v2_mask uhm_v1_mask  Count
# <chr> <lgl>       <lgl>        <int>
# 1 v1    FALSE       FALSE       646725  => Pass
# 2 v1    FALSE       TRUE         14252  => Fail
# 3 v1    TRUE        FALSE        23373  => Pass
# 4 v1    TRUE        TRUE         69984  => Maybe
# 5 NA    FALSE       NA          173392  => Pass
# 6 NA    TRUE        NA           42078  => Maybe
#

#
# Load Previous Summary Files::
#
sum_dir <- safe_mkdir( file.path( opt$dat_path, "summary" ) )

uhm_all_csv <- file.path( sum_dir, "uhm_all_sum.csv.gz" )
rep_all_csv <- file.path( sum_dir, "rep_all_sum.csv.gz" )
ssh_all_csv <- file.path( sum_dir, "ssh_all_sum.csv.gz" )

uhm_all_sum <- NULL
rep_all_sum <- NULL
ssh_all_tib <- NULL

uhm_all_sum <- readr::read_csv( file = uhm_all_csv, show_col_types = FALSE )
rep_all_sum <- readr::read_csv( file = rep_all_csv, show_col_types = FALSE )
ssh_all_tib <- readr::read_csv( file = ssh_all_csv, show_col_types = FALSE )

#
# Quick Check on Sample Sheets::
#  To Be the code below demonstrates that the NA12873 was bad for EPIC v2...
#
# Example 1:: NA12873 was bad for EPIC v2
# ssh_all_tib %>% 
#   # tibble::as_tibble() %>%
#   # dplyr::ungroup() %>%
#   dplyr::arrange( Probe_PdP ) %>% 
#   # dplyr::filter( Poob != "poob-1" ) %>% 
#   # dplyr::filter( Sample_Base == "NA12873" ) %>%
#   # dplyr::filter( Sample_Base == "HeLa" ) %>%
#   dplyr::filter( Work == "iCDNPBV" ) %>%
#   dplyr::group_by( Manifest_Key,Poob,Work,Sample_Base ) %>% 
#   dplyr::summarise( Avg=mean( Probe_PdP, na.rm=TRUE), 
#                     Med=median( Probe_PdP, na.rm=TRUE), 
#                     Min=min( Probe_PdP, na.rm=TRUE ),
#                     Sds=sd( Probe_PdP, na.rm=TRUE ),
#                     Mad=mad( Probe_PdP, na.rm=TRUE ),
#                     .groups = "drop" ) %>%
#   dplyr::arrange( Med ) %>% split(.$Sample_Base) # print(n=1000)

#
# Example 2:: NA12873 was bad for EPIC v2
#
rep_NA2_sum_list <- NULL
rep_NA2_sum_list <- rep_all_sum %>% 
  dplyr::arrange( -Per ) %>% 
  dplyr::filter( Poob == "poob-1" ) %>% 
  # dplyr::filter( Manifest_Key != "EPIC_v1" ) %>% 
  dplyr::group_by( Manifest_Key,Poob,Work,Sample_Base ) %>% 
  dplyr::summarise( Avg=mean(Per, na.rm=TRUE), 
                    Med=median(Per, na.rm=TRUE), 
                    Min=min(Per, na.rm=TRUE ),
                    Sds=sd(Per, na.rm=TRUE ),
                    Mad=mad(Per, na.rm=TRUE ),
                    .groups = "drop" ) %>%
  dplyr::arrange( Med ) %>% split(.$Sample_Base)

#
# Load Both:: iDNPBV & iCDNPBV [except NA12873]
#
rep_top_files <- NULL
rep_top_files <- file_list( 
  path = file.path( opt$dat_path, "mask/iCDNPBV.negs-1.poob-0.1" ), 
  pattern = "CellLine.*.rep.r2.mask.dat.csv.gz$", 
  suffix = ".rep.r2.mask.dat.csv.gz$", 
  prefix = "iCDNPBV.negs-1.poob-01.Lightning.Auto.EPIC.[12].EPIC_" )

rep_top_tab <- NULL
rep_top_tab <- rep_top_files %>% 
  lapply( readr::read_csv, show_col_types = FALSE ) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate( Chip_VerStr = paste0("v",Chip_Version) )

rep_top_sum <- NULL
rep_top_sum <- rep_top_tab %>% 
  dplyr::group_by( Chip_VerStr,Sample_Base ) %>% 
  dplyr::summarise( Tot = n(), 
                    Nan = sum( is.na(mask) ), 
                    Mis = sum(  mask, na.rm=TRUE ), 
                    Pas = sum( !mask, na.rm=TRUE ),
                    Per = round( 100*Pas/Tot, 3 ),
                    .groups = "drop")

#
# Get Failure Counts::
#   - v1
#   - v2
#   - v1 + v2
#
rep_cnt_12_tib <- NULL
rep_cnt_12_tib <- rep_top_tab %>% 
  # dplyr::select( Probe_ID,Chip_VerStr,Sample_Base,mask ) %>%
  dplyr::select( Probe_ID,Chip_VerStr,mask ) %>%
  dplyr::group_by( Probe_ID,Chip_VerStr, mask ) %>%
  dplyr::summarise( Count=n(),
                    .groups = "drop" )

# QC Sanity Check::
#  rep_cnt_12_tib %>% dplyr::filter( Chip_VerStr == "v1" ) %>% dplyr::distinct( Probe_ID )
#  rep_cnt_12_tib %>% dplyr::filter( Chip_VerStr == "v2" ) %>% dplyr::distinct( Probe_ID )
rep_cnt_12_tab <- NULL
rep_cnt_12_tab <- rep_cnt_12_tib %>% 
  tidyr::pivot_wider( id_cols = c(Probe_ID), 
                      names_from = c(Chip_VerStr,mask), 
                      values_from = c(Count), 
                      names_sep = "_", values_fill = 0 ) %>% 
  dplyr::mutate( v1_Tot = v1_FALSE+v1_TRUE, 
                 v2_Tot = v2_FALSE+v2_TRUE,
                 v1_Per = round( 100*v1_FALSE/v1_Tot, 3),
                 v2_Per = round( 100*v2_FALSE/v2_Tot, 3) )

# QC Sanity Check::
#  rep_cnt_12_tab %>% dplyr::filter( is.na(v1_Per) )
#  rep_cnt_12_tab %>% dplyr::filter( is.na(v2_Per) )

#
#
# Put it all together::
#
#

# A tibble: 6 × 4
# v1    uhm_v2_mask uhm_v1_mask  Count
# <chr> <lgl>       <lgl>        <int>
# 1 v1    FALSE       FALSE       646725  => Pass
# 2 v1    FALSE       TRUE         14252  => Fail
# 3 v1    TRUE        FALSE        23373  => Pass
# 4 v1    TRUE        TRUE         69984  => Maybe
# 5 NA    FALSE       NA          173392  => Pass
# 6 NA    TRUE        NA           42078  => Maybe

fin_scr_tib <- NULL
fin_scr_tib <- epic_uhm_man_tib %>%
  dplyr::left_join( rep_cnt_12_tab, by="Probe_ID" ) %>%
  dplyr::mutate( 
    dplyr::across( v1, ~tidyr::replace_na(.x, "v2" ) ),
    uhm_score = dplyr::case_when(
      !uhm_v2_mask & !uhm_v1_mask ~  1,
      !uhm_v2_mask &  uhm_v1_mask ~ -1,
      uhm_v2_mask &  !uhm_v1_mask ~  1,
      uhm_v2_mask &   uhm_v1_mask ~  0,
      !uhm_v2_mask & is.na(uhm_v1_mask) ~ 1,
      uhm_v2_mask &  is.na(uhm_v1_mask) ~ 0,
      TRUE ~ NA_real_ ),
    rep_score = dplyr::case_when(
      v2_Per <   75 ~ -1,
      v2_Per <  100 ~  0,
      v2_Per == 100 ~  1,
      TRUE ~ NA_real_ ),
    fin_score = uhm_score + rep_score,
    mask = dplyr::case_when(
      v1 == "v1" & fin_score >= 0 ~ "FALSE",
      v1 == "v2" & fin_score >= 1 ~ "FALSE",
      TRUE ~ "TRUE" )
  )


if ( FALSE ) {
  fin_scr_tib %>%
    dplyr::summarise( Tot = n(),
                      v1_Avg = mean( v1_Per, na.rm=TRUE ),
                      v2_Avg = mean( v2_Per, na.rm=TRUE ),
                      v1_Med = median( v1_Per, na.rm=TRUE ),
                      v2_Med = median( v2_Per, na.rm=TRUE ) )
  
  fin_scr_tib %>% dplyr::group_by( v2_Per ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  fin_scr_tib %>% dplyr::group_by( uhm_score ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  fin_scr_tib %>% dplyr::group_by( fin_score ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  
  fin_scr_tib %>% dplyr::group_by( v1,fin_score ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  fin_scr_tib %>% dplyr::group_by( mask ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  
  fin_scr_tib %>% dplyr::group_by( uhm_score,rep_score,fin_score ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  #  920012 = 39438+148992+62747+87246+581589
}

fin_man_tib <- NULL
fin_man_csv <- file.path( opt$out_path, "EPIC_v2.masked.manifest.v0.csv.gz" )
fin_man_tib <- fin_scr_tib %>% dplyr::select( Probe_ID:v1, mask )
# readr::write_csv( x = fin_man_tib, file = fin_man_csv )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                 Start Previously Processed Data Loading
#                                   END
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                   BEG
#                               Build Plots
#
#  1. Summary:: uhm (v1/v2)
#  2. Summary:: rep (v1/v2)
#
#  3. Summary:: beta
#     - Reload uhm/rep beta matrix...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

uhm_top_v1_tab <- NULL
uhm_top_v1_tab <- uhm_top_v1_tib %>% 
  dplyr::left_join( rep_cnt_12_tab %>% dplyr::select(Probe_ID,v1_Per), by=c("Probe_ID") ) %>% 
  dplyr::rename( rpp =  v1_Per ) %>%
  dplyr::mutate( rpp = rpp / 100 ) %>%
  tidyr::pivot_longer( cols = c(dHU,dMH,dMU,rpp), 
                       names_to = "Metric_Key", 
                       values_to = "Metric_Val", 
                       values_drop_na = TRUE )

uhm_top_v2_tab <- NULL
uhm_top_v2_tab <- uhm_top_v2_tib %>% 
  dplyr::left_join( rep_cnt_12_tab %>% dplyr::select(Probe_ID,v2_Per), by=c("Probe_ID") ) %>% 
  dplyr::rename( rpp =  v2_Per ) %>%
  dplyr::mutate( rpp = rpp / 100 ) %>%
  tidyr::pivot_longer( cols = c(dHU,dMH,dMU,rpp), 
                       names_to = "Metric_Key", 
                       values_to = "Metric_Val",
                       values_drop_na = TRUE )


uhm_top_tab <- NULL
uhm_top_tab <- dplyr::bind_rows( uhm_top_v1_tab,uhm_top_v2_tab ) %>%
  dplyr::filter( !is.na(Metric_Val) )

# uhm_fin_tab %>% dplyr::filter( Chip_Version == 1 )
# rep_cnt_12_tab

uhm_fin_tab <- NULL
uhm_fin_tab <- fin_scr_tib %>% 
  dplyr::select( Probe_ID, uhm_score,rep_score,fin_score ) %>% 
  dplyr::left_join( dplyr::select( ses_rec_masked_tib, Probe_ID,Sesame_Mask ), by=c("Probe_ID") ) %>%
  dplyr::inner_join( uhm_top_tab, by=c("Probe_ID") )

# Better be zero::
# uhm_fin_tab %>% dplyr::filter( is.na(fin_score) ) %>% dplyr::group_by( Manifest_Key ) %>% dplyr::summarise( Count=n(), .groups = "drop" )

#
# TBD::
#   - Add Rep Passing Percent
#   - Add Sesame Mask Column
#
plot_dir <- safe_mkdir( dir = file.path( opt$out_path, "plots" ) )
dpi_val <- 320
plot_type <- "pdf"
build_plots <- TRUE
if ( build_plots ) {

  uhm_rep_den_gg <- NULL
  uhm_rep_den_gg <- uhm_fin_tab %>% 
    ggplot2::ggplot( aes( x=Metric_Val, color=Manifest_Key, fill=Manifest_Key ) ) + 
    ggplot2::geom_density( alpha = 0.2 ) + 
    ggplot2::facet_grid( rows = vars(fin_score,Sesame_Mask), 
                         cols = vars(Metric_Key), 
                         scales = "free_y" ) +
    ggplot2::theme( legend.position = "bottom" )
  
  # uhm_rep_den_gg
  uhm_rep_den_pdf <- file.path( plot_dir, paste("uhm-replicate-percent-passing.density",plot_type, sep=".") )
  if ( p0 ) cat(glue::glue("{pmssg} Plotting[{plot_type}] = '{uhm_rep_den_pdf}'...{RET2}"))
  ggplot2::ggsave( filename = uhm_rep_den_pdf, 
                   plot = uhm_rep_den_gg, 
                   device = plot_type, 
                   dpi = dpi_val, 
                   width = 7, height = 7 )
  
  uhm_den_gg <- NULL
  uhm_den_gg <- uhm_fin_tab %>% 
    dplyr::filter( Metric_Key != "rpp" ) %>%
    ggplot2::ggplot( aes( x=Metric_Val, color=Manifest_Key, fill=Manifest_Key ) ) + 
    ggplot2::geom_density( alpha = 0.2 ) + 
    ggplot2::facet_grid( rows = vars(fin_score,Sesame_Mask), 
                         cols = vars(Metric_Key), 
                         scales = "free_y" ) +
    ggplot2::theme( legend.position = "bottom" )
  
  # title = "Methylation Titration (All)", 

  # uhm_den_gg
  uhm_den_pdf <- file.path( plot_dir, paste("uhm.density",plot_type, sep=".") )
  if ( p0 ) cat(glue::glue("{pmssg} Plotting[{plot_type}] = '{uhm_den_pdf}'...{RET2}"))
  ggplot2::ggsave( filename = uhm_den_pdf, 
                   plot = uhm_den_gg, 
                   device = plot_type, 
                   dpi = dpi_val, 
                   width = 7, height = 7 )
  
  rep_den_gg <- NULL
  rep_den_gg <- uhm_fin_tab %>% 
    dplyr::filter( Metric_Key == "rpp" ) %>%
    ggplot2::ggplot( aes( x=Metric_Val, color=Manifest_Key, fill=Manifest_Key ) ) + 
    ggplot2::geom_density( alpha = 0.2 ) + 
    ggplot2::facet_grid( rows = vars(fin_score,Sesame_Mask), 
                         cols = vars(Metric_Key), 
                         scales = "free_y" ) +
    ggplot2::theme( legend.position = "bottom" )
  
  # rep_den_gg
  rep_den_pdf <- file.path( plot_dir, paste("replicate-percent-passing.density",plot_type, sep=".") )
  if ( p0 ) cat(glue::glue("{pmssg} Plotting[{plot_type}] = '{rep_den_pdf}'...{RET2}"))
  ggplot2::ggsave( filename = rep_den_pdf, plot = rep_den_gg, 
                   device = plot_type, 
                   dpi = dpi_val, 
                   width = 7, height = 7 )
  
}

# epic2_beta_rds <- file.path( opt$dat_path, "sdf/work/EPIC_v2/mutate_sdfs/GSIBIOINFO-638-v0/work-iCDNPBV/negs-1/poob-0.1/mutate_sdfs.GSIBIOINFO-638-v0.work-iCDNPBV.negs-1.poob-0.1.rds" )


# sample_sheet_tib
epic1_beta_mat <- NULL
epic1_beta_rds <- file.path( opt$dat_path, "sdf/work/EPIC_v1/mutate_sdfs/GSIBIOINFO-638-v0/work-iDNPBV/negs-1/poob-1/mutate_sdfs.GSIBIOINFO-638-v0.work-iDNPBV.negs-1.poob-1.rds" )
epic1_beta_mat <- readr::read_rds( file = epic1_beta_rds )

epic2_beta_mat <- NULL
epic2_beta_rds <- file.path( opt$dat_path, "sdf/work/EPIC_v2/mutate_sdfs/GSIBIOINFO-638-v0/work-iDNPBV/negs-1/poob-1/mutate_sdfs.GSIBIOINFO-638-v0.work-iDNPBV.negs-1.poob-1.rds" )
epic2_beta_mat <- readr::read_rds( file = epic2_beta_rds )

epic12_cgn <- NULL
epic12_cgn <- intersect( rownames(epic1_beta_mat), rownames(epic2_beta_mat) )

betaI_mat <- NULL
betaI_mat <- merge( epic1_beta_mat[ epic12_cgn, ], epic2_beta_mat[ epic12_cgn, ],
                    by = 'row.names', all = TRUE)
rownames( betaI_mat ) <- epic12_cgn
# betaI_mat %>% colnames()

#
# Loop Over Plots::
#
sample_names <- c( "Epigen","HeLa","Jurkat","MCF7","Raji","Zymo")

for ( sample_name in sample_names ) {
  # sample_name <- "Jurkat"
  sample_pref <- sample_name %>% stringr::str_sub( 1,1 )
  sheet_prep <- "Lightning"
  sheet_proc <- "Auto"
  
  # Add print()
  
  tar_ssh_tib <- NULL
  tar_ssh_tib <- sample_sheet_tib %>% 
    dplyr::filter( Sheet_Prep == sheet_prep & Sheet_Proc == sheet_proc ) %>%
    dplyr::filter( Sample_Base == sample_name ) %>% 
    # dplyr::filter( Sample_Group == "MeTritration" ) %>% 
    dplyr::arrange( Manifest_Key,Sample_Base,Concentration ) %>% 
    dplyr::group_by( Chip_Version,Sample_Name ) %>% 
    dplyr::mutate( Rep_Idx = row_number() ) %>% dplyr::ungroup() %>% 
    dplyr::mutate( Plot_Name = paste0(stringr::str_sub( Sample_Name,1,1 ), 
                                      Chip_Version,"_",Concentration,"_",
                                      Rep_Idx ) ) %>% 
    dplyr::arrange( Concentration,Chip_Version,Manifest_Key ) %>% 
    dplyr::distinct( Manifest_Key,Sample_Base,Concentration, .keep_all = TRUE )
  
  #
  # Plotting::
  #
  min_pval <- 0.05
  grp_key <- "col"
  cond_grp_vec <- NULL
  
  beta_dat_tib <- NULL
  beta_dat_tib <- betaI_mat[ , tar_ssh_tib$Sentrix_Name ] %>%
    as.data.frame() %>%
    tibble::as_tibble( rownames = "Probe_ID" ) %>%
    purrr::set_names( c( "Probe_ID", tar_ssh_tib$Plot_Name ) )
  
  # Need to add Manifest Columns back...
  pval_all_tib <- NULL
  beta_all_tib <- NULL
  beta_all_tib <- dplyr::select( manifest_tibs[[1]], Probe_ID,col ) %>% 
    dplyr::inner_join(beta_dat_tib, by=c("Probe_ID") )
  pval_all_tib <- beta_all_tib
  
  beta_sub_tib <- beta_all_tib %>% 
    dplyr::select( Probe_ID,col, dplyr::starts_with( sample_pref ) ) %>% 
    # dplyr::select( Probe_ID,col, dplyr::starts_with("E") ) %>% 
    # dplyr::select( Probe_ID,col, dplyr::starts_with("H") ) %>% 
    head( n=10000 )
  
  # Leaving the col column is actually kind of cool...
  rank_cols <- c(2:base::ncol(beta_sub_tib) )
  # Remove the col column::
  rank_cols <- c(3:base::ncol(beta_sub_tib) )
  
  
  #
  # Need to update beta_panel_gg() for when pval_data == NULL
  #
  pairs_gg <- NULL
  pairs_gg <- GGally::ggpairs( 
    data = beta_sub_tib,
    mapping = ggplot2::aes( color = grp_key,
                            fill  = grp_key,
                            alpha = 0.2 ),
    columns = rank_cols,
    
    upper = list(
      combo = "box_no_facet",
      continuous = GGally::wrap( beta_panel_gg,
                                 plot_type = "dB",
                                 pval_min = min_pval,
                                 beta_min = 0.3,
                                 beta_max = 0.7,
                                 dB_min = 0.2,
                                 
                                 beta_data = beta_all_tib,
                                 pval_data = pval_all_tib,
                                 pval_cols = rank_cols,
                                 full_plot = FALSE,
                                 grp_key = grp_key,
                                 grp_vec = cond_grp_vec,
                                 
                                 out_dir = plot_dir,
                                 run_tag = "Upper", 
                                 reload = opt$reload,
                                 reload_min = 10, 
                                 ret_data = FALSE,
                                 parallel = opt$parallel,
                                 
                                 vb=vb,vt=vt+1,tc=tc+1,tt=tt ) ),
    
    diag = list(
      continuous = GGally::wrap( beta_panel_gg,
                                 plot_type = "dg",
                                 pval_min = min_pval,
                                 beta_min = 0.3,
                                 beta_max = 0.7,
                                 dB_min = 0.2,
                                 
                                 beta_data = beta_all_tib,
                                 pval_data = pval_all_tib,
                                 pval_cols = rank_cols,
                                 full_plot = FALSE,
                                 grp_key = grp_key,
                                 grp_vec = cond_grp_vec,
                                 
                                 out_dir = plot_dir,
                                 run_tag = "Diag", 
                                 reload = opt$reload,
                                 reload_min = 10, 
                                 ret_data = FALSE,
                                 parallel = opt$parallel,
                                 
                                 vb=vb,vt=vt+1,tc=tc+1,tt=tt ) ),
    
    lower = list(
      combo = "box_no_facet",
      continuous = GGally::wrap( beta_panel_gg,
                                 plot_type = "r2",
                                 pval_min = min_pval,
                                 beta_min = 0.3,
                                 beta_max = 0.7,
                                 dB_min = 0.2,
                                 
                                 beta_data = beta_all_tib,
                                 pval_data = pval_all_tib,
                                 pval_cols = rank_cols,
                                 full_plot = FALSE,
                                 grp_key = grp_key,
                                 grp_vec = cond_grp_vec,
                                 
                                 out_dir = plot_dir,
                                 run_tag = "Lower", 
                                 reload = opt$reload,
                                 reload_min = 10, 
                                 ret_data = FALSE,
                                 parallel = opt$parallel,
                                 
                                 vb=vb,vt=vt+1,tc=tc+1,tt=tt ) )
  ) + labs( title = paste0( sample_name ) )
  #
  # Need to update this...
  #
  # pairs_gg <- pairs_gg + 
  #   labs( title=top_tag, subtitle=sub_tag, caption=sub_per_str )
  
  pairs_pdf <- NULL
  pairs_pdf <- file.path( plot_dir, paste(sample_name,"pairs",plot_type, sep=".") )
  if ( p0 ) cat(glue::glue("{pmssg} Plotting[{plot_type}] = '{pairs_pdf}'...{RET2}"))
  ggplot2::ggsave( filename = pairs_pdf, 
                   plot = pairs_gg, 
                   device = plot_type, 
                   dpi = dpi_val, 
                   width = 7, height = 7 )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                               Build Plots
#                                   END
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt,tc=tc,tt=tt )

# End of file
