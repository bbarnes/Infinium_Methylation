
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Stable Docker Manifest:: 
#               Scratch Space to Fix EPIC v2 Docker Manifest
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                Set Run Environment:: RStudio/Command-Line
#                            Source All Functions
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par  <- list()
args <- commandArgs(trailingOnly = FALSE)

par$src_path <- NULL
par$run_mode <- args[1]
par$date_str <- Sys.Date() %>% as.character()
par$prgm_dir <- 'stable'
par$prgm_tag <- 'stable_docker_manifest'
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

par$run_name <- "EPICv2"

# Multiple copies of v0 (bk and bk24)
# par$version <- 0
# par$version <- 1
# par$version <- 2
# par$version <- 3
# par$version <- 4
par$version <- 5

opt <- NULL
# opt <- swifthoof_sesame_options( pars = par, args = args, vb = par$verbose )
opt <- stable_options( pars = par, args = args, vb = par$verbose )
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

success = TRUE

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Workflow::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: Record Workflow...

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Pre-processing:: Load Sesame CGN's
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sesameData::sesameDataCacheAll()

ses_list_tib <- NULL
ses_list_tib <- sesameData::sesameDataList() %>% 
  as.data.frame() %>% tibble::as_tibble()

ses_inf_list <- NULL
ses_inf_list <- ses_list_tib %>% 
  dplyr::filter( Title %>% stringr::str_ends("inference") ) %>%
  split( .$Title )

#
# Target Inference Fields::
#
# Age: "age.inference"
#   - "SkinBlood"
# "ethnicity.inference"
#.  - Only one option
# "sex.inference"
#.  - This should be all the x/y chromosomes...
#.  - Should be inferred from overlapping EPIC v1 chrX/chrY
#
# NOTE: GCT Score might already be calculated...
#

#
# Target:: age.inference + SkinBlood
#
age_dat <- NULL
age_dat <- sesameData::sesameDataGet( title = "age.inference" )

age_skin_tib <- NULL
age_skin_tib <- age_dat$SkinBlood$CpGmarker %>% 
  # as.data.frame() %>%
  tibble::as_tibble() %>% 
  dplyr::rename( Probe_ID = value ) %>%
  dplyr::filter( Probe_ID != "intercept" ) %>%
  dplyr::mutate( Rm_Suffix_Age = "skin_blood_age" ) %>%
  dplyr::arrange( Probe_ID )

#
# Target:: ethnicity.inference
#
eth_dat <- NULL
eth_dat <- sesameData::sesameDataGet( title = "ethnicity.inference" )

eth_cgn_tib <- NULL
eth_cgn_tib <- tibble::tibble(
  Probe_ID = c( eth_dat$ccs.probes, eth_dat$rs.probes ),
  Rm_Suffix_Ethnicity = "ethnicity"
) %>%
  dplyr::arrange( Probe_ID )


#
# Probably don't need this loop...
#
if ( FALSE ) {
  for ( title in names(ses_inf_list) ) {
    tar_title <- ses_inf_list[[title]]$Title
    if ( p1 ) cat(glue::glue("{pmssg} Title='{ses_inf_list[[title]]$Title}', Target='{tar_title}'.{RET}"))
    
    cur_inf_dat <- sesameData::sesameDataGet( title = tar_title )
    
    cur_inf_dat %>% names() %>% print()
    
    if ( opt$single ) break
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Pre-processing:: Load Manifest Mask
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mask_tib <- NULL
mask_csv <- file.path( opt$top_path, "data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1-FlaggedProbes/EPIC-8v2-0_A1-FlaggedProbes.csv.gz" )
mask_tib <- readr::read_csv( file = mask_csv, show_col_types = FALSE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Pre-processing:: Load Manifest v1
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_v1_tib <- NULL
man_v1_csv <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_manifests.sub8.csv.gz" )
man_v1_tib <- readr::read_csv( file = man_v1_csv, show_col_types = FALSE ) %>% 
  dplyr::filter( Species == "homo_sapiens", 
                 Manifest == "EPIC", 
                 Manifest_Version == "B4" )

tar_sex_tib <- NULL
tar_sex_tib <- man_v1_tib %>% 
  dplyr::filter( Chromosome == "chrX" | Chromosome == "chrY" ) %>%
  dplyr::select( Probe_ID, Chromosome ) %>%
  dplyr::rename( Rm_Suffix_Chromosome = Chromosome ) %>%
  dplyr::arrange( Probe_ID )

tar_sex_sum <- NULL
tar_sex_sum <- tar_sex_tib %>% 
  dplyr::group_by( Rm_Suffix_Chromosome ) %>%
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p1 ) tar_sex_sum %>% print( n=base::nrow(tar_sex_sum) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                  Pre-processing:: Join All Suffix Data
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_suf_tib <- NULL
all_suf_tib <- age_skin_tib %>% 
  dplyr::full_join( eth_cgn_tib, by=c("Probe_ID") ) %>% 
  dplyr::full_join( tar_sex_tib,  by=c("Probe_ID") ) %>%
  dplyr::arrange( Probe_ID )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Pre-processing:: Load Manifest v2
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_v2_tib <- NULL
man_v2_ver <- "raw"
man_v2_ver <- "unq"
man_v2_dir <- file.path( opt$top_path, "data/manifests/methylation/stable_docker_manifest/core6/manifests", man_v2_ver  )
man_v2_csv <- file.path( man_v2_dir, "EPICv2-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
man_v2_tib <- readr::read_csv( file = man_v2_csv, show_col_types = FALSE ) %>% 
  dplyr::mutate(
    # This needs to happen first before Probe_ID's have suffix removed...
    mask = dplyr::case_when(
      Probe_ID %in% mask_tib$IlmnID ~ TRUE,
      TRUE ~ FALSE
    ),
    Full_ID = Probe_ID,
    Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"),
    Probe_ID = dplyr::case_when(
      Loci_ID %in% age_skin_tib$Probe_ID ~ Loci_ID,
      TRUE ~ Probe_ID
    ),
    Probe_ID = dplyr::case_when(
      Loci_ID %in% eth_cgn_tib$Probe_ID ~ Loci_ID,
      TRUE ~ Probe_ID
    ),
    Probe_ID = dplyr::case_when(
      Loci_ID %in% tar_sex_tib$Probe_ID ~ Loci_ID,
      TRUE ~ Probe_ID
    )
  ) %>% 
  dplyr::left_join( all_suf_tib, by=c("Loci_ID"="Probe_ID") ) %>% 
  dplyr::distinct( Probe_ID, .keep_all = TRUE )

age_mat_cnt1 <- man_v2_tib %>% dplyr::filter( Probe_ID %in% age_skin_tib$Probe_ID ) %>% base::nrow()
eth_mat_cnt1 <- man_v2_tib %>% dplyr::filter( Probe_ID %in% eth_cgn_tib$Probe_ID ) %>% base::nrow()
sex_mat_cnt1 <- man_v2_tib %>% dplyr::filter( Probe_ID %in% tar_sex_tib$Probe_ID ) %>% base::nrow()

age_mat_cnt2 <- man_v2_tib %>% dplyr::filter( !is.na(Rm_Suffix_Age) ) %>% base::nrow()
eth_mat_cnt2 <- man_v2_tib %>% dplyr::filter( !is.na(Rm_Suffix_Ethnicity) ) %>% base::nrow()
sex_mat_cnt2 <- man_v2_tib %>% dplyr::filter( !is.na(Rm_Suffix_Chromosome) ) %>% base::nrow()

if ( age_mat_cnt1 != age_mat_cnt2 ||
     eth_mat_cnt1 != eth_mat_cnt2 ||
     sex_mat_cnt1 != sex_mat_cnt2 ) {
  
  cat(glue::glue("{pmssg} Failed Count Metrics::{RET2}",
                 "{pmssg} Age: {age_mat_cnt1} != {age_mat_cnt2}.{RET}",
                 "{pmssg} Eth: {eth_mat_cnt1} != {eth_mat_cnt2}.{RET}",
                 "{pmssg} Sex: {sex_mat_cnt1} != {sex_mat_cnt2}.{RET2}") )
                 
} else {
  
  if ( p1 ) 
    cat(glue::glue("{pmssg} Passed Count Metrics::{RET2}",
                   "{pmssg} Age: {age_mat_cnt1} != {age_mat_cnt2}.{RET}",
                   "{pmssg} Eth: {eth_mat_cnt1} != {eth_mat_cnt2}.{RET}",
                   "{pmssg} Sex: {sex_mat_cnt1} != {sex_mat_cnt2}.{RET2}") )

}

#
# Validate Suffix Clean Probes::
#
# This should equal: age_mat_cnt2 + eth_mat_cnt2 + sex_mat_cnt2
man_suf_v2_tib <- man_v2_tib %>% 
  dplyr::filter( !Probe_ID %>% stringr::str_detect("_") )

#
# Masked Summary::
#
man_v2_mask_sum <- NULL
man_v2_mask_sum <- man_v2_tib %>% 
  dplyr::group_by( mask ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p1 ) man_v2_mask_sum %>% print( n=base::nrow(man_v2_mask_sum) )

#
# Probe_Type Summary::
#
man_v2_prbs_sum <- man_v2_tib %>% 
  dplyr::group_by( Probe_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p1 ) man_v2_prbs_sum %>% print( n=base::nrow(man_v2_prbs_sum) )

# Inspect non-cpg probes::
#
# man_v2_tib %>% dplyr::filter( Probe_Type == "rs" )
# man_v2_tib %>% dplyr::filter( Probe_Type == "ch" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                               Write Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Subset without extra annotation columns::
sub_v2_tib <- NULL
sub_v2_tib <- man_v2_tib %>% 
  dplyr::select( !mask ) %>%
  dplyr::select( !Loci_ID ) %>%
  dplyr::select( !Full_ID ) %>%
  dplyr::select( !Rm_Suffix_Age ) %>%
  dplyr::select( !Rm_Suffix_Ethnicity ) %>%
  dplyr::select( !Rm_Suffix_Chromosome )

out_core_ver <- "core7"
out_core_ver <- "core8"

out_man_path <- file.path( opt$out_path, "manifests", out_core_ver, man_v2_ver )
out_man_path <- safe_mkdir( dir = out_man_path )

sub_man_csv <- file.path( out_man_path, "EPICv2-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
all_man_csv <- file.path( out_man_path, "EPICv2-A1.manifest.sesame-base.annotated.cpg-sorted.csv.gz" )

readr::write_csv( x = man_v2_tib, file = all_man_csv )
readr::write_csv( x = sub_v2_tib, file = sub_man_csv )

#
# Build EPICv1 with all non-inference	probes with EPICv2 names
#
sel_sub_v1_vec <- NULL
sel_sub_v1_vec <- c("Probe_ID", "M", "U", "DESIGN", "COLOR_CHANNEL", "col", 
                    "Probe_Type", "Probe_Source", "Next_Base", "Probe_Design" ) # "mask")

# v1 needs controls::
# ctl_v1_tib <- man_v1_tib %>% dplyr::filter( Probe_Type != "cg" ) %>% dplyr::filter( Probe_Type != "ch" ) %>% dplyr::filter( Probe_Type != "rs" )
ctl_v2_tib <- NULL
ctl_v2_tib <- sub_v2_tib %>% 
  dplyr::filter( Probe_Type != "cg" ) %>% 
  dplyr::filter( Probe_Type != "ch" ) %>% 
  dplyr::filter( Probe_Type != "rs" ) %>%
  dplyr::filter( Probe_Type != "nv" )

ctl_v2_sum <- NULL
ctl_v2_sum <- ctl_v2_tib %>% 
  dplyr::group_by( Probe_Type ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

sub_v1_tib <- NULL
sub_v1_tib <- man_v2_tib %>% 
  dplyr::select( Probe_ID, Loci_ID, DESIGN, COLOR_CHANNEL, 
                 # col,
                 # Probe_Type, 
                 Probe_Source, Next_Base, Probe_Design ) %>% 
  # dplyr::right_join( man_v1_tib, by=c("Loci_ID"="Probe_ID") ) %>% 
  dplyr::inner_join( man_v1_tib, by=c("Loci_ID"="Probe_ID") ) %>% 
  dplyr::select( dplyr::all_of(sel_sub_v1_vec) ) %>% 
  dplyr::bind_rows( ctl_v2_tib ) %>%
  dplyr::distinct( Probe_ID, .keep_all = TRUE )

# Make sure the controls are present i.e. non-empty value below::
# sub_v1_tib %>% dplyr::filter( Probe_Type != "cg" ) %>% dplyr::filter( Probe_Type != "ch" ) %>% dplyr::filter( Probe_Type != "rs" )

man_suf_v1_tib <- NULL
man_suf_v1_tib <- sub_v1_tib %>% 
  dplyr::filter( !Probe_ID %>% stringr::str_detect("_") )

out_man_path <- file.path( opt$out_path, "manifests", out_core_ver, man_v2_ver )
out_man_path <- safe_mkdir( dir = out_man_path )

sub_v1_csv <- file.path( out_man_path, "EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz" )
# all_man_csv <- file.path( out_man_path, "EPIC-B4.manifest.sesame-base.annotated.cpg-sorted.csv.gz" )

readr::write_csv( x = sub_v1_tib, file = sub_v1_csv )
# readr::write_csv( x = man_v2_tib, file = all_man_csv )


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
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt+3,tc=tc,tt=tt )

sysTime <- Sys.time()
if ( p0 ) cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))


# End of file
