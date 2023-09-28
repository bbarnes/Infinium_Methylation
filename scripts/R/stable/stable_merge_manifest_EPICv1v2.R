
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
par$prgm_dir <- 'stable'
par$prgm_tag <- 'stable_merge_manifest_EPICv1v2'
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
par$version <- 0
# par$version <- 1
# par$version <- 2
# par$version <- 3
# par$version <- 4

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
#                     Pre-processing:: Load Manifests
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_dir <- file.path( opt$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/full" )

man_v1_csv <- file.path( man_dir, "EPIC_v1.gs_to_sesame.csv.gz" )
man_v2_csv <- file.path( man_dir, "EPIC_v2.gs_to_sesame.csv.gz" )

man_v1_tib <- NULL
man_v2_tib <- NULL

man_v1_tib <- readr::read_csv( file = man_v1_csv, show_col_types = FALSE )
man_v2_tib <- readr::read_csv( file = man_v2_csv, show_col_types = FALSE )

#
# Investigation of fields to join
#
if ( FALSE ) {
  man_tib1 <- dplyr::inner_join( man_v1_tib, man_v2_tib, by=c("AlleleA_ProbeSeq"), suffix=c("_v1", "_v2") )
  man_tib2 <- dplyr::inner_join( man_v1_tib, man_v2_tib, by=c("Probe_ID","AlleleA_ProbeSeq"), suffix=c("_v1", "_v2") )
  
  man_tib1 %>% dplyr::anti_join( man_tib2, by=c("AlleleA_ProbeSeq") )
  man_tib1 %>% dplyr::anti_join( man_tib2, by=c("Probe_ID_v1"="Probe_ID") )
  man_tib1 %>% dplyr::anti_join( man_tib2, by=c("Probe_ID_v2"="Probe_ID") )
}

man_tib <- NULL
man_tib <- dplyr::inner_join( 
  man_v1_tib, man_v2_tib, 
  by=c("Probe_ID","AlleleA_ProbeSeq","AlleleB_ProbeSeq","col","mask"), 
  suffix=c("_v1", "_v2") )

man_csv <- file.path( man_dir, "EPIC_v1_v2.gs_to_sesame.csv.gz" )
readr::write_csv( x = man_tib, file = man_csv )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt+3,tc=tc,tt=tt )

sysTime <- Sys.time()
if ( p0 ) cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))


# End of file
