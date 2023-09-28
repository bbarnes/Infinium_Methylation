
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Script for Parsing Genomes
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                              Source Packages::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rm(list=ls(all=TRUE))

# Options, Tidy Practices and Parallel Computing Packages::
suppressWarnings(suppressPackageStartupMessages( 
  base::require("optparse",   quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages( 
  base::require("tidyverse",  quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("doParallel", quietly = TRUE) ) )

# suppressWarnings(suppressPackageStartupMessages( 
#   base::require("microbenchmark", quietly = TRUE) ) )

# Load Plotting Packages::
# suppressWarnings(suppressPackageStartupMessages(
#   base::require("ggplot2", quietly = TRUE) ) )
# suppressWarnings(suppressPackageStartupMessages(
#   base::require("GGally", quietly = TRUE) ) )

# BiocManager::install("sesame")

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
par$prgm_tag <- 'stable_evonik_control_search'
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
par <- source_functions( pars = par, rcpp = FALSE, vb = par$verbose )
par <- params_check( pars = par, args = args, 
                     prgm_aux_check = FALSE, vb = par$verbose )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                           Get Program Options::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par$version <- 0
# par$version <- 1
# par$version <- 2
# par$version <- 3
# par$version <- 4
# par$version <- 5
# par$version <- 6

par$run_name <- "evonik"

opt <- NULL
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
opt_reqs <- c( 'out_path', 
               'Rscript', 'verbose' )

opt$rcpp <- 2
opt$rcpp <- 0
prgm_dat <- program_init( name = par$prgm_tag,
                          opts = opt, opt_reqs = opt_reqs,
                          pars = par, par_reqs = par_reqs, 
                          rcpp = opt$rcpp,
                          vb = opt$verbose, vt=3, tc=0 )

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

success = TRUE;

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                          Load Controls Tangos::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cntl_cols <- NULL
cntl_cols <- readr::cols(
  Address = readr::col_integer(),
  Type = readr::col_character(),
  Name = readr::col_character()
)

cntl_data_tib <- NULL
cntl_data_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/infinium-methylationepic-v-1-0-b5-manifest-file.controls-only-3.csv.gz" )
cntl_data_tib <- readr::read_csv( file = cntl_data_csv,
                                  col_names = names(cntl_cols$cols), 
                                  col_types = cntl_cols )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                Load iDats::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Load iDats:: Evonik
evon_idat_list <- NULL
evon_idat_path <- file.path( opt$top_path, "Projects.new/Evonik/Validation" )
evon_idat_list <- sesame::searchIDATprefixes( dir.name = evon_idat_path )
evon_idat_size <- evon_idat_list %>% length()
if ( p0 ) cat(glue::glue("{pmssg} Evonk Idats Size: {evon_idat_size}.{RET}"))


# Load iDats:: EPIC v2
epic_idat_list <- NULL
epic_idat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/idats" )
epic_idat_list <- sesame::searchIDATprefixes( dir.name = epic_idat_path )
epic_idat_size <- epic_idat_list %>% length()
if ( p0 ) cat(glue::glue("{pmssg} EPIC Idats Size: {epic_idat_size}.{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Search for Control Tangos::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

cntl_cnts_tib <- NULL

idat_src_key <- "Evonik"
for ( prefix in names(evon_idat_list) ) {
  if ( p0 ) cat(glue::glue("{pmssg} Evonk prefix: {prefix}.{RET}"))
  
  grn_idat_path <- NULL
  grn_idat_path <- file.path( paste0(evon_idat_list[[prefix]], "_Grn.idat" ) )
  
  grn_idat_dat <- NULL
  grn_idat_dat <- illuminaio::readIDAT( file = grn_idat_path )
  grn_idat_tib <- NULL
  grn_idat_tib <- grn_idat_dat$Quants %>% as.data.frame() %>% 
    tibble::rownames_to_column( var = "Address" ) %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate( Address = Address %>% as.integer() )
  
  grn_join_tib <- NULL
  grn_join_tib <- cntl_data_tib %>% 
    dplyr::inner_join( grn_idat_tib, by = "Address" )
  
  grn_cntl_cnt <- NULL
  grn_cntl_cnt <- grn_join_tib %>% base::nrow()
  
  grn_zero_cnt <- NULL
  grn_zero_cnt <- grn_join_tib %>% dplyr::filter( Mean == 0 ) %>% base::nrow()
  
  #
  # Red Idat Analysis::
  #
  red_idat_path <- NULL
  red_idat_path <- file.path( paste0(evon_idat_list[[prefix]], "_Red.idat" ) )
  
  red_idat_dat <- NULL
  red_idat_dat <- illuminaio::readIDAT( file = red_idat_path )
  red_idat_tib <- NULL
  red_idat_tib <- red_idat_dat$Quants %>% as.data.frame() %>% 
    tibble::rownames_to_column( var = "Address" ) %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate( Address = Address %>% as.integer() )
  
  red_join_tib <- NULL
  red_join_tib <- cntl_data_tib %>% 
    dplyr::inner_join( red_idat_tib, by = "Address" )
  
  red_cntl_cnt <- NULL
  red_cntl_cnt <- red_join_tib %>% base::nrow()
  
  red_zero_cnt <- NULL
  red_zero_cnt <- red_join_tib %>% dplyr::filter( Mean == 0 ) %>% base::nrow()
  
  #
  # Join Data::
  #
  cnt_tib <- NULL
  cnt_tib <- tibble::tibble(
    Source = idat_src_key,
    Grn_Cnt = grn_cntl_cnt,
    Red_Cnt = red_cntl_cnt
  )
  
  break
  
  cntl_cnts_tib <- cntl_cnts_tib %>%
    dplyr::bind_rows( cnt_tib )
}
cntl_cnts_tib %>% print( n=base::nrow(cntl_cnts_tib) )

idat_src_key <- "EPIC"
for ( prefix in names(epic_idat_list) ) {
  if ( p0 ) cat(glue::glue("{pmssg} EPIC prefix: {prefix}.{RET}"))
  
  grn_idat_path <- NULL
  grn_idat_path <- file.path( paste0(epic_idat_list[[prefix]], "_Grn.idat.gz" ) )
  
  grn_idat_dat <- NULL
  grn_idat_dat <- illuminaio::readIDAT( file = grn_idat_path )
  grn_idat_tib <- NULL
  grn_idat_tib <- grn_idat_dat$Quants %>% as.data.frame() %>% 
    tibble::rownames_to_column( var = "Address" ) %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate( Address = Address %>% as.integer() )
  
  grn_cntl_cnt <- NULL
  grn_cntl_cnt <- cntl_data_tib %>% 
    dplyr::inner_join( grn_idat_tib, by = "Address" ) %>%
    base::nrow()
  
  red_idat_path <- NULL
  red_idat_path <- file.path( paste0(epic_idat_list[[prefix]], "_Red.idat.gz" ) )
  
  red_idat_dat <- NULL
  red_idat_dat <- illuminaio::readIDAT( file = red_idat_path )
  red_idat_tib <- NULL
  red_idat_tib <- red_idat_dat$Quants %>% as.data.frame() %>% 
    tibble::rownames_to_column( var = "Address" ) %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate( Address = Address %>% as.integer() )
  
  red_cntl_cnt <- NULL
  red_cntl_cnt <- cntl_data_tib %>% 
    dplyr::inner_join( red_idat_tib, by = "Address" ) %>%
    base::nrow()
  
  cnt_tib <- NULL
  cnt_tib <- tibble::tibble(
    Source = idat_src_key,
    Grn_Cnt = grn_cntl_cnt,
    Red_Cnt = red_cntl_cnt
  )
  
  cntl_cnts_tib <- cntl_cnts_tib %>%
    dplyr::bind_rows( cnt_tib )
}
cntl_cnts_tib %>% print( n=base::nrow(cntl_cnts_tib) )

cntl_cnts_tib %>% dplyr::group_by( Source,Grn_Cnt,Red_Cnt ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
# A tibble: 3 × 4
# Source Grn_Cnt Red_Cnt Count
# <chr>    <int>   <int> <int>
# 1 EPIC       633     633    40
# 2 EPIC       635     635    16
# 3 Evonik     623     623  1140

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
