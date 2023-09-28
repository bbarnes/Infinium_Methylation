
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                Swifthoof:: 
#                Analyze idats <= Sesame => to pval/beta
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
par$prgm_tag <- 'swifthoof'
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

par$local_run_type <- "AKE-MVP-Failed-v1"
par$local_run_type <- "EX_CAM1"
par$local_run_type <- "NA12878"

opt <- NULL
opt <- swifthoof_local_defaults( pars = par, args = args, vb = par$verbose )
vb  <- opt$verbose

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                         Program Initialization::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

par_reqs <- c( 'run_mode',
               'src_path', 'scr_path', 'exe_path', 'prgm_dir', 'prgm_tag' )
opt_reqs <- c( 'out_path', # 'ref_path', 'ref_file', 'ref_build', 'ref_species',
               'Rscript', 'verbose' )

#
# TBD:: Update docker defaults after building new docker branch::
# TBD:: Add auxilary files to to program_init check::
#
prgm_dat <- program_init( name = par$prgm_tag,
                          opts = opt, opt_reqs = opt_reqs,
                          pars = par, par_reqs = par_reqs, 
                          rcpp = 0,
                          vb = vb, vt=3, tc=0 )

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

pval_vec     <- split_str_to_vec(opt$pval)
min_pval_vec <- split_str_to_vec(opt$minPval)
min_perc_vec <- split_str_to_vec(opt$minPerc)
workflow_vec <- split_str_to_vec(opt$workflow)

# Remove "r/raw" and force "raw" to be first::
workflow_vec <- c("raw",workflow_vec[!workflow_vec %in% c("r","raw")])

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Done. Parsing Inputs.{RET}{RET}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Select Chips from idats and/or Target Manifest::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

chipPrefixes <- NULL
chipPrefixes <- sesame::searchIDATprefixes( opt$idat_path )
sampleCounts <- chipPrefixes %>% names() %>% length()

if (is.null(chipPrefixes) || length(chipPrefixes)==0)
  stop(glue::glue("{RET}[{par$prgmTag}]: chipPrefixes is null or length=0!!!{RET2}"))

if (opt$verbose>0)
  cat(glue::glue("[{par$prgmTag}]: Found sample counts={sampleCounts}!{RET2}"))

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                  Main::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                             Load Manifest(s)::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( is.null(opt$manifest) ) {
  if ( is.null( opt$aux_man_path) ) {
    
  }
}

sesameData::sesameDataList()
ses_450k_grs <- sesameData::sesameDataGet("HM450.hg19.manifest")
ses_450k_tib <- ses_450k_grs %>% as.data.frame() %>% tibble::as_tibble( rownames = "Loci_ID" )
ses_450k_tib %>% dplyr::group_by( probeType ) %>% dplyr::summarise( Count = n() )

ses_epic_grs <- sesameData::sesameDataGet("EPIC.hg19.manifest")
ses_epic_tib <- ses_epic_grs %>% as.data.frame() %>% tibble::as_tibble( rownames = "Probe_ID" )

ses_epic_tib %>% dplyr::select(Probe_ID, address_A, address_B, designType, channel, probeType, nextBase )
# Probe_ID,M,U,DESIGN,COLOR_CHANNEL,col,Probe_Type,Probe_Source,Next_Base,Probe_Design

ses_epic_tib %>% dplyr::group_by( probeType ) %>% dplyr::summarise( Count = n() )

ses_epic_probeInfo <- sesameData::sesameDataGet("EPIC.probeInfo")

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("[{par$prgmTag}]: Finished(time={sysTime}){RET}{RET}"))

# End of file
