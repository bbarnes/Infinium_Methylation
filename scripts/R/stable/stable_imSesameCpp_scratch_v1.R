
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                Script for Screen Probe Analytical Performance
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
par$prgm_tag <- 'stable_imSesameCpp'
par$verbose  <- 3
local_paths  <- c( 
  "/Users/bbarnes/Documents/tools/imSuite/scripts/R"
)

# local_paths  <- c( 
#   "/Users/bbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
#   "/Users/bretbarnes/Documents/tools/Workhorse-Unstained/scripts/R",
#   "/illumina/scratch/darkmatter/tools/Workhorse-Unstained/scripts/R",
#   "/repo/Workhorse-Unstained/scripts/R" )

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

par$run_name <- "EPICv1"

# run_previous <- TRUE
# run_previous <- FALSE
# if ( run_previous ) {
#   par$run_name <- paste0( par$run_name, ".prev" )
# } else {
#   par$run_name <- paste0( par$run_name, ".bacr" )
# }

opt <- NULL
opt <- imProbeQC_options( pars = par, args = args, vb = par$verbose )
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

opt$rcpp <- 0
opt$rcpp <- 2
opt$rcpp <- 3
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
#                             Pre-processing::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_sub_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_manifests.sub8.rds")
rev_sub_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/rev_manifests.sub8.rds")
neg_ctl_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_negative_ctls.rds" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Negative Controls
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
neg_ctl_tib  <- NULL
neg_ctl_tib  <- readr::read_rds( neg_ctl_rds )
print( neg_ctl_tib )

epic_ctl_tib <- NULL
epic_ctl_rds <- file.path( opt$top_path, "data/manifests/methylation/bgz/epic_ctls.rds" )
epic_ctl_tib <- readr::read_rds( file = epic_ctl_rds )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Sample Sheet
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ssh_tib <- NULL
ssh_csv <- file.path( opt$top_path, "data/sampleSheets/EPIC/EPIC_Auto_SampleSheets_2019-12-29.csv.gz" )
ssh_tib <- readr::read_csv( file = ssh_csv, show_col_types = FALSE ) %>%
  dplyr::filter( PassDetpPoob_Percent_CG >= 90 ) %>%
  dplyr::mutate( Sentrix_Path = paste0( opt$top_path,"/data/idats/EPIC/",Sentrix_Barcode,"/",Sentrix_Name ) )

ssh_vec <- NULL
ssh_vec <- ssh_tib %>% dplyr::pull( Sentrix_Path )  %>% as.vector()

# opt$prefix <- file.path( opt$top_path, "data/idats/EPIC/202296710014/202296710014_R01C01" )
# opt$prefix <- file.path( opt$top_path, "data/idats/EPIC/203806760067/203806760067_R01C01" )
# pre_vec <- c( opt$prefix )

pre_vec <- head( ssh_vec, n=3 )
pre_vec %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Manifest
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
all_man_tib  <- NULL
all_man_tib  <- readr::read_rds( all_sub_rds )
print( all_man_tib )
# all_man_tib %>% dplyr::group_by( Species, Manifest, Manifest_Version ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
# all_man_tib %>% dplyr::filter( Chromosome == "0" ) %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Run-Time Variables
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tc <- 0
vt <- 0

# vb <- 5 # Two level deep (p3)
# vb <- 4 # One level deep (p2)
# vb <- 3 # Standard
vb <- 2 # Light
# vb <- 1 # Min
# vb <- 0 # None

success <- FALSE

opt$run_pair <- FALSE
opt$run_pair <- TRUE

opt$min_pval <- 0.05
opt$min_beta <- 0.3
opt$max_beta <- 0.7
opt$min_perO <- 0.75
opt$min_perI <- 0.05

# Mask Methods:: csJjKk
# workflow_vec = c( "dc",
#                    "ds",
#                    "dJ",
#                    "dj",
#                    "dK",
#                    "dk" )
workflow_vec = c( "cd",
                  "sd",
                  "Jd",
                  "jd",
                  "Kd",
                  "kd" )
workflow_vec = c( "id" )

workflow_vec = c( "i" )

opt$return_df <- 0
opt$return_df <- 1
opt$return_df <- 2
opt$return_df <- 3
opt$return_df <- 4

for ( prefix in pre_vec ) {
  
  pair_df <- NULL
  if ( opt$run_pair )
    pair_df <- read_idat_pair_rcpp( prefix_path  = prefix, 
                                    # prefix_path  = opt$prefix, 
                                    output_path  = opt$out_path,
                                    workflow_vec = workflow_vec %>% as.vector(),
                                    
                                    pval_add_vec = neg_ctl_tib$Address %>% as.vector(), 
                                    addU_man_vec = all_man_tib$U %>% as.vector(),
                                    addM_man_vec = all_man_tib$M %>% as.vector(),
                                    
                                    cgns_man_vec = all_man_tib$Probe_ID %>% as.vector(), 
                                    cols_man_vec = all_man_tib$col %>% as.vector(),
                                    keys_man_vec = all_man_tib$Manifest %>% as.vector(),
                                    anns_man_vec = all_man_tib$Annotation %>% as.vector(), 
                                    chrs_man_vec = all_man_tib$Chromosome %>% as.vector(),
                                    
                                    min_pval     = opt$min_pval,
                                    min_beta     = opt$min_beta,
                                    max_beta     = opt$max_beta,
                                    min_perO     = opt$min_perO, 
                                    min_perI     = opt$min_perI, 
                                    read_bgz     = FALSE,
                                    write_bgz    = FALSE,
                                    return_df    = opt$return_df,
                                    vb=vb,vt=vt,tc=tc ) %>% 
      tibble::as_tibble() %>% 
      dplyr::select( Probe_ID, UG, UR, MG, MR, col, mask_0 ) %>% 
      dplyr::rename( mask = mask_0 ) %>% 
      utils::type.convert(as.is=TRUE) %>% 
      dplyr::mutate(across(where(is.factor), as.character) )
  
  
  if ( !is.null(pair_df) ) {
    cat( glue::glue("\n\npair_df::\n") )
    pair_df %>% print()
    
    #
    # Test Sesame SigDF functionality::
    #
    # sesame::SigDF()
    
    # The check below is only to validate that Idx=0 gives the same values. It does
    #  so this can be removed...
    if ( FALSE ) {
      beta_dif_cnt <- pair_df %>% dplyr::filter( Beta_0 != Beta_0.1 ) %>% base::nrow()
      pval_dif_cnt <- pair_df %>% dplyr::filter( Pval_0 != Pval_0.1 ) %>% base::nrow()
      
      cat("\n")
      cat( glue::glue("{TAB}beta_dif_cnt = {beta_dif_cnt}{RET}") )
      cat( glue::glue("{TAB}pval_dif_cnt = {pval_dif_cnt}{RET}") )
    }
    
    # pair_df %>% dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) ) %>%
    #   dplyr::group_by( Probe_Type ) %>%
    #   dplyr::summarise( Count=n(), .groups = "drop" ) %>% print()
    
    
    # sesame::SigDF( )
    rcpp_tib <- NULL
    rcpp_tib <- pair_df %>% 
      dplyr::select( Probe_ID, UG, UR, MG, MR, col, mask_0 ) %>% 
      dplyr::rename( mask = mask_0 ) %>% 
      utils::type.convert(as.is=TRUE) %>% 
      dplyr::mutate(across(where(is.factor), as.character) )
    
    #
    # Split Data::
    #
    data_tib <- NULL
    data_tib <- rcpp_tib %>% 
      dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") )
    data_sum <- data_tib %>% 
      dplyr::group_by( mask ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    
    #
    # Split Controls::
    #
    ctls_tib <- NULL
    ctls_tib <- rcpp_tib %>% 
      dplyr::filter(  stringr::str_starts( Probe_ID, pattern = "ct") )
    ctls_sum <- ctls_tib %>% 
      dplyr::group_by( mask ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )

    #
    # Build Full SDF::
    #
    full_sdf <- NULL
    full_sdf <- data_tib %>%
      as.data.frame() %>%
      sesame::SigDF( platform = "EPIC", 
                     ctl = ctls_tib
      )
    
    
    #
    # Control Selection::
    #
    ctls_ses_tib <- NULL
    ctls_ses_tib <- rcpp_tib %>%
      sesame::SigDF( platform = "EPIC" ) %>% 
      sesame::controls() %>%
      tibble::as_tibble()
    
    full_ses_sdf <- NULL
    full_ses_sdf <- data_tib %>%
      as.data.frame() %>%
      sesame::SigDF( platform = "EPIC", 
                     ctl = ctls_ses_tib
      )
    
    if ( FALSE ) {
      
      work_ses_tib <- NULL
      work_ses_tib <- mutate_sdf_simple( 
        sdf = full_ses_sdf,
        steps = "DB", 
        negs_min = 1.0, 
        poob_min = 1.0, 
        vb=vb, vt=vt, tc=tc ) %>% 
        tibble::as_tibble()
      
      beta_ses_vec <- NULL
      beta_ses_vec <- mutate_sdf_simple( 
        sdf = work_ses_tib, 
        steps = "v", 
        negs_min = 1.0, 
        poob_min = 1.0, 
        vb=vb, vt=vt, tc=tc )
      beta_ses_vec %>% head() %>% print()
      
      
      detectionPnegEcdf2( sdf = full_ses_sdf, return.pval = TRUE, pval.threshold = 1.0, use_type = FALSE )
      
      funcG <- stats::ecdf( negctls$UG )
      funcR <- stats::ecdf( negctls$UR )
      
      negs_ses_vec <- NULL
      negs_ses_vec <- stats::setNames(BiocGenerics::pmin(
        1-funcR(BiocGenerics::pmax(full_ses_sdf$MR, full_ses_sdf$UR, na.rm=TRUE)),
        1-funcG(BiocGenerics::pmax(full_ses_sdf$MG, full_ses_sdf$UG, na.rm=TRUE))), full_ses_sdf$Probe_ID)
      
      
      negs_ses_vec <- NULL
      negs_ses_vec <- mutate_sdf_simple(
        # sdf = data_tib,
        # sdf = full_ses_sdf,
        # sdf = work_ses_tib,
        steps = "n", 
        negs_min = 1.0, 
        poob_min = 1.0, 
        vb=vb, vt=vt, tc=tc )
      
      
      
      
      #
      #
      # Below is old code...
      #
      #
      #
      
      #
      # Build Sesame SDF::
      #
      full_sdf <- NULL
      full_sdf <- rcpp_tib %>%
        sesame::SigDF( platform = "EPIC", 
                       ctl = epic_ctl_tib
                       # ctl =  ctls_tib
        )
      
      full_tib <- NULL
      full_tib <- full_sdf %>% 
        tibble::as_tibble()
      
    }
    
    if ( FALSE ) {
      
      # TBD:: Run workflow::
      # mutate_sdf( ) or mutate_sdf_simple()
      # mutate_sdf( sdf = pair_sdf, detp = 1.0, work = "XYZ" )
      work_tib <- NULL
      work_tib <- mutate_sdf_simple( 
        sdf = full_sdf,
        steps = "DB", 
        negs_min = 1.0, 
        poob_min = 1.0, 
        vb=vb, vt=vt, tc=tc ) %>% 
        tibble::as_tibble()
      
      beta_vec <- NULL
      beta_vec <- mutate_sdf_simple( 
        sdf = work_tib, 
        steps = "v", 
        negs_min = 1.0, 
        poob_min = 1.0, 
        vb=vb, vt=vt, tc=tc )
      # beta_vec %>% head()
      
      negs_vec <- NULL
      negs_vec <- mutate_sdf_simple( 
        sdf = full_sdf,
        steps = "n", 
        negs_min = 1.0, 
        poob_min = 1.0, 
        vb=vb, vt=vt, tc=tc )
      
      negctls <- NULL
      # negctls <- pair_df %>% sesame::controls() %>% 
      negctls <- full_sdf %>% sesame::controls() %>% 
        tibble::as_tibble() %>% 
        dplyr::mutate( Type = Probe_ID %>% stringr::str_to_upper() ) %>%
        dplyr::filter( Type %>% stringr::str_detect("NEGATIVE") )
      
      negs_vec <- NULL
      negs_vec <- stats::setNames(BiocGenerics::pmin(
        1-funcR(BiocGenerics::pmax(full_sdf$MR, full_sdf$UR, na.rm=TRUE)),
        1-funcG(BiocGenerics::pmax(full_sdf$MG, full_sdf$UG, na.rm=TRUE))), sdf$Probe_ID)
      
      
      # negs_vec <- NULL
      # negs_vec <- mutate_sdf_simple( 
      #   sdf = work_tib %>% dplyr::mutate(
      #     MG = MG %>% as.integer(),
      #     MR = MR %>% as.integer(),
      #     UG = UG %>% as.integer(),
      #     UR = UR %>% as.integer()
      #   ) %>% as.data.frame(),
      #   steps = "n", 
      #   negs_min = 1.0, poob_min = 1.0, 
      #   vb=vb, vt=vt, tc=tc )
      
    }
    
    # detectionPnegEcdf2(            sdf = sdf, return.pval = TRUE,  pval.threshold = negs_min, use_type = FALSE )
  }
  
  break
}
cat( glue::glue("\n\nDONE(test_read_idats): success='{success}'\n\n") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
