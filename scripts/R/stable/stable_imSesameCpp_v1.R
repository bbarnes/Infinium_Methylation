
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

par$run_name <- "Embarkv1"
par$run_name <- "EPICv1"
par$run_name <- "COREv1"

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
# print( neg_ctl_tib )

epic_ctl_tib <- NULL
epic_ctl_rds <- file.path( opt$top_path, "data/manifests/methylation/bgz/epic_ctls.rds" )
epic_ctl_tib <- readr::read_rds( file = epic_ctl_rds )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Manifest
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
all_man_tib  <- NULL
all_man_tib  <- readr::read_rds( all_sub_rds )
# print( all_man_tib )
# all_man_tib %>% dplyr::group_by( Species, Manifest, Manifest_Version ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
# all_man_tib %>% dplyr::filter( Chromosome == "0" ) %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Sample Sheet
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ssh_tib <- NULL
sel_tib <- NULL

if ( par$run_name == "EPICv1" || par$run_name == "COREv1" ) {
  
  # TBD:: Unsupervised or "Suggested Suppervised" Sample Cluster...
  #  - Set Sample_Name from predicted Auto_Sample_BETA_bDB_Key for now...
  #
  ssh_csv <- file.path( opt$top_path, "data/sampleSheets/EPIC/EPIC_Auto_SampleSheets_2019-12-29.csv.gz" )
  ssh_tib <- readr::read_csv( file = ssh_csv, show_col_types = FALSE ) %>%
    dplyr::mutate( Sentrix_Path = paste0( opt$top_path,"/data/idats/EPIC/",Sentrix_Barcode,"/",Sentrix_Name ) ) %>%
    dplyr::mutate( Sample_Name = Auto_Sample_BETA_bDB_Key )
  
  # NOTE: Only Using Replicate Samples for now...
  #
  sel_tib <- ssh_tib %>% 
    dplyr::filter( Sample_Name == "HELA" |
                     Sample_Name == "JURKAT" |
                     Sample_Name == "MCF7" |
                     Sample_Name == "RAJI" ) %>%
    dplyr::filter( PassDetpPoob_Percent_CG >= 90 )
  
  if ( par$run_name == "COREv1" ) {
    #
    # Get all EPIC v1 DELTA Core/Fail Samples Seach::
    #
    epic_core_path <- file.path( opt$top_path, "data/idats/idats_EPIC-8x1-DELTA-Core" )
    epic_fail_path <- file.path( opt$top_path, "data/idats/idats_EPIC-8x1-DELTA-Fail" )
    
    epic_core_list <- NULL
    epic_core_list <- sesame::searchIDATprefixes( dir.name = epic_core_path, recursive = TRUE )
    epic_core_cnts <- epic_core_list %>% length()
    
    epic_fail_list <- NULL
    epic_fail_list <- sesame::searchIDATprefixes( dir.name = epic_fail_path, recursive = TRUE )
    epic_fail_cnts <- epic_fail_list %>% length()
    
    sel_tib <- dplyr::bind_rows(
      ssh_tib %>% dplyr::filter( Sentrix_Name %in% names(epic_core_list) ) %>%
        dplyr::mutate( Select_Group = "Pass" ),
      ssh_tib %>% dplyr::filter( Sentrix_Name %in% names(epic_fail_list) ) %>%
        dplyr::mutate( Select_Group = "Fail" ),
    )
    
  }

  # opt$prefix <- file.path( opt$top_path, "data/idats/EPIC/202296710014/202296710014_R01C01" )
  # opt$prefix <- file.path( opt$top_path, "data/idats/EPIC/203806760067/203806760067_R01C01" )
  # pre_vec <- c( opt$prefix )
  
  #
  # TBD:: Load Sesame Masked Probes
  #
  # epic_ses_dat <- sesameData::sesameDataGet( title = "EPIC.probeInfo" )
  
} else if ( par$run_name == "Embarkv1" ) {
  
  #
  # Embark::
  #
  #  "Projects.new/Embark/data/raw/3-Experiment_TechnicalReplicates"
  #  "Projects.new/Embark/data/raw/7-Experiment_MethylationTitration"
  #
  
  # ssh_tib <- NULL
  # ssh_csv <- file.path( opt$top_path, "Projects.new/Embark/Manifest/Manifest_v3_140920222/Embark-UCSC-v4/auto_sample_sheets/mask0/mask0.auto_sample_sheet.csv.gz" )
  # ssh_tib <- readr::read_csv( file = ssh_csv, show_col_types = FALSE ) %>%
  #   dplyr::filter( Pass_Perc >= 80 ) %>%
  #   dplyr::mutate(
  #     Sentrix_Barcode = Sentrix_Name %>% stringr::str_remove( "_.*$"),
  #     Sentrix_Path = paste0( opt$top_path,"/data/idats/EPIC/",Sentrix_Barcode,"/",Sentrix_Name )
  #   )
  
  #
  # TBD:: Split Replicate and Titration Data or label them...
  #  - Filtering out MRC5
  #  - CLeaning up Sample_Names
  #
  ssh_tib <- NULL
  ssh_csv <- file.path( opt$top_path, "Projects.new/Embark/data/sampleSheets/sampleSheet.human.v1.csv" )
  ssh_tib <- readr::read_csv( file = ssh_csv, show_col_types = FALSE ) %>% 
    dplyr::filter( Experiment_Folder %>% stringr::str_starts("3-") ) %>% 
    # dplyr::filter( Experiment_Folder %>% stringr::str_starts("3-") | 
    #                  Experiment_Folder %>% stringr::str_starts("7-") ) %>% 
    dplyr::mutate( Sentrix_Name = paste(Sentrix_ID,Sentrix_Position, sep="_") ) %>%
    dplyr::filter( Sample_Name %>% stringr::str_detect( "_pooled_") ) %>% 
    dplyr::mutate( Sample_Name = Sample_Name %>% stringr::str_remove("_pooled_.*$") )
  
  #
  # TBD:: Rename sel_tib to something else... And retain sample information...
  # TBD:: Only using Roco for testing...
  #
  sel_tib <- NULL
  sel_dir <- file.path( opt$top_path, "data/idats/idats_Embark_v1" )
  sel_tib <- sesame::searchIDATprefixes( dir.name = sel_dir ) %>% 
    cbind() %>% as.data.frame() %>% 
    tibble::rownames_to_column( var = "Sentrix_Name" ) %>% 
    tibble::as_tibble() %>% 
    magrittr::set_names( c("Sentrix_Name","Sentrix_Path") ) %>%
    dplyr::inner_join( ssh_tib, by = c("Sentrix_Name") )
  
  # Single Sample Testing...
  # sel_tib <- sel_tib %>% 
  #   dplyr::filter( Sample_Name == "Roco" )
  
} else {
  stop( glue::glue("{perrs} Failed to find idats! Exiting...{RET}") )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Sample Date Summary
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
date_decode_ssh_sum <- NULL
date_decode_ssh_sum <- ssh_tib %>% 
  dplyr::group_by( iscan_Decoding_Year ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p4 ) date_decode_ssh_sum %>% print( n=base::nrow(date_decode_ssh_sum) )

date_decode_sel_sum <- NULL
date_decode_sel_sum <- sel_tib %>% 
  dplyr::group_by( iscan_Decoding_Year ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p4 ) date_decode_sel_sum %>% print( n=base::nrow(date_decode_sel_sum) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Sample Count Summary
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ssh_sum <- NULL
ssh_sum <- ssh_tib %>%
  dplyr::group_by( Sample_Name ) %>% 
  dplyr::summarise( Tot_Poob=n(),
                    Min_Poob=min( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                    Max_Poob=max( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                    Avg_Poob=mean( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                    Sds_Poob=sd( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                    Med_Poob=median( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                    Mad_Poob=mad( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                    .groups = "drop" )
if ( p1 ) ssh_sum %>% print( n=base::nrow(ssh_sum) )

#
# TBD:: Use:: "PassDetpNegs_Percent_CG"
#
# sel_sum <- NULL
# sel_sum <- sel_tib %>%
#   dplyr::filter( iscan_Decoding_Year > 2017 ) %>%
#   dplyr::group_by( Sample_Name, iscan_Decoding_Year ) %>% 
#   dplyr::summarise( Tot_Poob=n(),
#                     Min_Poob=min( PassDetpPoob_Percent_CG, na.rm = TRUE ),
#                     Max_Poob=max( PassDetpPoob_Percent_CG, na.rm = TRUE ),
#                     Avg_Poob=mean( PassDetpPoob_Percent_CG, na.rm = TRUE ),
#                     Sds_Poob=sd( PassDetpPoob_Percent_CG, na.rm = TRUE ),
#                     Med_Poob=median( PassDetpPoob_Percent_CG, na.rm = TRUE ),
#                     Mad_Poob=mad( PassDetpPoob_Percent_CG, na.rm = TRUE ),
#                     .groups = "drop" )
# if ( p1 ) sel_sum %>% print( n=base::nrow(sel_sum) )

if ( par$run_name == "COREv1" ) {
  
  sel_sum <- NULL
  sel_sum <- sel_tib %>%
    dplyr::group_by( Sample_Name,Select_Group ) %>% 
    dplyr::summarise( Tot_Poob=n(),
                      Min_Poob=min( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Max_Poob=max( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Avg_Poob=mean( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Sds_Poob=sd( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Med_Poob=median( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Mad_Poob=mad( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      .groups = "drop" )
  if ( p1 ) sel_sum %>% print( n=base::nrow(sel_sum) )

}


#
# TBD:: LEFT OFF HERE
#  - Point is we need to select only the DELTA EPIC GOOD Sentrix_Names...
#  - All of them should be found in the EPIC directory, also upload this to BOX or share its...
#
# sel_tib %>% dplyr::mutate( Sentrix_Path_Group = Sentrix_Path %>% stringr::str_remove("^/Users/bbarnes/Documents/data/idats/") %>% stringr::str_remove("\\/.*$") ) %>% dplyr::group_by( Sentrix_Path_Group ) %>% dplyr::summarise( Count=n(), .groups = "drop" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Run-Time Variables
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

tc <- 0
vt <- 0

# vb <- 5 # Two level deep (p3)
# vb <- 4 # One level deep (p2)
# vb <- 3 # Standard
# vb <- 2 # Light
# vb <- 1 # Min
vb <- 0 # None

success <- FALSE

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
# opt$return_df <- 4

opt$max_sam <- 4
opt$max_rep <- 8

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Data-processing:: Replicate Sample List
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
rep_data_tib <- NULL

rep_sel_list <- NULL
rep_sel_list <- sel_tib %>% split( .$Sample_Name ) %>%
  head( n=opt$max_sam )
for ( sample_name in names(rep_sel_list) ) {
  all_beta_tib <- NULL
  all_poob_tib <- NULL
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Data-processing:: Sentrix Sample List
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  ssh_sel_list <- NULL
  ssh_sel_list <- rep_sel_list[[sample_name]] %>% 
    dplyr::mutate( Rep_Idx = dplyr::row_number(),
                   Rep_Key = paste( Sample_Name,Rep_Idx, sep="_") ) %>% 
    split( .$Sentrix_Name ) %>% head( n=opt$max_rep )

  for ( sentrix_name in names(ssh_sel_list) ) {
    prefix <- ssh_sel_list[[sentrix_name]]$Sentrix_Path
    if ( p1 ) cat(glue::glue("{pmssg} Current prefix = '{prefix}'{RET}"))
    
    # }
    # for ( prefix in pre_vec ) {
    #   if ( p1 ) cat(glue::glue("{pmssg} Current prefix = '{prefix}'{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Data-processing:: Loading
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    idat_tib <- NULL
    idat_tib <- read_idat_pair_rcpp( prefix_path  = prefix,
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
      # dplyr::select( Probe_ID, UG, UR, MG, MR, col, mask_0 ) %>% 
      # dplyr::rename( mask = mask_0 ) %>% 
      utils::type.convert(as.is=TRUE) %>% 
      dplyr::mutate(across(where(is.factor), as.character) ) %>%
      dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )
    
    #
    # TBD:: Return select columns with appropriate data types in Rcpp...
    #   - This will allow the code above to be removed!
    #
    
    if ( is.null(idat_tib) ) {
      stop( glue::glue("{perrs} Failed to load idats: prefix='{prefix}'{RET}") )
      break
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Data-processing:: Summary
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    idat_sum <- NULL
    idat_sum <- idat_tib %>% 
      dplyr::group_by( Probe_Type ) %>%
      dplyr::summarise( Count=n(), .groups = "drop" )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Data-processing:: Parse Analytical Data
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    data_tib <- NULL
    data_tib <- idat_tib %>% 
      dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") )
    data_sum <- data_tib %>% 
      dplyr::group_by( mask ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Data-processing:: Parse Control Data
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    ctls_tib <- NULL
    ctls_tib <- idat_tib %>% 
      dplyr::filter(  stringr::str_starts( Probe_ID, pattern = "ct") )
    ctls_sum <- ctls_tib %>% 
      dplyr::group_by( mask ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Data-processing:: Workflow Pipeline
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    work_sdf <- NULL
    work_sdf <- data_tib %>%
      as.data.frame() %>%
      sesame::SigDF( platform = "EPIC", 
                     ctl = ctls_tib
      ) %>%
      mutate_sdf_simple( 
        steps = "DB", 
        negs_min = 1.0, 
        poob_min = 1.0, 
        vb=vb, vt=vt, tc=tc )

    rep_key <- ssh_sel_list[[sentrix_name]]$Rep_Key
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Data-processing:: Extract Beta Values
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    beta_tib <- NULL
    beta_tib <- mutate_sdf_simple( 
      sdf = work_sdf,
      steps = "v", 
      negs_min = 1.0, 
      poob_min = 1.0, 
      vb=vb, vt=vt, tc=tc ) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column( var = "Probe_ID") %>% 
      tibble::as_tibble() %>% 
      magrittr::set_names( value = c("Probe_ID",rep_key) )
      # magrittr::set_names( value = c("Probe_ID","Beta") )
    
    if ( is.null( all_beta_tib ) ) {
      all_beta_tib <- beta_tib
    } else {
      all_beta_tib <- all_beta_tib %>% 
        inner_join( beta_tib, by = c("Probe_ID") )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Data-processing:: Extract pooBAH P-Values
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    poob_tib <- NULL
    poob_tib <- mutate_sdf_simple( 
      sdf = work_sdf,
      steps = "O", 
      negs_min = 1.0, 
      poob_min = 1.0, 
      vb=vb, vt=vt, tc=tc ) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column( var = "Probe_ID") %>% 
      tibble::as_tibble() %>% 
      magrittr::set_names( value = c("Probe_ID",rep_key) )
      # magrittr::set_names( value = c("Probe_ID","Poob") )
    
    if ( is.null( all_poob_tib ) ) {
      all_poob_tib <- poob_tib
    } else {
      all_poob_tib <- all_poob_tib %>% 
        inner_join( poob_tib, by = c("Probe_ID") )
    }
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #
    #                   Data-processing:: Join Data & Output
    #
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    
    #
    # [Done Above]:: Join betas and pvals as a table
    #
    #  dplyr::inner_join( beta_tib, poob_tib, by=c("Probe_ID") )
    
  }
  if ( p0 ) cat(glue::glue("{pmssg} Finished Processing Idat Pairs sample='{sample_name}'.{RET2}"))

  #
  # Simple dB Score Method
  #

  prd_mat <- NULL
  min_mat <- NULL
  dbs_mat <- NULL
  
  beta_mat <- NULL
  beta_mat <- all_beta_tib %>% 
    tibble::column_to_rownames( var = "Probe_ID" ) %>% 
    as.matrix()  

  poob_mat <- NULL
  poob_mat <- all_poob_tib %>% 
    tibble::column_to_rownames( var = "Probe_ID" ) %>% 
    as.matrix()
  
  #
  # LEFT OFF HERE::
  #
  #. - Need to implement dual input matricies
  #
  
  ncol_cnt <- ncol( beta_mat )
  for ( ii in c(1:ncol_cnt) ) {
    for ( jj in c(1:ncol_cnt) ) {
      if ( ii >= jj ) next
      
      if ( p1 ) cat(glue::glue("{pmssg} {sample_name}: {ii} x {jj}.{RET}"))
      
      dB_vec <- 1 - base::abs( beta_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() )
      dP_vec <- 1 - poob_mat[ ,c(ii,jj) ] %>% matrixStats::rowMaxs()
      
      # cbind( dB_vec, dP_vec ) %>% head() 
      
      #
      # 1. Product vector to be rowSum()
      # 2. dP Min vector to be rowSum()
      # 3. cbind( 1,2 ) and get the quotient
      #
      
      # Product Vector::
      prd_mat <- cbind( prd_mat,  cbind( dB_vec, dP_vec ) %>% matrixStats::rowProds() )
      min_mat <- cbind( min_mat, dP_vec )
      dbs_mat <- cbind( dbs_mat, dB_vec )
      
      # dB_pass_mat <- dB_pass_mat &
      #   beta_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() < min_dB
      
      # break
    }
    
    # break
  }
  
  dat_tib <- NULL
  dat_tib <- tibble::tibble(
    Sample = sample_name,
    Probe_ID = all_beta_tib$Probe_ID,
    sum_scr = prd_mat %>% matrixStats::rowSums2(),
    sum_det = min_mat %>% matrixStats::rowSums2(),
    avg_dbs = dbs_mat %>% matrixStats::rowMeans2(),
    med_dbs = dbs_mat %>% matrixStats::rowMedians()
  ) %>% dplyr::mutate(
    fin_scr = sum_scr / sum_det
  )

  if ( FALSE ) {
    dat_tib %>% ggplot2::ggplot( aes( x=fin_scr ) ) +
      ggplot2::geom_density()
    
    dat_tib %>% ggplot2::ggplot( aes( x=fin_scr, y=avg_dbs ) ) +
      ggplot2::geom_point( )
    
    dat_tib %>% ggplot2::ggplot( aes( x=fin_scr, y=avg_dbs ) ) +
      ggplot2::geom_density2d()
    
    dat_tib %>% ggplot2::ggplot( aes( x=fin_scr, y=med_dbs ) ) +
      ggplot2::geom_density2d()
  }
  
  rep_data_tib <- rep_data_tib %>% dplyr::bind_rows( dat_tib )
  
  if ( opt$single ) break
}

rep_data_sum <- rep_data_tib %>% 
  dplyr::group_by( Sample ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

# Looks good::
rep_den_pdf <- file.path( opt$out_path, paste(opt$run_name,"rep_density.pdf", sep=".") )
rep_den_gg <- rep_data_tib %>% 
  ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
  ggplot2::geom_density( alpha = 0.2 ) +
  ggplot2::facet_grid( rows = vars(Sample) )

ggplot2::ggsave( filename = rep_den_pdf, 
                 device = "pdf", 
                 width = 7, 
                 height = 7, 
                 dpi = 320 )

# More or less worthless::
rep_2den_gg <- rep_data_tib %>% 
  ggplot2::ggplot( aes( x=fin_scr, y=med_dbs ) ) +
  ggplot2::geom_density2d() +
  ggplot2::facet_grid( rows = vars(Sample) )

if ( par$run_name == "COREv1" ) {
  
  # Need to split COREv1 into Select_Group above...
  
  # rep_mask_tib <- rep_data_tib %>% 
  #   dplyr::left_join( epic_mask_tib, by=c("Probe_ID") ) %>% 
  #   dplyr::mutate(
  #     Masked = dplyr::case_when(
  #       is.na(Masked) ~ FALSE,
  #       TRUE ~ TRUE )
  #   )
  
}

if ( par$run_name == "EPICv1" ) {
  
  # rep_data_tib %>% dplyr::filter( Probe_ID %in% epic_ses_dat$mask )
  
  epic_mask_tib <- tibble::tibble( Probe_ID = epic_ses_dat$mask, Masked = TRUE )
  rep_mask_tib <- rep_data_tib %>% 
    dplyr::left_join( epic_mask_tib, by=c("Probe_ID") ) %>% 
    dplyr::mutate(
      Masked = dplyr::case_when(
        is.na(Masked) ~ FALSE,
        TRUE ~ TRUE )
    )
  
  rep_mask_sum <- rep_mask_tib %>% 
    dplyr::group_by( Masked ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )

  # Quick Validation::
  #
  # rep_mask_tib %>% 
  #   dplyr::group_by( Masked ) %>% 
  #   dplyr::summarise( Count=n(), .groups = "drop" )
  
  sn <- 0.000001
  
  # Apply Cutoffs/Log
  rep_den_mask_gg <- rep_mask_tib %>% 
    dplyr::filter( fin_scr >= 0.9 ) %>%
    # ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
    ggplot2::ggplot( aes( x=log(fin_scr) + sn, color = Sample, fill = Sample ) ) +
    ggplot2::geom_density( alpha = 0.2 )  +
    ggplot2::facet_grid( rows = vars(Masked) )
  
  # Conclusion:: Score does show a difference for Replicate score vs. third
  #  party analysis...
  
  rep_den_mask_pdf <- file.path( opt$out_path, paste(opt$run_name,"rep_den_mask.pdf", sep=".") )
  rep_den_mask_gg <- rep_mask_tib %>% 
    dplyr::filter( fin_scr >= 0.9 ) %>%
    ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
    # ggplot2::ggplot( aes( x=log(fin_scr) + sn, color = Sample, fill = Sample ) ) +
    ggplot2::geom_density( alpha = 0.2 )  +
    ggplot2::facet_grid( rows = vars(Masked),
                         cols = vars(Sample) )
  
  ggplot2::ggsave( filename = rep_den_mask_pdf, device = "pdf", width = 7, height = 7, dpi = 320 )

}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
