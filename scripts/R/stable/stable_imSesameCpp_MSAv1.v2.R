
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
par$prgm_tag <- 'stable_imSesameCpp_MSAv1'
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
par$version <- 1
par$version <- 2
par$version <- 3
par$version <- 4
# par$version <- 5
# par$version <- 6

# par$run_name <- "EPICv1"
# par$run_name <- "Embarkv1"
# par$run_name <- "FAILv1"
# par$run_name <- "COREv1"

par$run_name <- "MSAv03"
par$run_name <- "MSAv10"

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
#             Pre-processing:: Search for Idats (MSA v.1.0)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

msa10_idat_list <- NULL
msa10_idat_path <- file.path( opt$top_path, "data/pre-idats/MSA" )
msa10_idat_list <- sesame::searchIDATprefixes( dir.name = msa10_idat_path )

epic2_idat_list <- NULL
epic2_idat_path <- file.path( opt$top_path, "data/idats/idats_EPIC_v2-20220912-Alpha_subset" )
epic2_idat_list <- sesame::searchIDATprefixes( dir.name = epic2_idat_path )

# A tibble: 655 × 2
idat_path_tib <- NULL
idat_path_tib <- c( msa10_idat_list,epic2_idat_list ) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column( var = "Sentrix_Name") %>% 
  magrittr::set_names( c("Sentrix_Name", "Sentrix_Path") ) %>%
  tibble::as_tibble() %>%
  dplyr::distinct( Sentrix_Name, .keep_all = TRUE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Pre-processing:: Sample Sheet (MSA.v.1.0)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

msa_ssh_tib <- NULL
msa_ssh_tib <- format_msa_sshs( 
  paths = c( file.path( opt$top_path, "data/pre-idats/MSA/sample_sheets/rep0"),
             file.path( opt$top_path, "data/pre-idats/MSA/sample_sheets/uhm0"),
             file.path( opt$top_path, "data/pre-idats/MSA/sample_sheets/uhm1") ),
  
  out_dir    = file.path( opt$out_path, "sample_sheets" ),
  run_tag    = opt$run_name,
  reload     = opt$reload,
  reload_min = 10,
  ret_data   = FALSE,
  parallel   = opt$parallel,
  write_out  = FALSE,
  
  vb=vb,vt=vt+1,tc=tc+1, tt=tt )

msa_ssh_sum <- NULL
msa_ssh_sum <- msa_ssh_tib %>%
  print_sum( vec = c("Sample_Group","Source_Key","Source_Name","Sample_Base",
                     "Sample_Input","Sample_Titration","User_Format"),
             vb=vb,vt=vt+1,tc=tc+1, tt=tt )

# [TBD]: Update Source_Key
# [TBD]: Add Platform (or can we just preeict it), look at Rccp output...

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Pre-processing:: Sample Sheets (EPICv1/v2)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

epi_ssh_tib <- NULL
epi_ssh_tib <- format_epic_ssh_alpha( 
  file = file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/SampleSheets/formatted/EPICv2-UCSC-v0.LightningAuto.select.sample_sheet.csv.gz" ),
  write_out  = FALSE,
  vb=vb,vt=vt+3,tc=tc+1, tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Pre-processing:: Sample Sheets (Combine All)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# [TBD]: Where did these go???
#
# A tibble: 25 × 8
ssh_mis_tib <- NULL
ssh_mis_tib <- dplyr::bind_rows( epi_ssh_tib,msa_ssh_tib ) %>% 
  dplyr::distinct( Sentrix_Name, .keep_all = TRUE ) %>%
  dplyr::anti_join( idat_path_tib, 
                    by=c("Sentrix_Name") ) # %>% as.data.frame()

# Join and Match::
ssh_tib <- NULL
ssh_tib <- dplyr::bind_rows( epi_ssh_tib,msa_ssh_tib ) %>% 
  dplyr::distinct( Sentrix_Name, .keep_all = TRUE ) %>%
  dplyr::inner_join( idat_path_tib, 
                     by=c("Sentrix_Name") ) # %>% as.data.frame()

ssh_sum0 <- NULL
ssh_sum0 <- ssh_tib %>%
  print_sum( vec = c("Sample_Base","Source_Key",
                     "Sample_Input","Sample_Titration") ,
             vb=vb,vt=vt+1,tc=tc+1, tt=tt )

ssh_sum1 <- NULL
ssh_sum1 <- ssh_tib %>%
  print_sum( vec = c("Sample_Group","Sample_Base",
                     "Sample_Input","Sample_Titration") ,
             vb=vb,vt=vt+3,tc=tc+1, tt=tt )

ssh_sum2 <- NULL
ssh_sum2 <- ssh_tib %>%
  print_sum( vec = c("Sample_Group","Sample_Base",
                     "Sample_Input","Sample_Titration",
                     "Source_Key","Source_Name" ),
             vb=vb,vt=vt+3,tc=tc+1, tt=tt )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#             [TBD] Pre-processing:: Build Sample Experiments
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Only T0,T1,R0,E2
# Only >= 50 ng
Sample_Input_min <- 50
Source_Key_vec <- c( "T0","T1","R0","E2" )
base_ssh_tib <- NULL
base_ssh_tib <- ssh_tib %>% 
  dplyr::filter( Source_Key %in% Source_Key_vec ) %>% 
  dplyr::filter( Sample_Input >= Sample_Input_min ) %>%
  dplyr::arrange( Sample_Titration,Sample_Input,Sample_Base )

# WHAT???
ssh_tib %>% dplyr::filter( Sample_Group == "MeTitration" )

#
# - Sample_Group = [  MeTitration ] Titration
uhm_base_ssh_tib <- NULL
uhm_base_ssh_tib <- base_ssh_tib %>% 
  dplyr::filter( Sample_Group == "MeTitration" ) %>%
  dplyr::filter( Source_Key != "E2" )

uhm_ssh_list <- NULL
uhm_ssh_list <- dplyr::bind_rows(
  uhm_base_ssh_tib %>% dplyr::filter( Sample_Titration ==  0 | Sample_Titration ==  50 ) %>% dplyr::mutate( Exp_Str = 1 ),
  uhm_base_ssh_tib %>% dplyr::filter( Sample_Titration ==  0 | Sample_Titration == 100 ) %>% dplyr::mutate( Exp_Str = 2 ),
  uhm_base_ssh_tib %>% dplyr::filter( Sample_Titration == 50 | Sample_Titration == 100 ) %>% dplyr::mutate( Exp_Str = 3 )
) %>% split( .$Exp_Str )

# uhm_ssh_list %>% lapply( function(x) { print( x, n=1000 ) })
  
#
# - Sample_Group = [ !MeTitration ] Replicate
base_ssh_tib %>% dplyr::filter( Sample_Group != "MeTitration" )


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Pre-processing:: Previous Products
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_man_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_manifests.sub8.rds")
rev_man_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/rev_manifests.sub8.rds")
neg_ctl_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_negative_ctls.rds" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Negative Controls
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

neg_ctl_tib  <- NULL
neg_ctl_tib  <- readr::read_rds( neg_ctl_rds ) %>%
  dplyr::mutate( Probe_ID = paste0(Probe_ID,"_0") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Manifest (ALL)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_man_tib  <- NULL
all_man_tib  <- readr::read_rds( all_man_rds )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Manifest (MSA10)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

msa_man_tib <- NULL
msa_man_tib <- format_msa_mans( 
  v03_man_csv = file.path( opt$top_path, "data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_SS_BP123_A1.csv.gz" ),
  v10_man_csv = file.path( opt$top_path, "data/pre-idats/MSA/MSA-48v1-0-Post-PQC_2B.body.csv.gz" ), 
  
  out_dir    = file.path( opt$out_path, "manifests/inputs" ),
  run_tag    = opt$run_name,
  reload     = opt$reload,
  reload_min = 10,
  ret_data   = FALSE,
  parallel   = opt$parallel,
  write_out  = FALSE,
  
  vb=vb,vt=vt+1,tc=tc+1, tt=tt )

#
# Interesting Summaries to be Reviewed...
#
if ( FALSE ) {

  msa_man_tib %>% 
    dplyr::mutate( Probe_Type = stringr::str_sub( Probe_ID, 1,2) ) %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  all_ctl_tib <- NULL
  all_ctl_tib <- all_man_tib %>% 
    dplyr::filter( Probe_Type != "cg" ) %>% 
    dplyr::filter( Probe_Type != "ch" ) %>% 
    dplyr::filter( Probe_Type != "rs" ) %>% 
    dplyr::filter( Probe_Type != "nv" )
  
  all_ctl_sum <- NULL
  all_ctl_sum <- all_ctl_tib %>%
    dplyr::group_by( Manifest,Manifest_Version,Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
}

# v1_body_csv <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv1/rank_historic_probes/EPICv2-UCSC-v5.rank_historic_probes.csv.gz" )
# v2_body_csv <- file.path( opt$top_path, "scratch/stable_scratch_EPICv2_docker_snps/EPICv2-UCSC-v5/EPICv2/rank_historic_probes/EPICv2-UCSC-v5.rank_historic_probes.csv.gz" )
# v1_body_tib <- NULL
# v1_body_tib <- readr::read_csv( file = v1_body_csv, show_col_types = FALSE )
# v2_body_tib <- NULL
# v2_body_tib <- readr::read_csv( file = v2_body_csv, show_col_types = FALSE )

v2_dock_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/docker/manifest/EPICv2/EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
v2_dock_tib <- NULL
v2_dock_tib <- readr::read_csv( file = v2_dock_csv, show_col_types = FALSE ) %>%
  clean_tib()
v2_dock_sum <- NULL
v2_dock_sum <- v2_dock_tib %>% print_sum( vec = c("Probe_Type"),
                                          vb=vb,vt=vt+1,tc=tc+1, tt=tt )

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

#
# TBD:: Scatter plots with random sub-selection::
#
opt$max_sam <- 10000
opt$max_rep <- 6

# opt$max_sam <- 3
# opt$max_rep <- 3
# opt$single <- TRUE

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Data-processing:: Replicate Sample List
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# LEFT OFF HERE!!
#
# [TBD]: Write VCF!
#
opt$single <- FALSE
opt$single <- TRUE
opt$min_dB <- 0.2

test_samples_vec <- NULL
test_samples_vec <- c( "CORIELL","HELA", "JURKAT","K562","MCF7",
                       "NA11922","NA12752","NA12877","NA12878","NA12882",
                       "R311861","R311874","R311882","R311886","R311895",
                       "RAJI" )
test_samples_vec <- c( "HELA","JURKAT","K562","MCF7",
                       "NA11922","NA12752","NA12877","NA12878","NA12882" )

org_sel_list <- NULL
org_sel_list <- ssh_all_tib %>% 
  # dplyr::filter( Source_Key == "R0" ) %>%
  dplyr::filter( Sample_Base != "EPIGENDX" ) %>%
  dplyr::filter( Sample_Base %in% test_samples_vec ) %>%
  # dplyr::filter( Sample_Input == 50 ) %>%
  dplyr::filter( Sample_Input >= 50 ) %>%
  dplyr::filter( User_Format != "EX24" ) %>%
  split( .$Sample_Base )

# Subset::
top_sel_list <- NULL
top_sel_list <- org_sel_list %>% lapply( head )
bot_sel_list <- NULL
bot_sel_list <- org_sel_list %>% lapply( tail )

rep_sel_list <- NULL
rep_sel_list <- dplyr::bind_rows(
  dplyr::bind_rows(top_sel_list, .id = "Sample_Base" ),
  dplyr::bind_rows(bot_sel_list, .id = "Sample_Base" )
) %>% 
  dplyr::distinct( Sentrix_Name, .keep_all = TRUE ) %>%
  dplyr::mutate( Sample_Name = paste(Sample_Base,Sample_Input, sep="_") ) %>%
  split( .$Sample_Base )

uhm_sel_list <- NULL
uhm_sel_list <- ssh_all_tib %>% 
  dplyr::filter( Sample_Base == "EPIGENDX" ) %>%
  # dplyr::filter( Sample_Input == 50 ) %>%
  dplyr::mutate( Sample_Name = paste(Sample_Base,Sample_Titration, sep="_") ) %>%
  split( .$Sample_Name )

#
# LEFT OFF HERE::
# [TBD]: Review other 'LEFT OFF HERE'
# [TBD]: Fix Negative Controls
# 
# Compare function:
#  [TBD]: Merge replicate/titration functions to single Rcpp function
#  [TBD]: Return a table instead of tibble
#  [TBD]: Run with both Poob/Negs/None(Pre-Masked)
#  [TBD]: Provide % Passing and Weighted Value
#

rep_dB_tib <- NULL
rep_data_tib <- NULL
negs_sdf_tab <- NULL
rm_outliers_vec <- c(TRUE, FALSE)
for ( rm_outliers in rm_outliers_vec ) {
  for ( sample_name in names(rep_sel_list) ) {
    all_beta_tib <- NULL
    all_poob_tib <- NULL
    all_negs_tib <- NULL
    
    rm_outliers_str <- "rm-outliers"
    if ( !rm_outliers ) rm_outliers_str <- "ok-outliers"
    
    cur_out_path <- NULL
    cur_out_path <- safe_mkdir( dir = file.path(opt$out_path, "samples",rm_outliers_str, paste0(sample_name,"_",rm_outliers_str) ) )
    
    # sample_name <- "HELA_200"
    # sample_name = "HELA"
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                    Data-processing:: Sentrix Sample List
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    
    ssh_sel_list <- NULL
    ssh_sel_list <- rep_sel_list[[sample_name]] %>% 
      dplyr::inner_join( idat_path_tib, by=c("Sentrix_Name") ) %>%
      dplyr::mutate( Rep_Idx = dplyr::row_number(),
                     Rep_Key = paste( Sample_Base,Rep_Idx, sep="_") ) %>% 
      split( .$Sentrix_Name ) %>% head( n=opt$max_rep )
    
    # sentrix_name <- names(ssh_sel_list)[1]
    # idat_path_tib
    
    for ( sentrix_name in names(ssh_sel_list) ) {
      prefix   <- ssh_sel_list[[sentrix_name]]$Sentrix_Path
      src_key  <- ssh_sel_list[[sentrix_name]]$Source_Key
      ng_input <- ssh_sel_list[[sentrix_name]]$Sample_Input
      if ( p1 ) cat(glue::glue("{pmssg} Current prefix({sample_name}) = '{prefix}'{RET}"))

      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                       Data-processing:: Loading
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      # prefix <- "/Users/bbarnes/Documents/data/idats/idats_MSAv03_48x1_Alpha/207545400008/207545400008_R04C01"
      # prefix <- "/Users/bbarnes/Documents/data/idats/EPIC/203010460029/203010460029_R12C01"
      # prefix <- "/Users/bbarnes/Documents/data/idats/EPIC/201502830033/201502830033_R02C01"
      # all_man_tib  <- NULL
      # all_man_tib  <- readr::read_rds( all_sub_rds ) # %>% dplyr::bind_rows( msa_man_tib )
      # all_man_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% dplyr::group_by( Manifest, Manifest_Version ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
      
      cur_out_name <- NULL
      cur_out_name <- paste0(sentrix_name,".",sample_name,"_",rm_outliers_str)
      cur_negs_csv <- NULL
      cur_negs_csv <- file.path( cur_out_path, paste0(cur_out_name,".negs.sdf.csv.gz") )
      
      idat_tib <- NULL
      idat_tib <- read_idat_pair_rcpp( 
        prefix_path      = prefix,
        output_path      = opt$out_path,
        workflow_vec     = workflow_vec %>% as.vector(),
        
        pval_add_vec     = all_neg_tib$Address %>% as.vector(), 
        addU_man_vec     = msa_man_tib$U %>% as.vector(),
        addM_man_vec     = msa_man_tib$M %>% as.vector(),
        
        cgns_man_vec     = msa_man_tib$Probe_ID %>% as.vector(), 
        cols_man_vec     = msa_man_tib$col %>% as.vector(),
        keys_man_vec     = msa_man_tib$Manifest %>% as.vector(),
        anns_man_vec     = msa_man_tib$Annotation %>% as.vector(), 
        chrs_man_vec     = msa_man_tib$Chromosome %>% as.vector(),
        
        min_pval         = opt$min_pval,
        min_beta         = opt$min_beta,
        max_beta         = opt$max_beta,
        min_perO         = opt$min_perO, 
        min_perI         = opt$min_perI, 
        read_bgz         = FALSE,
        write_bgz        = FALSE, 
        rm_pval_outliers = rm_outliers,
        return_df        = opt$return_df,
        vb=vb,vt=vt,tc=tc ) %>% 
        tibble::as_tibble() %>% clean_tib() %>%
        utils::type.convert(as.is=TRUE) %>% 
        dplyr::mutate(across(where(is.factor), as.character) ) %>%
        dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )
      
      # How to measure outlier removal::
      # 
      negs_sdf <- NULL
      negs_sdf <- idat_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% 
        dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative_") )
      readr::write_csv( x = negs_sdf, file = cur_negs_csv )
      
      negs_sdf_tab <- negs_sdf_tab %>% 
        dplyr::bind_rows(
          dplyr::mutate( negs_sdf,
                         Rm_Outliers = rm_outliers,
                         Sample = sample_name ) %>%
            dplyr::select( Rm_Outliers,Sample, dplyr::everything() )
        )
      
      # negs_sdf %>%
      #   dplyr::group_by( mask ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
      
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
      idat_sum <- print_sum( tib = idat_tib, vec = c("mask","Probe_Type","col"),
                             vb=vb,vt=vt,tc=tc, tt=tt )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Data-processing:: Parse Control Data
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      # data_tib <- NULL
      # data_tib <- idat_tib
      # data_tib <- idat_tib %>% 
      #   dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") )
      
      ctls_tib <- NULL
      ctls_tib <- idat_tib %>% 
        dplyr::filter(  stringr::str_starts( Probe_ID, pattern = "ct") )
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Data-processing:: Workflow Pipeline
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      work_sdf <- NULL
      # work_sdf <- data_tib %>%
      work_sdf <- idat_tib %>%
        # dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") ) %>%
        as.data.frame() %>%
        sesame::SigDF( platform = "EPIC", 
                       ctl = ctls_tib
                       # ctl = NULL
        ) %>%
        mutate_sdf_simple( 
          steps = "DB", 
          negs_min = 1.0, 
          poob_min = 1.0, 
          vb=vb+3,vt=vt,tc=tc )
      
      rep_key <- NULL
      rep_key <- ssh_sel_list[[sentrix_name]]$Rep_Key
      
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      #                   Data-processing:: Extract Beta Values
      # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
      
      beta_tib <- NULL
      beta_tib <- mutate_sdf_simple( 
        sdf = work_sdf,
        # steps = "v", # No Masking...
        steps = "V", # With masking...
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
      
      if ( FALSE ) {
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        #           Data-processing:: Extract detectionPnegEcdf2 P-Values
        # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
        
        negs_tib <- NULL
        negs_tib <- mutate_sdf_simple( 
          sdf = work_sdf,
          steps = "N", 
          negs_min = 1.0, 
          poob_min = 1.0, 
          vb=vb, vt=vt, tc=tc ) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column( var = "Probe_ID") %>% 
          tibble::as_tibble() %>% 
          magrittr::set_names( value = c("Probe_ID",rep_key) )
        # magrittr::set_names( value = c("Probe_ID","Poob") )
        
        if ( is.null( all_negs_tib ) ) {
          all_negs_tib <- negs_tib
        } else {
          all_negs_tib <- all_negs_tib %>% 
            inner_join( negs_tib, by = c("Probe_ID") )
        }
      }
      
    }
    if ( p0 ) cat(glue::glue("{pmssg} Finished Processing Idat Pairs sample='{sample_name}'.{RET2}"))
    
    #
    # LEFT OFF HERE::
    #  - [TBD]: Incorporate: {ng_input,rep_key}, could just pull from the sample sheet...
    #
    #  - [TBD]: Need to implement dual input matrices
    #     - This will allow Titration efforts
    #
    # [TBD]: Incorporate weight score
    # [TBD]: Add sample names
    # [TBD]: Reduce output
    # [TBD]: Add EPICv2
    # [TBD]: Create Comparison Groups Sheet
    #
    # Add More Stats...
    cur_dB_tib <- NULL
    cur_dB_tib <- calc_dBs( beta_tibA = all_beta_tib, 
                            pval_tibA = all_poob_tib,
                            name_strA = sample_name,
                            
                            min_dB = opt$min_dB, cmp_str = "lt", is_abs = TRUE,
                            
                            out_dir    = file.path( opt$out_path, "replicate" ),
                            run_tag    = sample_name, 
                            reload     = opt$reload,
                            reload_min = 10, 
                            ret_data   = FALSE,
                            parallel   = opt$parallel,
                            
                            vb=vb+2,vt=vt+1,tc=tc+1,tt=tt )
    
    # [TBD]: Need to plot EPIC correlation against db_pas_per & wS_per, which is better
    # cur_dat$ret_tib %>% dplyr::filter( is.na(wS_per) )


    #
    # General Test Case Plotting::
    #
    if ( FALSE ) {
      box_pas_wS_ggg <- NULL
      box_pas_wS_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes(x=db_pas_per,y=wS_per, group = db_pas_per) ) + 
        ggplot2::geom_boxplot( varwidth = TRUE )
      box_pas_wS_ggg
      
      pnt_pas_wS_ggg <- NULL
      pnt_pas_wS_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes(x=db_pas_per,y=wS_per) ) + 
        ggplot2::geom_point()
      pnt_pas_wS_ggg
      
      den_wS_ggg <- NULL
      den_wS_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per ) ) +
        ggplot2::geom_density()
      den_wS_ggg
      
      pnt_wS_avg_ggg <- NULL
      pnt_wS_avg_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_avg ) ) +
        ggplot2::geom_point()
      pnt_wS_avg_ggg
      
      d2d_wS_avg_ggg <- NULL
      d2d_wS_avg_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_avg ) ) +
        ggplot2::geom_density2d()
      d2d_wS_avg_ggg
      
      pnt_wS_med_ggg <- NULL
      pnt_wS_med_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_med ) ) +
        ggplot2::geom_point()
      pnt_wS_avg_ggg

      d2d_wS_med_ggg <- NULL
      d2d_wS_med_ggg <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_med ) ) +
        ggplot2::geom_density2d()
      d2d_wS_med_ggg
      
      d2d_wS_med_ggg2 <- NULL
      d2d_wS_med_ggg2 <- cur_dB_tib %>% 
        ggplot2::ggplot( aes( x=wS_per, y=dB_med ) ) +
        ggplot2::geom_density2d() +
        ggplot2::facet_grid(
          rows = vars(db_pas_per)
        )
      d2d_wS_med_ggg2
      

    }
    
    rep_data_tib <- rep_data_tib %>% 
      dplyr::bind_rows( cur_dB_tib )
    
    if ( opt$single ) break
  }
  if ( opt$single ) break
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                    Summarize and Plot Results::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

rep_data_sum <- NULL
rep_data_sum <- rep_data_tib %>% 
  print_sum( vec = c("Sample","Rm_Outliers"), vb=vb+3,vt=vt,tc=tc )

if ( par$run_name == "EPICv1" ||
     par$run_name == "COREv1" ||
     par$run_name == "FAILv1" ) {
  
  # rep_data_tib %>% dplyr::filter( Probe_ID %in% epic_ses_dat$mask )
  
  epic_mask_tib <- NULL
  epic_mask_tib <- tibble::tibble( Probe_ID = epic_ses_dat$mask, Masked = TRUE )
  
  rep_mask_tib <- NULL
  rep_mask_tib <- rep_data_tib %>% 
    dplyr::left_join( epic_mask_tib, by=c("Probe_ID") ) %>% 
    dplyr::mutate(
      Masked = dplyr::case_when(
        is.na(Masked) ~ FALSE,
        TRUE ~ TRUE )
    )
  
  rep_mask_sum <- NULL
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
  rep_den_mask_gg <- NULL
  rep_den_mask_gg <- rep_mask_tib %>% 
    dplyr::filter( fin_scr >= 0.9 ) %>%
    # ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
    ggplot2::ggplot( aes( x=log(fin_scr) + sn, color = Sample, fill = Sample ) ) +
    ggplot2::geom_density( alpha = 0.2 )  +
    ggplot2::facet_grid( rows = vars(Masked) )
  
  # Conclusion:: Score does show a difference for Replicate score vs. third
  #  party analysis...
  #
  rep_den_mask_pdf <- file.path( opt$out_path, paste(opt$run_name,"rep_den_mask.pdf", sep=".") )
  rep_den_mask_gg <- NULL
  rep_den_mask_gg <- rep_mask_tib %>% 
    dplyr::filter( fin_scr >= 0.9 ) %>%
    ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
    # ggplot2::ggplot( aes( x=log(fin_scr) + sn, color = Sample, fill = Sample ) ) +
    ggplot2::geom_density( alpha = 0.2 )  +
    ggplot2::facet_grid( rows = vars(Masked),
                         cols = vars(Sample) )
  ggplot2::ggsave( filename = rep_den_mask_pdf, 
                   device = "pdf", width = 7, height = 7, dpi = 320 )
  
} else {
  
  # Looks good::
  rep_den_pdf <- file.path( opt$out_path, paste(opt$run_name,"rep_density.pdf", sep=".") )
  rep_den_gg <- NULL
  rep_den_gg <- rep_data_tib %>% 
    ggplot2::ggplot( aes( x=fin_scr, color = Sample, fill = Sample ) ) +
    ggplot2::geom_density( alpha = 0.2 ) +
    ggplot2::facet_grid( cols = vars(Sample),
                         rows = vars(Rm_Outliers) )
  
  ggplot2::ggsave( filename = rep_den_pdf, 
                   device = "pdf", 
                   width = 7, 
                   height = 7, 
                   dpi = 320 )
  
  #
  # Plot Negative w/w0 Outlier Removal
  #
  negs_sdf_ggg <- negs_sdf_tab %>% 
    ggplot2::ggplot( aes(x=UG, fill=Rm_Outliers) ) + 
    ggplot2::geom_density( alpha=0.2 ) +
    ggplot2::facet_grid( cols = vars(Sample),
                         rows = vars(Rm_Outliers) )
  
  # More or less worthless::
  rep_2den_gg <- NULL
  rep_2den_gg <- rep_data_tib %>% 
    ggplot2::ggplot( aes( x=fin_scr, y=med_dbs ) ) +
    ggplot2::geom_density2d() +
    ggplot2::facet_grid( rows = vars(Sample) )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
