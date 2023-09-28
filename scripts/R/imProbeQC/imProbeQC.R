
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
par$prgm_tag <- 'imProbeQC'
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
# par$version <- 3
# par$version <- 4
# par$version <- 5
# par$version <- 6

par$run_name <- "EPICv2"
# par$run_name <- "MSA"

run_previous <- TRUE
run_previous <- FALSE

if ( run_previous ) {
  par$run_name <- paste0( par$run_name, ".prev" )
} else {
  par$run_name <- paste0( par$run_name, ".bacr" )
}

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

success = TRUE

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                             Load Manifest::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_tib <- NULL
# man_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/EPIC-8v2-0_A1.body.csv.gz" )
man_csv <- file.path( opt$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/full/EPIC_v2.gs_to_sesame.csv.gz" )
man_tib <- readr::read_csv( file = man_csv )

man_fail_tib <- NULL
man_fail_csv <- file.path( opt$top_path, "data/manifests/methylation/MethylationEPIC_v2-2/EPIC-8v2-0_A1-FlaggedProbes/EPIC-8v2-0_A1-FlaggedProbes.csv.gz" )
man_fail_tib <- readr::read_csv( file = man_fail_csv, show_col_types = FALSE ) %>% 
  dplyr::select( IlmnID ) %>% 
  dplyr::rename( Probe_ID = IlmnID) %>% 
  dplyr::mutate( Masked = TRUE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Load Sample Sheet::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssh_all_tib <- NULL
ssh_all_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/imProbeQC_sampleSheet/EPICv2-UCSC-v0/EPICv2-UCSC-v0.LightningAuto.select.sample_sheet.csv.gz" )
ssh_all_tib <- readr::read_csv( file = ssh_all_csv, show_col_types = FALSE )

ssh_rep_tib <- NULL
ssh_rep_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/imProbeQC_sampleSheet/EPICv2-UCSC-v0/EPICv2-UCSC-v0.LightningAuto.replicate.sample_sheet.csv.gz" )
ssh_rep_tib <- readr::read_csv( file = ssh_rep_csv, 
                                skip = 10,
                                show_col_types = FALSE )

ssh_tit_tib <- NULL
ssh_tit_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/imProbeQC_sampleSheet/EPICv2-UCSC-v0/EPICv2-UCSC-v0.LightningAuto.titration.sample_sheet.csv.gz" )
ssh_tit_tib <- readr::read_csv( file = ssh_tit_csv,
                                skip = 10,
                                show_col_types = FALSE )

# Should be zero::
#  base::nrow(ssh_rep_tib) + base::nrow(ssh_tit_tib) - base::nrow(ssh_all_tib)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Load Beta Values::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( run_previous ) {
  # Previous Version:: 20022023
  rep_beta_tsv <- file.path( opt$top_path, "Projects.new/imProbeQC/output/methylation_pipeline_EPICv2_LightingAuto_Replicate/post_processing/sesame_betas_filtered_pOOBAH.txt.gz" )
  tit_beta_tsv <- file.path( opt$top_path, "Projects.new/imProbeQC/output/methylation_pipeline_EPICv2_LightingAuto_Titration/post_processing/sesame_betas_filtered_pOOBAH.txt.gz" )
} else {
  # Latest Version:: 18042023
  rep_beta_tsv <- file.path( opt$top_path, "Projects.new/imProbeQC/output/bacr/methylation_pipeline_EPICv2_LightingAuto_Replicate/post_processing/betas_filtered_pOOBAH.txt.gz" )
  tit_beta_tsv <- file.path( opt$top_path, "Projects.new/imProbeQC/output/bacr/methylation_pipeline_EPICv2_LightingAuto_Titration/post_processing/betas_filtered_pOOBAH.txt.gz" )
}

rep_beta_tib <- NULL
tit_beta_tib <- NULL

rep_beta_tib <- readr::read_tsv( file = rep_beta_tsv, show_col_types = FALSE )
tit_beta_tib <- readr::read_tsv( file = tit_beta_tsv, show_col_types = FALSE )

names(rep_beta_tib)[1] <- "Probe_ID"
names(tit_beta_tib)[1] <- "Probe_ID"

rep_name_vec <- NULL
tit_name_vec <- NULL

rep_name_vec <- rep_beta_tib %>% names()
tit_name_vec <- tit_beta_tib %>% names()

rep_name_vec <- rep_name_vec[2:length(rep_name_vec)]
tit_name_vec <- tit_name_vec[2:length(tit_name_vec)]

#
# Split Replicate and Titration Sample Sheets::
#
rep_ssh_tib <- NULL
rep_ssh_tib <- ssh_all_tib %>% dplyr::filter( Sentrix_Name %in% rep_name_vec ) %>% 
  dplyr::filter( Manifest_Key == "EPIC_v2") %>% 
  dplyr::filter( Sample_Base != "NA12873" ) %>% 
  dplyr::filter( Sample_Base != "PromegaMouse" ) %>% 
  dplyr::filter( Concentration >= 250 )

tit_ssh_tib <- NULL
tit_ssh_tib <- ssh_all_tib %>% 
  dplyr::filter( Sentrix_Name %in% tit_name_vec ) %>%
  dplyr::filter( Manifest_Key == "EPIC_v2") %>%
  dplyr::filter( Sample_Base == "Epigen" )

#
# Convert Sample Sheets to Lists
#
rep_ssh_lst <- NULL
rep_ssh_lst <- rep_ssh_tib %>% split(.$Sample_Base)

tit_ssh_lst <- NULL
tit_ssh_lst <- tit_ssh_tib %>% split(.$Sample_Name)

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                      Analyze Samples:: Replicates
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$min_dB <- 0.2

# dB_tib <- NULL
# dB_tib <- tibble::tibble( Probe_ID = rep_beta_tib$Probe_ID )

rep_tib <- NULL
rep_tib <- tibble::tibble( Probe_ID = rep_beta_tib$Probe_ID )

rep_cnt <- rep_ssh_lst %>% length()
for ( sample in names(rep_ssh_lst) ) {
  cat(glue::glue("{pmssg} Replicate Analysis: {sample}.{RET}"))

  samp_col <- NULL
  beta_tib <- NULL
  # beta_mat <- NULL
  
  if ( rep_ssh_lst[[sample]] %>% base::nrow() < 3 ) next
  
  samp_col <- c("Probe_ID", rep_ssh_lst[[sample]]$Sentrix_Name )
  beta_tib <- rep_beta_tib %>% dplyr::select( dplyr::all_of( samp_col ) )
  # beta_mat <- beta_tib %>% tibble::column_to_rownames( var = "Probe_ID" ) %>% 
  #   as.matrix()
  
  #
  # Calculate dB Failures::
  #
  cur_tib <- NULL
  cur_tib <- probe_replicate_screen( tib = beta_tib, 
                                     sample = sample, 
                                     min_dB = opt$min_dB, 
                                     
                                     out_dir = opt$out_path,
                                     run_tag = "Lower", 
                                     reload = opt$reload,
                                     reload_min = 10, 
                                     ret_data = FALSE,
                                     parallel = opt$parallel,
                                     
                                     vb=vb,vt=vt+1,tc=tc+1,tt=tt )
  
  rep_tib <- rep_tib %>% dplyr::bind_cols( cur_tib )
  
  if ( opt$single ) break
}

# Example Summary::
#  rep_tib %>% dplyr::group_by( dB_Pass_HeLa ) %>% dplyr::summarise( Count=n(), .groups = "drop" )

prb_rep_tab <- NULL
prb_rep_tab <- rep_tib %>% 
  tidyr::pivot_longer( cols = - Probe_ID, 
                       names_to = "Sample", 
                       values_to = "Pass" )

prb_rep_sum <- NULL
prb_rep_sum <- prb_rep_tab %>% 
  dplyr::group_by( Probe_ID, Pass ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% 
  tidyr::pivot_wider( id_cols = Probe_ID, 
                      names_from = Pass, 
                      values_from = Count, 
                      values_fill = 0 ) %>% 
  dplyr::rename( Pass_Cnt = 'TRUE', Unk_Cnt = 'NA', Fail_Cnt = 'FALSE' ) %>% 
  dplyr::mutate( Pass_Per = Pass_Cnt / ( Pass_Cnt + Unk_Cnt + Fail_Cnt ) )

# Summary of Summary::
prb_rep_sum_sum <- NULL
prb_rep_sum_sum <- prb_rep_sum %>% 
  dplyr::group_by( Pass_Per ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

#
# TBD:: Build Sample Based Pass Percent
#
sam_rep_sum <- NULL
sam_rep_sum <- prb_rep_tab %>% 
  dplyr::group_by( Sample, Pass ) %>% 
  dplyr::summarise( Count = n(), .groups = "drop" ) %>% 
  tidyr::pivot_wider( id_cols = Sample, 
                      names_from = Pass, 
                      values_from = Count ) %>% 
  dplyr::rename( Pass_Cnt = 'TRUE', Unk_Cnt = 'NA', Fail_Cnt = 'FALSE' ) %>% 
  dplyr::mutate( Pass_Per = Pass_Cnt / ( Pass_Cnt + Unk_Cnt + Fail_Cnt ) )

if ( p1 ) sam_rep_sum %>% print( n=base::nrow(sam_rep_sum) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#.                      Analyze Samples:: Titration
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$run_titration <- FALSE
opt$run_titration <- TRUE

if ( opt$run_titration ) {
 
  tit_tib <- NULL
  tit_tib <- tibble::tibble( Probe_ID = tit_beta_tib$Probe_ID )
  
  name_vec <- names(tit_ssh_lst)
  name_cnt <- name_vec %>% length()

  for ( ii in c(1:name_cnt) ) {
    for ( jj in c(1:name_cnt) ) {
      if ( ii >= jj ) next
      
      samp_colA <- NULL
      beta_tibA <- NULL
      samp_colB <- NULL
      beta_tibB <- NULL
      
      samp_colA <- c("Probe_ID", tit_ssh_lst[[name_vec[ii]]]$Sentrix_Name )
      beta_tibA <- tit_beta_tib %>% dplyr::select( dplyr::all_of( samp_colA ) )
      
      samp_colB <- c("Probe_ID", tit_ssh_lst[[name_vec[jj]]]$Sentrix_Name )
      beta_tibB <- tit_beta_tib %>% dplyr::select( dplyr::all_of( samp_colB ) )
      
      #
      # Update this later
      #
      
      min_dB <- 0.2
      less_than <- FALSE
      if ( name_vec[ii] %>% stringr::str_ends("_0") &&
           name_vec[jj] %>% stringr::str_ends("_100") ) min_dB <- 0.5
      
      if ( name_vec[ii] %>% stringr::str_ends("_0") &&
           name_vec[jj] %>% stringr::str_ends("_50") ) min_dB <- 0.2
      
      if ( name_vec[ii] %>% stringr::str_ends("_50") &&
           name_vec[jj] %>% stringr::str_ends("_100") ) min_dB <- 0.1

      if ( name_vec[ii] %>% stringr::str_ends("_100") &&
           name_vec[jj] %>% stringr::str_ends("_0") ) {
        min_dB <- 0.5
        less_than <- TRUE
      }
      if ( name_vec[ii] %>% stringr::str_ends("_50") &&
           name_vec[jj] %>% stringr::str_ends("_0") ) {
        min_dB <- 0.2
        less_than <- TRUE
      }
      if ( name_vec[ii] %>% stringr::str_ends("_100") &&
           name_vec[jj] %>% stringr::str_ends("_50") ) {
        min_dB <- 0.1
        less_than <- TRUE
      }
      
      cur_tib <- NULL
      cur_tib <- probe_replicate_screen2( tibA = beta_tibA,
                                          tibB = beta_tibB,
                                          sampleA = name_vec[ii], 
                                          sampleB = name_vec[jj], 
                                          min_dB = min_dB, 
                                          less_than = less_than,
                                          
                                          out_dir = opt$out_path,
                                          run_tag = "Lower", 
                                          reload = opt$reload,
                                          reload_min = 10, 
                                          ret_data = FALSE,
                                          parallel = opt$parallel,
                                          
                                          vb=vb,vt=vt+1,tc=tc+1,tt=tt )
      
      tit_tib <- tit_tib %>% dplyr::bind_cols( cur_tib )
      
      if ( opt$single ) break
      break
    }
    break
  }
  
  prb_tit_tab <- NULL
  prb_tit_tab <- tit_tib %>%
    tidyr::pivot_longer( cols = - Probe_ID, 
                         names_to = "Sample", 
                         values_to = "Pass" )
  
  prb_tit_sum <- NULL
  prb_tit_sum <- prb_tit_tab %>% 
    dplyr::group_by( Probe_ID, Pass ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ) %>% 
    tidyr::pivot_wider( id_cols = Probe_ID, 
                        names_from = Pass, 
                        values_from = Count, 
                        values_fill = 0 ) %>% 
    dplyr::rename( Pass_Cnt = 'TRUE', Unk_Cnt = 'NA', Fail_Cnt = 'FALSE' ) %>% 
    dplyr::mutate( Pass_Per = Pass_Cnt / ( Pass_Cnt + Unk_Cnt + Fail_Cnt ) )
  
  # Summary of Summary::
  prb_tit_sum_sum <- NULL
  prb_tit_sum_sum <- prb_tit_sum %>% 
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
}

prb_all_tab <- NULL
prb_all_tab <- rep_tib %>% 
  dplyr::left_join( tit_tib, by=c("Probe_ID") ) %>% 
  tidyr::pivot_longer( cols = - Probe_ID, 
                       names_to = "Sample", 
                       values_to = "Pass" )

prb_all_sum <- NULL
prb_all_sum <- prb_all_tab %>% 
  dplyr::group_by( Probe_ID, Pass ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" ) %>% 
  tidyr::pivot_wider( id_cols = Probe_ID, 
                      names_from = Pass, 
                      values_from = Count, 
                      values_fill = 0 ) %>% 
  dplyr::rename( Pass_Cnt = 'TRUE', Unk_Cnt = 'NA', Fail_Cnt = 'FALSE' ) %>% 
  dplyr::mutate( Pass_Per = Pass_Cnt / ( Pass_Cnt + Unk_Cnt + Fail_Cnt ) )

# Summary of Summary::
prb_all_sum_sum <- NULL
prb_all_sum_sum <- prb_all_sum %>% 
  dplyr::group_by( Pass_Per ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                                 Summary::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

mask_sum_tib <- NULL
mask_sum_tib <- prb_rep_sum %>% 
  dplyr::left_join( prb_tit_sum, by=c("Probe_ID"), suffix=c("_Rep","_UHM") ) %>%
  dplyr::left_join( man_fail_tib, by=c("Probe_ID") ) %>%
  dplyr::mutate(
    Masked = dplyr::case_when(
      is.na( Masked ) ~ FALSE,
      TRUE ~ TRUE
    )
  ) %>%
  dplyr::group_by( Pass_Per_Rep, Pass_Per_UHM, Masked ) %>%
  dplyr::summarise( Count=n(), .groups = "drop" )

mask_sum_tib %>% print( n=1000 )

if ( FALSE ) {
  mask_sum_tib %>% 
    # dplyr::filter( Masked == TRUE ) %>%
    ggplot2::ggplot( aes(x=Count, color=Masked ) ) + 
    ggplot2::geom_density() + 
    ggplot2::facet_grid( rows = vars(Pass_Per_Rep), 
                         cols = vars(Pass_Per_UHM), 
                         scales = "free" )
  
}

# All Summary:
full_sum_tib <- NULL
full_sum_tib <- dplyr::left_join( prb_rep_sum, 
                                  prb_tit_sum, by=c("Probe_ID"), 
                                  suffix=c("_Rep","_UHM") ) %>%
  dplyr::group_by( Pass_Per_Rep, Pass_Per_UHM ) %>%
  dplyr::summarise( Count=n(), .groups = "drop" )

#
# CGN Only Summary::
#
cgn_all_sum <- NULL
cgn_all_sum <- dplyr::full_join(
  prb_all_sum %>% dplyr::filter( Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("cg") ) %>%
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  prb_all_sum %>% dplyr::filter( !Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("cg") ) %>%
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  by=c("Pass_Per"), suffix=c("_Fail", "_Pass")
)

cgn_rep_sum <- NULL
cgn_rep_sum <- dplyr::full_join(
  prb_rep_sum %>% dplyr::filter( Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("cg") ) %>%
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  prb_rep_sum %>% dplyr::filter( !Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("cg") ) %>%
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  by=c("Pass_Per"), suffix=c("_Fail", "_Pass")
)  

cgn_tit_sum <- NULL
cgn_tit_sum <- dplyr::full_join(
  prb_tit_sum %>% dplyr::filter( Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("cg") ) %>%
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  prb_tit_sum %>% dplyr::filter( !Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("cg") ) %>%
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  by=c("Pass_Per"), suffix=c("_Fail", "_Pass")
)  

#
# All Probe Summary::
#
all_sum <- NULL
all_sum <- dplyr::full_join(
  prb_all_sum %>% dplyr::filter( Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  prb_all_sum %>% dplyr::filter( !Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  by=c("Pass_Per"), suffix=c("_Fail", "_Pass")
)

rep_sum <- NULL
rep_sum <- dplyr::full_join(
  prb_rep_sum %>% dplyr::filter( Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  prb_rep_sum %>% dplyr::filter( !Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  by=c("Pass_Per"), suffix=c("_Fail", "_Pass")
)  

tit_sum <- NULL
tit_sum <- dplyr::full_join(
  prb_tit_sum %>% dplyr::filter( Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  prb_tit_sum %>% dplyr::filter( !Probe_ID %in% man_fail_tib$Probe_ID ) %>% 
    dplyr::group_by( Pass_Per ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ),
  by=c("Pass_Per"), suffix=c("_Fail", "_Pass")
)  

fail_mask_sum <- NULL
fail_mask_sum <- prb_all_tab %>% 
  dplyr::filter( Pass == FALSE ) %>% 
  dplyr::distinct( Probe_ID ) %>%
  dplyr::left_join( man_tib, by=c("Probe_ID") ) %>%
  dplyr::group_by(mask) %>%
  dplyr::summarise( Count=n(), .groups = "drop" )

man_mask_sum <- NULL
man_mask_sum <- man_tib %>%
  dplyr::group_by(mask) %>%
  dplyr::summarise( Count=n(), .groups = "drop" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
