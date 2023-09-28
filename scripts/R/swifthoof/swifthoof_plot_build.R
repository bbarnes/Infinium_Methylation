
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

# Load Plotting Packages::
suppressWarnings(suppressPackageStartupMessages(
  base::require("ggplot2", quietly = TRUE) ) )
suppressWarnings(suppressPackageStartupMessages(
  base::require("GGally", quietly = TRUE) ) )

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
par$prgm_tag <- 'swifthoof_plot_build'
par$verbose  <- 3
local_paths  <- c( 
  "/Users/bbarnes/Documents/tools/imSuite/scripts/R"
)

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

par$run_name <- "EPICv2"

opt <- NULL
# opt <- swifthoof_sesame_options( pars = par, args = args, vb = par$verbose )
opt <- imProbeQC_options( pars = par, args = args, vb = par$verbose )
# opt <- stable_options( pars = par, args = args, vb = par$verbose )
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
#                            Local Run Parameters::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# opt$build_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/docker-1.11.15.1.p.0.4/EPICv2.1.11.15.1.p.0.4/swifthoof_main" )
opt$build_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/docker-1.11.15.1.p.0.4" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Load Sample Sheets::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# Data Files: Beta/Pvals::
dat_source_lst <- NULL
dat_source_lst <- file_list( path    = opt$build_path, 
                             prefix  = opt$build_path, 
                             suffix  = "_EPIC_.*_raw.call.dat.csv.gz",
                             pattern = "_raw.call.dat.csv.gz$", 
                             recursive = TRUE )


# Auto Sample Sheet::
if ( FALSE ) {
  ssh_auto_tib <- NULL
  ssh_auto_tib <- file_list( path = opt$build_path, 
                             prefix = opt$build_path, 
                             suffix = "_EPIC_A1_AutoSampleSheet.csv.gz",
                             pattern = "_EPIC_A1_AutoSampleSheet.csv.gz$" ) %>% 
    lapply( read_csv, show_col_types = FALSE ) %>%
    dplyr::bind_rows()
  
  # AutoSample_R2_Key_1
  # AutoSample_dB_Key_1
  ssh_auto_tib %>% 
    dplyr::filter( AutoSample_R2_Key_1 != AutoSample_dB_Key_1 )
}

# Raw Human Provided Sample Sheet::
# ssh_source_tib <- NULL
# ssh_source_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/SampleSheets/SampleSheet-LightningAuto.csv.gz" )
# ssh_source_tib <- readr::read_csv( file = ssh_source_csv, skip = 7,
#                                    show_col_types = FALSE )

# Formatted Human Provided Sample Sheet::
ssh_source_tib <- NULL
ssh_source_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/GSIBIOINFO-638/SampleSheets/formatted/EPICv2-UCSC-v0.LightningAuto.select.sample_sheet.csv.gz" )
ssh_source_tib <- readr::read_csv( file = ssh_source_csv,
                                   show_col_types = FALSE ) %>%
  dplyr::filter( Sentrix_Name %in% names(dat_source_lst) )

ssh_source_sum <- NULL
ssh_source_sum <- ssh_source_tib %>% 
  dplyr::group_by( Manifest_Key,Sample_Base ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
ssh_source_sum %>% print( n=base::nrow(ssh_source_sum) )

ssh_source_lst <- NULL
ssh_source_lst <- ssh_source_tib %>% split( .$Sample_Base )
# ssh_source_lst <- ssh_source_tib %>% split( .$Sample_Name )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Process Each Sample::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

for ( sample_name in names(ssh_source_lst) ) {
  # sample_name <- "HeLa"
  
  cur_ssh_tib  <- NULL
  cur_ssh_tib  <- ssh_source_lst[[sample_name]] %>%
    dplyr::filter( Concentration == 500 ) %>%
    dplyr::mutate( Plot_Name = paste0( "v",Chip_Version,"_",Sample_Name ) ) %>%
    dplyr::group_by( Plot_Name ) %>%
    dplyr::mutate( Rep_Idx = dplyr::row_number(),
                   Plot_Name = paste( Plot_Name,Rep_Idx, sep="_") )
  
  sample_group <- NULL
  sample_group <- cur_ssh_tib %>% dplyr::distinct(Sample_Group) %>%
    dplyr::pull( Sample_Group )
  
  if ( sample_group[1] != "CellLine" ) next

  # Load All Data
  cur_dat_lst <- NULL
  cur_dat_lst <- dat_source_lst[cur_ssh_tib$Sentrix_Name] %>% 
    lapply( read_csv, show_col_types = FALSE )

  # Update Names with Plot_Names and split data by beta/pval
  cur_dat_tib <- NULL
  cur_dat_tib <- cur_dat_lst %>%
    dplyr::bind_rows( .id = "Sentrix_Name" ) 
  
  cur_dat_tab <- NULL
  cur_dat_tab <- cur_ssh_tib %>% 
    dplyr::select( Sentrix_Name, Plot_Name ) %>% 
    dplyr::inner_join( cur_dat_tib, multiple = "all", by=c("Sentrix_Name") )
  
  #
  # QUICK CHECK::
  #
  raw_beta_tab <- NULL
  raw_beta_tab <- cur_dat_tab %>% 
    dplyr::select( -Sentrix_Name ) %>% 
    dplyr::select( Probe_ID,Plot_Name,betas ) %>%
    tidyr::pivot_wider( id_cols = c(Probe_ID), 
                        names_from = c(Plot_Name), 
                        values_from = c("betas") ) %>% 
    dplyr::filter( !is.na(v1_HeLa_500_1) ) %>% dplyr::filter( !is.na(v2_HeLa_500_1) )

  raw_pval_tab <- NULL
  raw_pval_tab <- cur_dat_tab %>% 
    dplyr::select( -Sentrix_Name ) %>% 
    dplyr::select( Probe_ID,Plot_Name,pvals_pOOBAH ) %>%
    tidyr::pivot_wider( id_cols = c(Probe_ID), 
                        names_from = c(Plot_Name), 
                        values_from = c("pvals_pOOBAH") ) %>% 
    dplyr::filter( !is.na(v1_HeLa_500_1) ) %>% dplyr::filter( !is.na(v2_HeLa_500_1) )
  
  #
  # Quick Pearson Check::
  #
  raw_r2_mat <- raw_beta_tab %>% 
    tibble::column_to_rownames( var = "Probe_ID" ) %>% 
    as.matrix() %>% 
    stats::cor( use = "pairwise.complete.obs", method = "pearson" )
  
  #
  # Full Check::
  #
  # cur_dat_tab %>% 
  #   dplyr::select( -Sentrix_Name ) %>% 
  #   tidyr::pivot_wider( id_cols = c(Probe_ID), 
  #                       names_from = c(Plot_Name), 
  #                       values_from = c("pvals_pOOBAH", "pvals_PnegEcdf", "betas") )
  # 
  # cur_dat_tab %>% dplyr::select( Probe_ID, Plot_Name, betas ) %>% dplyr::filter( Plot_Name == "" )
  
  min_pval <- 0.1
  grp_key <- NULL
  cond_grp_vec <- NULL
  rank_cols <- NULL
  plot_dir <- safe_mkdir( dir = file.path( opt$out_path, "plots") )
  
  # Plot Over All
  pairs_gg <- NULL
  pairs_gg <- GGally::ggpairs( 
    data = raw_beta_tab %>% head(n=1000),
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
                                 
                                 beta_data = raw_beta_tab,
                                 pval_data = raw_pval_tab,
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
                                 
                                 beta_data = raw_beta_tab,
                                 pval_data = raw_pval_tab,
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
                                 
                                 beta_data = raw_beta_tab,
                                 pval_data = raw_pval_tab,
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
  
  
  
  
  # Plot Over Manifest_Key
  
  if ( opt$single ) break
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
