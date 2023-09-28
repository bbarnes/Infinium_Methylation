
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
par$prgm_tag <- 'stable_docker_join_manifest'
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
par$version <- 1
par$version <- 2
# par$version <- 3
# par$version <- 4
# par$version <- 5

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

#
# Manifest v1 -> v2 Workflow::
#
#  1. Load EPIC v1 Genome Studio Manifest
#  2. Load EPIC v2 Genome Studio Manifest
#  3. Inspect manifests
#
#  4. Intersect v1/v2 by {Name,U,M}
#  5. Inspect results
#
#  6. Convert to Workhorse Docker Manifest format
#  7  Download v.1.11.15 core/EPIC-B4 manifest
#  8  Compare against original v.1.11.15 format
#
#  9. Clean Up Files
#

# gzip -dc /Users/bbarnes/Documents/data/manifests/methylation/GenomeStudio/infinium-methylationepic-v-1-0-b5-manifest-file.csv.gz | head -n 865927 | tail -n 865920 > /Users/bbarnes/Documents/data/manifests/methylation/GenomeStudio/EPIC-8v1-0_B4.body.csv

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#         Pre-processing:: Load Orignal Docker EPICv1 Manifest
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

org_v1_prbs_tib <- NULL
org_v1_prbs_csv <- file.path( opt$top_path, "data/manifests/methylation/workhorse_docker/EPIC-B4.manifest.sesame-base.cpg-sorted.csv.gz" )
org_v1_prbs_tib <- readr::read_csv( file = org_v1_prbs_csv, show_col_types = FALSE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#       Pre-processing:: Load Genome Studio EPIC v1/v2 Manifests
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_v1_full_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/infinium-methylationepic-v-1-0-b5-manifest-file.csv.gz" )
man_v2_full_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/EPIC-8v2-0_A1.csv.gz" )
man_v2_mask_csv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/EPIC-8v2-0_A1.mask.csv.gz" )

# man_epic_path <- safe_mkdir( dir = file.path( opt$top_path, "data/manifests/methylation/workhorse_docker/epic" ) )

man_epic_path <- safe_mkdir( dir = file.path( opt$out_path, "manifests" ) )

man_v1_ctls_csv <- file.path( man_epic_path, "EPIC-8v1-0_A1.ctl.csv.gz" )
man_v1_prbs_csv <- file.path( man_epic_path, "EPIC-8v1-0_A1.prb.csv.gz" )

man_v2_ctls_csv <- file.path( man_epic_path, "EPIC-8v2-0_A1.ctl.csv.gz" )
man_v2_prbs_csv <- file.path( man_epic_path, "EPIC-8v2-0_A1.prb.csv.gz" )
man_v3_prbs_csv <- file.path( man_epic_path, "EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz" )

man_v1_ctls_tib <- NULL
man_v1_prbs_tib <- NULL

man_v2_ctls_tib <- NULL
man_v2_prbs_tib <- NULL

probe_source <- "EPIC"

#
# Define Columns::
#
sel_ctls_vec <- NULL
sel_ctls_vec <- c("Probe_ID", "M", "U", "mask", "DESIGN", "COLOR_CHANNEL", "col", 
                  "Probe_Type", "Probe_Source", "Next_Base", "Probe_Design" )

sel_prbs_vec <- NULL
sel_prbs_vec <- c("Probe_ID", "M", "U", "mask", "DESIGN", "COLOR_CHANNEL", "col", 
                  "Probe_Type", "Probe_Source", "Next_Base", "Probe_Design",
                  "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")

if ( file.exists(man_v1_ctls_csv) &&
     file.exists(man_v1_prbs_csv) &&
     file.exists(man_v2_ctls_csv) &&
     file.exists(man_v2_prbs_csv) ) {
  
  
} else {
  #
  # Load EPICv2 Controls::
  #
  man_v0_ctls_tib <- NULL
  man_v0_ctls_tib <- org_v1_prbs_tib %>% 
    dplyr::filter( Probe_Type != "cg" ) %>% 
    dplyr::filter( Probe_Type != "ch" ) %>% 
    dplyr::filter( Probe_Type != "rs" ) %>%
    dplyr::mutate( mask = FALSE ) %>%
    dplyr::select( dplyr::all_of(sel_ctls_vec) ) %>%
    dplyr::mutate( Probe_Design = Probe_Design %>% as.integer() )

  readr::write_csv( x = man_v0_ctls_tib, file = man_v1_ctls_csv )
  readr::write_csv( x = man_v0_ctls_tib, file = man_v2_ctls_csv )
  
  #
  # Load EPICv1::
  #
  man_v1_full_dat <- NULL
  man_v1_full_dat <- load_genome_studio_manifest( 
    file = man_v1_full_csv, 
    load_clean = TRUE,
    load_controls = TRUE, 
    cols_convert = FALSE,
    write_clean = FALSE,
    overwrite = FALSE,
    ret_data = TRUE,
    vb=vb,vt=vt,tc=tc,tt=tt )
  
  # man_v1_ctls_tib <- man_v1_full_dat$controls
  man_v1_prbs_tib <- man_v1_full_dat$probes
  
  #
  # Load EPICv2::
  #
  man_v2_full_dat <- NULL
  man_v2_full_dat <- load_genome_studio_manifest( 
    file = man_v2_full_csv, 
    load_clean = TRUE,
    load_controls = TRUE, 
    cols_convert = FALSE,
    write_clean = FALSE,
    overwrite = FALSE,
    ret_data = TRUE,
    vb=vb,vt=vt,tc=tc,tt=tt )
  
  # man_v2_ctls_tib <- man_v2_full_dat$controls
  man_v2_prbs_tib <- man_v2_full_dat$probes
  
  #
  # Load EPICv2 Maked Probes::
  #
  man_v2_mask_tib <- NULL
  man_v2_mask_tib <- readr::read_csv( file = man_v2_mask_csv, show_col_types = FALSE )
  
  #
  # Format Datasets::
  #
  man_v1_prbs_tib2 <- NULL
  man_v1_prbs_tib2 <- man_v1_prbs_tib %>%
    dplyr::rename( Probe_ID = IlmnID, U = AddressA_ID, M = AddressB_ID,
                   DESIGN = Infinium_Design_Type, COLOR_CHANNEL = Color_Channel
    ) %>%
    dplyr::mutate(
      Probe_Source = probe_source,
      Probe_Type = Probe_ID %>% stringr::str_sub( 1,2 ),
      Probe_Design = dplyr::case_when(
        DESIGN == "I"  ~ "1",
        DESIGN == "II" ~ "2",
        TRUE ~ "0"
      ),
      # Color_Channel
      col = stringr::str_sub( COLOR_CHANNEL, 1,1 ),
      col = dplyr::case_when(
        is.na(COLOR_CHANNEL) ~ "2",
        TRUE ~ col
      ),
      COLOR_CHANNEL = dplyr::case_when(
        is.na(COLOR_CHANNEL) ~ "Both",
        TRUE ~ COLOR_CHANNEL
      ),
      mask = FALSE
    ) %>%
    dplyr::select( dplyr::all_of(sel_prbs_vec) ) %>%
    dplyr::arrange( Probe_ID ) %>% 
    dplyr::rename( Full_ID = Probe_ID ) %>% 
    dplyr::left_join( org_v1_prbs_tib %>% dplyr::select( Probe_ID, M, U ), 
                      by=c("M","U") ) %>% 
    dplyr::select( dplyr::all_of( c(sel_prbs_vec, "Full_ID")) ) %>% 
    dplyr::arrange( Probe_ID ) %>%
    dplyr::mutate( Probe_Design = Probe_Design %>% as.integer() ) %>%
    dplyr::bind_rows( man_v0_ctls_tib )
  
  # Sanity Check::
  #  org_v1_prbs_tib %>% dplyr::filter( !U %in% man_v1_prbs_tib2$U )
  #  man_v1_prbs_tib2 %>% dplyr::filter( !U %in% org_v1_prbs_tib$U )
  #
  # tmp_tib <- man_v1_prbs_tib2 %>% 
  #   dplyr::left_join( org_v1_prbs_tib %>% dplyr::select( Probe_ID, M, U ), 
  #                     by=c("M","U") )
  # tmp_tib %>% dplyr::filter( Probe_ID != Full_ID )
  
  readr::write_csv( x = man_v1_prbs_tib2, file = man_v1_prbs_csv )
  
  man_v2_prbs_tib2 <- NULL
  man_v2_prbs_tib2 <- man_v2_prbs_tib %>% 
    dplyr::rename( Probe_ID = IlmnID, U = AddressA_ID, M = AddressB_ID,
                   DESIGN = Infinium_Design_Type, COLOR_CHANNEL = Color_Channel,
                   Probe_Design = Infinium_Design
    ) %>%
    dplyr::mutate(
      Probe_Source = probe_source,
      # col = dplyr::case_when(
      #   is.na(col) ~ "2",
      #   TRUE ~ col
      # ),
      COLOR_CHANNEL = dplyr::case_when(
        is.na(COLOR_CHANNEL) ~ "Both",
        TRUE ~ COLOR_CHANNEL
      ),
      mask = FALSE
    ) %>%
    dplyr::select( dplyr::all_of(sel_prbs_vec) ) %>%
    dplyr::arrange( Probe_ID )
  
  man_v2_mask_tib2 <- NULL
  man_v2_mask_tib2 <- man_v2_mask_tib %>% 
    dplyr::rename( Probe_ID = IlmnID, U = AddressA_ID, M = AddressB_ID,
                   DESIGN = Infinium_Design_Type, COLOR_CHANNEL = Color_Channel,
                   Probe_Design = Infinium_Design
    ) %>%
    dplyr::mutate(
      Probe_Source = probe_source,
      # col = dplyr::case_when(
      #   is.na(col) ~ "2",
      #   TRUE ~ col
      # ),
      COLOR_CHANNEL = dplyr::case_when(
        is.na(COLOR_CHANNEL) ~ "Both",
        TRUE ~ COLOR_CHANNEL
      ),
      mask = TRUE
    ) %>%
    dplyr::select( dplyr::all_of(sel_prbs_vec) ) %>%
    dplyr::arrange( Probe_ID )
  
  epic_v2_all_tib2 <- NULL
  epic_v2_all_tib2 <- man_v2_mask_tib2 %>%
    dplyr::bind_rows( man_v2_prbs_tib2 ) %>% 
    dplyr::distinct( Probe_ID, .keep_all = TRUE ) %>%
    dplyr::mutate( 
      Full_ID = Probe_ID, 
      Loci_ID = Probe_ID %>% stringr::str_remove("_.*$"),
      Rep_Num = Full_ID %>% stringr::str_remove( "^.*_[TB][CO][12]") %>% as.integer() ) %>% 
    dplyr::arrange( Probe_ID, -Rep_Num ) %>% 
    dplyr::select( -Probe_ID ) %>%
    dplyr::left_join( dplyr::select( man_v1_prbs_tib2,Probe_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq ), 
                      by=c("AlleleA_ProbeSeq","AlleleB_ProbeSeq")
    ) %>%
    dplyr::mutate(
      Probe_ID = dplyr::case_when(
        !is.na( Probe_ID ) & Rep_Num == 1 ~ Loci_ID,
        TRUE ~ Full_ID
      ),
      Probe_ID = dplyr::case_when(
        U == 71611891 ~ Full_ID,
        TRUE ~ Probe_ID
      ),
      Under_ID = Probe_ID %>% stringr::str_replace_all("\\.", "_"), 
      Probe_ID = dplyr::case_when( 
        Probe_Type == "ch" ~ Under_ID, 
        TRUE ~ Probe_ID )
    ) %>%
    dplyr::select( dplyr::all_of( c(sel_prbs_vec, "Full_ID", "Loci_ID", "Rep_Num")) ) %>%
    dplyr::arrange( Probe_ID ) %>% 
    dplyr::mutate( Probe_Design = Probe_Design %>% as.integer() ) %>%
    dplyr::bind_rows( man_v0_ctls_tib )
  
  # %>% dplyr::mutate( Under_ID = Probe_ID %>% stringr::str_replace_all("\\.", "_"), Probe_ID = dplyr::case_when( Probe_Type == "ch" ~ Under_ID, TRUE ~ Probe_ID ) )
  # dplyr::filter( man_v2_prbs_tib, Probe_Type == "ch" ) %>% dplyr::mutate( Under_ID = Probe_ID %>% stringr::str_replace_all("\\.", "_"), Probe_ID = dplyr::case_when( Probe_Type == "ch" ~ Under_ID, TRUE ~ Probe_ID ) ) %>% head() %>% as.data.frame()
  
  # Sanity Checks::
  #  v2_cnt_tib <- epic_v2_all_tib2 %>% dplyr::add_count( Probe_ID, name="Rep_Cnt" )
  #  v2_cnt_tib %>% dplyr::filter( Rep_Cnt != 1 ) %>% head() %>% as.data.frame()
  #  epic_v2_all_tib2 %>% dplyr::filter( Loci_ID == "cg00002033" ) %>% as.data.frame()
  
  readr::write_csv( x = epic_v2_all_tib2, file = man_v2_prbs_csv )
  
  #
  # Short verions as previously used...
  #
  
  epic_v2_all_tib3 <- NULL
  epic_v2_all_tib3 <- epic_v2_all_tib2 %>% 
    dplyr::select( Probe_ID,M,U,DESIGN,COLOR_CHANNEL,col,Probe_Type,
                   Probe_Source,Next_Base,Probe_Design,mask ) %>%
    dplyr::mutate( Probe_Source = "EPIC-A1" )
  
  readr::write_csv( x = epic_v2_all_tib3, file = man_v3_prbs_csv )
  
  if ( p1 ) cat(glue::glue("{pmssg} Done writing clean files!{RET2}") )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#     Pre-processing:: Re-Load Genome Studio EPIC v1/v2 Manifests
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

man_v1_ctls_tib <- NULL
man_v1_prbs_tib <- NULL

man_v2_ctls_tib <- NULL
man_v2_prbs_tib <- NULL

man_v1_prbs_tib <- readr::read_csv( file = man_v1_prbs_csv, show_col_types = FALSE )
man_v2_prbs_tib <- readr::read_csv( file = man_v2_prbs_csv, show_col_types = FALSE )

man_v1_ctls_tib <- readr::read_csv( file = man_v1_ctls_csv, show_col_types = FALSE )
man_v2_ctls_tib <- readr::read_csv( file = man_v2_ctls_csv, show_col_types = FALSE )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#     Pre-processing:: Validate Original Docker vs. New Manifests
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

org_v1_prbs_tib %>% dplyr::anti_join( man_v1_prbs_tib, by=c("Probe_ID") )
man_v1_prbs_tib %>% dplyr::anti_join( org_v1_prbs_tib, by=c("Probe_ID") )

org_v1_prbs_tib %>% dplyr::anti_join( man_v1_prbs_tib, by=c("Probe_ID") ) %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
man_v1_prbs_tib %>% dplyr::anti_join( org_v1_prbs_tib, by=c("Probe_ID") ) %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )

org_v1_prbs_tib %>% dplyr::anti_join( man_v2_prbs_tib, by=c("Probe_ID") ) %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
man_v2_prbs_tib %>% dplyr::anti_join( org_v1_prbs_tib, by=c("Probe_ID") ) %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )


org_v1_prbs_tib %>% dplyr::filter( Probe_Type == "rs" )
man_v1_prbs_tib %>% dplyr::filter( Probe_Type == "rs" )
man_v2_prbs_tib %>% dplyr::filter( Probe_Type == "rs" )

org_v1_prbs_tib %>% dplyr::filter( Probe_Type == "ch" )
man_v1_prbs_tib %>% dplyr::filter( Probe_Type == "ch" )
man_v2_prbs_tib %>% dplyr::filter( Probe_Type == "ch" )



#
# Merge EPICv1 + EPICv2 
#

# rs/ch need to be handled differently

# dplyr::inner_join( man_v2_prb_tib, man_v2_prb_tib, by=c()


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                     Post-processing:: Write Manifest
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

sel_sub_vec <- NULL
sel_sub_vec <- c("Probe_ID", "M", "U", "DESIGN", "COLOR_CHANNEL", "col", 
                 "Probe_Type", "Probe_Source", "Next_Base", "Probe_Design" ) # "mask")

man_out_dir <- safe_mkdir( dir = file.path( opt$out_path, "manifest" ) )

man_v2_sub_csv <- file.path( man_out_dir, "EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
man_v2_sub_tib <- NULL
man_v2_sub_tib <- man_v2_prbs_tib %>% 
  dplyr::mutate( Probe_Source = "EPIC-A1" ) %>% 
  dplyr::select( dplyr::all_of(sel_sub_vec) )

readr::write_csv( x = man_v2_sub_tib, file = man_v2_sub_csv )

# man_v2_sub_tib %>% dplyr::filter( Probe_ID %>% stringr::str_detect("_") )






# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Post-analysis:: Load Sample-Sheets/Beta-Pvals
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# TBD:: Use standard pair plotting::
#

if ( FALSE ) {
  
  epic_v12_ssh_tib <- NULL
  epic_v12_ssh_csv <- file.path( opt$top_path, "scratch/stable_epic_v2_manifest_plots/GSIBIOINFO-638-v1/build_epicv2_sample_sheet/GSIBIOINFO-638-v1.build_epicv2_sample_sheet.csv.gz" )
  epic_v12_ssh_tib <- readr::read_csv( file = epic_v12_ssh_csv, show_col_types = FALSE )
  epic_v12_ssh_list <- epic_v12_ssh_tib %>% split(.$Sample_Base)
  
  dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/docker-1.11.15.1.p.0.4/EPICv2.1.11.15.1.p.0.4/swifthoof_main" )
  ssh_list <- NULL
  ssh_list <- file_list( path = dat_path, # prefix = "",
                          suffix  = "_EPIC_.*_AutoSampleSheet.csv.gz", 
                          pattern = "_AutoSampleSheet.csv.gz$" )
  
  beta_list <- NULL
  beta_list <- file_list( path = dat_path, # prefix = "",
                          suffix  = "_EPIC_.*_raw.call.dat.csv.gz", 
                          pattern = "_raw.call.dat.csv.gz$" )
  
  
  tar_sam_vec <- c("HeLa", "Jurkat", "MCF7")
  for ( sample in tar_sam_vec ) {
    cur_tib <- epic_v12_ssh_list[[ sample ]] %>% dplyr::filter( Sheet_Prep == "Lightning" & Sheet_Proc == "Auto" )
    
    cur_ssh_list <- ssh_list[ names(ssh_list) %>% intersect( cur_tib$Sentrix_Name ) ]
    cur_beta_list <- beta_list[ names(ssh_list) %>% intersect( cur_tib$Sentrix_Name ) ]
    
    #
    # Need to make the beta_* cur_beta_*
    #
    beta_dat <- NULL
    beta_dat <- beta_list %>% # head() %>% 
      lapply( readr::read_csv, show_col_types = FALSE )
    
    beta_tab <- NULL
    beta_tab <- beta_dat %>%
      dplyr::bind_rows( .id = "Sentrix_ID" )
    
    beta_tib <- NULL
    beta_tib <- beta_tab %>%
      dplyr::select( Sentrix_ID, Probe_ID, betas ) %>%
      tidyr::pivot_wider( id_cols = c(Probe_ID), 
                          names_from = c(Sentrix_ID), 
                          values_from = c( betas ) )
    
    
    break
  }

  
  
  
  
  
  
  
  #
  # OLD CODE::
  #
  dat_path <- file.path( opt$top_path, "Projects.new/EPIC_v2/scratch/docker-1.11.15.1.p.0.4/EPICv2.1.11.15.1.p.0.4/swifthoof_main" )
  
  ssh_list <- NULL
  ssh_list <- file_list( path = dat_path, # prefix = "",
                         suffix  = "_AutoSampleSheet.csv.gz", 
                         pattern = "_AutoSampleSheet.csv.gz$" )
  
  ssh_list %>% length()
  
  ssh_tab <- NULL
  ssh_tab <- ssh_list %>% # head() %>% 
    lapply( readr::read_csv, show_col_types = FALSE ) %>%
    dplyr::bind_rows( )
  
  ssh_tab_sum2 <- NULL
  ssh_tab_sum2 <- ssh_tab %>% 
    dplyr::filter( AutoSample_R2_2_Key_1 == "HELA" |
                     AutoSample_R2_2_Key_1 == "JURKAT" |
                     AutoSample_R2_2_Key_1 == "RAJI" ) %>%
    dplyr::group_by( AutoSample_R2_2_Key_1, detect_manifest ) %>% 
    dplyr::summarise( Avg_Pass2 = median(cg_calls_pass_perc_1),
                      Med_Pass2 = median(cg_calls_pass_perc_1),
                      Cnt = n(),
                      .groups = "drop" ) %>%
    dplyr::filter( Cnt > 3 )

  ssh_tab_sum1 <- NULL
  ssh_tab_sum1 <- ssh_tab %>% 
    dplyr::filter( AutoSample_dB_1_Key_1 == "HELA" |
                     AutoSample_dB_1_Key_1 == "JURKAT" |
                     AutoSample_dB_1_Key_1 == "RAJI" ) %>%
    dplyr::group_by( AutoSample_dB_1_Key_1, detect_manifest ) %>% 
    dplyr::summarise( Avg_Pass1 = median(cg_calls_pass_perc_1),
                      Med_Pass1 = median(cg_calls_pass_perc_1),
                      Cnt = n(),
                      .groups = "drop" ) %>%
    dplyr::filter( Cnt > 3 )
  
  # ssh_tab %>% 
  #   dplyr::group_by( AutoSample_R2_2_Key_2, detect_manifest ) %>% 
  #   dplyr::summarise( Avg_Pass = median(cg_calls_pass_perc_2),
  #                     Med_Pass = median(cg_calls_pass_perc_2),
  #                     Cnt = n(),
  #                     .groups = "drop" ) %>%
  #   dplyr::filter( Cnt > 3 )
  

  beta_list <- NULL
  beta_list <- file_list( path = dat_path, # prefix = "",
                          suffix  = "_raw.call.dat.csv.gz", 
                          pattern = "_raw.call.dat.csv.gz$" )
  
  beta_dat <- NULL
  beta_dat <- beta_list %>% # head() %>% 
    lapply( readr::read_csv, show_col_types = FALSE )
  
  beta_tab <- NULL
  beta_tab <- beta_dat %>%
    dplyr::bind_rows( .id = "Sentrix_ID" )
  
  beta_tib <- NULL
  beta_tib <- beta_tab %>%
    dplyr::select( Sentrix_ID, Probe_ID, betas ) %>%
    tidyr::pivot_wider( id_cols = c(Probe_ID), 
                        names_from = c(Sentrix_ID), 
                        values_from = c( betas ) )
  
}







# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                      Pre-analysis:: OLD-CODE-BELOW
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {
  
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
  #                  Pre-processing:: Load Join Manifest v1_v2
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  src_12_dir <- file.path( opt$top_path, "data/manifests/methylation/Sesame/genome_studio_conversion/full" )
  
  src_v1_tib <- NULL
  src_v2_tib <- NULL
  src_12_tib <- NULL
  
  src_v1_csv <- file.path( src_12_dir, "EPIC_v1.gs_to_sesame.csv.gz" )
  src_v2_csv <- file.path( src_12_dir, "EPIC_v2.gs_to_sesame.csv.gz" )
  src_12_csv <- file.path( src_12_dir, "EPIC_v1.gs_to_sesame.csv.gz" )
  
  src_v1_tib <- readr::read_csv( file = src_v1_csv, show_col_types = FALSE )
  src_v2_tib <- readr::read_csv( file = src_v2_csv, show_col_types = FALSE )
  src_12_tib <- readr::read_csv( file = src_12_csv, show_col_types = FALSE)
  
  src_v2_tib
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Pre-processing:: Load Manifest v2
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  raw_v2_tib <- NULL
  man_v2_ver <- "unq"
  man_v2_ver <- "raw"
  
  man_v2_dir <- file.path( opt$top_path, "data/manifests/methylation/stable_docker_manifest/core6/manifests", man_v2_ver  )
  man_v2_csv <- file.path( man_v2_dir, "EPICv2-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
  
  raw_v2_tib <- readr::read_csv( file = man_v2_csv, show_col_types = FALSE )
  
  src_raw_anti_tib <- src_v2_tib %>% dplyr::anti_join( raw_v2_tib, by=c("Probe_ID") )
  raw_src_anti_tib <- raw_v2_tib %>% dplyr::anti_join( src_v2_tib, by=c("Probe_ID") )
  
  sNr_v2_tib <- NULL
  sNr_v2_tib <- dplyr::bind_rows( 
    src_12_tib %>% dplyr::select( Probe_ID ),
    raw_v2_tib %>% dplyr::select( Probe_ID )
  ) %>% 
    dplyr::distinct( Probe_ID )
  
  raw_v2_cnt <- raw_v2_tib %>% base::nrow()
  src_v1_cnt <- src_v1_tib %>% base::nrow()
  src_v2_cnt <- src_v2_tib %>% base::nrow()
  src_12_cnt <- src_12_tib %>% base::nrow()
  sNr_v2_cnt <- sNr_v2_tib %>% base::nrow()
  
  cat(glue::glue("{pmssg} Counts By Manifest::{RET2}",
                 "{pmssg} RAW: {raw_v2_cnt}.{RET}",
                 "{pmssg} Sv1: {src_v1_cnt}.{RET}",
                 "{pmssg} Sv2: {src_v2_cnt}.{RET}",
                 "{pmssg} SRC: {src_12_cnt}.{RET}",
                 "{pmssg} sNr: {sNr_v2_cnt}.{RET2}") )
  
  all_v2_tib <- raw_v2_tib %>% 
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
    ) %>% 
    dplyr::left_join( all_suf_tib, by=c("Loci_ID"="Probe_ID") ) %>% 
    dplyr::distinct( Probe_ID, .keep_all = TRUE )
  
  raw_mis_tib <- NULL
  raw_mis_tib <- src_v2_tib %>% 
    dplyr::anti_join( raw_v2_tib, by=c("Probe_ID", "U","M") ) %>% 
    dplyr::inner_join( raw_v2_tib, by=("U"), suffix=c("_SRC","_RAW") )
  
  # raw_mis_tib %>% head() %>% as.data.frame()
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                 Pre-processing:: Make Suffix Adjustment
  #                               Manifest v2
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  man_v2_tib <- NULL
  man_v2_tib <- raw_v2_tib %>% 
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
  out_core_ver <- "core9"
  
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
  #                         Scratch Compare Results::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( FALSE ) {
    
    top_dat_dir <- "/Users/bbarnes/Documents/Projects.new/EPIC_v2"
    
    v1_file_list <- NULL
    v1_file_path <- file.path( top_dat_dir, "scratch/docker/EPICv1.206203800149.1.11.15.1.p.1/swifthoof_main" )
    v1_file_path <- file.path( top_dat_dir, "scratch/docker/EPICv1.206203800149.1.11.15.1.p.1.auto2/swifthoof_main" )
    v1_file_list <- file_list( path = v1_file_path, 
                               suffix  = "_EPIC_B4_raw.call.dat.csv.gz",
                               pattern = "_EPIC_B4_raw.call.dat.csv.gz$"
    )
    
    v2_file_list <- NULL
    v2_file_path <- file.path( top_dat_dir, "scratch/docker/EPICv2.206891110001.1.11.15.1.p.15.auto2/swifthoof_main" )
    # v2_file_list <- list.files( path = v2_file_path, pattern = "EPIC_A1_raw.call.dat.csv.gz", full.names = TRUE )
    v2_file_list <- file_list( path = v2_file_path, 
                               suffix  = "_EPIC_A1_raw.call.dat.csv.gz",
                               pattern = "_EPIC_A1_raw.call.dat.csv.gz$"
    )
    
    #
    # Load Data::
    #
    v1_dat_list <- NULL
    v1_dat_list <- lapply( v1_file_list, readr::read_csv, show_col_types = FALSE ) %>% 
      dplyr::bind_rows( .id = "Sentrix_ID" ) %>% 
      dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_remove("_.*$") )
    
    v2_dat_list <- NULL
    v2_dat_list <- lapply( v2_file_list, readr::read_csv, show_col_types = FALSE ) %>% 
      dplyr::bind_rows( .id = "Sentrix_ID" ) %>% 
      dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_remove("_.*$") )
    
    betas_tib <- NULL
    betas_tib <- dplyr::bind_rows( v1_dat_list,v2_dat_list ) %>% 
      dplyr::select( Sentrix_ID, Probe_ID, betas ) %>%
      tidyr::pivot_wider( id_cols = c(Probe_ID), 
                          names_from = c(Sentrix_ID), 
                          values_from = c( betas ) )
    # values_from = c(pvals_pOOBAH, pvals_PnegEcdf, betas) )
    
    #
    # Single Example::
    #
    v1_tib <- v1_file_list[[1]] %>% readr::read_csv( show_col_types = FALSE )
    v2_tib <- v2_file_list[[1]] %>% readr::read_csv( show_col_types = FALSE ) %>%
      dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_remove("_.*$") ) %>%
      dplyr::distinct( Probe_ID, .keep_all = TRUE )
    
    beta_tib1 <- v1_tib %>% dplyr::inner_join( v2_tib, by=c("Probe_ID"), suffix = c("_v1", "_v2") )
    
    beta_tib1 %>% head( n=10000) %>% 
      ggplot2::ggplot( aes(x=betas_v1, y=betas_v2) ) + 
      ggplot2::geom_point()
  }
  
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
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt+3,tc=tc,tt=tt )

sysTime <- Sys.time()
if ( p0 ) cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))


# End of file
