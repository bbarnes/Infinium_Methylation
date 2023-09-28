
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
par$version <- 3
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
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt+3,tc=tc,tt=tt )

sysTime <- Sys.time()
if ( p0 ) cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
