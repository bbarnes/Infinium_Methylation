
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
par$prgm_tag <- 'stable_format_controls_human'
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

opt <- NULL
# opt <- stable_options( pars = par, args = args, vb = par$verbose )
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

# rcppp
# opt$rcpp <- 0
# opt$rcpp <- 2
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
#                 Pre-processing:: Reference Genomes
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

gen_path <- NULL
gen_path <-
  file_list( path = opt$ref_path,
             subs = opt$ref_species,
             unique = FALSE,
             dir_only = TRUE,
             ret_type = COM,
             subs_exists = FALSE,
             vb = vb ) %>%
  file_list( subs = opt$ref_source,
             unique = FALSE,
             dir_only = TRUE,
             ret_type = COM,
             subs_exists = FALSE,
             vb = vb ) %>%
  file_list( subs = opt$ref_build,
             unique = FALSE,
             dir_only = TRUE,
             ret_type = COM,
             subs_exists = FALSE,
             paths_exists = FALSE,
             vb = vb )

dir_add_str <- NULL
if ( opt$ref_source == "NCBI") dir_add_str <- "Fasta"

ref_seqs <- NULL
ref_seqs <-
  file_list( path = gen_path,
             subs = paste0("Sequence/WholeGenome",dir_add_str),
             file = opt$ref_file,
             unique = FALSE,
             suffix = c("\\.gz", "\\.fa$", "\\.genome$"),
             ret_type = "list",
             subs_exists = FALSE,
             paths_exists = FALSE,
             files_exists = FALSE,
             vb = vb )

chr_dirs <- NULL
chr_dirs <-
  file_list( path = gen_path,
             subs = "Sequence/Chromosomes",
             names = names(ref_seqs),
             unique = FALSE,
             dir_only = TRUE,
             ret_type = "list",
             subs_exists = FALSE,
             paths_exists = FALSE,
             vb = vb )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                     Pre-processing:: Sample Sheets
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

ssh_rep_tib <- NULL
ssh_rep_csv <- file.path( opt$top_path, "Projects.new/MSA/MethylationQC/MSAv03/post_processing/samplesheet.csv.gz" )
ssh_rep_tib <- readr::read_csv( file = ssh_rep_csv, show_col_types = FALSE ) %>%
  dplyr::mutate( Sample_Concentration = Sample_Group %>% 
                   stringr::str_remove(" .*$") %>% as.integer(),
                 Base_Sample_Name = Sample_Well %>% stringr::str_to_upper(),
                 Sample_Name = paste( Base_Sample_Name,Sample_Concentration, sep="_"),
                 Sentrix_Name = Sample_ID, # paste(Sentrix_ID,Sentrix_Position, sep="_"),
                 Sentrix_ID = Sentrix_Name %>% stringr::str_remove("_.*$"),
                 Sentrix_Path = paste0( opt$top_path,"/data/idats/idats_MSAv03_48x1_Alpha/",Sentrix_ID,"/",Sentrix_Name ) ) %>%
  dplyr::select( Sentrix_Name, Sample_Name,Base_Sample_Name,Sample_Concentration, Sentrix_Path ) %>% 
  dplyr::mutate( Select_Group := par$run_name )

ssh_rep_sum <- NULL
ssh_rep_sum <- ssh_rep_tib %>% 
  dplyr::group_by( Base_Sample_Name, Sample_Concentration ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p1 ) ssh_rep_sum %>% print( n=base::nrow(ssh_rep_sum) )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                       Pre-processing:: Manifests (BGZ)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

all_sub_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_manifests.sub8.rds")
rev_sub_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/rev_manifests.sub8.rds")
neg_ctl_rds  <- file.path( opt$top_path, "data/manifests/methylation/bgz/all_negative_ctls.rds" )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#              Pre-processing:: Negative Controls (BGZ)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

neg_ctl_tib  <- NULL
neg_ctl_tib  <- readr::read_rds( neg_ctl_rds )
neg_org_tib  <- readr::read_rds( neg_ctl_rds )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: EPIC Controls (BGZ)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

epic_ctl_tib <- NULL
epic_ctl_rds <- file.path( opt$top_path, "data/manifests/methylation/bgz/epic_ctls.rds" )
epic_ctl_tib <- readr::read_rds( file = epic_ctl_rds ) %>%
  dplyr::select( Probe_ID,U,Type ) %>% 
  dplyr::distinct( U, .keep_all = TRUE ) %>% 
  dplyr::rename( Probe_ID_epic = Probe_ID, 
                 Probe_Type_epic = Type ) %>%
  clean_tib()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#             Pre-processing:: Example MSA GS Controls (ACA)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

example_gs_ctl_tib <- NULL
example_gs_ctl_csv <- file.path( opt$top_path, "Projects.new/MSA/MethylationQC/MSAv03/manifest/MSAv03_controls.csv.gz" )
example_gs_ctl_tib <- readr::read_csv( file = example_gs_ctl_csv, 
                                       show_col_types = FALSE ) %>% clean_tib()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#         Pre-processing:: MSA GS Manifest Partitions (ACA)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

msa_gs_ctls_src_col <- NULL
msa_gs_ctls_src_col <- readr::cols(
  Address = readr::col_integer(),
  Type = readr::col_character(),
  Color_Channel = readr::col_character(),
  Name = readr::col_character()
  # Probe_ID = readr::col_character(), 
)

msa_gs_full_src_csv <- "/Users/bbarnes/Documents/data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_BP123_A1.csv.gz"
msa_gs_head_src_csv <- "/Users/bbarnes/Documents/data/manifests/methylation/EX/MSAEX03/genome-studio-partions/MSA-Interm-48v0-3_BP123_A1.gs-header.csv"
msa_gs_ctls_src_csv <- "/Users/bbarnes/Documents/data/manifests/methylation/EX/MSAEX03/genome-studio-partions/MSA-Interm-48v0-3_BP123_A1.gs-controls.csv"
msa_gs_body_src_csv <- "/Users/bbarnes/Documents/data/manifests/methylation/EX/MSAEX03/genome-studio-partions/MSA-Interm-48v0-3_BP123_A1.gs-analytical.csv"

msa_gs_head_src_line <- NULL
if ( !file.exists(msa_gs_head_src_csv) ) {
  msa_gs_head_parse_cmd <- paste0("gzip -dc ",msa_gs_full_src_csv," | head -n 7 > ",msa_gs_head_src_csv)
  base::system( msa_gs_head_parse_cmd )
}
if ( file.exists(msa_gs_head_src_csv) ) msa_gs_head_src_line <- 
  readr::read_lines( file = msa_gs_head_src_csv )

msa_gs_ctls_src_tib <- NULL
msa_gs_ctls_src_line <- NULL
if ( !file.exists(msa_gs_ctls_src_csv) ) {
  msa_gs_ctls_parse_cmd <- paste0("gzip -dc ",msa_gs_full_src_csv," | tail -n 642 > ",msa_gs_ctls_src_csv)
  base::system( msa_gs_ctls_parse_cmd )
}
if ( file.exists(msa_gs_ctls_src_csv) ) msa_gs_ctls_src_tib  <- 
  readr::read_csv( file = msa_gs_ctls_src_csv, 
                   col_names = names(msa_gs_ctls_src_col$cols), 
                   col_types = msa_gs_ctls_src_col, # lazy = TRUE,
                   skip = 1 ) %>% 
  dplyr::select( dplyr::all_of( names(msa_gs_ctls_src_col$cols) ) ) %>%
  clean_tib()
if ( file.exists(msa_gs_ctls_src_csv) ) msa_gs_ctls_src_line <- 
  readr::read_lines( file = msa_gs_ctls_src_csv )

#
# TBD:: Compare [ example_gs_ctl_tib vs. msa_gs_ctls_src_tib ]
#

msa_gs_body_src_tib <- NULL
if ( !file.exists(msa_gs_body_src_csv) ) {
  msa_gs_body_parse_cmd <- paste0("gzip -dc ",msa_gs_full_src_csv," | tail -n 77851 | head -n 77209 > ",msa_gs_body_src_csv)
  base::system( msa_gs_body_parse_cmd )
}
if ( file.exists(msa_gs_body_src_csv) ) msa_gs_body_src_tib <- 
  readr::read_csv( file = msa_gs_body_src_csv, show_col_types = FALSE ) %>%
  clean_tib()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                                   BEG
#                            Controls Format::
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

format_ctls <- FALSE
format_ctls <- TRUE

if ( format_ctls ) {
  if ( p1 ) cat(glue::glue("{pmssg} BEG: Formatting Raw Controls...{RET}"))
  
  #
  # This is the correct stuff!!!
  #  [Done] Need to copy this function:: Workhorse-Unstained/scripts/R/trifecta/functions/control_functions.R
  #  - Need to create new addControlType() that follows this naming scheme full_ctl_tib
  #
  full_ctl_tib <- NULL
  full_ctl_csv <- file.path( opt$top_path, "tools/docker/Infinium_Methylation_Workhorse/dat/manifest/controls/Infinium_Methylation_Controls_1983_full.csv.gz" )
  full_ctl_tib <- readr::read_csv( file = full_ctl_csv, 
                                   show_col_types = FALSE ) %>% 
    dplyr::distinct( U,Probe_Seq, .keep_all = TRUE ) %>% 
    dplyr::select( Probe_ID,U,Probe_Seq,
                   Control_Group,Control_Name,Sample_Dependent,
                   Order_ID,PIDX,Control_Group_Str,Probe_Type ) %>%
    dplyr::rename( Probe_ID_full = Probe_ID, 
                   Control_Type_full = Control_Group, 
                   Control_Name_full = Control_Name, 
                   Order_ID_full = Order_ID, 
                   Control_Group_Str_full = Control_Group_Str, 
                   PIDX_full = PIDX, 
                   Probe_Type_full = Probe_Type, 
                   Sample_Dependent_full = Sample_Dependent ) %>% 
    dplyr::filter( !( U == 22711390 | U == 46651360 ) ) %>% clean_tib()
  
  #
  # NOTE: Not real confident why these tango addresses show up with
  #.  Sample Dependent Controls, so removing them above...
  #
  # A tibble: 2 × 10
  # Probe_ID_full                 U Probe_Seq Control_Type_full      Control_Name_full  Sample_Dependent_full Order_ID_full PIDX_full Control_Group_Str_full Probe_Type_full
  # <chr>                     <dbl> <chr>     <chr>                  <chr>              <chr>                 <chr>         <chr>     <chr>                  <chr>          
  # 1 ctl_BS_Conversion_I_C1 22711390 NA        BISULFITE CONVERSION I BS Conversion I C1 Dependent             NA            1         BISULFITE_CONVERSION_I BS             
  # 2 ctl_BS_Conversion_I_U1 46651360 NA        BISULFITE CONVERSION I BS Conversion I U1 Dependent             NA            1         BISULFITE_CONVERSION_I BS     
  
  # Sanity Sequence Duplication Check::
  #. full_ctl_tib %>% dplyr::add_count( Probe_Seq, name="Dup_Seq_Cnt" ) %>% dplyr::filter( Dup_Seq_Cnt != 1 ) %>% as.data.frame()
  # Sanity Tango Address Duplication Check::
  #. full_ctl_tib %>% dplyr::add_count( U, name="Dup_U_Cnt" ) %>% dplyr::filter( Dup_U_Cnt != 1 ) %>% as.data.frame()
  #
  
  full_ctl_sum <- NULL
  full_ctl_sum <- full_ctl_tib %>% 
    dplyr::group_by( Control_Type_full ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) cat(glue::glue("{pmssg} Full Control Summary::{RET}"))
  if ( p4 ) full_ctl_sum %>% print( n=nrow(full_ctl_sum) )
  if ( p4 ) cat(glue::glue("{RET}"))
  
  #
  # Original Match File::
  #
  match_ctl_tib <- NULL
  match_ctl_tsv <- file.path( opt$top_path, "data/Controls/01152015_DarkMatterControls.unique.probe.match.tsv.gz" )
  match_ctl_tib <- readr::read_tsv( file = match_ctl_tsv, show_col_types = FALSE ) %>%
    dplyr::rename( Probe_ID = probe_id,
                   Probe_Seq = sequence ) %>%
    dplyr::mutate( U = address_name %>% stringr::str_remove("^1") %>%
                     stringr::str_remove("^0+") %>%
                     as.integer(),
                   Probe_ID = Probe_ID %>% 
                     stringr::str_remove("_1_A") %>%
                     stringr::str_remove("_1_B")
    ) %>% 
    dplyr::select( Probe_ID, U, Probe_Seq ) %>%
    dplyr::distinct( U, .keep_all = TRUE ) %>%
    dplyr::distinct( Probe_Seq, .keep_all = TRUE ) %>%
    addControlType() %>%
    addControlTypeSlim() %>%
    dplyr::mutate( 
      Is_EPIC = dplyr::case_when(
        U %in% epic_ctl_tib$U ~ TRUE,
        TRUE ~ FALSE ),
      Last_Base = Probe_Seq %>% stringr::str_sub( 50,50 )
    ) %>% 
    dplyr::rename( Probe_ID_match = Probe_ID, 
                   Control_Type_match = Control_Type,
                   Control_Class_match = Control_Class,
                   Is_EPIC_match = Is_EPIC ) %>%
    dplyr::select( Probe_ID_match,U,Probe_Seq,Last_Base, 
                   dplyr::everything() ) %>% clean_tib()
  
  # SNP Check::
  match_ctl_snp_tib <- NULL
  match_ctl_snp_tib <- match_ctl_tib %>% 
    dplyr::filter( Probe_ID_match %>% stringr::str_starts("rs") )
  match_ctl_snp_cnt <- base::nrow(match_ctl_snp_tib)
  if ( p2 ) cat(glue::glue("{pmssg} Control SNP Count='{match_ctl_snp_cnt}'{RET2}"))
  
  match_ctl_sum <- NULL
  match_ctl_sum <- match_ctl_tib %>% 
    dplyr::group_by( Control_Type_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) cat(glue::glue("{pmssg} Control Type::{RET}"))
  if ( p4 ) match_ctl_sum %>% print( n=nrow(match_ctl_sum) )
  if ( p4 ) cat(glue::glue("{RET}"))
  
  match_cls_sum <- NULL
  match_cls_sum <- match_ctl_tib %>% 
    dplyr::group_by( Control_Class_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) cat(glue::glue("{pmssg} Control Class::{RET}"))
  if ( p4 ) match_cls_sum %>% print( n=nrow(match_cls_sum) )
  if ( p4 ) cat(glue::glue("{RET}"))
  
  epic_ctl_sum <- NULL
  epic_ctl_sum <- match_ctl_tib %>% 
    dplyr::filter( Is_EPIC_match ) %>%
    dplyr::group_by( Control_Type_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) cat(glue::glue("{pmssg} EPIC Control Type::{RET}"))
  if ( p4 ) epic_ctl_sum %>% print( n=nrow(epic_ctl_sum) )
  if ( p4 ) cat(glue::glue("{RET}"))
  
  # Just Some Printing Sanity Checking Stuff
  #
  # match_ctl_tib %>% head() %>% as.data.frame()
  # full_ctl_tib  %>% head() %>% as.data.frame()
  # epic_ctl_tib  %>% head() %>% as.data.frame()
  #
  # NOTE:: Nomenclature for sufficies in file::
  #.   - _match => from original match file.
  #.   - _epic => from EPIC Sesame manifest. 
  #.   - _full => from Full list in Infinium Methylation Workhorse docker image
  #
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            Join all Datasets::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  join_ctl_dir <- safe_mkdir( file.path( opt$out_path, "join") )
  join_ctl_csv <- file.path( join_ctl_dir, "Infinium_Methylation_Controls_1971_join-names.raw.csv.gz" )
  join_ctl_csv <- file.path( join_ctl_dir, "Infinium_Methylation_Controls_1971_join-names.raw.csv" )
  
  join_ctl_tib <- NULL
  join_ctl_tib <- match_ctl_tib %>% 
    dplyr::left_join( epic_ctl_tib, by=c("U") ) %>%
    dplyr::full_join( full_ctl_tib, by=c("U","Probe_Seq") ) %>%
    dplyr::arrange( Control_Class_match ) %>% clean_tib()
  readr::write_csv( x = join_ctl_tib, file = join_ctl_csv )
  
  # NOTE:: Sample_Dependent_full is not trustworthy...
  #
  # join_ctl_tib %>% dplyr::group_by( Sample_Dependent_full, Control_Type_match ) %>%
  #   dplyr::summarise( Count=n(), .groups = "drop" )
  
  seqs_ctl_tib <- NULL
  seqs_ctl_tib <- join_ctl_tib %>% 
    dplyr::filter( !is.na(Probe_Seq) ) %>%
    dplyr::mutate( Probe_ID_bsp = paste0( Probe_Seq,"_",U ) )
  
  seqs_ctl_sum <- NULL
  seqs_ctl_sum <- seqs_ctl_tib %>% 
    dplyr::group_by( Control_Class_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) cat(glue::glue("{pmssg} Seqs Control Type Summary::{RET}"))
  if ( p4 ) seqs_ctl_sum %>% print( n=base::nrow(seqs_ctl_sum) )
  if ( p4 ) cat(glue::glue("{RET}"))
  
  #
  # Build New Categorie: "Control_Type" and compare against "Probe_Type_full"
  #. - Don't think this works because it base on previous EPIC data...
  type_ctl_sum <- NULL
  type_ctl_sum <- seqs_ctl_tib %>% 
    dplyr::group_by( Probe_Type_full ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) cat(glue::glue("{pmssg} Seqs Control Class Summary::{RET}"))
  if ( p4 ) type_ctl_sum %>% print( n=base::nrow(type_ctl_sum) )
  if ( p4 ) cat(glue::glue("{RET}"))
  
  if ( p1 ) cat(glue::glue("{pmssg} END: Formatting Raw Controls.{RET2}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Split/Join Infinium I Probes::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# NOTE:: This is currently only a helper function...
#

split_infs <- FALSE
split_infs <- TRUE

if ( split_infs ) {
  if ( p1 ) cat(glue::glue("{pmssg} BEG: Splitting Infinium I Controls...{RET}"))
  
  # Data Structures::
  #
  #.  join_ctl_tib  - Includes everything (Independent Controls as well)
  #.  seqs_ctl_tib  - join_ctl_tib with only Dependent
  #   all_class_tib - seqs_ctl_tib with reduced columns
  #
  
  # Conclusion: I don't think the Independent/Dependent means anything...
  join_dep_sum <- NULL
  join_dep_sum <- join_ctl_tib %>% 
    dplyr::group_by( Sample_Dependent_full, Control_Class_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) cat(glue::glue("{pmssg} (In)Dependent Control Type Summary::{RET}"))
  if ( p4 ) join_dep_sum %>% print( n=base::nrow(join_dep_sum) )
  if ( p4 ) cat(glue::glue("{RET}"))
  
  #
  # Next Steps::
  #   1. Combine Infinium I probes
  #.     - TBD:: Investigate Count from seqs_ctl_tib for InfI probes...
  #.     - Might require Address/Sequence flipping...
  #.  2. Filter for Unique Matches
  #
  #.  seqs_ctl_tib %>% dplyr::filter( Probe_ID_match %>% stringr::str_detect("_I_") ) %>% as.data.frame()
  
  all_class_tib <- NULL
  all_class_tib <- seqs_ctl_tib %>% 
    dplyr::select( -Probe_ID_bsp, -Is_EPIC_match, -Sample_Dependent_full, 
                   -Control_Type_full, -Control_Name_full, 
                   -Control_Group_Str_full, -Order_ID_full, 
                   -Probe_Type_epic ) %>% 
    dplyr::arrange(Probe_ID_full)
  
  # Look Over Categories::
  # A tibble: 29 × 1
  all_catg_tib <- all_class_tib %>% dplyr::distinct( Control_Class_match )
  
  # Look Over Categories:: Used in EPIC
  # A tibble: 21 × 1
  epi_catg_tib <- all_class_tib %>% dplyr::filter( !is.na(Probe_ID_epic) ) %>% 
    dplyr::distinct( Control_Class_match )
  
  # Look Over Categories:: NOT Used in EPIC
  # A tibble: 21 × 1
  non_catg_tib <- all_catg_tib %>% 
    dplyr::filter( !Control_Class_match %in% epi_catg_tib$Control_Class_match )
  
  # Investigate EPIC Naming Match/Mismatch::
  #.  Conclusion only a handful of BS controls have name mismatches...
  #
  all_epic_mis_tib <- NULL
  all_epic_mis_tib <- all_class_tib %>% 
    dplyr::filter( !is.na(Probe_ID_epic) ) %>% 
    dplyr::filter( Probe_ID_epic != Probe_ID_full )
  
  all_epic_mat_tib <- NULL
  all_epic_mat_tib <- all_class_tib %>% 
    dplyr::filter( !is.na(Probe_ID_epic) ) %>% 
    dplyr::filter( Probe_ID_epic == Probe_ID_full )
  
  # All General Summaries::
  #
  all_class_sum <- NULL
  all_class_sum <- all_class_tib %>% 
    dplyr::group_by( Control_Class_match ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) cat(glue::glue("{pmssg} All Control Class::{RET}"))
  if ( p4 ) all_class_sum %>% print( n=base::nrow(all_class_sum) )
  if ( p4 ) cat(glue::glue("{RET}"))
  
  all_cross_sum <- NULL
  all_cross_sum <- all_class_tib %>% 
    dplyr::group_by( Control_Class_match,Control_Type_match ) %>%
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) cat(glue::glue("{pmssg} All Control Cross::{RET}"))
  if ( p4 ) all_cross_sum %>% print( n=base::nrow(all_cross_sum) )
  if ( p4 ) cat(glue::glue("{RET}"))
  
  # B1: Need to figure out pairing
  #.  - Add last base for Infinium U/M = C/U (which is backwards)
  #
  # Example::
  #.  Address Type      Color_Channel Name          Probe_ID
  #  12798255 BISULFITE CONVERSION I  Magenta       BS Conversion I-C1 ctl_12798255
  #
  all_B1_tib <- NULL
  all_B1_tib <- all_class_tib %>% dplyr::filter( Control_Class_match == "B1" )
  # all_B1_tib %>% print(n=1000)
  # all_B1_tib %>% dplyr::filter( !is.na(Probe_ID_epic) ) %>% print(n=1000)
  
  all_B1U_tib <- NULL
  all_B1U_tib <- all_B1_tib %>% 
    dplyr::filter( Last_Base == "A" ) %>% 
    dplyr::select( -PIDX_full ) %>%
    dplyr::rename( Probe_Seq_U = Probe_Seq, 
                   Last_Base_U = Last_Base, 
                   Probe_ID_epic_U = Probe_ID_epic, 
                   Probe_ID_full_U = Probe_ID_full,
                   Probe_Type_full_U = Probe_Type_full )
  
  all_B1M_tib <- NULL
  all_B1M_tib <- all_B1_tib %>% 
    dplyr::filter( Last_Base == "G" ) %>% 
    dplyr::select( -PIDX_full ) %>%
    dplyr::rename( Probe_Seq_M = Probe_Seq, 
                   Last_Base_M = Last_Base, 
                   Probe_ID_epic_M = Probe_ID_epic, 
                   Probe_ID_full_M = Probe_ID_full,
                   Probe_Type_full_M = Probe_Type_full, M = U )
  
  # NOTE:: 
  #. - col is made up! Should be extracted from 122mer
  #  - Need to now split this back out into single columns
  #
  join_B1_tib <- NULL
  join_B1_tib <- dplyr::inner_join( 
    all_B1U_tib, all_B1M_tib,
    by = join_by(Probe_ID_match, Control_Type_match, Control_Class_match)
  ) %>% 
    dplyr::mutate( col = "R",
                   Prb_Idx = dplyr::row_number(),
                   Probe_ID = paste0("ctl_BS_Conversion_I_",Prb_Idx) ) %>%
    dplyr::select( Probe_ID,U,M,Probe_Seq_U,Probe_Seq_M,col,
                   Control_Type_match,Control_Class_match,Prb_Idx,
                   Last_Base_U,Last_Base_M,
                   Probe_ID_epic_U,Probe_ID_epic_M,
                   Probe_ID_full_U,Probe_ID_full_M,
                   Probe_Type_full_U,Probe_Type_full_M,
                   Probe_ID_match, dplyr::everything() )
  # join_B1_tib %>% head() %>% as.data.frame()
  
  #. - Need to match Genome Studio Output...
  #
  # Example::
  #.  Address Type      Color_Channel Name          Probe_ID
  #  12798255 BISULFITE CONVERSION I  Magenta       BS Conversion I-C1 ctl_12798255
  #
  gs_B1_tib <- NULL
  gs_B1_tib <- dplyr::bind_rows(
    join_B1_tib %>% 
      dplyr::mutate( 
        Address = U,
        Type = paste0("BS Conversion I"),
        Color_Channel = "Magenta",
        Probe_ID = paste0("BS Conversion I-C",Prb_Idx," ctl_",U),
        Probe_Seq = Probe_Seq_U
      ) %>% dplyr::select( Address,Type,Color_Channel,Probe_ID,Probe_Seq ),
    join_B1_tib %>% 
      dplyr::mutate( 
        Address = M,
        Type = paste0("BS Conversion I"),
        Color_Channel = "Magenta",
        Probe_ID = paste0("BS Conversion I-U",Prb_Idx," ctl_",U),
        Probe_Seq = Probe_Seq_M
      ) %>% dplyr::select( Address,Type,Color_Channel,Probe_ID,Probe_Seq )
  )
  
  #
  # B2:
  #
  
  
  if ( p1 ) cat(glue::glue("{pmssg} END: Splitting Infinium I Controls.{RET2}"))
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Pre-processing:: Manifest (MSA03)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_msa <- FALSE
run_msa <- TRUE

if ( run_msa && par$run_name == "MSAv03" ) {
  #
  # Adding MSA03
  #
  msa_man_csv <- file.path( opt$top_path, "data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_SS_BP123_A1.csv.gz" )
  msa_raw_tib <- NULL
  msa_raw_tib <- readr::read_csv( file = msa_man_csv, show_col_types = FALSE )
  msa_man_tib <- NULL
  msa_man_tib <- msa_raw_tib %>% 
    dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace("^ZZ", "ctl_"),
                   # Probe_ID = Probe_ID %>% stringr::str_replace("^ctl_neg", "ctl_Negative_"),
                   Probe_ID = dplyr::case_when(
                     Probe_Type == "NEGATIVE" ~ Probe_ID %>% stringr::str_replace("^ctl_neg", "ctl_Negative"),
                     TRUE ~ Probe_ID
                   ),
                   Name = Probe_ID %>% stringr::str_remove("_[^_]+$"),
                   M = dplyr::case_when(
                     is.na(M) ~ 0.0,
                     TRUE ~ M
                   ),
                   U = U %>% as.integer(),
                   M = M %>% as.integer(),
                   col = dplyr::case_when(
                     is.na(col) ~ "2",
                     TRUE ~ col
                   ),
                   Manifest = "MSA03", 
                   Manifest_Version = "V3", 
                   Annotation = BP,
                   Locus_Name = Probe_ID %>% stringr::str_remove("_.*$"),
                   # Probe_Type = Probe_ID %>% stringr::str_sub(1,2),
                   SNP_Probe_ID = NA ) %>% 
    dplyr::rename( Chromosome = CHR,
                   Coordinate = MAPINFO ) %>%
    # Probe_ID_msa = Probe_ID ) %>%
    dplyr::select( Probe_ID, U, M, col, Name, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                   Species, Manifest, Manifest_Version, Annotation, 
                   Chromosome, Coordinate, Probe_Type, Locus_Name, 
                   SNP_Probe_ID ) %>% clean_tib()
  
  #
  # Investigate Match All by Sequence::
  #
  # NOTES::
  #.  - msa_join_tib should have all the controls with match names
  #.  - Need to format into msa_join_GS_tib
  #.    Use: "data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_BP123_A1.csv.gz" to build GS
  #.  - Need to Sesame directly with msa_man_tib
  #.  - Need to extract Intensities for each control type directly
  msa_join_tib <- NULL
  msa_join_tib <- dplyr::inner_join(
    dplyr::rename( msa_man_tib, Probe_ID_msa = Probe_ID ),
    dplyr::rename( all_class_tib, Address = U ),
    by=c("AlleleA_ProbeSeq"="Probe_Seq"),
    suffix=c("_msa","_ctl")
  )
  # msa_join_tib %>% dplyr::arrange(Probe_ID_match) %>% head(n=20) %>% as.data.frame()  
  
  # These need to match
  gst_type_vec <- example_gs_ctl_tib %>% dplyr::distinct( Type ) %>% arrange( Type )
  msa_type_vec <- msa_join_tib %>% dplyr::distinct( Probe_Type ) %>% arrange( Probe_Type )
  mat_type_vec <- msa_man_tib %>% dplyr::distinct( Probe_Type ) %>% arrange( Probe_Type )
  
  # Order of Probes::
  # > example_gs_ctl_tib %>% dplyr::distinct( Type )
  # A tibble: 17 × 1
  # Type                   
  # <chr>                  
  # +  1 STAINING               
  # +  2 EXTENSION
  # +  3 HYBRIDIZATION
  #
  # +  4 TARGET REMOVAL         
  # +  5 BISULFITE CONVERSION I 
  # +  6 BISULFITE CONVERSION II
  # +  7 SPECIFICITY I
  # +  8 SPECIFICITY II
  # +  9 NON-POLYMORPHIC        
  # + 10 NEGATIVE               
  # + 11 RESTORATION            
  # + 12 NORM_T                 
  #   13 NORM_G                 
  #   14 NORM_C                 
  #   15 NORM_A
  #
  # + 16 Stringency
  # + 17 Non-Specific Binding
  
  # Need to Include Infinium Controls:
  #
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               1 STAINING::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Probe_Type => STAINING [COPY]
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "STAINING" )
  # A tibble: 6 × 5
  # Address Type     Color_Channel Name          Probe_ID    
  # <dbl> <chr>    <chr>         <chr>         <chr>       
  # 1 21630339 STAINING -99           DNP(20K)      ctl_21630339
  # 2 24669308 STAINING -99           Biotin(5K)    ctl_24669308
  # 3 27630314 STAINING Red           DNP (High)    ctl_27630314
  # 4 34648333 STAINING Blue          Biotin (Bkg)  ctl_34648333
  # 5 41666334 STAINING Green         Biotin (High) ctl_41666334
  # 6 43603326 STAINING Purple        DNP (Bkg)     ctl_43603326
  # msa_man_tib %>% dplyr::filter( U==21630339 | M == 21630339 ) %>% as.data.frame()
  # msa_man_tib %>% dplyr::filter( Probe_Type %>% stringr::str_starts("STAINING") )
  
  # OLD Direct Method:: Addresses are already in the Pool for Infinium Controls!
  # ctl_tib1 <- NULL
  # ctl_tib1 <- example_gs_ctl_tib %>% dplyr::filter( Type == "STAINING" ) %>% 
  #   dplyr::mutate( Address = Address %>% as.integer(), 
  #                  AlleleA_ProbeSeq = NA_character_, 
  #                  AlleleB_ProbeSeq = NA_character_ )
  
  gs_ctl_tib1 <- NULL
  gs_ctl_tib1 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type %>% stringr::str_starts("STAINING") ) %>% 
    dplyr::mutate( Type = "STAINING" ) %>% 
    dplyr::select( Probe_ID,U,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Color_Channel = dplyr::case_when(
        Probe_ID == "ctl_DNP_20K"     ~ "-99",
        Probe_ID == "ctl_Biotin_5K"   ~ "-99",
        Probe_ID == "ctl_DNP_High"    ~ "Red",
        Probe_ID == "ctl_Biotin_Bkg"  ~ "Blue",
        Probe_ID == "ctl_Biotin_High" ~ "Green",
        Probe_ID == "ctl_DNP_Bkg"     ~ "Purple",
        TRUE ~ NA_character_
      ),
      Name = dplyr::case_when(
        Probe_ID == "ctl_DNP_20K"     ~ "DNP(20K)",
        Probe_ID == "ctl_Biotin_5K"   ~ "Biotin(5K)",
        Probe_ID == "ctl_DNP_High"    ~ "DNP (High)",
        Probe_ID == "ctl_Biotin_Bkg"  ~ "Biotin (Bkg)",
        Probe_ID == "ctl_Biotin_High" ~ "Biotin (High)",
        Probe_ID == "ctl_DNP_Bkg"     ~ "DNP (Bkg)",
        TRUE ~ NA_character_
      ),
      Address_Mate = as.integer(0)
    ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               2 EXTENSION::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Probe_Type => EXTENSION [UPDATE]
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "EXTENSION" )
  # A tibble: 4 × 5
  #  Address Type      Color_Channel Name          Probe_ID    
  # <dbl> <chr>     <chr>         <chr>         <chr>       
  # 11603365 EXTENSION Blue          Extension (G) ctl_11603365
  # 14607337 EXTENSION Purple        Extension (T) ctl_14607337
  # 17616306 EXTENSION Red           Extension (A) ctl_17616306
  # 12613307 EXTENSION Green         Extension (C) ctl_12613307
  # msa_man_tib %>% dplyr::filter( Probe_Type %>% stringr::str_starts("EXTENSION") )
  
  gs_ctl_tib2 <- NULL
  gs_ctl_tib2 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type %>% stringr::str_starts("EXTENSION") ) %>% 
    dplyr::mutate( Type = "EXTENSION" ) %>% 
    dplyr::select( Probe_ID,U,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Color_Channel = dplyr::case_when(
        Name == "EXTENSION (G)" ~ "Blue",
        Name == "EXTENSION (T)" ~ "Purple",
        Name == "EXTENSION (A)" ~ "Red",
        Name == "EXTENSION (C)" ~ "Green",
        TRUE ~ NA_character_
      ),
      Name = dplyr::case_when(
        Name == "EXTENSION (G)" ~ "Extension (G)",
        Name == "EXTENSION (T)" ~ "Extension (T)",
        Name == "EXTENSION (A)" ~ "Extension (A)",
        Name == "EXTENSION (C)" ~ "Extension (C)",
        TRUE ~ NA_character_
      ),
      Address_Mate = as.integer(0)
    ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            3 HYBRIDIZATION::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Probe_Type => HYBRIDIZATION [UPDATE]
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "HYBRIDIZATION" )
  # A tibble: 3 × 5
  # Address Type          Color_Channel Name         Probe_ID    
  # <dbl> <chr>         <chr>         <chr>        <chr>       
  # 1 23617335 HYBRIDIZATION Black         Hyb (Low)    ctl_23617335
  # 2 20636378 HYBRIDIZATION Blue          Hyb (Medium) ctl_20636378
  # 3 19612319 HYBRIDIZATION Green         Hyb (High)   ctl_19612319
  # msa_man_tib %>% dplyr::filter( Probe_Type %>% stringr::str_starts("HYBRIDIZATION") )
  
  gs_ctl_tib3 <- NULL
  gs_ctl_tib3 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type %>% stringr::str_starts("HYBRIDIZATION") ) %>% 
    dplyr::mutate( Type = "HYBRIDIZATION" ) %>% 
    dplyr::select( Probe_ID,U,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Color_Channel = dplyr::case_when(
        Name == "HYBRIDIZATION-LOW"    ~ "Black",
        Name == "HYBRIDIZATION-MEDIUM" ~ "Blue",
        Name == "HYBRIDIZATION - HIGH" ~ "Green",
        Name == "HYBRIDIZATION-HIGH"   ~ "Green",
        TRUE ~ NA_character_
      ),
      Name = dplyr::case_when(
        Name == "HYBRIDIZATION-LOW"    ~ "Hyb (Low)",
        Name == "HYBRIDIZATION-MEDIUM" ~ "Hyb (Medium)",
        Name == "HYBRIDIZATION - HIGH" ~ "Hyb (High)",
        Name == "HYBRIDIZATION-HIGH"   ~ "Hyb (High)",
        TRUE ~ NA_character_
      ),
      Address_Mate = as.integer(0)
    ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             16 Stringency::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # STRINGENCY I  => STRINGENCY [COPY]
  # STRINGENCY II => STRINGENCY [COPY]
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "Stringency" )
  # A tibble: 2 × 5
  # Address Type       Color_Channel Name        Probe_ID    
  # <dbl> <chr>      <chr>         <chr>       <chr>       
  # 1 32629312 Stringency Red           String (PM) ctl_32629312
  # 2 33668307 Stringency Purple        String (MM) ctl_33668307
  # msa_man_tib %>% dplyr::filter( Probe_Type %>% stringr::str_starts("STRINGENCY") )
  
  gs_ctl_tib16 <- NULL
  gs_ctl_tib16 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type %>% stringr::str_starts("STRINGENCY") ) %>% 
    dplyr::mutate( Type = "Stringency" ) %>% 
    dplyr::select( Probe_ID,U,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Color_Channel = dplyr::case_when(
        Probe_ID == "ctl_Stringency_1" ~ "Red",
        Probe_ID == "ctl_Stringency_2" ~ "Purple",
        TRUE ~ NA_character_
      ),
      Name = dplyr::case_when(
        Probe_ID == "ctl_Stringency_1" ~ "String (PM)",
        Probe_ID == "ctl_Stringency_2" ~ "String (MM)",
        TRUE ~ NA_character_
      ),
      Address_Mate = as.integer(0)
    ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                        17 Non-Specific Binding::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Probe_ID: Non-Specific_Binding => "Non-Specific Binding" [COPY]
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "Non-Specific Binding" )
  # A tibble: 4 × 5
  # Address Type                 Color_Channel Name        Probe_ID    
  # <dbl> <chr>                <chr>         <chr>       <chr>       
  # 1 26619332 Non-Specific Binding Red           NSB (Bgnd)1 ctl_26619332
  # 2 27624356 Non-Specific Binding Purple        NSB (Bgnd)2 ctl_27624356
  # 3 25617343 Non-Specific Binding Blue          NSB (Bgnd)3 ctl_25617343
  # 4 24616350 Non-Specific Binding Green         NSB (Bgnd)4 ctl_24616350
  # msa_man_tib %>% dplyr::filter( U==26619332 | M == 26619332 ) %>% as.data.frame()
  # msa_man_tib %>% dplyr::filter( Probe_Type %>% stringr::str_starts("Non-Specific_Binding") )
  
  gs_ctl_tib17 <- NULL
  gs_ctl_tib17 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type %>% stringr::str_starts("Non-Specific_Binding") ) %>% 
    dplyr::mutate( Type = "Non-Specific Binding" ) %>% 
    dplyr::select( Probe_ID,U,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Color_Channel = dplyr::case_when(
        Probe_ID == "ctl_Non-Specific_Binding_1" ~ "Red",
        Probe_ID == "ctl_Non-Specific_Binding_2" ~ "Purple",
        Probe_ID == "ctl_Non-Specific_Binding_3" ~ "Blue",
        Probe_ID == "ctl_Non-Specific_Binding_4" ~ "Green",
        TRUE ~ NA_character_
      ),
      Name = dplyr::case_when(
        Probe_ID == "ctl_Non-Specific_Binding_1" ~ "NSB (Bgnd)1",
        Probe_ID == "ctl_Non-Specific_Binding_2" ~ "NSB (Bgnd)2",
        Probe_ID == "ctl_Non-Specific_Binding_3" ~ "NSB (Bgnd)3",
        Probe_ID == "ctl_Non-Specific_Binding_4" ~ "NSB (Bgnd)4",
        TRUE ~ NA_character_
      )
    ) %>% dplyr::select( Address,Type,Color_Channel,Name,
                         Probe_ID, AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                     Infinium Methylation Controls::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # Previous Work::
  #
  # all_B1U_tib
  # all_B1M_tib
  # join_B1_tib
  # gs_B1_tib
  
  # +  4 TARGET REMOVAL
  # +  5 BISULFITE CONVERSION I 
  #      + Normal
  #      + No Bias
  # +  6 BISULFITE CONVERSION II
  # +  7 SPECIFICITY I
  # +  8 SPECIFICITY II
  # +  9 NON-POLYMORPHIC
  # - 10 NEGATIVE
  # - 11 RESTORATION
  # - 12 NORM_T
  # - 13 NORM_G
  # - 14 NORM_C
  # - 15 NORM_A
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            4 TARGET REMOVAL::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #  
  # > example_gs_ctl_tib %>% dplyr::filter( Type %>% stringr::str_starts("TARGET") )
  # A tibble: 2 × 5
  # Address Type           Color_Channel Name             Probe_ID    
  # <dbl> <chr>          <chr>         <chr>            <chr>       
  # 1 31623323 TARGET REMOVAL Green         Target Removal 1 ctl_31623323
  # 2 96694467 TARGET REMOVAL Magenta       Target Removal 2 ctl_96694467
  #
  # NOTE: 96694467 was completely resynthesised...
  gs_ctl_tib4 <- NULL
  gs_ctl_tib4 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type %>% stringr::str_starts("Target") ) %>% 
    dplyr::mutate( Type = "TARGET REMOVAL" ) %>% 
    dplyr::select( Probe_ID,U,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Name = dplyr::case_when(
        Address == 96694467 ~ "Target Removal 2",
        TRUE ~ Name
      ),
      Color_Channel = dplyr::case_when(
        Name == "Target Removal 1" ~ "Green",
        Name == "Target Removal 2" ~ "Magenta",
        TRUE ~ NA_character_
      )
    ) %>% dplyr::select( Address,Type,Color_Channel,Name,
                         Probe_ID, AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    5 BISULFITE CONVERSION I:: Normal
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "BISULFITE CONVERSION I" )
  # A tibble: 10 × 5
  #     Address Type                   Color_Channel Name               Probe_ID    
  #       <dbl> <chr>                  <chr>         <chr>              <chr>       
  #  1 12798255 BISULFITE CONVERSION I Magenta       BS Conversion I-C1 ctl_12798255
  #  2 52757471 BISULFITE CONVERSION I Magenta       BS Conversion I-C2 ctl_52757471
  #  3 53659162 BISULFITE CONVERSION I Magenta       BS Conversion I-C3 ctl_53659162
  #  4 71719211 BISULFITE CONVERSION I Magenta       BS Conversion I-C4 ctl_71719211
  #  5 80665422 BISULFITE CONVERSION I Magenta       BS Conversion I-C5 ctl_80665422
  #  6 95745587 BISULFITE CONVERSION I Magenta       BS Conversion I-U1 ctl_95745587
  #  7 10793931 BISULFITE CONVERSION I Magenta       BS Conversion I-U2 ctl_10793931
  #  8  5640503 BISULFITE CONVERSION I Magenta       BS Conversion I-U3 ctl_5640503 
  #  9  3784326 BISULFITE CONVERSION I Magenta       BS Conversion I-U4 ctl_3784326 
  # 10 28717193 BISULFITE CONVERSION I Magenta       BS Conversion I-U5 ctl_28717193
  
  name_chr <- "C"
  gs_ctl_tib5C <- NULL
  gs_ctl_tib5C <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "BISULFITE CONVERSION I" ) %>% 
    dplyr::filter( !Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_2_0.1_") ) %>% 
    dplyr::mutate( Type = "BISULFITE CONVERSION I" ) %>% 
    dplyr::select( Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = M, 
                   Address_Mate = U,
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Prb_Idx = Probe_ID %>% 
        stringr::str_remove("^ctl_BS_Conversion_I_") %>%
        stringr::str_remove("^ctl_BS_conversion_NoBiasHIII_ASPE_2_0.1_") %>%
        stringr::str_remove("_.*$") %>%
        as.integer(),
      Rep_Idx = Probe_ID %>% 
        stringr::str_remove("^.*_") %>%
        as.integer(),
      Name = paste0("BS Conversion I-",name_chr,Prb_Idx),
      Probe_ID = paste0("ctl_BS_Conversion_I_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Magenta"
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  name_chr <- "U"
  gs_ctl_tib5U <- NULL
  gs_ctl_tib5U <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "BISULFITE CONVERSION I" ) %>% 
    dplyr::filter( !Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_2_0.1_") ) %>% 
    dplyr::mutate( Type = "BISULFITE CONVERSION I" ) %>% 
    dplyr::select( Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U,
                   Address_Mate = M,
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Prb_Idx = Probe_ID %>% 
        stringr::str_remove("^ctl_BS_Conversion_I_") %>%
        stringr::str_remove("^ctl_BS_conversion_NoBiasHIII_ASPE_2_0.1_") %>%
        stringr::str_remove("_.*$") %>%
        as.integer(),
      Rep_Idx = Probe_ID %>% 
        stringr::str_remove("^.*_") %>%
        as.integer(),
      Name = paste0("BS Conversion I-",name_chr,Prb_Idx),
      Probe_ID = paste0("ctl_BS_Conversion_I_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Magenta"
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                  5 BISULFITE CONVERSION I:: NoBiasHIII
  #
  #                  ctl_BS_conversion_NoBiasHIII_ASPE_2_0.1
  #                 ctl_BS_conversion_NoBiasHIII_ASPE_12_0.1
  #                  ctl_BS_conversion_NoBiasHIII_ASPE_9_0.1
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  name_chr <- "C"
  gs_ctl_tib5CH <- NULL
  gs_ctl_tib5CH <- msa_man_tib %>% 
    dplyr::filter( Probe_Type %>% stringr::str_starts("BISULFITE CONVERSION I") ) %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_2_0.1_") |
                     Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_12_0.1_") |
                     Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_9_0.1_") ) %>% 
    dplyr::arrange(AlleleA_ProbeSeq) %>% 
    dplyr::group_by( AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::mutate( Type = "NOBIAS BISULFITE CONVERSION I",
                   Prb_Idx = dplyr::cur_group_id(), 
                   Rep_Idx = dplyr::row_number() ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select( Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                   Prb_Idx,Rep_Idx ) %>%
    dplyr::rename( Address = M, 
                   Address_Mate = U,
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Name = paste0("NoBias BS Conversion I-",name_chr,Prb_Idx),
      Probe_ID = paste0(stringr::str_remove(Probe_ID,"_[0-9]-[0-9]$"),"_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Pink"
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  name_chr <- "U"
  gs_ctl_tib5UH <- NULL
  gs_ctl_tib5UH <- msa_man_tib %>% 
    dplyr::filter( Probe_Type %>% stringr::str_starts("BISULFITE CONVERSION I") ) %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_2_0.1_") |
                     Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_12_0.1_") |
                     Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_9_0.1_") ) %>% 
    dplyr::arrange(AlleleA_ProbeSeq) %>% 
    dplyr::group_by( AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::mutate( Type = "NOBIAS BISULFITE CONVERSION I",
                   Prb_Idx = dplyr::cur_group_id(), 
                   Rep_Idx = dplyr::row_number() ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select( Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq,
                   Prb_Idx,Rep_Idx ) %>%
    dplyr::rename( Address = U, 
                   Address_Mate = M,
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Name = paste0("NoBias BS Conversion I-",name_chr,Prb_Idx),
      Probe_ID = paste0(stringr::str_remove(Probe_ID,"_[0-9]-[0-9]$"),"_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Pink"
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    6 BISULFITE CONVERSION II:: Normal
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "BISULFITE CONVERSION II" )
  # A tibble: 4 × 5
  # Address Type                    Color_Channel Name               Probe_ID    
  # <dbl> <chr>                   <chr>         <chr>              <chr>       
  # 1 79804842 BISULFITE CONVERSION II Magenta       BS Conversion II-1 ctl_79804842
  # 2 83810932 BISULFITE CONVERSION II Magenta       BS Conversion II-2 ctl_83810932
  # 3 35736226 BISULFITE CONVERSION II Magenta       BS Conversion II-3 ctl_35736226
  # 4 67713990 BISULFITE CONVERSION II Magenta       BS Conversion II-4 ctl_67713990
  
  gs_ctl_tib6 <- NULL
  gs_ctl_tib6 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "BISULFITE CONVERSION II" ) %>% 
    dplyr::filter( !Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_2_0.1_") ) %>% 
    dplyr::filter( !Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_12_0.1_") ) %>% 
    dplyr::filter( !Probe_ID %>% stringr::str_starts("ctl_BS_conversion_NoBiasHIII_ASPE_9_0.1_") ) %>% 
    dplyr::mutate( Type = "BISULFITE CONVERSION II" ) %>% 
    dplyr::select( Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Address_Mate = M,
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Prb_Idx = Probe_ID %>% 
        stringr::str_remove("^ctl_BS_Conversion_II_") %>%
        stringr::str_remove("_.*$") %>%
        as.integer(),
      Rep_Idx = Probe_ID %>% 
        stringr::str_remove("^.*_") %>%
        as.integer(),
      Name = paste0("BS Conversion II-",Prb_Idx),
      # Probe_ID = paste0("ctl_BS_Conversion_II_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Magenta"
    ) %>% # dplyr::filter( is.na(Prb_Idx) )
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             7 SPECIFICITY I::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "SPECIFICITY I" )
  # A tibble: 12 × 5
  #     Address Type          Color_Channel Name               Probe_ID    
  # <dbl> <chr>         <chr>         <chr>              <chr>       
  #  1 18668925 SPECIFICITY I Magenta       GT Mismatch 1 (PM) ctl_18668925
  #  2 63647415 SPECIFICITY I Magenta       GT Mismatch 2 (PM) ctl_63647415
  #  3 38600876 SPECIFICITY I Magenta       GT Mismatch 3 (PM) ctl_38600876
  #  4  7688128 SPECIFICITY I Magenta       GT Mismatch 1 (MM) ctl_7688128 
  #  5 10741434 SPECIFICITY I Magenta       GT Mismatch 2 (MM) ctl_10741434
  #  6 91670295 SPECIFICITY I Magenta       GT Mismatch 3 (MM) ctl_91670295
  #  7 33652355 SPECIFICITY I Magenta       GT Mismatch 4 (PM) ctl_33652355
  #  8 94661891 SPECIFICITY I Magenta       GT Mismatch 5 (PM) ctl_94661891
  #  9 21773943 SPECIFICITY I Magenta       GT Mismatch 6 (PM) ctl_21773943
  # 10 88646924 SPECIFICITY I Magenta       GT Mismatch 4 (MM) ctl_88646924
  # 11 31773497 SPECIFICITY I Magenta       GT Mismatch 5 (MM) ctl_31773497
  # 12 85637245 SPECIFICITY I Magenta       GT Mismatch 6 (MM) ctl_85637245
  #
  # Investigation::
  #.  TBD:: Need to understand the GT_Mismatch_1_MM/PM naming convention
  #
  # GT Mismatch 1 (PM)
  # msa_man_tib %>% dplyr::filter( U == 18668925 ) %>% as.data.frame()
  # msa_man_tib %>% dplyr::filter( AlleleA_ProbeSeq == "ATAATACAATAAACCAATCTTACTACTTAACTACTCTACTCCAAACTAAA" ) %>% as.data.frame()
  #
  # GT Mismatch 1 (MM)
  # msa_man_tib %>% dplyr::filter( U == 7688128 ) %>% as.data.frame()
  # msa_man_tib %>% dplyr::filter( AlleleB_ProbeSeq == "ATAATACAATAAACCAATCTTACTACTTAACTACTCTACTCCAAACTAAG" ) %>% as.data.frame()
  
  name_chr <- "PM"
  gs_ctl_tib7P <- NULL
  gs_ctl_tib7P <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "Specificity I" ) %>%
    dplyr::mutate( Type = "SPECIFICITY I" ) %>% 
    dplyr::select( Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Address_Mate = M,
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Prb_Idx = Probe_ID %>% 
        stringr::str_remove("^ctl_Non_Specific_I_") %>%
        stringr::str_remove("_.*$") %>%
        as.integer(),
      Rep_Idx = Probe_ID %>% 
        stringr::str_remove("^.*_") %>%
        as.integer(),
      Name = paste0("GT Mismatch ",Prb_Idx," (",name_chr,")"),
      Probe_ID = paste0("ctl_Non_Specific_I_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Magenta"
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  name_chr <- "MM"
  gs_ctl_tib7M <- NULL
  gs_ctl_tib7M <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "Specificity I" ) %>%
    dplyr::mutate( Type = "SPECIFICITY I" ) %>% 
    dplyr::select( Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = M, 
                   Address_Mate = U,
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Prb_Idx = Probe_ID %>% 
        stringr::str_remove("^ctl_Non_Specific_I_") %>%
        stringr::str_remove("_.*$") %>%
        as.integer(),
      Rep_Idx = Probe_ID %>% 
        stringr::str_remove("^.*_") %>%
        as.integer(),
      Name = paste0("GT Mismatch ",Prb_Idx," (",name_chr,")"),
      Probe_ID = paste0("ctl_Non_Specific_I_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Magenta"
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                             8 SPECIFICITY II::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "SPECIFICITY II" )
  # A tibble: 3 × 5
  # Address Type           Color_Channel Name          Probe_ID    
  # <dbl> <chr>          <chr>         <chr>         <chr>       
  # 1 60799910 SPECIFICITY II Magenta       Specificity 1 ctl_60799910
  # 2 58794576 SPECIFICITY II Magenta       Specificity 2 ctl_58794576
  # 3 38789133 SPECIFICITY II Magenta       Specificity 3 ctl_38789133
  # msa_man_tib %>% dplyr::filter( Probe_Type == "Specificity II" ) %>% as.data.frame()
  
  gs_ctl_tib8 <- NULL
  gs_ctl_tib8 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "Specificity II" ) %>%
    dplyr::mutate( Type = "SPECIFICITY II" ) %>% 
    dplyr::select( Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Address_Mate = M,
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Prb_Idx = Probe_ID %>% 
        stringr::str_remove("^ctl_Non_Specific_II_") %>%
        stringr::str_remove("_.*$") %>%
        as.integer(),
      Rep_Idx = Probe_ID %>% 
        stringr::str_remove("^.*_") %>%
        as.integer(),
      Name = paste0("Specificity ",Prb_Idx),
      Probe_ID = paste0("ctl_Non_Specific_II_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Magenta"
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           9 NON-POLYMORPHIC::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "NON-POLYMORPHIC" )
  # A tibble: 9 × 5
  #    Address Type            Color_Channel Name     Probe_ID    
  # <dbl> <chr>           <chr>         <chr>    <chr>       
  # 1 34633358 NON-POLYMORPHIC Red           NP (A)   ctl_34633358
  # 2 16648324 NON-POLYMORPHIC Purple        NP (T)   ctl_16648324
  # 3 43641328 NON-POLYMORPHIC Green         NP (C)   ctl_43641328
  # 4 13642359 NON-POLYMORPHIC Blue          NP (G)   ctl_13642359
  # 5 58619842 NON-POLYMORPHIC Magenta       NP (G) 1 ctl_58619842
  # 6  6722343 NON-POLYMORPHIC Magenta       NP (G) 2 ctl_6722343 
  # 7 31622874 NON-POLYMORPHIC Magenta       NP (G) 3 ctl_31622874
  # 8 98683256 NON-POLYMORPHIC Magenta       NP (G) 4 ctl_98683256
  # 9 21644544 NON-POLYMORPHIC Magenta       NP (G) 5 ctl_21644544
  #
  # msa_man_tib %>% dplyr::filter( Probe_Type == "NON-POLYMORPHIC" ) %>% as.data.frame()
  # all_class_tib %>% dplyr::filter( Control_Class_match == "NG" )
  # all_class_tib %>% dplyr::filter( Probe_Seq == "TTTTTAGAGGATTATATTGAAGGTAAGGGGTTATTGAGGATGTTGAGGCC" ) %>% as.data.frame()
  # all_class_tib %>% dplyr::filter( Probe_Seq == "TTTTTGGAGAATAAGAAGGATTTAAAGGAAGGTAGAGAGTAGTTTGTTTT" ) %>% as.data.frame()
  # all_class_tib %>% dplyr::filter( Probe_Seq == "TAGGTGTATTTGGTTAGGAAGAGGGTATAGGTAGGGGTATGAGTTTAGCC" ) %>% as.data.frame()
  # all_class_tib %>% dplyr::filter( Probe_Seq == "TATTAGTTATGGGAGTTATGAGATTGTTTAAGGAGTATATGTAGTGTGTC" ) %>% as.data.frame()
  # all_class_tib %>% dplyr::filter( Probe_Seq == "TTTAGGGGGATTATTGAATGGAATTGTAGTTTGAGAGTGTTGTTATAAAA" ) %>% as.data.frame()
  
  gs_ctl_tib9 <- NULL
  gs_ctl_tib9 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "NON-POLYMORPHIC" ) %>%
    dplyr::mutate( Type = "NON-POLYMORPHIC" ) %>% 
    dplyr::select( Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Address_Mate = M,
                   Name = Probe_Type ) %>% 
    dplyr::mutate(
      Color_Channel = dplyr::case_when(
        Address == 34633358 ~ "Red",
        Address == 16648324 ~ "Purple",
        Address == 43641328 ~ "Green",
        Address == 13642359 ~ "Blue",
        TRUE ~ "Magenta"
      ),
      Name = dplyr::case_when(
        Address == 34633358 ~ "NP (A)",
        Address == 16648324 ~ "NP (T)",
        Address == 43641328 ~ "NP (C)",
        Address == 13642359 ~ "NP (G)",
        
        Address == 58619842 ~ "NP (G) 1",
        Address ==  6722343 ~ "NP (G) 2",
        Address == 31622874 ~ "NP (G) 3",
        Address == 98683256 ~ "NP (G) 4",
        Address == 21644544 ~ "NP (G) 5",
        TRUE ~ NA_character_
      ),
      Probe_ID = dplyr::case_when(
        Address == 34633358 ~ paste0(Probe_ID),
        Address == 16648324 ~ paste0(Probe_ID),
        Address == 43641328 ~ paste0(Probe_ID),
        Address == 13642359 ~ paste0(Probe_ID),
        
        Address == 58619842 ~ paste0(Probe_ID,"_1-1"),
        Address ==  6722343 ~ paste0(Probe_ID,"_2-1"),
        Address == 31622874 ~ paste0(Probe_ID,"_3-1"),
        Address == 98683256 ~ paste0(Probe_ID,"_4-1"),
        Address == 21644544 ~ paste0(Probe_ID,"_5-1"),
        TRUE ~ NA_character_
      ),
      Prb_Idx = dplyr::case_when(
        Address == 34633358 ~ 0,
        Address == 16648324 ~ 0,
        Address == 43641328 ~ 0,
        Address == 13642359 ~ 0,
        
        Address == 58619842 ~ 1,
        Address ==  6722343 ~ 2,
        Address == 31622874 ~ 3,
        Address == 98683256 ~ 4,
        Address == 21644544 ~ 5,
        TRUE ~ NA_integer_
      ),
      Rep_Idx = dplyr::case_when(
        Address == 34633358 ~ 1,
        Address == 16648324 ~ 2,
        Address == 43641328 ~ 3,
        Address == 13642359 ~ 4,
        
        Address == 58619842 ~ 1,
        Address ==  6722343 ~ 1,
        Address == 31622874 ~ 1,
        Address == 98683256 ~ 1,
        Address == 21644544 ~ 1,
        TRUE ~ NA_integer_
      )
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                              10 NEGATIVE::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "NEGATIVE" )
  # A tibble: 411 × 5
  # Address Type     Color_Channel Name          Probe_ID    
  # <dbl> <chr>    <chr>         <chr>         <chr>       
  # 1 99791588 NEGATIVE Magenta       Negative 1286 ctl_99791588
  # 2 33760143 NEGATIVE Magenta       Negative 265  ctl_33760143
  # 3 19612924 NEGATIVE Magenta       Negative 1284 ctl_19612924
  # 4 47791180 NEGATIVE Magenta       Negative 583  ctl_47791180
  # 5 57653116 NEGATIVE Magenta       Negative 318  ctl_57653116
  # msa_man_tib %>% dplyr::filter( Probe_Type == "NEGATIVE" ) %>% head() %>% as.data.frame()
  #
  # NOTE:: 
  #.  Rep_Idx: last number on Probe_ID
  #.  Prb_Idx: Strip Rep Group Sequence,Probe_ID = group_id_cnt...
  #.  Check Unique Final Probe_ID
  #
  
  gs_ctl_tib10 <- NULL
  gs_ctl_tib10 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "NEGATIVE" ) %>%
    dplyr::filter( Probe_ID %>% stringr::str_starts("ctl_Negative") ) %>%
    dplyr::mutate( Type = "NEGATIVE",
                   Target_ID = Probe_ID %>% stringr::str_remove("_[0-9]+$"),
                   Rep_Idx = Probe_ID %>% stringr::str_remove("^.*_") %>% 
                     as.integer()
    ) %>%
    dplyr::select( Target_ID,Rep_Idx,Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Address_Mate = M,
                   Name = Probe_Type ) %>% 
    dplyr::arrange( Probe_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::group_by( Target_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::mutate( Prb_Idx = dplyr::cur_group_id(), 
                   Rep_Idx2 = dplyr::row_number() ) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(
      Name = paste0("Negative ",Prb_Idx),
      Probe_ID = paste0(Target_ID,"_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Magenta"
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                            11 RESTORATION::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "RESTORATION" )
  # A tibble: 1 × 5
  #    Address Type        Color_Channel Name    Probe_ID    
  #      <dbl> <chr>       <chr>         <chr>   <chr>       
  # 1 28637363 RESTORATION Green         Restore ctl_28637363
  # msa_man_tib %>% dplyr::filter( Probe_Type == "RESTORATION" ) %>% head() %>% as.data.frame()
  
  gs_ctl_tib11 <- NULL
  gs_ctl_tib11 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "RESTORATION" ) %>%
    dplyr::filter( U != 28637363 ) %>%
    dplyr::mutate( Type = "RESTORATION",
                   Target_ID = Probe_ID %>% stringr::str_remove("_[0-9]+$"),
                   Rep_Idx = Probe_ID %>% stringr::str_remove("^.*_") %>% 
                     as.integer()
    ) %>%
    dplyr::select( Target_ID,Rep_Idx,Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Address_Mate = M,
                   Name = Probe_Type ) %>% 
    dplyr::arrange( Probe_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::group_by( Target_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::mutate( Prb_Idx = dplyr::cur_group_id(), 
                   Rep_Idx2 = dplyr::row_number() ) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(
      Name = paste0("Restore ",Prb_Idx),
      Probe_ID = paste0(Target_ID,"_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Lime"
    ) %>%
    dplyr::arrange( Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # NOTE: Remove 28637363 and add back in manually to the top of the tib...
  #  Address Type        Color_Channel Name    Probe_ID    
  # 28637363 RESTORATION Green         Restore ctl_28637363
  
  gs_ctl_tib11 <- dplyr::bind_rows( 
    tibble::tibble(
      Address = as.integer(28637363),
      Type = "RESTORATION",
      Color_Channel = "Green",
      Name = "Restore",
      Probe_ID = "ctl_Restore",
      Address_Mate = 0,
      AlleleA_ProbeSeq = NA_character_,
      AlleleB_ProbeSeq = NA_character_
    ), 
    gs_ctl_tib11 )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                               12 NORM_T::
  #                               13 NORM_G::
  #                               14 NORM_C::
  #                               15 NORM_A::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # 12 NORM_T::
  #
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "NORM_T" )
  # A tibble: 58 × 5
  #    Address Type   Color_Channel Name     Probe_ID    
  #      <dbl> <chr>  <chr>         <chr>    <chr>       
  # 1 46644213 NORM_T Magenta       Norm_T46 ctl_46644213
  # 2 65661414 NORM_T Magenta       Norm_T2  ctl_65661414
  # 3 64627961 NORM_T Magenta       Norm_T8  ctl_64627961
  # msa_man_tib %>% dplyr::filter( Probe_Type == "NORM_T" )
  #
  # msa_norm_ctl_sum <- msa_man_tib %>% 
  #   dplyr::filter( Probe_Type %>% stringr::str_starts("NORM_") ) %>% 
  #   dplyr::group_by( Probe_Type ) %>% 
  #   dplyr::summarise( Count=n(), .groups = "drop" )
  #
  # 13 NORM_G::
  #
  # > example_gs_ctl_tib %>% dplyr::filter( Type == "NORM_G" )
  # A tibble: 27 × 5
  #    Address Type   Color_Channel Name     Probe_ID    
  #      <dbl> <chr>  <chr>         <chr>    <chr>       
  # 1 55654548 NORM_G Magenta       Norm_G22 ctl_55654548
  # 2  5676968 NORM_G Magenta       Norm_G24 ctl_5676968 
  # 3 66644179 NORM_G Magenta       Norm_G16 ctl_66644179
  
  gs_ctl_tib12 <- NULL
  gs_ctl_tib12 <- msa_man_tib %>% 
    dplyr::filter( Probe_Type %>% stringr::str_starts("NORM_") ) %>%
    dplyr::mutate( Type = Probe_Type %>% stringr::str_remove("ed$"),
                   Target_ID = Probe_ID %>% stringr::str_remove("_[0-9]+$"),
                   Rep_Idx = Probe_ID %>% stringr::str_remove("^.*_") %>% 
                     as.integer(),
                   Prb_Grp = Type %>% stringr::str_remove("^.*_")
    ) %>%
    dplyr::select( Target_ID,Rep_Idx,Prb_Grp,Probe_ID,U,M,Probe_Type,Type, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::rename( Address = U, 
                   Address_Mate = M,
                   Name = Probe_Type ) %>% 
    dplyr::arrange( Prb_Grp,Probe_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::group_by( Target_ID,AlleleA_ProbeSeq,AlleleB_ProbeSeq ) %>%
    dplyr::mutate( Prb_Idx = dplyr::cur_group_id(), 
                   Rep_Idx2 = dplyr::row_number() ) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(
      Name = paste0("Norm_",Prb_Grp,Prb_Idx),
      Probe_ID = paste0(Target_ID,"_",Prb_Idx,"-",Rep_Idx),
      Color_Channel = "Magenta"
    ) %>%
    dplyr::arrange( Prb_Grp,Prb_Idx,Rep_Idx ) %>% 
    dplyr::select( Address,Type,Color_Channel,Name,
                   Probe_ID,Address_Mate, 
                   AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                      Remaining Controls Probes::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #             Join Data & Compare Against Master Match File::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  # A tibble: 8,370 × 8
  use_gs_ctl_tib <- NULL
  use_gs_ctl_tib <- dplyr::bind_rows(
    gs_ctl_tib1,
    gs_ctl_tib2,
    gs_ctl_tib3,
    gs_ctl_tib4,
    gs_ctl_tib5C,
    gs_ctl_tib5CH,
    gs_ctl_tib5U,
    gs_ctl_tib5UH,
    gs_ctl_tib6,
    gs_ctl_tib7M,
    gs_ctl_tib7P,
    gs_ctl_tib8,
    gs_ctl_tib9,
    gs_ctl_tib10,
    gs_ctl_tib11,
    gs_ctl_tib12
  ) %>% clean_tib()
  
  msa_ctl_tab <- NULL
  msa_ctl_tab <- msa_man_tib %>% 
    dplyr::filter( Probe_ID %>% stringr::str_starts("ct") ) %>% 
    dplyr::select( Probe_ID, U,M ) %>% 
    tidyr::pivot_longer( cols = c(U,M), 
                         names_to = "UM", 
                         values_to = "Address" ) %>% 
    clean_tib() %>% dplyr::filter( Address != 0 )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Pre-processing:: Sample Sheet
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  if ( TRUE || par$run_name == "MSAv03" ) {
    ssh_rep_tib <- NULL
    ssh_rep_csv <- file.path( opt$top_path, "Projects.new/MSA/MethylationQC/MSAv03/post_processing/samplesheet.csv.gz" )
    ssh_rep_tib <- readr::read_csv( file = ssh_rep_csv, show_col_types = FALSE ) %>%
      dplyr::mutate( Sample_Concentration = Sample_Group %>% 
                       stringr::str_remove(" .*$") %>% as.integer(),
                     Base_Sample_Name = Sample_Well %>% stringr::str_to_upper(),
                     Sample_Name = paste( Base_Sample_Name,Sample_Concentration, sep="_"),
                     Sentrix_Name = Sample_ID, # paste(Sentrix_ID,Sentrix_Position, sep="_"),
                     Sentrix_ID = Sentrix_Name %>% stringr::str_remove("_.*$"),
                     Sentrix_Path = paste0( opt$top_path,"/data/idats/idats_MSAv03_48x1_Alpha/",Sentrix_ID,"/",Sentrix_Name ) ) %>%
      dplyr::select( Sentrix_Name, Sample_Name,Base_Sample_Name,Sample_Concentration, Sentrix_Path )
    
    ssh_rep_sum <- ssh_rep_tib %>% 
      dplyr::group_by( Base_Sample_Name, Sample_Concentration ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    if ( p1 ) ssh_rep_sum %>% print( n=base::nrow(ssh_rep_sum) )
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                        Load Manifest:: ACA A2
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  ACA_cols <- NULL
  ACA_cols <- readr::cols(
    Address       = readr::col_integer(),
    Type          = readr::col_character(),
    Color_Channel = readr::col_character(),
    Name          = readr::col_character(),
    Probe_ID      = readr::col_character()
  )
  
  # Update:: Use A2::
  #  CMD = tail -n 8801 Projects.new/MSA/data/MSA-Interm-48v0-3/MSA-Interm-48v0-3_BP123_A2.csv | perl -pe 's/,+\r+$//gi' > Projects.new/MSA/data/MSA-Interm-48v0-3/MSA-Interm-48v0-3_BP123_A2.ctl-only.csv
  man_A2_ctl_tib <- NULL
  man_A2_ctl_csv <- file.path( opt$top_path, "Projects.new/MSA/data/MSA-Interm-48v0-3/MSA-Interm-48v0-3_BP123_A3.ctl-only.csv" )
  man_A2_ctl_tib <- readr::read_csv( file = man_A2_ctl_csv, 
                                     col_names = names(ACA_cols$cols), 
                                     col_types = ACA_cols ) %>% clean_tib()
  
  # Update:: Use A3::
  #  CMD = tail -n 8801 Projects.new/MSA/data/MSA-Interm-48v0-3/MSA-Interm-48v0-3_BP123_A3.csv | perl -pe 's/,+\r+$//gi' > Projects.new/MSA/data/MSA-Interm-48v0-3/MSA-Interm-48v0-3_BP123_A3.ctl-only.csv
  man_A3_ctl_tib <- NULL
  man_A3_ctl_csv <- file.path( opt$top_path, "Projects.new/MSA/data/MSA-Interm-48v0-3/MSA-Interm-48v0-3_BP123_A3.ctl-only.csv" )
  man_A3_ctl_tib <- readr::read_csv( file = man_A3_ctl_csv, 
                                     col_names = names(ACA_cols$cols), 
                                     col_types = ACA_cols ) %>% clean_tib()
  
  #
  # Compare ACA: A2 vs. A3
  #
  #  man_A2_ctl_tib %>% dplyr::inner_join( man_A3_ctl_tib, by=c("Address"), suffix=c("_A2","_A3") ) %>% dplyr::filter( Type_A2 != Type_A3 )
  #
  # Here is the real difference::
  # CMD:: diff Projects.new/MSA/data/MSA-Interm-48v0-3/MSA-Interm-48v0-3_BP123_A2.csv  Projects.new/MSA/data/MSA-Interm-48v0-3/MSA-Interm-48v0-3_BP123_A3.csv
  #
  
  tar_A2_ctl_tib <- NULL
  tar_A2_ctl_csv <- file.path( opt$top_path, "Projects.new/MSA/data/from-Jenn-MSA3/MSAV03_A2_v2/ctrl_probes_errors.csv" )
  tar_A2_ctl_tib <- readr::read_csv( file = tar_A2_ctl_csv, 
                                     show_col_types = FALSE ) %>% clean_tib()
  
  #
  # TBD:: Compare Control differences:: A3 vs. Bret
  #
  #  use_gs_ctl_tib %>% dplyr::filter( !Address %in% man_A3_ctl_tib$Address )
  #  man_A3_ctl_tib %>% dplyr::filter( !Address %in% use_gs_ctl_tib$Address )
  #  use_gs_ctl_tib %>% dplyr::inner_join( man_A3_ctl_tib, by=c("Address"), suffix = c("_B1","_A3") ) %>% dplyr::filter( Type_B1 != Type_A3 )
  
  man_A2_ctl_sum <- man_A2_ctl_tib %>% dplyr::group_by( Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  man_A3_ctl_sum <- man_A3_ctl_tib %>% dplyr::group_by( Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  use_gs_ctl_sum <- use_gs_ctl_tib %>% dplyr::group_by( Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  
  type_ctl_dif_tib <- NULL
  type_ctl_dif_tib <- man_A3_ctl_tib %>% 
    dplyr::filter( Type != "STAINING" ) %>% 
    dplyr::filter( Type != "EXTENSION" ) %>% 
    dplyr::filter( Type != "HYBRIDIZATION" ) %>% 
    dplyr::filter( Type != "TARGET REMOVAL" ) %>% 
    dplyr::inner_join( 
      use_gs_ctl_tib,
      by=c("Address"), 
      suffix = c("_A3","_B1") ) %>% 
    dplyr::filter( Type_B1 != Type_A3 )
  
  type_ctl_dif_sum <- NULL
  type_ctl_dif_sum <- type_ctl_dif_tib %>% 
    dplyr::group_by( Type_A3,Type_B1) %>%
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  #
  # Conclusion from above::
  #.  - ACA uses best probe and calls the rest some sort repeat.
  #.    - This is what causes the better performance
  #.  - Bret provides all of the controls. 
  #.    - These probes need to be screened... This should provide better performance
  #.  - This was missed by Bret due to not understanding the naming schemes.
  #   - One Restoration and 8 SPECIFICY I probes are named incorrectly in ACA
  #
  #
  # Next Step::
  #
  #  - TBD: Box Plots
  #  - TBD: Compare Intensities against EPIC v1 for the same controls...
  #
  
  neg_ctl_tib <- NULL
  neg_ctl_tib <- use_gs_ctl_tib %>% dplyr::filter( Type == "NEGATIVE" )
  
  msa_prefixes <- NULL
  msa_pre_path <- file.path( opt$top_path, "data/idats/idats_MSAv03_48x1_Alpha" )
  msa_prefixes <- sesame::searchIDATprefixes( dir.name = msa_pre_path )
  
  workflow_vec <- NULL
  workflow_vec <- c( "i" )
  
  opt$return_df <- 0
  opt$return_df <- 1
  opt$return_df <- 2
  opt$return_df <- 3
  # opt$return_df <- 4
  
  opt$min_pval <- 0.05
  opt$min_beta <- 0.3
  opt$max_beta <- 0.7
  opt$min_perO <- 0.75
  opt$min_perI <- 0.05
  
  ssh_rep_lst <- NULL
  ssh_rep_lst <- ssh_rep_tib %>% split(.$Base_Sample_Name)
  
  ictl_tib <- NULL
  sample_list <- NULL
  sample_list <- names(ssh_rep_lst)
  # sample_list <- c( "HELA", "RAJI" )
  for ( sample_name in sample_list ) {
    # sample_name <- "HELA"
    
    cur_prefixes <- NULL
    cur_prefixes <- ssh_rep_lst[[sample_name]] %>% 
      dplyr::filter( Sentrix_Path %in% msa_prefixes )
    
    #for ( prefix in cur_prefixes$Sentrix_Path ) {
    for ( s_idx in c(1:base::nrow(cur_prefixes) ) ) {
      
      prefix <- cur_prefixes$Sentrix_Path[s_idx]
      sentrix_name <- cur_prefixes$Sentrix_Name[s_idx]
      
      cur_idat_tib <- NULL
      cur_idat_tib <- read_idat_pair_rcpp( 
        prefix_path  = prefix,
        output_path  = opt$out_path,
        workflow_vec = workflow_vec %>% as.vector(),
        
        pval_add_vec = neg_ctl_tib$Address %>% as.vector(), 
        addU_man_vec = msa_man_tib$U %>% as.vector(),
        addM_man_vec = msa_man_tib$M %>% as.vector(),
        
        cgns_man_vec = msa_man_tib$Probe_ID %>% as.vector(), 
        cols_man_vec = msa_man_tib$col %>% as.vector(),
        keys_man_vec = msa_man_tib$Manifest %>% as.vector(),
        anns_man_vec = msa_man_tib$Annotation %>% as.vector(), 
        chrs_man_vec = msa_man_tib$Chromosome %>% as.vector(),
        
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
        clean_tib() %>%
        dplyr::mutate(across(where(is.factor), as.character) ) %>%
        dplyr::mutate( Probe_Type = Probe_ID %>% stringr::str_sub(1,2) )
      
      # cur_idat_tib %>% dplyr::filter( Probe_Type == "ct" )
      # cur_idat_tib %>% dplyr::filter( Probe_ID %in% use_gs_ctl_tib$Probe_ID )
      # cur_idat_tib %>% dplyr::inner_join( use_gs_ctl_tib, by=c("Probe_ID"), multiple = "all" )
      
      cur_ictl_tib <- NULL
      cur_ictl_tib <- cur_idat_tib %>% 
        dplyr::inner_join( msa_ctl_tab, by=c("Probe_ID"), multiple = "all" ) %>% 
        dplyr::rename( Probe_ID_Src = Probe_ID ) %>% 
        dplyr::inner_join( use_gs_ctl_tib, by=c("Address"), multiple = "all" ) %>% 
        dplyr::mutate( Sample_Name  = sample_name, 
                       Sentrix_Name = sentrix_name )
      
      cur_ictl_sum <- NULL
      cur_ictl_sum <- cur_ictl_tib %>% 
        dplyr::group_by( Type ) %>% 
        dplyr::summarise( Count=n(), .groups = "drop" )
      
      ictl_tib <- ictl_tib %>% dplyr::bind_rows( cur_ictl_tib )
      
      if ( FALSE ) {
        grn_idat <- file.path( paste0( prefix, "_Grn.idat.gz" ) )
        red_idat <- file.path( paste0( prefix, "_Red.idat.gz" ) )
        
        red_dat <- NULL
        red_dat <- illuminaio::readIDAT( file = red_idat )
        red_tib <- NULL
        red_tib <- red_dat$Quants %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column( var = "Address" ) %>% 
          tibble::as_tibble() %>% clean_tib()
        
        grn_dat <- NULL
        grn_dat <- illuminaio::readIDAT( file = grn_idat )
        grn_tib <- NULL
        grn_tib <- grn_dat$Quants %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column( var = "Address" ) %>% 
          tibble::as_tibble() %>% clean_tib()
        
        idat_tib <- NULL
        idat_tib <- dplyr::inner_join( grn_tib,red_tib, 
                                       by=c("Address"), suffix=c("_Grn","_Red") )
        
        b2_ctl_tib <- NULL
        b2_ctl_tib <- use_gs_ctl_tib %>% dplyr::filter( Type %in% tar_A2_ctl_tib$Type )
        
        b2_dat_tib <- NULL
        b2_dat_tib <- idat_tib %>% 
          dplyr::inner_join( b2_ctl_tib, by=c("Address") )
        
        #
        # Points::
        #
        plots_dir <- safe_mkdir( dir = file.path( opt$out_path, "plots") )
        gg_pnts_pdf <- file.path( plots_dir, "BISULFITE_CONVERSION_II.pdf" )
        
        gg_pnts <- NULL
        gg_pnts <- b2_dat_tib %>% 
          dplyr::filter( Type == "BISULFITE CONVERSION II" ) %>%
          dplyr::mutate(
            Flagged = dplyr::case_when(
              Address %in% tar_A2_ctl_tib$Address ~ TRUE,
              TRUE ~ FALSE
            )
          ) %>%
          dplyr::group_by( Name ) %>% 
          ggplot2::ggplot( aes( x=Mean_Grn, y=Mean_Red, group=Name, color=Name) ) + 
          ggplot2::geom_point() + 
          ggplot2::facet_grid( rows = vars(Flagged), 
                               cols = vars(Type) ) + 
          ggplot2::theme( legend.position='none' )
        
        ggplot2::ggsave( filename = gg_pnts_pdf, plot = gg_pnts, 
                         device = "pdf", 
                         # width = 680, height = 680, 
                         dpi = 320 )
      }
      
      # if ( TRUE ) break
    }
    
    
    # if ( opt$single ) break
    # if ( TRUE ) break
  }
  
  # Summary Stats::
  # ictl_tib %>% dplyr::group_by( Type,Sample_Name,Sentrix_Name ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  #
  plot_path <- safe_mkdir( dir = file.path( opt$out_path, "plots" ) )
  
  #
  # ACA3 and Core are the same for Negative Controls...
  #
  ictl_neg_tib <- NULL
  ictl_neg_tib <- ictl_tib %>% 
    dplyr::filter( Type == "NEGATIVE" ) %>%
    dplyr::filter( mask ) %>% 
    dplyr::select( Address, UG,UR, Probe_ID,Name,Sentrix_Name,Sample_Name ) %>% 
    dplyr::rename( Grn=UG, Red=UR) %>% 
    tidyr::pivot_longer( cols = c(Grn,Red), 
                         names_to = "Color", 
                         values_to = "Intentisty" )
  
  ictl_neg_den_ggg <- NULL
  ictl_neg_den_pdf <- file.path( plot_path, "NEGATIVE.both.den.pdf" )
  ictl_neg_den_ggg <- ictl_neg_tib %>% 
    dplyr::group_by( Sentrix_Name ) %>%
    dplyr::mutate( Sentrix_Name = Sentrix_Name %>% as.factor() ) %>%
    # dplyr::mutate( Color = Color %>% as.factor() ) %>%
    ggplot2::ggplot( aes( x=Intentisty, color=Color, group=Sentrix_Name ) ) +
    ggplot2::geom_density( alpha = 0.4 ) + 
    ggplot2::facet_grid( rows = vars(Sample_Name),
                         cols = vars(Color) ) +
    ggplot2::theme( legend.position='none' )
  
  ggplot2::ggsave( filename = ictl_neg_den_pdf, plot = ictl_neg_den_ggg, 
                   device = "pdf", dpi = 320, width = 7, height = 7 )
  
  # ictl_neg3_den_ggg <- NULL
  # ictl_neg3_den_pdf <- file.path( plot_path, "NEGATIVE.aca3.den.pdf" )
  # ictl_neg3_den_ggg <- ictl_neg_tib %>% 
  #   dplyr::group_by( Sentrix_Name ) %>%
  #   dplyr::mutate( Sentrix_Name = Sentrix_Name %>% as.factor() ) %>%
  #   # dplyr::mutate( Color = Color %>% as.factor() ) %>%
  #   ggplot2::ggplot( aes( x=Intentisty, color=Color, group=Sentrix_Name ) ) +
  #   ggplot2::geom_density( alpha = 0.4 ) + 
  #   ggplot2::facet_grid( rows = vars(Sample_Name),
  #                        cols = vars(Color) ) +
  #   ggplot2::theme( legend.position='none' )
  # 
  # ggplot2::ggsave( filename = ictl_neg3_den_pdf, plot = ictl_neg3_den_ggg, 
  #                  device = "pdf", dpi = 320, width = 7, height = 7 )
  
  #
  # BS II:: TBD: Plot Grn vs. Red...
  #
  ictl_bs2_tib <- NULL
  ictl_bs2_tib <- ictl_tib %>% 
    dplyr::filter( Type == "BISULFITE CONVERSION II" ) %>%
    # dplyr::filter( mask ) %>% 
    dplyr::select( Address, UG,UR, Probe_ID,Name,Sentrix_Name,Sample_Name ) %>% 
    dplyr::rename( Grn=UG, Red=UR) # %>% 
  # tidyr::pivot_longer( cols = c(Grn,Red), 
  #                      names_to = "Color", 
  #                      values_to = "Intentisty" )
  
  ictl_bs2_den_ggg <- NULL
  ictl_bs2_den_pdf <- file.path( plot_path, "NEGATIVE.umsk.pnt.pdf" )
  ictl_bs2_den_ggg <- ictl_bs2_tib %>% 
    dplyr::group_by( Sentrix_Name ) %>%
    dplyr::mutate( Sentrix_Name = Sentrix_Name %>% as.factor() ) %>%
    # dplyr::mutate( Color = Color %>% as.factor() ) %>%
    ggplot2::ggplot( aes( x=Grn, y=Red, color=Sentrix_Name, group=Sentrix_Name ) ) +
    ggplot2::geom_point(  ) + 
    # ggplot2::geom_density( alpha = 0.4 ) + 
    ggplot2::facet_grid( rows = vars(Sample_Name) ) +
    # cols = vars(Color) ) +
    ggplot2::theme( legend.position='none' )
  
  ggplot2::ggsave( filename = ictl_bs2_den_pdf, plot = ictl_bs2_den_ggg, 
                   device = "pdf", dpi = 320, width = 7, height = 7 )
  
  
  
  
  if ( FALSE ) {
    ictl_neg_den_gg <- ictl_neg_tib %>% 
      dplyr::group_by( Name ) %>%
      dplyr::mutate( Name = Name %>% as.factor() ) %>%
      # dplyr::mutate( Color = Color %>% as.factor() ) %>%
      ggplot2::ggplot( aes( x=Intentisty, color=Color, group=Name ) ) +
      ggplot2::geom_density( alpha = 0.4 ) + 
      ggplot2::facet_grid( rows = vars(Color) ) +
      ggplot2::theme( legend.position='none' )
    
    ictl_neg_box_gg <- NULL
    ictl_neg_box_gg <- ictl_neg_tib %>% 
      dplyr::group_by( Name,Color ) %>% 
      dplyr::arrange( Name,Color ) %>% 
      # dplyr::filter( Color == "Red" ) %>%
      dplyr::add_tally( wt = Intentisty, name = "Sum_Int" ) %>% 
      dplyr::add_count( Name,Color, name="Count") %>% 
      dplyr::mutate( Int_Avg = Sum_Int / Count ) %>% 
      dplyr::arrange( Int_Avg ) %>%
      dplyr::mutate( Row = dplyr::row_number() ) %>%
      # dplyr::filter( Row %% 177 == 0 ) %>%
      dplyr::arrange( Int_Avg ) %>%
      head(n=2000) %>%
      ggplot2::ggplot( aes( Name, Intentisty ) ) +
      # ggplot2::ggplot( aes( Intentisty, Name ) ) +
      ggplot2::geom_boxplot(  ) + 
      ggplot2::facet_grid( rows = vars(Color) ) +
      ggplot2::theme( legend.position='none' )
    
    
    ictl_neg_box_gg <- ictl_neg_tib %>% 
      dplyr::arrange( Name )
    # dplyr::group_by( Name ) %>%
    # dplyr::mutate( Name = Name %>% as.factor() ) %>%
    # dplyr::mutate( Color = Color %>% as.factor() ) %>%
    ggplot2::ggplot( aes( Name, Intentisty ) ) +
      ggplot2::geom_boxplot(  ) + 
      ggplot2::facet_grid( rows = vars(Color) ) +
      ggplot2::theme( legend.position='none' )
  }
  
}






if ( FALSE ) {
  
  # illuminaio::readIDAT()
  red_idat <- file.path( opt$top_path, "data/idats/idats_MSAv03_48x1_Alpha/207545400007/207545400007_R01C01_Red.idat.gz" )
  grn_idat <- file.path( opt$top_path, "data/idats/idats_MSAv03_48x1_Alpha/207545400007/207545400007_R01C01_Grn.idat.gz" )
  
  red_dat <- NULL
  red_dat <- illuminaio::readIDAT( file = red_idat )
  red_tib <- NULL
  red_tib <- red_dat$Quants %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column( var = "Address" ) %>% 
    tibble::as_tibble() %>% clean_tib()
  
  grn_dat <- NULL
  grn_dat <- illuminaio::readIDAT( file = grn_idat )
  grn_tib <- NULL
  grn_tib <- grn_dat$Quants %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column( var = "Address" ) %>% 
    tibble::as_tibble() %>% clean_tib()
  
  idat_tib <- NULL
  idat_tib <- dplyr::inner_join( grn_tib,red_tib, 
                                 by=c("Address"), suffix=c("_Grn","_Red") )
  
  b2_ctl_tib <- NULL
  b2_ctl_tib <- use_gs_ctl_tib %>% dplyr::filter( Type %in% tar_A2_ctl_tib$Type )
  # b2_ctl_tib <- use_gs_ctl_tib %>% dplyr::filter( Type == "BISULFITE CONVERSION II" )
  #  man_A2_ctl_tib %>% dplyr::filter( Type == "BISULFITE CONVERSION II" )
  #  use_gs_ctl_tib %>% dplyr::filter( Type == "BISULFITE CONVERSION II" )
  
  # idat_tib %>% dplyr::filter( Address %in% b2_ctl_tib$Address )
  
  b2_dat_tib <- NULL
  b2_dat_tib <- idat_tib %>% dplyr::inner_join( b2_ctl_tib, by=c("Address") )
  # b2_dat_tib <- idat_tib %>% dplyr::inner_join( use_gs_ctl_tib, by=c("Address") )
  
  
  #
  # Points::
  #
  plots_dir <- safe_mkdir( dir = file.path( opt$top_path, "plots") )
  gg_pnts_pdf <- file.path( plots_dir, "BISULFITE_CONVERSION_II.pdf" )
  
  gg_pnts <- NULL
  gg_pnts <- b2_dat_tib %>% 
    dplyr::filter( Type == "BISULFITE CONVERSION II" ) %>%
    dplyr::mutate(
      Flagged = dplyr::case_when(
        Address %in% tar_A2_ctl_tib$Address ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    dplyr::group_by( Name ) %>% 
    ggplot2::ggplot( aes( x=Mean_Grn, y=Mean_Red, group=Name, color=Name) ) + 
    ggplot2::geom_point() + 
    ggplot2::facet_grid( rows = vars(Flagged), 
                         cols = vars(Type) ) + 
    ggplot2::theme( legend.position='none' )
  
  ggplot2::ggsave( filename = gg_pnts_pdf, plot = gg_pnts, 
                   device = "pdf", 
                   # width = 680, height = 680, 
                   dpi = 320 )
  
  #
  # Density::
  #
  gg_dens <- NULL
  gg_dens <- dplyr::bind_rows( dplyr::mutate( grn_tib, Color="Grn" ), 
                               dplyr::mutate( red_tib, Color="Red" ) ) %>% 
    dplyr::inner_join( b2_ctl_tib, by=c("Address") ) %>% 
    dplyr::mutate(
      Flagged = dplyr::case_when(
        Address %in% tar_A2_ctl_tib$Address ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    dplyr::group_by( Name ) %>% 
    ggplot2::ggplot( aes( x=Mean, group=Color, color=Name) ) + 
    ggplot2::geom_density(  ) + 
    ggplot2::facet_grid( rows = vars(Flagged), 
                         cols = vars(Type), scales = "free" ) + 
    ggplot2::theme(legend.position='none')
  
  # gg_dens
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                          Other Scratch Code::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  #
  # What's missing::
  #
  mis_class_tib <- NULL
  mis_class_tib <- all_class_tib %>% 
    dplyr::filter( !Probe_Seq %in% use_gs_ctl_tib$AlleleA_ProbeSeq )
  
  mis_class_sum <- NULL
  mis_class_sum <- mis_class_tib %>% 
    dplyr::group_by( Control_Type_match,Control_Class_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  
  if ( p2 ) cat(glue::glue("{pmssg} mis_class_tib Summary::{RET}"))
  if ( p2 ) mis_class_sum %>% print( n=base::nrow(mis_class_sum) )
  if ( p2 ) cat(glue::glue("{RET}"))
  
  #
  # Conclusion:: We need to build the proper Master_Match_Human_Controls
  #.  - Single
  #.  - Paired
  #
  
  # This is the single table version, showing what's missing
  #  - Group by Control_Class_match and evaluate each...
  #
  exp_sig_class_tib <- NULL
  exp_sig_class_tib <- use_gs_ctl_tib %>% 
    dplyr::filter( !is.na(AlleleA_ProbeSeq) ) %>% 
    dplyr::distinct( AlleleA_ProbeSeq, .keep_all = TRUE ) %>% 
    dplyr::right_join( all_class_tib, by=c("AlleleA_ProbeSeq"="Probe_Seq") )
  
  exp_sig_class_lst <- NULL
  exp_sig_class_lst <- exp_sig_class_tib %>% split( .$Control_Class_match )
  
  # TBD:: Apply the same rules above to all the missing categories::
  #  Example: exp_sig_class_lst[["B1"]] %>% dplyr::filter( is.na(Address) )
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Post Processing:: Write Manifests
#
# 1. Join all Genome Studio Controls::
# 2. Remove Probes from manifest::
# 3. Load or split: "data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_BP123_A1.csv.gz"
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

write_manifests <- TRUE
write_manifests <- FALSE

if ( write_manifests ) {
  
  all_gs_ctl_tib <- NULL
  all_gs_ctl_tib <- msa_man_tib
  
  if ( FALSE ) {
    # all_gs_ctl_tib <- msa_man_tib %>% head() %>% as.data.frame()
    # all_gs_ctl_tib <- NULL
    # all_gs_ctl_tib <- msa_man_tib %>% head() %>% as.data.frame()
    # msa_man_tib %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
    
    # A tibble: 1,598 × 2
    unq_gs_ctl_seq_tib <- NULL
    unq_gs_ctl_seq_tib <- all_gs_ctl_tib %>% 
      dplyr::filter( !is.na(AlleleA_ProbeSeq) ) %>% 
      dplyr::distinct( AlleleA_ProbeSeq,AlleleB_ProbeSeq )
    
    # A tibble: 8,351 × 1
    unq_gs_ctl_add_tib <- NULL
    unq_gs_ctl_add_tib <- all_gs_ctl_tib %>% 
      dplyr::filter( !is.na(AlleleA_ProbeSeq) ) %>% 
      dplyr::distinct( Address )
    
    # A tibble: 8,351 × 3
    unq_gs_ctl_asq_tib <- NULL
    unq_gs_ctl_asq_tib <- all_gs_ctl_tib %>% 
      dplyr::filter( !is.na(AlleleA_ProbeSeq) ) %>% 
      dplyr::distinct( Address,AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  }
  
  #
  # 2. Remove Probes from manifest::
  #
  msa_gs_body_cln_tib <- NULL
  msa_gs_body_cln_tib <- msa_gs_body_src_tib %>% 
    dplyr::filter( !AddressA_ID %in% use_gs_ctl_tib$Address ) # %>% 
  # dplyr::filter( !is.na(AlleleB_ProbeSeq) & !AddressB_ID %in% use_gs_ctl_tib$Address ) %>%
  # dplyr::filter( IlmnID %>% stringr::str_starts("ZZ") )
  
  # msa_gs_body_cln_tib <- msa_gs_body_src_tib %>% 
  #   dplyr::filter( 
  #     ( !is.na(AlleleA_ProbeSeq) & !AddressA_ID %in% all_gs_ctl_tib$Address ) &
  #       ( !is.na(AlleleB_ProbeSeq) & !AddressB_ID %in% all_gs_ctl_tib$Address ) )
  # msa_gs_body_cln_tib %>% dplyr::filter( IlmnID %>% stringr::str_starts("ZZ") )
  # msa_gs_body_cln_tib %>% dplyr::filter( IlmnID %>% stringr::str_starts("ZZ") ) %>% as.data.frame()
  # msa_gs_body_cln_tib %>% dplyr::filter( IlmnID %>% stringr::str_starts("ZZ") ) %>% dplyr::group_by( Probe_Type ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  
  #
  # Thought:: To remove or not to remove: that's the analytical categorization...
  #
  #. - Provide both all[ removed and not removed ]
  #. - Provide both top[ removed and not removed ]
  #
  #  TBD: Full version of standard 2020 controls 
  #. - 122mers and all designs
  #. - improbe thermodynamic scores
  #. - Anntotation for all 2020 probes...
  #
  
  
  #
  # TBD:: Update msa_gs_body_src_tib like msa_man_tib!!!
  #
  
  man_out_path <- safe_mkdir( dir = file.path( opt$out_path, "manifests") )
  ctl_gs_line_src_line <- "[Controls],,,,,,,,,,,,,,,,,,,,,,,,,,"
  # msa_gs_ctls_src_tib  %>% readr::write_csv(   file = msa_gs_ctlA_all_csv, append = TRUE )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Write Manifest:: + Analytical Controls: All
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  man_cntA_all_dir <- safe_mkdir( dir = file.path( man_out_path, "cntlA-all" ) )
  msa_gs_ctlA_all_csv  <- file.path( man_cntA_all_dir, "MSA-Interm-48v0-3_BP123_A1.genome-studio.cntlA-all.csv.gz" )
  
  # Write Data::
  msa_gs_head_src_line %>% readr::write_lines( file = msa_gs_ctlA_all_csv, append = FALSE )
  msa_gs_body_src_tib  %>% readr::write_csv(   file = msa_gs_ctlA_all_csv, append = TRUE, col_names = TRUE )
  ctl_gs_line_src_line %>% readr::write_lines( file = msa_gs_ctlA_all_csv, append = TRUE )
  use_gs_ctl_tib  %>% readr::write_csv( file = msa_gs_ctlA_all_csv, append = TRUE )
  
  # Reload data::
  dbl_gs_ctlA_all_tib <- NULL
  dbl_gs_ctlA_all_tib <- load_genome_studio_manifest( 
    file = msa_gs_ctlA_all_csv, load_controls = TRUE, load_clean = TRUE, 
    write_clean = FALSE, overwrite = FALSE, ret_data = TRUE )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Write Manifest:: - Analytical Controls: All
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  man_cnt0_all_dir <- safe_mkdir( dir = file.path( man_out_path, "cntl0-all" ) )
  msa_gs_ctl0_all_csv  <- file.path( man_cnt0_all_dir, "MSA-Interm-48v0-3_BP123_A1.genome-studio.cntl0-all.csv.gz" )
  
  msa_gs_head_src_line %>% readr::write_lines( file = msa_gs_ctl0_all_csv, append = FALSE )
  msa_gs_body_cln_tib  %>% readr::write_csv(   file = msa_gs_ctl0_all_csv, append = TRUE, col_names = TRUE )
  ctl_gs_line_src_line %>% readr::write_lines( file = msa_gs_ctl0_all_csv, append = TRUE )
  use_gs_ctl_tib  %>% readr::write_csv( file = msa_gs_ctl0_all_csv, append = TRUE )
  
  # Reload data::
  dbl_gs_ctl0_all_tib <- NULL
  dbl_gs_ctl0_all_tib <- load_genome_studio_manifest( 
    file = msa_gs_ctl0_all_csv, load_controls = TRUE, load_clean = TRUE, 
    write_clean = FALSE, overwrite = FALSE, ret_data = TRUE )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Write Manifest:: + No Controls
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  man_cnt_NULL_dir <- safe_mkdir( dir = file.path( man_out_path, "cntl-NULL" ) )
  msa_gs_ctl_NULL_csv  <- file.path( man_cnt_NULL_dir, "MSA-Interm-48v0-3_BP123_A1.genome-studio.cntl-NULL.csv.gz" )
  
  msa_gs_head_src_line %>% readr::write_lines( file = msa_gs_ctl_NULL_csv, append = FALSE )
  msa_gs_body_cln_tib  %>% readr::write_csv(   file = msa_gs_ctl_NULL_csv, append = TRUE, col_names = TRUE )
  ctl_gs_line_src_line %>% readr::write_lines( file = msa_gs_ctl_NULL_csv, append = TRUE )
  # use_gs_ctl_tib  %>% readr::write_csv( file = msa_gs_ctl_NULL_csv, append = TRUE )
  
  # Reload data::
  dbl_gs_ctl0_all_tib <- NULL
  dbl_gs_ctl0_all_tib <- load_genome_studio_manifest( 
    file = msa_gs_ctl_NULL_csv, load_controls = TRUE, load_clean = TRUE, 
    write_clean = FALSE, overwrite = FALSE, ret_data = TRUE )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                      Write Manifest:: + No Controls
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  man_cnt_EPIC_dir <- safe_mkdir( dir = file.path( man_out_path, "cntl-EPIC" ) )
  msa_gs_ctl_EPIC_csv  <- file.path( man_cnt_EPIC_dir, "MSA-Interm-48v0-3_BP123_A1.genome-studio.cntl-EPIC.csv.gz" )
  
  msa_gs_head_src_line %>% readr::write_lines( file = msa_gs_ctl_EPIC_csv, append = FALSE )
  msa_gs_body_cln_tib  %>% readr::write_csv(   file = msa_gs_ctl_EPIC_csv, append = TRUE, col_names = TRUE )
  ctl_gs_line_src_line %>% readr::write_lines( file = msa_gs_ctl_EPIC_csv, append = TRUE )
  # use_gs_ctl_tib %>% readr::write_csv( file = msa_gs_ctl_EPIC_csv, append = TRUE )
  use_gs_epi_tib <- use_gs_ctl_tib %>% dplyr::filter( !is.na(AlleleA_ProbeSeq) ) %>% dplyr::inner_join( all_class_tib %>% dplyr::filter( !is.na(Probe_ID_epic) ), by=c("AlleleA_ProbeSeq"="Probe_Seq") ) %>% dplyr::select( Address:AlleleB_ProbeSeq )
  use_gs_epi_tib %>% readr::write_csv( file = msa_gs_ctl_EPIC_csv, append = TRUE )

  # Reload data::
  dbl_gs_ctl0_all_tib <- NULL
  dbl_gs_ctl0_all_tib <- load_genome_studio_manifest( 
    file = msa_gs_ctl_EPIC_csv, load_controls = TRUE, load_clean = TRUE, 
    write_clean = FALSE, overwrite = FALSE, ret_data = TRUE )
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Write Manifest:: + Analytical Controls: Top
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                   Write Manifest:: - Analytical Controls: Top
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
}


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                            Controls Format::
#                                   END
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #









# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                                   BEG
#                       Controls:: Probe Alignment
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Next Steps::
#.  - Align all controls sequences using U+Sequence
#.  - Sub-string all Negative controls using Probe_ID_match => bisulfite convert => count ACTG's
#

# TBD:: These Parameters Should be Moved to Program Inputs::
#
ref_tar <- "hg19"
ref_tar <- "GRCh37"
bsp_exe <- file.path( opt$top_path, "tools/programs/BSMAPz/bsmapz" )
opt$bsmap_exe <- bsp_exe

align_ctls <- TRUE
align_ctls <- FALSE

if ( align_ctls ) {
  
  run_bsp <- FALSE
  run_bsp <- TRUE
  
  if ( run_bsp ) {
    if ( p1 ) cat(glue::glue("{pmssg} BEG: Launching BSMAP Alignment...{RET}"))
    
    fas_path <- NULL
    bed_path <- NULL
    bsp_path <- NULL
    
    fas_path <- safe_mkdir( dir = file.path( opt$out_path, "fas" ) )
    bed_path <- safe_mkdir( dir = file.path( opt$out_path, "bed" ) )
    bsp_path <- safe_mkdir( dir = file.path( opt$out_path, "bsp" ) )
    
    prb_fas <- NULL
    prb_fas <- file.path( fas_path, "Infinium_Methylation_Controls_1971.seq-tan.fa.gz" )
    seqs_ctl_tib %>% 
      dplyr::mutate( FAS = paste0(">",Probe_Seq,"_",U,"\n",Probe_Seq) ) %>%
      dplyr::pull( FAS ) %>% 
      as.vector() %>% 
      readr::write_lines( file = prb_fas )
    
    prb_ids_vec <- NULL
    prb_ids_vec <- seqs_ctl_tib %>% 
      dplyr::mutate( bsp_id = paste0(Probe_Seq,"_",U) ) %>% 
      dplyr::pull( bsp_id )
    
    #
    # Run BSMAP:: ids
    #
    bsp_tsv <- NULL
    bsp_tsv <- run_bsmap( ref_fas = ref_seqs[[ref_tar]],
                          can_fas = prb_fas,
                          bsp_exe = opt$bsmap_exe,
                          bsp_dir = opt$bsmap_dir,
                          slim = TRUE,
                          out_dir = bsp_path,
                          run_tag = opt$run_name,
                          reload = opt$reload,
                          reload_min = 10,
                          vb=vb,vt=vt,tc=tc,tt=tt )
    
    bsp_tib <- NULL
    bsp_tib <- load_bsmap( file = bsp_tsv,
                           sort = FALSE,
                           add_cnt = TRUE,
                           parse_id = FALSE,
                           out_dir = bsp_path,
                           run_tag = opt$run_name,
                           reload = opt$reload,
                           reload_min = 10,
                           vb=vb,vt=vt,tc=tc,tt=tt )
    
    # QUESTION: Why is there only 1,071 unique seqs?
    #
    # bsp_tib %>% dplyr::distinct( Probe_ID )
    
    #
    # Write BED::
    #
    bsp_bed <- NULL
    bsp_bed <- file.path( bed_path, "Infinium_Methylation_Controls_1971.seq-tan.bed" )
    bsp_tib %>% dplyr::select( Chr,Beg,Probe_ID ) %>% 
      dplyr::mutate( Beg = Beg -60, End = Beg + 122 ) %>% 
      dplyr::select( Chr,Beg,End, Probe_ID ) %>%
      dplyr::arrange( Chr,Beg,End ) %>%
      readr::write_tsv( file = bsp_bed, col_names = FALSE )
    
    #
    # Extract 122mer::
    #
    fwd_122_fas <- file.path( fas_path, "Infinium_Methylation_Controls_1971.seq-tan.122mer.fa.gz" )
    twoBitToFa_exe <- "/Users/bbarnes/Documents/tools/ucsc/twoBitToFa"
    twoBitToFa_cmd <- paste0( twoBitToFa_exe,
                              " -bed=",bsp_bed,
                              paste0(" ",ref_seqs[[ref_tar]] %>% stringr::str_remove(".gz$"),".2bit"),
                              " ",fwd_122_fas )
    base::system( command = twoBitToFa_cmd )
    
    ids_122_dat <- NULL
    ids_122_dat <- Biostrings::readDNAStringSet( filepath = fwd_122_fas, format = "fasta" )
    # ids_122_dat@ranges
    
    dat_122_tib <- NULL
    dat_122_tib <- dplyr::bind_cols(
      ids_122_dat@ranges %>% as.data.frame() %>% tibble::as_tibble(),
      ids_122_dat %>% as.data.frame() %>% 
        tibble::as_tibble() %>% 
        magrittr::set_names( value = c("Forward_Sequence") )    
    )
    
    if ( p1 ) cat(glue::glue("{pmssg} END: Launching BSMAP Alignment.{RET2}"))
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
