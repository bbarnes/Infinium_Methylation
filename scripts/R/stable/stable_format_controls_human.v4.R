
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
# par$version <- 2
# par$version <- 3
# par$version <- 4
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
#.                               Workflow::
#
# 1. Write all HSA 2020 controls
# 2. Write probe_to_122mer() in function/trifecta_functions.R
#    - Test on HSA 2020 controls and { cg, ch, rs }
# 3. Update MSA Manifest
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

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
#                       Pre-processing:: Manifests
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
neg_org_tib  <- readr::read_rds( neg_ctl_rds )

# print( neg_ctl_tib )

epic_ctl_tib <- NULL
epic_ctl_rds <- file.path( opt$top_path, "data/manifests/methylation/bgz/epic_ctls.rds" )
epic_ctl_tib <- readr::read_rds( file = epic_ctl_rds ) %>%
  dplyr::select( Probe_ID,U,Type ) %>% 
  dplyr::distinct( U, .keep_all = TRUE ) %>% 
  dplyr::rename( Probe_ID_epic = Probe_ID, 
                 Probe_Type_epic = Type )

#
# Example Genome Studio Format Below::
#
example_gs_ctl_tib <- NULL
example_gs_ctl_csv <- file.path( opt$top_path, "Projects.new/MSA/MethylationQC/MSAv03/manifest/MSAv03_controls.csv.gz" )
example_gs_ctl_tib <- readr::read_csv( file = example_gs_ctl_csv, show_col_types = FALSE )

#
# NOTE:: This is pretty much garabage at this point. I'll keep in the stable
#. scratch as a record...
#
# Original Formatted File with Names::
#
if ( FALSE ) {
  names_ctl_tib <- NULL
  names_ctl_csv <- file.path( opt$top_path, "data/Controls/Infinium_Methylation_Controls_15_1983_sesame.csv.gz" )
  names_ctl_tib <- readr::read_csv( file = names_ctl_csv, show_col_types = FALSE )
  
  names_tmp_csv <- names_ctl_csv <- file.path( opt$top_path, "data/Controls/tmp/Infinium_Methylation_Controls_15_1983_sesame.tmp.csv" )
  names_tmp_tib <- names_ctl_tib %>% 
    dplyr::mutate( Probe_ID = Probe_ID %>% 
                     stringr::str_remove("^ctl_") %>%
                     stringr::str_remove("_[^_]+$") ) %>%
    dplyr::select( U,Probe_ID )
  readr::write_csv( x = names_tmp_tib, file = names_tmp_csv )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                                   BEG
#                    Controls:: Correct Proccess Below
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
    dplyr::filter( !( U == 22711390 | U == 46651360 ) )
  
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
  if ( p2 ) cat(glue::glue("{pmssg} Full Control Summary::{RET}"))
  if ( p2 ) full_ctl_sum %>% print( n=nrow(full_ctl_sum) )
  
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
    dplyr::select( Probe_ID_match,U,Probe_Seq,Last_Base, dplyr::everything() )
  
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
  if ( p2 ) cat(glue::glue("{pmssg} Control Type::{RET}"))
  if ( p2 ) match_ctl_sum %>% print( n=nrow(match_ctl_sum) )
  if ( p2 ) cat(glue::glue("{RET}"))
  
  match_cls_sum <- NULL
  match_cls_sum <- match_ctl_tib %>% 
    dplyr::group_by( Control_Class_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) cat(glue::glue("{pmssg} Control Class::{RET}"))
  if ( p2 ) match_cls_sum %>% print( n=nrow(match_cls_sum) )
  if ( p2 ) cat(glue::glue("{RET}"))
  
  epic_ctl_sum <- NULL
  epic_ctl_sum <- match_ctl_tib %>% 
    dplyr::filter( Is_EPIC_match ) %>%
    dplyr::group_by( Control_Type_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) cat(glue::glue("{pmssg} EPIC Control Type::{RET}"))
  if ( p2 ) epic_ctl_sum %>% print( n=nrow(epic_ctl_sum) )
  if ( p2 ) cat(glue::glue("{RET}"))
  
  if ( p1 ) cat(glue::glue("{pmssg} END: Formatting Raw Controls.{RET2}"))
}

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
# TBD:: 
#.  - Review the three data input sources: { match, epic, full } and add a new 
#.    column with a proper name (i.e. ctl_***)
#.  - Add missing data via bsmap alignment (122mer, etc.)
#.  - Clean columns to match standard manifest output. 
#

join_ctl_dir <- safe_mkdir( file.path( opt$out_path, "join") )
join_ctl_csv <- file.path( join_ctl_dir, "Infinium_Methylation_Controls_1971_join-names.raw.csv.gz" )
join_ctl_csv <- file.path( join_ctl_dir, "Infinium_Methylation_Controls_1971_join-names.raw.csv" )

join_ctl_tib <- NULL
join_ctl_tib <- match_ctl_tib %>% 
  dplyr::left_join( epic_ctl_tib, by=c("U") ) %>%
  dplyr::full_join( full_ctl_tib, by=c("U","Probe_Seq") ) %>%
  dplyr::arrange( Control_Class_match )
  # dplyr::arrange( Control_Type_match )
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
  # dplyr::group_by( Control_Type_match ) %>% 
  dplyr::group_by( Control_Class_match ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) cat(glue::glue("{pmssg} Seqs Control Type Summary::{RET}"))
if ( p2 ) seqs_ctl_sum %>% print( n=base::nrow(seqs_ctl_sum) )
if ( p2 ) cat(glue::glue("{RET}"))

#
# Build New Categorie: "Control_Type" and compare against "Probe_Type_full"
#. - Don't think this works because it base on previous EPIC data...
type_ctl_sum <- NULL
type_ctl_sum <- seqs_ctl_tib %>% 
  dplyr::group_by( Probe_Type_full ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )
if ( p2 ) cat(glue::glue("{pmssg} Seqs Control Class Summary::{RET}"))
if ( p2 ) type_ctl_sum %>% print( n=base::nrow(type_ctl_sum) )
if ( p2 ) cat(glue::glue("{RET}"))

#
# Rebuild the Control_Type method...
#

#
# Next Steps::
#.  - Align all controls sequences using U+Sequence
#.  - Sub-string all Negative controls using Probe_ID_match => bisulfite convert => count ACTG's
#
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                    Controls:: Correct Proccess Below
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

# TBD:: These Parameters Should be Moved to Program Inputs::
#
ref_tar <- "hg19"
ref_tar <- "GRCh37"
bsp_exe <- file.path( opt$top_path, "tools/programs/BSMAPz/bsmapz" )
opt$bsmap_exe <- bsp_exe

align_ctls <- FALSE
align_ctls <- TRUE

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
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                      Split/Join Infinium I Probes::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
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
    if ( p2 ) cat(glue::glue("{pmssg} (In)Dependent Control Type Summary::{RET}"))
    if ( p2 ) join_dep_sum %>% print( n=base::nrow(join_dep_sum) )
    if ( p2 ) cat(glue::glue("{RET}"))
    
    

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
    if ( p2 ) cat(glue::glue("{pmssg} All Control Class::{RET}"))
    if ( p2 ) all_class_sum %>% print( n=base::nrow(all_class_sum) )
    if ( p2 ) cat(glue::glue("{RET}"))

    all_cross_sum <- NULL
    all_cross_sum <- all_class_tib %>% 
      dplyr::group_by( Control_Class_match,Control_Type_match ) %>%
      dplyr::summarise( Count=n(), .groups = "drop" )
    if ( p2 ) cat(glue::glue("{pmssg} All Control Cross::{RET}"))
    if ( p2 ) all_cross_sum %>% print( n=base::nrow(all_cross_sum) )
    if ( p2 ) cat(glue::glue("{RET}"))
    
    #
    # LEFT OFF HERE!
    #
    #. TBD:: Loop over each category and make new name assignments::
    #
    # example_gs_ctl_tib %>% print(n=1000)
    
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
    
    if ( FALSE ) {
      viz_inf1_tib <- NULL
      viz_inf1_tib <- seqs_ctl_tib %>% 
        dplyr::filter( Probe_ID_match %>% stringr::str_detect("_I_") ) %>% 
        dplyr::select( -Probe_ID_bsp, -Is_EPIC_match, -Sample_Dependent_full, 
                       -Control_Type_full, -Control_Name_full, 
                       -Control_Group_Str_full, -Order_ID_full, 
                       -Probe_Type_epic ) %>% 
        dplyr::arrange(Probe_ID_full)
      
      viz_inf1_sum <- viz_inf1_tib %>% 
        dplyr::group_by( Control_Type_match ) %>%
        dplyr::summarise( Count=n(), .groups = "drop" )
      if ( p2 ) viz_inf1_sum %>% print( n=base::nrow(viz_inf1_sum) )
      
      # Data Frame View::
      viz_inf1_dtf <- NULL
      viz_inf1_dtf <- viz_inf1_tib %>% as.data.frame()
      
      #
      # Non-Polymorphic Search::
      #
      viz_nonp_tib <- NULL
      viz_nonp_tib <- seqs_ctl_tib %>% 
        dplyr::filter( Control_Type_match %>% stringr::str_starts("NonPoly") ) %>% 
        dplyr::select( -Probe_ID_bsp, -Is_EPIC_match, -Sample_Dependent_full, 
                       -Control_Type_full, -Control_Name_full, 
                       -Control_Group_Str_full, -Order_ID_full, 
                       -Probe_Type_epic ) %>% 
        dplyr::arrange(Probe_ID_full)
      if ( p2 ) viz_nonp_tib %>% print(n=base::nrow(viz_nonp_tib) )
      
      #
      # Infinium I Normalization Search::
      #
      viz_norm_tib <- NULL
      viz_norm_tib <- seqs_ctl_tib %>% 
        dplyr::filter( Control_Type_match %>% stringr::str_starts("Normalization") ) %>% 
        dplyr::select( -Probe_ID_bsp, -Is_EPIC_match, -Sample_Dependent_full, 
                       -Control_Type_full, -Control_Name_full, 
                       -Control_Group_Str_full, -Order_ID_full, 
                       -Probe_Type_epic ) %>% 
        dplyr::arrange(Probe_ID_full)
      if ( p2 ) viz_norm_tib %>% print(n=base::nrow(viz_norm_tib) )
      
      #
      # Infinium I SNP Search::
      #
      viz_snps_tib <- NULL
      viz_snps_tib <- seqs_ctl_tib %>% 
        dplyr::filter( Control_Type_match %>% stringr::str_starts("SNP") ) %>% 
        dplyr::select( -Probe_ID_bsp, -Is_EPIC_match, -Sample_Dependent_full, 
                       -Control_Type_full, -Control_Name_full, 
                       -Control_Group_Str_full, -Order_ID_full, 
                       -Probe_Type_epic ) %>% 
        dplyr::arrange(Probe_ID_full)
      if ( p2 ) viz_snps_tib %>% print(n=base::nrow(viz_snps_tib) )
    }

    if ( p1 ) cat(glue::glue("{pmssg} END: Splitting Infinium I Controls.{RET2}"))
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #
  #                            Join Probe Names::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  join_names <- TRUE
  join_names <- FALSE
  
  if ( join_names ) {

    #
    # TBD:: Merge names back using: seqs_ctl_tib
    #.  - Need to make Merge_Key (Probe_Seq_U)
    #
    join_dat_tib <- NULL
    join_dat_tib <- seqs_ctl_tib %>% 
      dplyr::full_join( bsp_tib, multiple = "all", 
                        by=c("Probe_ID_bsp"="Probe_ID") )
    
    # > join_dat_tib %>% dplyr::distinct( U,Probe_Seq )
    # A tibble: 1,971 × 2
    #
    # > join_dat_tib %>% dplyr::filter( Tag == "UM" )
    # A tibble: 1,023 × 26
    #
    # > join_dat_tib %>% dplyr::filter( Tag == "UM" ) %>% dplyr::distinct( U,Probe_Seq )
    # A tibble: 1,023 × 2
    #
    # > join_dat_tib %>% dplyr::filter( Tag != "UM" ) %>% dplyr::distinct( U,Probe_Seq )
    # A tibble: 48 × 2
    #
    # 1023 + 48 = 1071
    # 1971 - 1071 = 900
    #
    
    #
    # Build a summary of un-aligned probes::
    #
    join_mis_aln_sum <- NULL
    join_mis_aln_sum <- join_dat_tib %>% 
      dplyr::filter( Tag != "UM" ) %>% 
      dplyr::distinct( U,Probe_Seq, .keep_all = TRUE ) %>% 
      dplyr::group_by( Control_Type_match ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    
    join_dat_tib %>% dplyr::filter( Tag == "UM" ) %>% head() %>% as.data.frame()
    
    #
    # TBD:: Investigate Core Controls (Specific,Non-Poly,BS) to Negative Controls
    #
    join_sel_aln_sum <- join_dat_tib %>% 
      dplyr::filter( Control_Type_match == "BsConversion_I" | 
                       Control_Type_match == "BsConversion_II" | 
                       Control_Type_match == "NonSpecific_I" | 
                       Control_Type_match == "NonSpecific_II" )

    #
    # TBD:: Compare non-Standard Controls (Specific,Non-Poly,BS) to Negative Controls
    #
    
    #
    # TBD:: Update This...
    #
    epic_ctl_tib %>% dplyr::filter( Annotation %>% stringr::str_starts("BISULFITE") )
    
    all_ctls_tib %>% dplyr::left_join( epic_ctl_tib, by=c("Address"="U") )
    
    # all_ctls_tib %>% dplyr::filter( epic_ctl_tib )
    # epic_ctl_tib %>% dplyr::filter( !Name %in% all_ctls_tib$Probe_ID )
  }
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                    Pre-processing:: Manifest (MSA03)
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

run_msa <- TRUE
run_msa <- FALSE

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
                   Probe_ID = Probe_ID %>% stringr::str_replace("^ctl_neg", "ctl_Negative_"),
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
                   SNP_Probe_ID )
  
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
  # *  2 EXTENSION              
  # *  3 HYBRIDIZATION
  #
  #  4 TARGET REMOVAL         
  #  5 BISULFITE CONVERSION I 
  #  6 BISULFITE CONVERSION II
  #  7 SPECIFICITY I          
  #  8 SPECIFICITY II         
  #  9 NON-POLYMORPHIC        
  # 10 NEGATIVE               
  # 11 RESTORATION            
  # 12 NORM_T                 
  # 13 NORM_G                 
  # 14 NORM_C                 
  # 15 NORM_A
  #
  # + 16 Stringency             
  # + 17 Non-Specific Binding   
  
  # Need to Include Infinium Controls:
  #
  
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
      )
    ) %>% dplyr::select( Address,Type,Color_Channel,Name,
                         Probe_ID, AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
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
      )
    ) %>% dplyr::select( Address,Type,Color_Channel,Name,
                         Probe_ID, AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  # ctl_tib2 %>% as.data.frame()
  
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
      )
    ) %>% dplyr::select( Address,Type,Color_Channel,Name,
                         Probe_ID, AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  ctl_tib3 %>% as.data.frame()
  
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
      )
    ) %>% dplyr::select( Address,Type,Color_Channel,Name,
                         Probe_ID, AlleleA_ProbeSeq,AlleleB_ProbeSeq )
  
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
  
  #  6 BISULFITE CONVERSION II
  #  7 SPECIFICITY I          
  #  8 SPECIFICITY II         
  #  9 NON-POLYMORPHIC        
  # 10 NEGATIVE               
  # 11 RESTORATION            
  # 12 NORM_T                 
  # 13 NORM_G                 
  # 14 NORM_C                 
  # 15 NORM_A
  
  #  4 TARGET REMOVAL
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
  
  #
  #  5 BISULFITE CONVERSION I 
  #
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
  
  #
  # This is where we need to do some splitting hacks...
  #
  gs_ctl_tib5 <- NULL
  gs_ctl_tib5 <- msa_man_tib %>% 
    # dplyr::filter( Probe_Type %>% stringr::str_starts("BISULFITE") ) %>% 
    dplyr::filter( Probe_Type == "BISULFITE CONVERSION I" ) %>% 
    dplyr::mutate( Type = "BISULFITE CONVERSION I" ) %>% 
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
  #
  #                      MORE OLD CODE SECTION BELOW::
  #
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  msa_gsU_tib <- NULL
  msa_gsU_tib <- example_gs_ctl_tib %>% 
    dplyr::inner_join( msa_join_tib, by=c("Address"="U"), 
                       suffix=c("_gst","_msa") )
  msa_gsM_tib <- NULL
  msa_gsM_tib <- example_gs_ctl_tib %>% 
    dplyr::inner_join( msa_join_tib, by=c("Address"="M"), 
                       suffix=c("_gst","_msa") )
  
  msa_gsU_tib %>% head() %>% as.data.frame()
  
  
  #
  # Name Mismatch::
  #
  mis_gsU_tib <- NULL
  mis_gsU_tib <- msa_gsU_tib %>% 
    dplyr::mutate(
      Probe_Type = Probe_Type %>% stringr::str_to_upper()
    ) %>%
    dplyr::select( Address,Type,Color_Channel,Name_gst, Probe_Type ) %>% 
    dplyr::filter( Type != Probe_Type )

  mis_gsM_tib <- NULL
  mis_gsM_tib <- msa_gsM_tib %>% 
    dplyr::mutate(
      Probe_Type = Probe_Type %>% stringr::str_to_upper()
    ) %>%
    dplyr::select( Address,Type,Color_Channel,Name_gst, Probe_Type ) %>% 
    dplyr::filter( Type != Probe_Type )
  
  # Target Removal Just needs the " [12]$" Removed...
  msa_gsU_tib %>% 
    dplyr::mutate(
      Probe_Type = Probe_Type %>% stringr::str_to_upper()
    ) %>%
    dplyr::filter( Probe_Type %>% stringr::str_starts("TARGET") )

  #
  # TBD:: Investigate Infinium Controls::
  #
  msa_man_tib %>% dplyr::filter( AlleleA_ProbeSeq == "GTTGACGCTATCAAGGCTGCTGGTCACGACGGTAAGGTCAAGATCGGAAA" )
  
  #
  # Lots of old investigation code below::
  #

  msa_both_tib <- NULL
  msa_both_tib <- dplyr::bind_rows(
    dplyr::inner_join(
      msa_man_tib, dplyr::rename( all_class_tib, Address = U ),
      by=c("AlleleA_ProbeSeq"="Probe_Seq"),
      suffix=c("_msa","_ctl")
    ),
    dplyr::inner_join(
      msa_man_tib %>% dplyr::filter( !is.na(AlleleB_ProbeSeq) ), 
      dplyr::rename( all_class_tib, Address = U ),
      by=c("AlleleB_ProbeSeq"="Probe_Seq"),
      suffix=c("_msa","_ctl")
    )
  )
  
  all_class_tib %>% dplyr::filter( !U %in% msa_join_tib$Address )
  all_class_tib %>% dplyr::filter( !U %in% msa_both_tib$Address )
  
  
  "ATATACACTCACACATACACACAAACACATATACATACACTCTATCTTTA"
  # 29            Non_Specific_I_30 10740401 AACTTAAATTTAATAAATAATATTTTAACATACATACTATCTACTCCATG         G      NonSpecific_I                  N1          TRUE       ctl_Negative_855
  msa_join_tib %>% dplyr::filter( AlleleB_ProbeSeq == "AACTTAAATTTAATAAATAATATTTTAACATACATACTATCTACTCCATG" ) %>% as.data.frame()
  
  msa_join_lst <- msa_join_tib %>% split( .$Control_Class_match )
  msa_join_lst[["B1"]] %>% as.data.frame()

  #
  # B1:: Bisulfite Conversion I
  #
  msa_B1_tib <- NULL
  msa_B1_tib <- dplyr::inner_join(
    msa_man_tib,
    dplyr::select(join_B1_tib, Probe_ID,Probe_Seq_U,Probe_Seq_M,Prb_Idx),
    by=c("AlleleA_ProbeSeq"="Probe_Seq_M", "AlleleB_ProbeSeq"="Probe_Seq_U" ),
    multiple = "all" ) %>% 
    dplyr::select( Probe_ID,Prb_Idx, dplyr::everything() ) %>% 
    dplyr::arrange( Prb_Idx )
  
  # Not all controls have 5x after AQP::
  msa_B1_cnt_tib <- msa_B1_tib %>% 
    dplyr::group_by( Prb_Idx ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )


  #
  # Older Stuff Below::
  #
  mat_ctl_tib <- NULL
  mat_ctl_tib <- msa_man_tib %>% 
    dplyr::inner_join( all_ctls_tib, 
                       by=c("AlleleA_ProbeSeq"="Sequence") )
  
  msa_ctl_sum <- NULL
  msa_ctl_sum <- msa_man_tib %>% 
    dplyr::filter( !Probe_Type == "cg" ) %>% 
    dplyr::filter( !Probe_Type == "ch" ) %>% 
    dplyr::filter( !Probe_Type == "rs" ) %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) msa_ctl_sum %>% print( n=base::nrow(msa_ctl_sum) )
  
  # CODE BELOW SEEMS TO BE OUT OF DATE::
  #
  if ( FALSE ) {
    all_ctl_sum <- NULL
    all_ctl_sum <- all_man_tib %>% 
      dplyr::filter( Manifest == "ctl" | Manifest == "ctlB" ) %>% 
      dplyr::group_by( Annotation ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    if ( p2 ) all_ctl_sum %>% print( n=base::nrow(all_ctl_sum) )
    
    unq_ctl_tib <- NULL
    unq_ctl_tib <- all_man_tib %>% 
      dplyr::filter( Manifest == "ctl" | Manifest == "ctlB" ) %>% 
      dplyr::distinct(U, .keep_all = TRUE ) 
  }
  
  # msa_man_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% dplyr::filter( Manifest != "Embark" ) %>% dplyr::filter( Manifest != "Evonik" ) 
  # all_man_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% dplyr::filter( Manifest != "Embark" ) %>% dplyr::filter( Manifest != "Evonik" ) %>% dplyr::filter( Annotation == "NEGATIVE" )
  # tmp_tib %>% dplyr::filter( Probe_Type == "NEGATIVE" ) %>% dplyr::filter( !Probe_ID %>% stringr::str_starts("ZZneg") ) %>% dplyr::pull( Probe_ID )
  # all_man_tib %>% dplyr::filter( Manifest == "ctl" & Manifest_Version == "B4" ) %>% pull(Name)
  #
  # tmp_tib %>% dplyr::filter( Probe_Type == "NEGATIVE" ) %>% dplyr::filter( !Probe_ID %>% stringr::str_starts("ZZneg") ) %>% dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_remove("^ZZ"), din=Probe_ID %>% stringr::str_sub(1,2) ) %>% dplyr::group_by( din ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
  
  
  msa_ann_sum <- NULL
  msa_ann_sum <- msa_man_tib %>% 
    dplyr::group_by( Annotation ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) msa_ann_sum %>% print( n=base::nrow(msa_ann_sum) )
  
  msa_din_sum <- NULL
  msa_din_sum <- msa_man_tib %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) msa_din_sum %>% print( n=base::nrow(msa_din_sum) )
  
  neg_msa_tib <- NULL
  neg_msa_tib <- msa_man_tib %>% 
    dplyr::filter( Probe_Type == "NEGATIVE" ) %>% 
    dplyr::mutate( Probe_ID = Probe_ID %>% stringr::str_replace("^ctl_Negative", "ctl_neg"),
                   Loci_ID = Probe_ID %>% stringr::str_remove("_[0-9]+$") ) %>%
    dplyr::arrange( Probe_ID ) %>% 
    dplyr::left_join( neg_ctl_tib %>% dplyr::select( Probe_ID,Sequence ), 
                      by=c("Loci_ID"="Probe_ID") ) %>%
    dplyr::rename( Address = U ) %>% 
    dplyr::select( Address, Probe_ID, Sequence )
  
  #
  # Make sure negatives are unique...
  #
  # neg_msa_tib %>% dplyr::filter( Address %in% neg_ctl_tib$Address )
  # neg_msa_tib %>% dplyr::filter( Address %in% all_man_tib$U )
  
  # TBD:: Temporary Fix for MSA03 Negative Controls::
  #
  neg_ctl_tib <- neg_msa_tib
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
