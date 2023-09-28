
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
# par$version <- 1
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
  full_ctl_tib <- readr::read_csv( file = full_ctl_csv, show_col_types = FALSE )%>% 
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
                   Sample_Dependent_full = Sample_Dependent )
  
  full_ctl_sum <- NULL
  full_ctl_sum <- full_ctl_tib %>% 
    dplyr::group_by( Control_Type_full ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
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
    addControlTypeSlim() %>%
    dplyr::mutate( Is_EPIC = dplyr::case_when(
      U %in% epic_ctl_tib$U ~ TRUE,
      TRUE ~ FALSE )
    ) %>% 
    dplyr::rename( Probe_ID_match = Probe_ID, 
                   Control_Type_match = Control_Type, 
                   Is_EPIC_match = Is_EPIC )
  
  # SNP Check::
  match_ctl_snp_tib <- NULL
  match_ctl_snp_tib <- match_ctl_tib %>% 
    dplyr::filter( Probe_ID_match %>% stringr::str_starts("rs") )
  # match_ctl_snp_tib %>% print( n=base::nrow(match_ctl_snp_tib) )
  
  match_ctl_sum <- NULL
  match_ctl_sum <- match_ctl_tib %>% 
    dplyr::group_by( Control_Type_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) match_ctl_sum %>% print( n=nrow(match_ctl_sum) )
  
  epic_ctl_sum <- NULL
  epic_ctl_sum <- match_ctl_tib %>% 
    dplyr::filter( Is_EPIC_match ) %>%
    dplyr::group_by( Control_Type_match ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p2 ) epic_ctl_sum %>% print( n=nrow(epic_ctl_sum) )
  
  if ( p1 ) cat(glue::glue("{pmssg} END: Formatting Raw Controls.{RET2}"))
}

# Just Some Printing Sanity Checking Stuff
#
# match_ctl_tib %>% head() %>% as.data.frame()
# full_ctl_tib  %>% head() %>% as.data.frame()
# epic_ctl_tib  %>% head() %>% as.data.frame()

join_ctl_dir <- safe_mkdir( file.path( opt$out_path, "join") )
join_ctl_csv <- file.path( join_ctl_dir, "Infinium_Methylation_Controls_1971_join-names.raw.csv.gz" )
join_ctl_csv <- file.path( join_ctl_dir, "Infinium_Methylation_Controls_1971_join-names.raw.csv" )

join_ctl_tib <- NULL
join_ctl_tib <- match_ctl_tib %>% 
  dplyr::left_join( epic_ctl_tib, by=c("U") ) %>%
  dplyr::full_join( full_ctl_tib, by=c("U","Probe_Seq") ) %>%
  dplyr::arrange( Control_Type_match )

readr::write_csv( x = join_ctl_tib, file = join_ctl_csv )

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

#
# Next Steps::
#.  - Align all controls sequences using U+Sequence
#.  - Substring all Negative controls using Probe_ID_match => bisulfite convert => count ACTG's
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
# NOTE:: Don't Run below yet, we have everything we need above. 
# NOTE:: Need to save this as scratch and build a clean version...
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

# TBD:: These Parameters Should be Moved to Program Inputs::
#
ref_tar <- "hg19"
ref_tar <- "GRCh37"
bsp_exe <- file.path( opt$top_path, "tools/programs/BSMAPz/bsmapz" )
opt$bsmap_exe <- bsp_exe

# NOTE: Below is all just validation that we only need match_ctl_tib
#
if ( FALSE ) {
  #
  # Load All Raw Control Sequences::
  #
  all_ctls_tib <- NULL
  all_ctls_tsv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/Controls/Control_Sequences.2020.tsv.gz" )
  all_ctls_tib <- readr::read_tsv( file = all_ctls_tsv, show_col_types = FALSE ) %>%
    dplyr::mutate( Address = Address %>% 
                     stringr::str_remove("^1") %>%
                     stringr::str_remove("^0+") %>%
                     as.integer() ) %>% 
    dplyr::rename( Control_Name = Probe_ID ) %>%
    dplyr::filter( !Control_Name %>% stringr::str_ends("\\.2$"))
  
  #
  # NOTE: Evidence for the filter line: dplyr::filter( !Control_Name %>% stringr::str_ends("\\.2$"))
  #
  # all_ctls_tib %>% dplyr::distinct( Sequence, Address, .keep_all = TRUE )
  # all_ctls_tib %>% dplyr::distinct( Sequence, .keep_all = TRUE )
  # all_ctls_tib %>% dplyr::distinct( Address, .keep_all = TRUE )
  # 
  # all_ctls_tib %>% 
  #   dplyr::add_count( Sequence, name="Dup_Cnt" ) %>%
  #   dplyr::filter( Dup_Cnt != 1 )
  
  ann_ctls_tib <- NULL
  ann_ctls_tib <- all_ctls_tib %>% 
    dplyr::inner_join( match_ctl_tib, 
                       by=c("Address"="U","Sequence"="Probe_Seq") )
}

align_ctls <- FALSE
align_ctls <- TRUE

if ( align_ctls ) {

  run_bsp <- FALSE
  run_bsp <- TRUE

  if ( run_bsp ) {
    if ( p1 ) cat(glue::glue("{pmssg} BEG: Launching BSMAP Alignment...{RET}"))
    
    fas_path <- NULL
    bed_path <- NULL
    man_path <- NULL
    
    fas_path <- safe_mkdir( dir = file.path( opt$out_path, "fas" ) )
    bed_path <- safe_mkdir( dir = file.path( opt$out_path, "bed" ) )
    man_path <- safe_mkdir( dir = file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/Controls/tmp" ) )
    
    #
    # TBD:: Loop over ids && seqs...
    #
    # OLD Locations from dry run::
    #  - epic_ctl_ids_fas <- file.path( opt$out_path, "data/manifests/methylation/GenomeStudio/Controls/EPIC_Control_Sequences.644.ids.fa.gz" )
    #  - epic_ctl_seq_fas <- file.path( opt$out_path, "data/manifests/methylation/GenomeStudio/Controls/EPIC_Control_Sequences.644.seq.fa.gz" )
    #
    
    #
    # TBD:: 
    #  - Match 1971 with 644
    #  - Resolve the rest manually...
    #
    all_ctl_ids_tsv <- NULL
    all_ctl_ids_tsv <- file.path( man_path, "EPIC_Control_Sequences.1971-644.tmp.tsv.gz" )
    
    # Need to sort this out...
    tmp_ctls_tib <- all_ctls_tib %>% 
      dplyr::left_join( epic_ctl_tib, by=c("Address"="U") ) %>%
      dplyr::select( -M ) %>%
      dplyr::mutate(
        col = dplyr::case_when(
          is.na(col) ~ "1",
          TRUE ~ col
        ) %>% as.integer()
      )
    
    
    epic_ctl_ids_fas <- NULL
    epic_ctl_ids_fas <- file.path( fas_path, "EPIC_Control_Sequences.644.ids.fa.gz" )
    all_ctls_tib %>% 
      dplyr::inner_join( epic_ctl_tib, by=c("Address"="U") ) %>% 
      dplyr::mutate( FAS = paste0(">",Probe_ID,"\n",Sequence) ) %>% 
      dplyr::pull(FAS) %>% as.vector() %>% 
      readr::write_lines( file = epic_ctl_ids_fas )
    
    epic_ctl_seq_fas <- NULL
    epic_ctl_seq_fas <- file.path( fas_path, "EPIC_Control_Sequences.644.seq.fa.gz" )
    all_ctls_tib %>% 
      dplyr::distinct( Sequence, .keep_all = TRUE ) %>%
      dplyr::inner_join( epic_ctl_tib, by=c("Address"="U") ) %>% 
      dplyr::mutate( FAS = paste0(">",Sequence,"\n",Sequence) ) %>% 
      dplyr::pull(FAS) %>% as.vector() %>% 
      readr::write_lines( file = epic_ctl_seq_fas )

    #
    # Run BSMAP:: ids
    #
    ids_bsp_tsv <- NULL
    ids_bsp_tsv <- run_bsmap( ref_fas = ref_seqs[[ref_tar]],
                              can_fas = epic_ctl_ids_fas,
                              bsp_exe = opt$bsmap_exe,
                              bsp_dir = opt$bsmap_dir,
                              slim = TRUE,
                              out_dir = file.path( opt$out_path, "ids" ),
                              run_tag = opt$run_name,
                              reload = opt$reload,
                              reload_min = 10,
                              vb=vb,vt=vt,tc=tc,tt=tt )
    
    ids_bsp_tib <- NULL
    ids_bsp_tib <- load_bsmap( file = ids_bsp_tsv,
                               sort = FALSE,
                               add_cnt = TRUE,
                               parse_id = FALSE,
                               out_dir = file.path( opt$out_path, "ids" ),
                               run_tag = opt$run_name,
                               reload = opt$reload,
                               reload_min = 10,
                               vb=vb,vt=vt,tc=tc,tt=tt )
    
    #
    # Write BED::
    #
    ids_bsp_bed <- file.path( bed_path, "PIC_Control_Sequences.644.GRCh37.ids.slim.bs-only.bed" )
    ids_122_fas <- file.path( bed_path, "PIC_Control_Sequences.644.GRCh37.ids.slim.bs-only.122mer.fa" )
    ids_bsp_tib %>% dplyr::select( Chr,Beg,Probe_ID ) %>% 
      dplyr::mutate( Beg = Beg -60, End = Beg + 122 ) %>% 
      dplyr::select( Chr,Beg,End, Probe_ID ) %>%
      dplyr::arrange( Chr,Beg,End ) %>%
      readr::write_tsv( file = ids_bsp_bed, col_names = FALSE )
    
    #
    # Extract 122mer::
    #
    twoBitToFa_exe <- "/Users/bbarnes/Documents/tools/ucsc/twoBitToFa"
    twoBitToFa_cmd <- paste0( twoBitToFa_exe,
                              " -bed=",ids_bsp_bed,
                              paste0(" ",ref_seqs[[ref_tar]] %>% stringr::str_remove(".gz$"),".2bit"),
                              " ",ids_122_fas )
    base::system( command = twoBitToFa_cmd )
    
    ids_122_dat <- NULL
    ids_122_dat <- Biostrings::readDNAStringSet( filepath = ids_122_fas, format = "fasta" )
    
    if ( p1 ) cat(glue::glue("{pmssg} END: Launching BSMAP Alignment.{RET2}"))
  }

  pre_bsp <- TRUE
  pre_bsp <- FALSE
  
  if ( pre_bsp ) {
    #
    # Load BSMAP::
    #
    # gzip -dc bsmap/EPIC_Control_Sequences.644.GRCh37.bsp.gz | grep BS_Conversion | grep UM | cut -f 1,2,4-11
    #
    bsmap_ids_tib <- NULL
    bsmap_ids_bsp <- file.path( opt$top_path, "Projects.new/Evonik/universal_controls/bsmap/EPIC_Control_Sequences.644.GRCh37.bsp.gz" )
    bsmap_ids_bsp <- file.path( opt$top_path, "Projects.new/Evonik/universal_controls/bsmap/EPIC_Control_Sequences.644.GRCh37.ids.slim.bs-only.tsv.gz" )
    bsmap_ids_tib <- load_bsmap( file = bsmap_ids_bsp, 
                                 slim = TRUE, 
                                 out_dir = file.path( opt$out_path, "bsmap/ids" ),
                                 run_tag = opt$run_name,
                                 vb=vb,vt=vt,tc=tc,tt=tt )
    
    bsmap_ids_bed <- file.path( opt$out_path, "bsmap/ids/EPIC_Control_Sequences.644.GRCh37.ids.slim.bs-only.bed" )
    bsmap_ids_tib %>% dplyr::select( Chr,Beg,Probe_ID ) %>% 
      dplyr::mutate( Beg = Beg -60, End = Beg + 122 ) %>% 
      dplyr::select( Chr,Beg,End, Probe_ID ) %>%
      readr::write_tsv( file = bsmap_ids_bed, col_names = FALSE )
    
    # /Users/bbarnes/Documents/tools/twoBitToFa -bed=/Users/bbarnes/Documents/scratch/stable_MSA_imSesameCpp/MSAv03-UCSC-v0/bsmap/ids/EPIC_Control_Sequences.644.GRCh37.ids.slim.bs-only.bed /Users/bbarnes/Documents/data/imGenomes/Homo_sapiens/NCBI/GRCh37/Sequence/WholeGenomeFasta/GRCh37.genome.fa.2bit /Users/bbarnes/Documents/scratch/stable_MSA_imSesameCpp/MSAv03-UCSC-v0/bsmap/ids/EPIC_Control_Sequences.644.GRCh37.ids.slim.bs-only.122mer.fa
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

if ( par$run_name == "MSAv03" ) {
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
    dplyr::select( Probe_ID, U, M, col, Name, 
                   Species, Manifest, Manifest_Version, Annotation, 
                   Chromosome, Coordinate, Probe_Type, Locus_Name, 
                   SNP_Probe_ID, AlleleA_ProbeSeq )
  
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
