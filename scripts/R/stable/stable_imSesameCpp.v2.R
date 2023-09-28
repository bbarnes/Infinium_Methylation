
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

# par$version <- 0
# par$version <- 1
# par$version <- 2
# par$version <- 3
# par$version <- 4
# par$version <- 5
# par$version <- 6

# par$run_name <- "EPICv1"

par$run_name <- "Embarkv1"
par$run_name <- "FAILv1"
par$run_name <- "COREv1"

par$version <- 0
par$run_name <- "MSAv03"

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
neg_org_tib  <- readr::read_rds( neg_ctl_rds )

# print( neg_ctl_tib )

epic_ctl_tib <- NULL
epic_ctl_rds <- file.path( opt$top_path, "data/manifests/methylation/bgz/epic_ctls.rds" )
epic_ctl_tib <- readr::read_rds( file = epic_ctl_rds )

all_ctls_tib <- NULL
all_ctls_tsv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/Controls/Control_Sequences.2020.tsv.gz" )
all_ctls_tib <- readr::read_tsv( file = all_ctls_tsv, show_col_types = FALSE ) %>%
  dplyr::mutate( Address = Address %>% stringr::str_remove("^1") %>%
                   stringr::str_remove("^0+") %>%
                   as.integer() ) %>% 
  dplyr::rename( Control_Name = Probe_ID )

epic_ctl_fas <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/Controls/EPIC_Control_Sequences.644.fas.gz" )
all_ctls_tib %>% 
  dplyr::inner_join( epic_ctl_tib, by=c("Address"="U") ) %>% 
  dplyr::mutate( FAS = paste0(">",Probe_ID,"\n",Sequence) ) %>% 
  dplyr::pull(FAS) %>% as.vector() %>% 
  readr::write_lines( file = epic_ctl_fas )


all_ctls_tib %>% dplyr::left_join( epic_ctl_tib, by=c("Address"="U") )

# all_ctls_tib %>% dplyr::filter( epic_ctl_tib )
# epic_ctl_tib %>% dplyr::filter( !Name %in% all_ctls_tib$Probe_ID )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Manifest (MSA03)
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

msa_man_tib <- NULL
if ( par$run_name == "MSAv03" ) {
  #
  # Adding MSA03
  #
  msa_man_tib <- NULL
  msa_man_csv <- file.path( opt$top_path, "data/manifests/methylation/EX/MSAEX03/MSA-Interm-48v0-3_SS_BP123_A1.csv.gz" )
  msa_man_tib <- readr::read_csv( file = msa_man_csv, show_col_types = FALSE ) %>% 
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
  
  tmp_tib <- NULL
  tmp_tib <- readr::read_csv( file = msa_man_csv, show_col_types = FALSE )
  
  mat_ctl_tib <- msa_man_tib %>% 
    dplyr::inner_join( all_ctls_tib, by=c("AlleleA_ProbeSeq"="Sequence") )
  
  msa_man_tib %>% 
    dplyr::filter( !Probe_Type == "cg" ) %>% 
    dplyr::filter( !Probe_Type == "ch" ) %>% 
    dplyr::filter( !Probe_Type == "rs" ) %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
  
  all_ctl_tib <- all_man_tib %>% 
    dplyr::filter( Manifest == "ctl" | Manifest == "ctlB" ) %>% 
    dplyr::group_by( Annotation ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" ) %>% print(n=1000)
  
  unq_ctl_tib <- all_man_tib %>% 
    dplyr::filter( Manifest == "ctl" | Manifest == "ctlB" ) %>% 
    dplyr::distinct(U, .keep_all = TRUE ) 
  
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
  msa_ann_sum %>% print( n=base::nrow(msa_ann_sum) )
  
  msa_din_sum <- NULL
  msa_din_sum <- msa_man_tib %>% 
    dplyr::group_by( Probe_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  msa_din_sum %>% print( n=base::nrow(msa_din_sum) )
  
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
#                    Pre-processing:: Manifest
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

#
# Combine manifests
#
all_man_tib  <- NULL
all_man_tib  <- readr::read_rds( all_sub_rds ) %>% 
  dplyr::bind_rows( msa_man_tib )

# A tibble: 940,443 × 14   [Without MSA03]
# A tibble: 1,018,281 × 14 [With MSA03]

# all_man_tib %>% dplyr::group_by( col ) %>% dplyr::summarise( Count=n() )

# print( all_man_tib )
# all_man_tib %>% dplyr::group_by( Species, Manifest, Manifest_Version ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
# all_man_tib %>% dplyr::filter( Chromosome == "0" ) %>% print()

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Sample Sheet
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
ssh_tib <- NULL
sel_tib <- NULL

#
# TBD:: Switch between COREv1 & FAILv1
#
if ( par$run_name == "MSAv03" ) {
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
  
  sin_tib <- NULL
  sin_tib <- ssh_rep_tib %>% 
    dplyr::mutate( Select_Group := par$run_name )
  
} else if ( par$run_name == "EPICv1" || 
     par$run_name == "COREv1" ||
     par$run_name == "FAILv1"
     ) {
  
  # TBD:: Unsupervised or "Suggested Suppervised" Sample Cluster...
  #  - Set Sample_Name from predicted Auto_Sample_BETA_bDB_Key for now...
  #
  ssh_tib <- NULL
  ssh_csv <- file.path( opt$top_path, "data/sampleSheets/EPIC/EPIC_Auto_SampleSheets_2019-12-29.csv.gz" )
  ssh_tib <- readr::read_csv( file = ssh_csv, show_col_types = FALSE ) %>%
    dplyr::mutate( Sentrix_Path = paste0( opt$top_path,"/data/idats/EPIC/",Sentrix_Barcode,"/",Sentrix_Name ) ) %>%
    dplyr::mutate( Sample_Name = Auto_Sample_BETA_bDB_Key,
                   Sample_ID = Sample_Name,
                   
                   SentrixBarcode_A = Sentrix_Barcode,
                   SentrixPosition_A = Sentrix_Poscode,
                   DNA = "CellLine",
                   Format = ChipFormat,
                   Test = "CellLine",
                   'Scan Setting' = "MethylationNXT",
                   LAST = "LAST"
                   )
  
  #
  # NOTE: Only Using Replicate Samples for now...
  #
  sel_tib <- NULL
  sel_tib <- ssh_tib %>% 
    dplyr::filter( Sample_Name == "HELA" |
                     Sample_Name == "JURKAT" |
                     Sample_Name == "MCF7" |
                     Sample_Name == "RAJI" )
  
  if ( par$run_name == "EPICv1" ) sel_tib <- sel_tib %>% 
    dplyr::filter( PassDetpPoob_Percent_CG >= 90 )
  
  # Sample_ID           SentrixBarcode_A SentrixPosition_A DNA           Format Test          `Scan Setting`
  # <chr>                          <dbl> <chr>             <chr>         <chr>  <chr>         <chr>         
  #   1 206203800149_R01C01     206203800149 R01C01            CellLine      8x1    CellLine      MethylationNXT
  sel_col_vec <- c( "Sample_ID","SentrixBarcode_A","SentrixPosition_A","DNA","Format","Test","Scan Setting" )
  
  #
  # LEFT OFF HERE!!!
  #
  # ssh_tib %>% head(n=2) %>% dplyr::select( "Sample_ID","SentrixBarcode_A","SentrixPosition_A","DNA","Format","Test","Scan Setting" )
  # sel_tib %>% head(n=2) %>% dplyr::select( "Sample_ID","SentrixBarcode_A","SentrixPosition_A","DNA","Format","Test","Scan Setting" )

  if ( par$run_name == "COREv1" || par$run_name == "FAILv1" ) {
    #
    # Get all EPIC v1 DELTA Core/Fail Samples Seach::
    #
    if ( par$run_name == "COREv1" ) 
      epic_path <- file.path( opt$top_path, "data/idats/idats_EPIC-8x1-DELTA-Core" )
    if ( par$run_name == "FAILv1" ) 
      epic_path <- file.path( opt$top_path, "data/idats/idats_EPIC-8x1-DELTA-Fail" )
    
    epic_list <- NULL
    epic_list <- sesame::searchIDATprefixes( dir.name = epic_path, recursive = TRUE )
    epic_cnts <- epic_list %>% length()
    if ( p1 ) cat(glue::glue("{pmssg} epic_cnts({par$run_name}) = {epic_cnts}.{RET}"))

    sin_tib <- NULL
    sin_tib <- sel_tib %>% 
      dplyr::filter( Sentrix_Name %in% names(epic_list) ) %>%
      dplyr::mutate( Select_Group := par$run_name )
  }

  #
  # Load Sesame Masked Probes
  #
  epic_ses_dat <- NULL
  epic_ses_dat <- sesameData::sesameDataGet( title = "EPIC.probeInfo" )
  
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
  
  sin_tib <- NULL
  sin_tib <- sel_tib %>% 
    # dplyr::filter( Sentrix_Name %in% names(epic_list) ) %>%
    dplyr::mutate( Select_Group := par$run_name )

} else {
  stop( glue::glue("{perrs} Failed to find Sample Sheets! Exiting...{RET}") )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Sample Date Summary
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
if ( par$run_name == "EPICv1" ||
     par$run_name == "COREv1" ||
     par$run_name == "FAILv1" ) {

  date_decode_ssh_sum <- NULL
  date_decode_ssh_sum <- ssh_tib %>% 
    dplyr::group_by( iscan_Decoding_Year ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) date_decode_ssh_sum %>% print( n=base::nrow(date_decode_ssh_sum) )
  
  date_decode_sin_sum <- NULL
  date_decode_sin_sum <- sin_tib %>% 
    dplyr::group_by( iscan_Decoding_Year ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  if ( p4 ) date_decode_sel_sum %>% print( n=base::nrow(date_decode_sel_sum) )
  
  # date_decode_sel_sum <- NULL
  # date_decode_sel_sum <- sel_tib %>% 
  #   dplyr::group_by( iscan_Decoding_Year ) %>% 
  #   dplyr::summarise( Count=n(), .groups = "drop" )
  # if ( p4 ) date_decode_sel_sum %>% print( n=base::nrow(date_decode_sel_sum) )
  
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Pre-processing:: Sample Count Summary
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( par$run_name != "MSAv03" ) {
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
}

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

if ( par$run_name == "COREv1" ||
     par$run_name == "FAILv1" ) {
  
  sin_sum <- NULL
  sin_sum <- sin_tib %>%
    dplyr::group_by( Sample_Name,Select_Group ) %>% 
    dplyr::summarise( Tot_Poob=n(),
                      Min_Poob=min( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Max_Poob=max( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Avg_Poob=mean( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Sds_Poob=sd( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Med_Poob=median( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      Mad_Poob=mad( PassDetpPoob_Percent_CG, na.rm = TRUE ),
                      .groups = "drop" )
  if ( p1 ) sin_sum %>% print( n=base::nrow(sin_sum) )
  
  # sel_sum <- NULL
  # sel_sum <- sel_tib %>%
  #   dplyr::group_by( Sample_Name,Select_Group ) %>% 
  #   dplyr::summarise( Tot_Poob=n(),
  #                     Min_Poob=min( PassDetpPoob_Percent_CG, na.rm = TRUE ),
  #                     Max_Poob=max( PassDetpPoob_Percent_CG, na.rm = TRUE ),
  #                     Avg_Poob=mean( PassDetpPoob_Percent_CG, na.rm = TRUE ),
  #                     Sds_Poob=sd( PassDetpPoob_Percent_CG, na.rm = TRUE ),
  #                     Med_Poob=median( PassDetpPoob_Percent_CG, na.rm = TRUE ),
  #                     Mad_Poob=mad( PassDetpPoob_Percent_CG, na.rm = TRUE ),
  #                     .groups = "drop" )
  # if ( p1 ) sel_sum %>% print( n=base::nrow(sel_sum) )
  
  #
  # TBD:: Write Summary File...
  #
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

#
# TBD:: Scatter plots with random sub-selection::
#
opt$max_sam <- 9
opt$max_rep <- 9

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                    Data-processing:: Replicate Sample List
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
rep_data_tib <- NULL
rep_sel_list <- NULL
rep_sel_list <- sin_tib %>% split( .$Sample_Name ) %>%
  head( n=opt$max_sam )

rep_sel_list <- NULL
rep_sel_list <- sin_tib %>% split( .$Sample_Name )

for ( sample_name in names(rep_sel_list) ) {
  all_beta_tib <- NULL
  all_poob_tib <- NULL
  
  sample_name <- "HELA_200"
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                    Data-processing:: Sentrix Sample List
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  ssh_sel_list <- NULL
  ssh_sel_list <- rep_sel_list[[sample_name]] %>% 
    dplyr::mutate( Rep_Idx = dplyr::row_number(),
                   Rep_Key = paste( Sample_Name,Rep_Idx, sep="_") ) %>% 
    split( .$Sentrix_Name ) %>% head( n=opt$max_rep )

  sentrix_name <- names(ssh_sel_list)[1]
  
  for ( sentrix_name in names(ssh_sel_list) ) {
    prefix <- ssh_sel_list[[sentrix_name]]$Sentrix_Path
    if ( p1 ) cat(glue::glue("{pmssg} Current prefix = '{prefix}'{RET}"))
    
    # }
    # for ( prefix in pre_vec ) {
    #   if ( p1 ) cat(glue::glue("{pmssg} Current prefix = '{prefix}'{RET}"))
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                       Data-processing:: Loading
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    # prefix <- "/Users/bbarnes/Documents/data/idats/idats_MSAv03_48x1_Alpha/207545400008/207545400008_R04C01"
    # prefix <- "/Users/bbarnes/Documents/data/idats/EPIC/203010460029/203010460029_R12C01"
    # prefix <- "/Users/bbarnes/Documents/data/idats/EPIC/201502830033/201502830033_R02C01"
    # all_man_tib  <- NULL
    # all_man_tib  <- readr::read_rds( all_sub_rds ) # %>% dplyr::bind_rows( msa_man_tib )
    # all_man_tib %>% dplyr::filter( Probe_Type == "ct" ) %>% dplyr::group_by( Manifest, Manifest_Version ) %>% dplyr::summarise( Count=n(), .groups = "drop" )
    
    # all_man_tib %>% dplyr::filter( U %in% neg_ctl_tib$Address ) %>% tail()
    
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
    if ( p1 ) idat_sum %>% print( n=base::nrow(idat_sum) )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Data-processing:: Parse Analytical Data
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    data_tib <- NULL
    data_tib <- idat_tib %>% 
      dplyr::filter( !stringr::str_starts( Probe_ID, pattern = "ct") )
    
    data_sum <- NULL
    data_sum <- data_tib %>% 
      dplyr::group_by( mask ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    if ( p1 ) data_sum %>% print( n=base::nrow(data_sum) )
    
    data_col_sum <- NULL
    data_col_sum <- data_tib %>% 
      dplyr::group_by( col ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    if ( p1 ) data_col_sum %>% print( n=base::nrow(data_col_sum) )
    
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    #                   Data-processing:: Parse Control Data
    # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
    ctls_tib <- NULL
    ctls_tib <- idat_tib %>% 
      dplyr::filter(  stringr::str_starts( Probe_ID, pattern = "ct") )
    
    ctls_sum <- NULL
    ctls_sum <- ctls_tib %>% 
      dplyr::group_by( mask ) %>% 
      dplyr::summarise( Count=n(), .groups = "drop" )
    if ( p1 ) ctls_sum %>% print( n=base::nrow(ctls_sum) )
    
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
  # Skip Samples with lesss than two samples...
  if ( ncol_cnt < 2 ) next
  
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
      prd_mat <- cbind( prd_mat, cbind( dB_vec, dP_vec ) %>% matrixStats::rowProds() )
      min_mat <- cbind( min_mat, dP_vec )
      dbs_mat <- cbind( dbs_mat, dB_vec )
      
      # dB_pass_mat <- dB_pass_mat &
      #   beta_mat[ ,c(ii,jj) ] %>% matrixStats::rowDiffs() < min_dB
      
      # break
    }
    
    # break
  }
  
  dat_tib <- NULL
  sam_tib <- tibble::tibble(
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
    sam_tib %>% ggplot2::ggplot( aes( x=fin_scr ) ) +
      ggplot2::geom_density()
    
    sam_tib %>% ggplot2::ggplot( aes( x=fin_scr, y=avg_dbs ) ) +
      ggplot2::geom_point( )
    
    sam_tib %>% ggplot2::ggplot( aes( x=fin_scr, y=avg_dbs ) ) +
      ggplot2::geom_density2d()
    
    sam_tib %>% ggplot2::ggplot( aes( x=fin_scr, y=med_dbs ) ) +
      ggplot2::geom_density2d()
  }
  
  rep_data_tib <- rep_data_tib %>% dplyr::bind_rows( sam_tib )
  
  if ( opt$single ) break
}

rep_data_sum <- NULL
rep_data_sum <- rep_data_tib %>% 
  dplyr::group_by( Sample ) %>% 
  dplyr::summarise( Count=n(), .groups = "drop" )

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
    ggplot2::facet_grid( rows = vars(Sample) )
  
  ggplot2::ggsave( filename = rep_den_pdf, 
                   device = "pdf", 
                   width = 7, 
                   height = 7, 
                   dpi = 320 )
  
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
