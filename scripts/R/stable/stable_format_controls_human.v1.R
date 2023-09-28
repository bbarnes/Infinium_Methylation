
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
par$prgm_tag <- 'stable_format_controls_human.R'
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
#.   - Test on HSA 2020 controls and { cg, ch, rs }
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #


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

#
# Original Formatted File with Names::
#
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

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                                   BEG
#                    Controls:: Correct Proccess Below
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

format_ctls <- FALSE
if ( format_ctls ) {
  #
  # This is the correct stuff!!!
  #  [Done] Need to copy this function:: Workhorse-Unstained/scripts/R/trifecta/functions/control_functions.R
  #  - Need to create new addControlType() that follows this naming scheme full_ctl_tib
  #
  full_ctl_tib <- NULL
  full_ctl_csv <- file.path( opt$top_path, "tools/docker/Infinium_Methylation_Workhorse/dat/manifest/controls/Infinium_Methylation_Controls_1983_full.csv.gz" )
  full_ctl_tib <- readr::read_csv( file = full_ctl_csv, show_col_types = FALSE )
  
  full_ctl_sum <- NULL
  full_ctl_sum <- full_ctl_tib %>% 
    dplyr::group_by(Control_Group) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  full_ctl_sum %>% print( n=nrow(full_ctl_sum) )
  
  # epic_ctl_tib %>% 
  #   dplyr::left_join( full_ctl_tib %>% dplyr::select(Probe_ID,U,Probe_Seq), by=c("U") ) %>% 
  #   as.data.frame()
  
  #
  # Original Match File::
  #
  match_ctl_tib <- NULL
  match_ctl_tsv <- file.path( opt$top_path, "data/Controls/01152015_DarkMatterControls.unique.probe.match.tsv.gz" )
  match_ctl_tib <- readr::read_tsv( file = match_ctl_tsv, show_col_types = FALSE ) %>%
    dplyr::rename( Control_ID = probe_id,
                   Probe_Seq = sequence ) %>%
    dplyr::mutate( U = address_name %>% stringr::str_remove("^1") %>%
                     stringr::str_remove("^0+") %>%
                     as.integer(),
                   Control_ID = Control_ID %>% 
                     stringr::str_remove("_1_A") %>%
                     stringr::str_remove("_1_B") 
    ) %>% 
    dplyr::select( Control_ID, U, Probe_Seq )
  
  # SNP Check::
  match_ctl_snp_tib <- NULL
  match_ctl_snp_tib <- match_ctl_tib %>% 
    dplyr::filter( Control_ID %>% stringr::str_starts("rs") )
  
  match_ctl_tib2 <- NULL
  match_ctl_tib2 <- readr::read_tsv( file = match_ctl_tsv, show_col_types = FALSE ) %>%
    dplyr::rename( Probe_ID = probe_id,
                   Probe_Seq = sequence ) %>%
    dplyr::mutate( U = address_name %>% stringr::str_remove("^1") %>%
                     stringr::str_remove("^0+") %>%
                     as.integer(),
                   Probe_ID = Probe_ID %>% 
                     stringr::str_remove("_1_A") %>%
                     stringr::str_remove("_1_B") 
    ) %>% 
    dplyr::select( Probe_ID, U, Probe_Seq )
  
  match_ctl_tib3 <- NULL
  match_ctl_tib3 <- addControlType( tib = match_ctl_tib2 ) %>%
    dplyr::mutate( Is_EPIC = dplyr::case_when(
      U %in% epic_ctl_tib$U ~ TRUE,
      TRUE ~ FALSE )
    )
  
  # match_ctl_tib3 %>% dplyr::add_count( U, name = "U_Dup" ) %>% dplyr::filter( U_Dup != 1 ) %>% dplyr::arrange( U )
  
  match_ctl_sum3 <- NULL
  match_ctl_sum3 <- match_ctl_tib3 %>% 
    dplyr::group_by( Control_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  match_ctl_sum3 %>% print( n=nrow(match_ctl_sum3) )
  
  epic_ctl_sum3 <- NULL
  epic_ctl_sum3 <- match_ctl_tib3 %>% 
    dplyr::filter( Is_EPIC ) %>%
    dplyr::group_by( Control_Type ) %>% 
    dplyr::summarise( Count=n(), .groups = "drop" )
  epic_ctl_sum3 %>% print( n=nrow(epic_ctl_sum3) )
  
  # Need to add these
  full_ctl_tib %>% dplyr::filter( !U %in% match_ctl_tib3$U )
}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#
#                    Controls:: Correct Proccess Below
#                                   END
#
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #





all_ctls_tib <- NULL
all_ctls_tsv <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/Controls/Control_Sequences.2020.tsv.gz" )
all_ctls_tib <- readr::read_tsv( file = all_ctls_tsv, show_col_types = FALSE ) %>%
  dplyr::mutate( Address = Address %>% stringr::str_remove("^1") %>%
                   stringr::str_remove("^0+") %>%
                   as.integer() ) %>% 
  dplyr::rename( Control_Name = Probe_ID )

write_fasta <- FALSE
if ( write_fasta ) {
  epic_ctl_fas <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/Controls/EPIC_Control_Sequences.644.fas.gz" )
  all_ctls_tib %>% 
    dplyr::inner_join( epic_ctl_tib, by=c("Address"="U") ) %>% 
    dplyr::mutate( FAS = paste0(">",Probe_ID,"\n",Sequence) ) %>% 
    dplyr::pull(FAS) %>% as.vector() %>% 
    readr::write_lines( file = epic_ctl_fas )
  
  epic_ctl_seq_fas <- file.path( opt$top_path, "data/manifests/methylation/GenomeStudio/Controls/EPIC_Control_Sequences.644.seq.fas.gz" )
  all_ctls_tib %>% 
    dplyr::distinct( Sequence, .keep_all = TRUE ) %>%
    dplyr::inner_join( epic_ctl_tib, by=c("Address"="U") ) %>% 
    dplyr::mutate( FAS = paste0(">",Sequence,"\n",Sequence) ) %>% 
    dplyr::pull(FAS) %>% as.vector() %>% 
    readr::write_lines( file = epic_ctl_seq_fas )
}

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

epic_ctl_tib %>% dplyr::filter( Annotation %>% stringr::str_starts("BISULFITE") )

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
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
