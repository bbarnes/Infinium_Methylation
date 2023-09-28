
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

success = TRUE;

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Load Sample Sheet::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$ssh_csv <- file.path( opt$top_path, "Projects.new/EPIC_v2/stable_epic_v2_manifest_plots/GSIBIOINFO-638-v1/build_epicv2_sample_sheet/GSIBIOINFO-638-v1.build_epicv2_sample_sheet.csv.gz" )
ssh_sel_csv <- file.path( opt$out_path, paste0( opt$run_name,".LightningAuto.select.sample_sheet.csv.gz" ) )
ssh_sel_tib <- NULL
ssh_sel_tib <- readr::read_csv( file = opt$ssh_csv, show_col_types = FALSE ) %>% 
  dplyr::filter( Sheet_Prep == "Lightning" & Sheet_Proc == "Auto" )
readr::write_csv( x = ssh_sel_tib, file = ssh_sel_csv )

# Methylation QC Pipeline Format::
#
# [Header],,,,,,
# Institute Name,ILMN,,,,,
# Investigator Name,SMG,,,,,
# Project Name,DemoDataEPICv2,,,,,
# Date,,,,,,
# ,,,,,,
# [Manifests],EPIC-8v2-0_A1.bpm,,,,,
# ,,,,,,
# ,,,,,,
# [Data],,,,,,
# Sample_ID,SentrixBarcode_A,SentrixPosition_A,DNA,Format,Test,Scan Setting
# 206891110004_R01C01,206891110004,R01C01,Cell Line,8x1,Control,MethylationNXT
# 206891110004_R07C01,206891110004,R07C01,Cell Line,8x1,Control,MethylationNXT
# 206891110004_R08C01,206891110004,R08C01,Coriell,8x1,Control,MethylationNXT
# ,,,,,,
# ,,,,,,
#

#
# Replicate Sample Sheet::
#
ssh_hed_vec <- NULL
ssh_hed_vec <- c( "[Header],,,,,,",
                  "Institute Name,ILMN,,,,,",
                  "Investigator Name,SMG,,,,,",
                  "Project Name,EPICv2_Replicate,,,,,",
                  "Date,,,,,,",
                  ",,,,,,",
                  "[Manifests],EPIC-8v2-0_A1.bpm,,,,,",
                  # "[Manifests],EPIC_v2_core.unmasked.manifest.v0.bpm,,,,,",
                  ",,,,,,",
                  ",,,,,,",
                  "[Data],,,,,,",
                  "Sample_ID,SentrixBarcode_A,SentrixPosition_A,DNA,Format,Test,Scan Setting" )

ssh_rep_csv <- file.path( opt$out_path, paste0( opt$run_name,".LightningAuto.replicate.sample_sheet.csv.gz" ) )
readr::write_lines( x = ssh_hed_vec, file = ssh_rep_csv )

ssh_rep_tib <- NULL
ssh_rep_tib <- ssh_sel_tib %>% 
  dplyr::filter( Sample_Group != "MeTritration" ) %>%
  dplyr::mutate( Sample_ID = Sentrix_Name, 
                 SentrixBarcode_A = Sample_ID %>% stringr::str_remove("_.*$"), 
                 SentrixPosition_A = Sample_ID %>% stringr::str_remove("^.*_"), 
                 DNA = Sample_Group, Format = "8x1", 
                 Test = Sample_Group,
                 # Test = dplyr::case_when(
                 #   Sample_Group == "MeTritration" ~ "Titration",
                 #   TRUE ~ "Replicate" ),
                 Scan_Setting = "MethylationNXT" ) %>%
  dplyr::select( Sample_ID,SentrixBarcode_A,SentrixPosition_A,
                 DNA,Format,Test,Scan_Setting )

readr::write_csv( x = ssh_rep_tib, file = ssh_rep_csv,
                  col_names = FALSE, 
                  append = TRUE )

#
# Titration Sample Sheet::
#
ssh_hed_vec <- NULL
ssh_hed_vec <- c( "[Header],,,,,,",
                  "Institute Name,ILMN,,,,,",
                  "Investigator Name,SMG,,,,,",
                  "Project Name,EPICv2_Titration,,,,,",
                  "Date,,,,,,",
                  ",,,,,,",
                  "[Manifests],EPIC-8v2-0_A1.bpm,,,,,",
                  # "[Manifests],EPIC_v2_core.unmasked.manifest.v0.bpm,,,,,",
                  ",,,,,,",
                  ",,,,,,",
                  "[Data],,,,,,",
                  "Sample_ID,SentrixBarcode_A,SentrixPosition_A,DNA,Format,Test,Scan Setting" )

ssh_tit_csv <- file.path( opt$out_path, paste0( opt$run_name,".LightningAuto.titration.sample_sheet.csv.gz" ) )
readr::write_lines( x = ssh_hed_vec, file = ssh_tit_csv )

ssh_tit_tib <- NULL
ssh_tit_tib <- ssh_sel_tib %>% 
  dplyr::filter( Sample_Group == "MeTritration" ) %>%
  dplyr::mutate( Sample_ID = Sentrix_Name, 
                 SentrixBarcode_A = Sample_ID %>% stringr::str_remove("_.*$"), 
                 SentrixPosition_A = Sample_ID %>% stringr::str_remove("^.*_"), 
                 DNA = Sample_Group, Format = "8x1", 
                 Test = Sample_Group,
                 # Test = dplyr::case_when(
                 #   Sample_Group == "MeTritration" ~ "Titration",
                 #   TRUE ~ "Replicate" ),
                 Scan_Setting = "MethylationNXT" ) %>%
  dplyr::select( Sample_ID,SentrixBarcode_A,SentrixPosition_A,
                 DNA,Format,Test,Scan_Setting )

readr::write_csv( x = ssh_tit_tib, file = ssh_tit_csv,
                  col_names = FALSE, 
                  append = TRUE )


# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                            Load Beta Values::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #



# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
