
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                      Script for Parsing Genomes
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
par$prgm_tag <- 'stable_msa_cph'
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

par$run_name <- "MSA"

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
#                              Load Raw Data::
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

opt$raw_cph_tsv <- file.path( opt$top_path, "Projects.new/MSA/cph/20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.tsv.gz" )
stopifnot( opt$raw_cph_tsv %>% file.exists() )

# Input format for cluster::
#  Locus_Name,Chromosome,Position_Start,Position_Stop,Reference_Allele,Alternate_Allele,Type

raw_cph_tib <- NULL
raw_cph_tib <- readr::read_tsv( file = opt$raw_cph_tsv, show_col_types = FALSE )

inp_cph_tib <- NULL
inp_cph_tib <- raw_cph_tib %>% 
  dplyr::rename( 
    Chromosome = chr,
    Position_Start = start,
    Position_Stop  = end ) %>%
  dplyr::mutate( 
    Position_Start = Position_Start %>% as.integer(),
    Position_Stop  = Position_Stop %>% as.integer(),
    Chromosome = Chromosome %>% stringr::str_remove("^chr"),
    Locus_Name=paste( "ch",Chromosome,Position_Start, sep="_" ),
    Reference_Allele = "",
    Alternate_Allele = "",
    Type = "ch" ) %>%
  dplyr::select( Locus_Name,
                 Chromosome,Position_Start,Position_Stop,
                 Reference_Allele,Alternate_Allele,Type )

opt$inp_cph_csv <- NULL
opt$inp_cph_csv <- file.path (opt$out_path, "ch.20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.v1.design_input.csv.gz" )
readr::write_csv( x = inp_cph_tib, file = opt$inp_cph_csv )

opt$inp_cph_csv <- NULL
opt$inp_cph_csv <- file.path (opt$out_path, "ch.20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.v1.design_input.csv" )
readr::write_csv( x = inp_cph_tib, file = opt$inp_cph_csv )

#
# Write BED File::
#

bed_cph_tib <- NULL
bed_cph_tib <- inp_cph_tib %>% 
  dplyr::mutate( chr = paste0("chr", Chromosome), 
                 beg = Position_Start - 60,
                 end = beg + 122 ) %>%
  dplyr::select( chr,beg,end,Locus_Name )

opt$inp_cph_bed <- NULL
opt$inp_cph_bed <- file.path (opt$out_path, "ch.20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.v1.design_input.bed" )
readr::write_tsv( x = bed_cph_tib, file = opt$inp_cph_bed, col_names = FALSE )

#
# Try a second more reduced format::
#
opt$inp_cph_csv2 <- NULL
opt$inp_cph_csv2 <- file.path (opt$out_path, "ch.20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.v2.design_input.csv.gz" )
inp_cph_tib2 <- NULL
inp_cph_tib2 <- inp_cph_tib %>% 
  dplyr::mutate( Locus_Name = "", Position_Stop = "" )
readr::write_csv( x = inp_cph_tib2, file = opt$inp_cph_csv2 )

opt$inp_cph_csv2 <- NULL
opt$inp_cph_csv2 <- file.path (opt$out_path, "ch.20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.v2.design_input.csv" )
readr::write_csv( x = inp_cph_tib2, file = opt$inp_cph_csv2 )

#
# May need a third try (v3)
#. ch.1.89724558F,1,89951970.0,89951970.0,,,ch
#
opt$inp_cph_csv3 <- NULL
opt$inp_cph_csv3 <- file.path (opt$out_path, "ch.20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.v3.design_input.csv.gz" )
inp_cph_tib3 <- NULL
inp_cph_tib3 <- inp_cph_tib %>% 
  dplyr::mutate( Locus_Name=paste( "ch",Chromosome,Position_Start, sep="." ),
                 Position_Stop = Position_Start )
readr::write_csv( x = inp_cph_tib3, file = opt$inp_cph_csv3 )

opt$inp_cph_csv3 <- NULL
opt$inp_cph_csv3 <- file.path (opt$out_path, "ch.20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.v3.design_input.csv" )
readr::write_csv( x = inp_cph_tib3, file = opt$inp_cph_csv3 )

# Cluster follow up rsynce command::
#. rsync -P -r scratch/stable_msa_cph/MSA-UCSC-v0 ussd-rnd:/illumina/scratch/darkmatter/Projects/MSA/scratch/stable_msa_cph/
#
# Silly follow up command; should just write the non-zipped version...
#. gzip -dc /illumina/scratch/darkmatter/Projects/MSA/scratch/stable_msa_cph/MSA-UCSC-v0/ch.20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.v2.design_input.csv.gz >  /illumina/scratch/darkmatter/Projects/MSA/scratch/stable_msa_cph/MSA-UCSC-v0/ch.20221221_Ecker_GeneBody_mCH_Markers_Top10Genes.v2.design_input.csv
#

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, vb=opt$verbose, tt=tt)

sysTime <- Sys.time()
cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
