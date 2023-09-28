
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                   Scratch Space for Docker Workhorse::
#                      Finger Printing Ivestigation
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
suppressWarnings(suppressPackageStartupMessages( 
  base::require("sesame",  quietly = TRUE) ) )

# install.packages("VennDiagram")
suppressWarnings(suppressPackageStartupMessages( 
  base::require("VennDiagram",  quietly = TRUE) ) )
# library("VennDiagram")

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
par$prgm_tag <- 'stable_forSteven_EPICv2_docker_conversions'
par$verbose  <- 3
local_paths  <- c(
  "/Users/bbarnes/Documents/tools/imSuite/scripts/R",
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

par$run_name <- "EPICv2_MVP"

# Multiple copies of v0 (bk and bk24)
par$version <- 0

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

# opt$pdf_path = safe_mkdir( dir = file.path( opt$out_path,"pdf") )
# opt$csv_path = safe_mkdir( dir = file.path( opt$out_path,"csv") )
# opt$rds_path = safe_mkdir( dir = file.path( opt$out_path,"rds") )
opt$snv_path = safe_mkdir( dir = file.path( opt$out_path,"snv") )

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#
#                           Quick RDS Conversion::
#
# This is to provide Steven and team with a text version of the annoI/S
#
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

if ( FALSE ) {

  convert_rds_to_csv <- TRUE
  convert_rds_to_csv <- FALSE
  
  if ( convert_rds_to_csv ) {
    dat_dir <- file.path( "/Users/bbarnes/Documents/Projects/EPIC_v2/docker/images/Infinium_Methylation_Workhorse_Centos.v.1.11.15.1.p.0.6.2/Infinium_Methylation_Workhorse/dat/dbSNP/b151" )
    
    annoI_b4_rds <- file.path( dat_dir, "EPIC-B4/EPIC-B4.annoI.rds" )
    annoS_b4_rds <- file.path( dat_dir, "EPIC-B4/EPIC-B4.annoS.rds" )
    annoI_b4_csv <- file.path( dat_dir, "EPIC-B4/EPIC-B4.annoI.csv" )
    annoS_b4_csv <- file.path( dat_dir, "EPIC-B4/EPIC-B4.annoS.csv" )
    
    annoI_b4_tib <- annoI_b4_rds %>% readr::read_rds() %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
    annoS_b4_tib <- annoS_b4_rds %>% readr::read_rds() %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
    readr::write_csv( x = annoI_b4_tib, file = annoI_b4_csv )
    readr::write_csv( x = annoS_b4_tib, file = annoS_b4_csv )
    
    annoI_a1_rds <- file.path( dat_dir, "EPIC-A1/EPIC-A1.annoI.rds" )
    annoS_a1_rds <- file.path( dat_dir, "EPIC-A1/EPIC-A1.annoS.rds" )
    annoI_a1_csv <- file.path( dat_dir, "EPIC-A1/EPIC-A1.annoI.csv" )
    annoS_a1_csv <- file.path( dat_dir, "EPIC-A1/EPIC-A1.annoS.csv" )
    
    annoI_a1_tib <- annoI_a1_rds %>% readr::read_rds() %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
    annoS_a1_tib <- annoS_a1_rds %>% readr::read_rds() %>% as.data.frame() %>% tibble::rownames_to_column( var = "Probe_ID" ) %>% tibble::as_tibble()
    readr::write_csv( x = annoI_a1_tib, file = annoI_a1_csv )
    readr::write_csv( x = annoS_a1_tib, file = annoS_a1_csv )
  }

  covnert_csv_to_rds <- FALSE
  covnert_csv_to_rds <- TRUE
  if ( covnert_csv_to_rds ) {
    #
    # Convert new annoI and annoS CSV to RDS
    #
    new_dat_dir <- file.path( opt$top_path, "Projects/EPIC_v2/docker/manifest/EPICv2/steven-test-swap" )
    
    annoI_a1_csv <- file.path( new_dat_dir, "EPIC-A1/EPIC-A1.annoI.csv.gz" )
    annoS_a1_csv <- file.path( new_dat_dir, "EPIC-A1/EPIC-A1.annoS.csv.gz" )
    
    annoI_a1_rds <- file.path( new_dat_dir, "EPIC-A1/EPIC-A1.annoI.rds" )
    annoS_a1_rds <- file.path( new_dat_dir, "EPIC-A1/EPIC-A1.annoS.rds" )

    # Load CSV::
    annoI_a1_tib <- NULL
    annoI_a1_tib <- readr::read_csv( file = annoI_a1_csv, show_col_types = FALSE )
    annoS_a1_tib <- NULL
    annoI_a1_tib <- readr::read_csv( file = annoS_a1_csv, show_col_types = FALSE )
    
    # Convert To Genomic Range::
    annoI_a1_grs <- NULL
    annoI_a1_grs <- GenomicRanges::GRanges( 
      seqnames = Rle( annoI_a1_tib$seqnames ),
      strand   = Rle( annoI_a1_tib$strand ),
      
      designType = annoI_a1_tib$designType,
      In.band    = annoI_a1_tib$In.band,
      REF = annoI_a1_tib$REF,
      ALT = annoI_a1_tib$ALT,
      rs  = annoI_a1_tib$rs,
      
      IRanges(start = annoI_a1_tib$start, 
              width = 1,
              # end   = annoI_a1_tib$end, 
              names = annoI_a1_tib$Probe_ID )
    )
    
    annoS_a1_grs <- NULL
    annoS_a1_grs <- GenomicRanges::GRanges( 
      seqnames = Rle( annoS_a1_tib$seqnames ),
      strand   = Rle( annoS_a1_tib$strand ),
      
      rs  = annoS_a1_tib$rs,
      designType = annoS_a1_tib$designType,
      U   = annoS_a1_tib$U,
      REF = annoS_a1_tib$REF,
      ALT = annoS_a1_tib$ALT,
      
      IRanges(start = annoS_a1_tib$start, 
              width = 1,
              # end   = annoS_a1_tib$end, 
              names = annoS_a1_tib$Probe_ID )
    )
    
    # Save RDS::
    readr::write_rds( x = annoI_a1_grs, file = annoI_a1_rds )
    readr::write_rds( x = annoS_a1_grs, file = annoS_a1_rds )
    
    # Copy to Docker::
    # docker cp /Users/bbarnes/Documents/Projects/EPIC_v2/docker/manifest/EPICv2/steven-test-swap/EPIC-A1/EPIC-A1.annoI.rds 04253273b6cf:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-A1/
    # docker cp /Users/bbarnes/Documents/Projects/EPIC_v2/docker/manifest/EPICv2/steven-test-swap/EPIC-A1/EPIC-A1.annoS.rds 04253273b6cf:/repo/Infinium_Methylation_Workhorse/dat/dbSNP/b151/EPIC-A1/

    # Test Docker::
    
  }
  
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  #                           Final Manifest Swap::
  # ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
  
  final_man_swap <- FALSE
  final_man_swap <- TRUE
  
  if ( final_man_swap ) {
    
    doc_dir <- "/Users/bbarnes/Documents/Projects/EPIC_v2/docker"

    new_man_dir <- file.path( doc_dir, "dat/EPICv2_VA_MVP-Data-2023_09_05" )
    new_man_csv <- file.path( new_man_dir, "EPICv2_VA_MVP_A1-Builder.body-sub.csv" )
    new_man_tib <- readr::read_csv( file = new_man_csv, show_col_types = FALSE )
    
    doc_man_dir <- file.path( doc_dir, "images/Infinium_Methylation_Workhorse_Centos.v.1.11.15.1.p.0.6.2/Infinium_Methylation_Workhorse/dat/manifest/core" )
    doc_man_csv <- file.path( doc_man_dir, "EPIC-A1.manifest.sesame-base.cpg-sorted.csv.gz" )
    doc_man_tib <- readr::read_csv( file = doc_man_csv, show_col_types = FALSE )
    
    # new_man_tib: A tibble: 936,382 × 21
    # doc_man_tib: A tibble: 937,690 × 11
    # 937690 - 936382 = 1308
    
    dif_man_tib <- NULL
    dif_man_tib <- new_man_tib %>% 
      dplyr::anti_join( doc_man_tib, c("IlmnID"="Full_ID") ) 
    
    add_man_tib <- NULL
    add_man_tib <- dif_man_tib %>% 
      dplyr::rename( 
        Full_ID = IlmnID, 
        Probe_ID = Name, 
        U = AddressA_ID, 
        M = AddressB_ID,
        COLOR_CHANNEL = Color_Channel ) %>%
      dplyr::mutate( 
        DESIGN = dplyr::case_when(
          !is.na(U) &  is.na(M) ~ "II",
          !is.na(U) & !is.na(M) ~ "I",
          TRUE ~ NA_character_ ),
        Probe_Design = dplyr::case_when(
          DESIGN ==  "I" ~ 1.0,
          DESIGN == "II" ~ 2.0,
          TRUE ~ NA_real_ ) %>% as.integer(),
        COLOR_CHANNEL = dplyr::case_when(
          DESIGN == "II" ~ "Both",
          TRUE ~ COLOR_CHANNEL ),
        Probe_Source = "EPIC-A1"
      ) %>%
      dplyr::select( Probe_ID, M, U, DESIGN, COLOR_CHANNEL, col,
                     Probe_Type, Probe_Source, Next_Base, Probe_Design, Full_ID )
    
    
    fin_man_csv <- file.path( doc_man_dir, "EPIC-A2.manifest.sesame-base.cpg-sorted.csv.gz" )
    fin_man_tib <- NULL
    fin_man_tib <- dplyr::bind_rows( doc_man_tib,add_man_tib ) %>% 
      dplyr::arrange( Probe_ID )
    readr::write_csv( x = fin_man_tib, file = fin_man_csv )
                                   
  }

}

# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #
#                                Finished::
# ----- ----- ----- ----- ----- -----|----- ----- ----- ----- ----- ----- #

prgm_ret_val <- 
  program_done(opts=opt, pars=par, precision=3, vb=vb,vt=vt+3,tc=tc,tt=tt )

sysTime <- Sys.time()
if ( p0 ) cat(glue::glue("{pmssg} Finished(time={sysTime}); Success={success}.{RET2}"))

# End of file
